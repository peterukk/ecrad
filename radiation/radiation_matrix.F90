! radiation_matrix.F90 - SPARTACUS matrix operations
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2018-10-15  R. Hogan    Added fast_expm_exchange_[23]
!   2020-12-xx  P. Ukkonen  Added an optimized expm routine for shortwave when nreg=3,
!                           and related kernels
!
! This module provides the neccessary mathematical functions for the
! SPARTACUS radiation scheme: matrix multiplication, matrix solvers
! and matrix exponentiation, but (a) multiple matrices are operated on
! at once with array access indended to facilitate vectorization, and
! (b) optimization for 2x2 and 3x3 matrices.  There is probably
! considerable scope for further optimization. Note that this module
! is not used by the McICA solver.

module radiation_matrix

  use parkind1, only : jprb, jprd, jpim, jprm
#ifdef USE_TIMING
  !
  ! Timing library
  !
  use gptl,                  only: gptlstart, gptlstop, gptlinitialize, gptlpr, gptlfinalize, gptlsetoption, &
                                   gptlpercent, gptloverhead
#endif
#ifdef USE_LIBXSMM
#ifdef SINGLE_PRECISION
  USE :: LIBXSMM, ONLY: LIBXSMM_BLASINT_KIND,                     &
  &                        LIBXSMM_MMFUNCTION => LIBXSMM_SMMFUNCTION,&
  &                        libxsmm_mmdispatch => libxsmm_smmdispatch,&
  &                        libxsmm_mmcall => libxsmm_smmcall
#else
  USE :: LIBXSMM, ONLY: LIBXSMM_BLASINT_KIND,                     &
  &                        LIBXSMM_MMFUNCTION => LIBXDMM_SMMFUNCTION,&
  &                        libxsmm_mmdispatch => libxdmm_smmdispatch,&
  &                        libxsmm_mmcall => libxsmm_dmmcall  
#endif
#endif
  implicit none
  public

  ! Codes to describe sparseness pattern, where the SHORTWAVE
  ! pattern is of the form:
  ! (x x x)
  ! (x x x)
  ! (0 0 x)
  ! where each element may itself be a square matrix.  
  integer, parameter :: IMatrixPatternDense     = 0
  integer, parameter :: IMatrixPatternShortwave = 1
#ifdef USE_TIMING
  integer :: ret
#endif
  public  :: mat_x_vec, singlemat_x_vec, mat_x_mat, &
       &     singlemat_x_mat, mat_x_singlemat, &
       &     identity_minus_mat_x_mat, solve_vec, solve_mat, expm, &
       &     fast_expm_exchange_2, fast_expm_exchange_3

  private :: solve_vec_2, solve_vec_3, solve_mat_2, &
       &     solve_mat_3, lu_factorization, lu_substitution, solve_mat_n, &
       &     diag_mat_right_divide_3

  interface fast_expm_exchange
    module procedure fast_expm_exchange_2, fast_expm_exchange_3
  end interface fast_expm_exchange

contains

  ! --- MATRIX-VECTOR MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A as n m-by-m square matrices (with the n dimension varying
  ! fastest) and b as n m-element vectors, and perform matrix-vector
  ! multiplications on first iend pairs
  function mat_x_vec(n,iend,m,A,b,do_top_left_only_in)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                   :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:) :: A
    real(jprb), intent(in), dimension(:,:)   :: b
    logical,    intent(in), optional         :: do_top_left_only_in
    real(jprb),             dimension(iend,m):: mat_x_vec

    integer :: j1, j2
    logical :: do_top_left_only

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',0,hook_handle)


    if (present(do_top_left_only_in)) then
      do_top_left_only = do_top_left_only_in
    else
      do_top_left_only = .false.
    end if

    ! Array-wise assignment
    mat_x_vec = 0.0_jprb

    if (do_top_left_only) then
      mat_x_vec(1:iend,1) = A(1:iend,1,1)*b(1:iend,1)
    else
      do j1 = 1,m
        do j2 = 1,m
          mat_x_vec(1:iend,j1) = mat_x_vec(1:iend,j1) &
               &               + A(1:iend,j1,j2)*b(1:iend,j2)
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',1,hook_handle)

  end function mat_x_vec

  pure function mat_x_vec_3(n,A,b)

    integer,    intent(in)                   :: n
    real(jprb), intent(in), dimension(n,3,3) :: A
    real(jprb), intent(in), dimension(n,3)   :: b
    real(jprb), dimension(n,3):: mat_x_vec_3
    integer :: j1

    ! Array-wise assignment
    mat_x_vec_3 = 0.0_jprb

    do j1 = 1,3
        mat_x_vec_3(:,j1) = mat_x_vec_3(:,j1) + A(:,j1,1)*b(:,1)
        mat_x_vec_3(:,j1) = mat_x_vec_3(:,j1) + A(:,j1,2)*b(:,2)
        mat_x_vec_3(:,j1) = mat_x_vec_3(:,j1) + A(:,j1,3)*b(:,3)
    end do

  end function mat_x_vec_3


  !---------------------------------------------------------------------
  ! Treat A as an m-by-m square matrix and b as n m-element vectors
  ! (with the n dimension varying fastest), and perform matrix-vector
  ! multiplications on first iend pairs
  function singlemat_x_vec(n,iend,m,A,b)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                    :: n, m, iend
    real(jprb), intent(in), dimension(m,m)    :: A
    real(jprb), intent(in), dimension(:,:)    :: b
    real(jprb),             dimension(iend,m) :: singlemat_x_vec

    integer    :: j1, j2
    real(jprb) :: hook_handle
    
    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',0,hook_handle)

    ! Array-wise assignment
    singlemat_x_vec = 0.0_jprb

    do j1 = 1,m
      do j2 = 1,m
        singlemat_x_vec(1:iend,j1) = singlemat_x_vec(1:iend,j1) &
             &                    + A(j1,j2)*b(1:iend,j2)
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',1,hook_handle)

  end function singlemat_x_vec


  ! --- SQUARE MATRIX-MATRIX MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A and B each as n m-by-m square matrices (with the n
  ! dimension varying fastest) and perform matrix multiplications on
  ! all n matrix pairs
  function mat_x_mat(n,iend,m,A,B,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                      :: n, m, iend
    integer,    intent(in), optional            :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:)    :: A, B

    real(jprb),             dimension(iend,m,m) :: mat_x_mat
    integer    :: j1, j2, j3
    integer    :: mblock, m2block
    integer    :: i_actual_matrix_pattern
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = IMatrixPatternDense
    end if

    ! Array-wise assignment
    mat_x_mat = 0.0_jprb

    if (i_actual_matrix_pattern == IMatrixPatternShortwave) then
      ! Matrix has a sparsity pattern
      !     (C D E)
      ! A = (F G H)
      !     (0 0 I)
      mblock = m/3
      m2block = 2*mblock 
      ! Do the top-left (C, D, F, G)
      do j2 = 1,m2block
        do j1 = 1,m2block
          do j3 = 1,m2block
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
      do j2 = m2block+1,m
        ! Do the top-right (E & H)
        do j1 = 1,m2block
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
        ! Do the bottom-right (I)
        do j1 = m2block+1,m
          do j3 = m2block+1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
    else
      ! Ordinary dense matrix
      do j2 = 1,m
        do j1 = 1,m
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',1,hook_handle)

  end function mat_x_mat

  pure function mat_x_mat_dense(n,m,A,B)

    integer,    intent(in)                      :: n, m
    real(jprb), intent(in), dimension(n,m,m)    :: A, B

    real(jprb),             dimension(n,m,m) :: mat_x_mat_dense
    integer    :: j1, j2, j3

    if (m==3) then
      do j2 = 1,m
        do j1 = 1,m
          mat_x_mat_dense(:,j1,j2) = A(:,j1,1)*B(:,1,j2) &
          &  + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
        end do
      end do

    else 
  
      ! Array-wise assignment
      mat_x_mat_dense = 0.0_jprb

      do j2 = 1,m
        do j1 = 1,m
          do j3 = 1,m
            mat_x_mat_dense(:,j1,j2) = mat_x_mat_dense(:,j1,j2) &
                  &                  + A(:,j1,j3)*B(:,j3,j2)
          end do
        end do
      end do

    end if

  end function mat_x_mat_dense
 
  !---------------------------------------------------------------------
  ! Treat A as an m-by-m matrix and B as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function singlemat_x_mat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(m,m)      :: A
    real(jprb), intent(in), dimension(:,:,:)    :: B
    real(jprb),             dimension(iend,m,m) :: singlemat_x_mat

    integer    :: j1, j2, j3
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',0,hook_handle)

    ! Array-wise assignment
    singlemat_x_mat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          singlemat_x_mat(1:iend,j1,j2) = singlemat_x_mat(1:iend,j1,j2) &
               &                        + A(j1,j3)*B(1:iend,j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',1,hook_handle)

  end function singlemat_x_mat


  !---------------------------------------------------------------------
  ! Treat B as an m-by-m matrix and A as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function mat_x_singlemat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:)    :: A
    real(jprb), intent(in), dimension(m,m)      :: B

    real(jprb),             dimension(iend,m,m) :: mat_x_singlemat
    integer    :: j1, j2, j3
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',0,hook_handle)

    ! Array-wise assignment
    mat_x_singlemat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          mat_x_singlemat(1:iend,j1,j2) = mat_x_singlemat(1:iend,j1,j2) &
               &                        + A(1:iend,j1,j3)*B(j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',1,hook_handle)

  end function mat_x_singlemat


  !---------------------------------------------------------------------
  ! Compute I-A*B where I is the identity matrix and A & B are n
  ! m-by-m square matrices
  function identity_minus_mat_x_mat(n,iend,m,A,B,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                   :: n, m, iend
    integer,    intent(in), optional         :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:) :: A, B
    real(jprb),             dimension(iend,m,m) :: identity_minus_mat_x_mat

    integer    :: j
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,A,B,i_matrix_pattern)
    else
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,A,B)
    end if

    identity_minus_mat_x_mat = - identity_minus_mat_x_mat
    do j = 1,m
      identity_minus_mat_x_mat(1:iend,j,j) &
           &     = 1.0_jprb + identity_minus_mat_x_mat(1:iend,j,j)
    end do

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',1,hook_handle)

  end function identity_minus_mat_x_mat


  ! --- REPEATEDLY SQUARE A MATRIX ---

  !---------------------------------------------------------------------
  ! Square m-by-m matrix "A" nrepeat times. A will be corrupted by
  ! this function.
  function repeated_square(m,A,nrepeat,i_matrix_pattern)
    integer,    intent(in)           :: m, nrepeat
    real(jprb), intent(inout)        :: A(m,m)
    integer,    intent(in), optional :: i_matrix_pattern
    real(jprb)                       :: repeated_square(m,m)

    integer :: j1, j2, j3, j4
    integer :: mblock, m2block
    integer :: i_actual_matrix_pattern

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = IMatrixPatternDense
    end if

    if (i_actual_matrix_pattern == IMatrixPatternShortwave) then
      ! Matrix has a sparsity pattern
      !     (C D E)
      ! A = (F G H)
      !     (0 0 I)
      mblock = m/3
      m2block = 2*mblock
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        ! Do the top-left (C, D, F & G)
        do j2 = 1,m2block
          do j1 = 1,m2block
            do j3 = 1,m2block
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        do j2 = m2block+1, m
          ! Do the top-right (E & H)
          do j1 = 1,m2block
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
          ! Do the bottom-right (I)
          do j1 = m2block+1, m
            do j3 = m2block+1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          A = repeated_square
        end if
      end do
    else
      ! Ordinary dense matrix
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        do j2 = 1,m
          do j1 = 1,m
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          A = repeated_square
        end if
      end do
    end if

  end function repeated_square

  pure subroutine repeated_square_sw_9(nrepeat,A,B)
    integer,    intent(in)            :: nrepeat
    real(jprb), intent(inout)         :: A(9,9)
    real(jprb), intent(out)           :: B(9,9)

    integer :: j1, j2, j3, j4

    do j4 = 1,nrepeat
      B = 0.0_jprb
      ! Do the top-left (C, D, F & G)
      do j2 = 1,6
          do j3 = 1,6
            B(1:6,j2) = B(1:6,j2) + A(1:6,j3)*A(j3,j2)
          end do
      end do
      do j2 = 7, 9
          do j3 = 1,9 ! Do the top-right (E & H)
            B(1:6,j2) = B(1:6,j2) + A(1:6,j3)*A(j3,j2)
          end do
          do j3 = 7,9  ! Do the bottom-right (I)
            B(7:9,j2) = B(7:9,j2) + A(7:9,j3)*A(j3,j2)
          end do
      end do
      if (j4 < nrepeat) then
          A = B
      end if
    end do
  end subroutine repeated_square_sw_9


  ! --- SOLVE LINEAR EQUATIONS ---

  !---------------------------------------------------------------------
  ! Solve Ax=b to obtain x.  Version optimized for 2x2 matrices using
  ! Cramer's method: "A" contains n 2x2 matrices and "b" contains n
  ! 2-element vectors; returns A^-1 b.
  pure subroutine solve_vec_2(n,iend,A,b,x)

    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  A(1:iend,1,1)*A(1:iend,2,2) &
         &                - A(1:iend,1,2)*A(1:iend,2,1))

    x(1:iend,1) = inv_det*(A(1:iend,2,2)*b(1:iend,1)-A(1:iend,1,2)*b(1:iend,2))
    x(1:iend,2) = inv_det*(A(1:iend,1,1)*b(1:iend,2)-A(1:iend,2,1)*b(1:iend,1))

  end subroutine solve_vec_2


  !---------------------------------------------------------------------
  ! Solve AX=B to obtain X, i.e. the matrix right-hand-side version of
  ! solve_vec_2, with A, X and B all containing n 2x2 matrices;
  ! returns A^-1 B using Cramer's method.
  pure subroutine solve_mat_2(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)
    real(jprb), intent(out) :: X(:,:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  A(1:iend,1,1)*A(1:iend,2,2) &
         &                - A(1:iend,1,2)*A(1:iend,2,1))

    X(1:iend,1,1) = inv_det*( A(1:iend,2,2)*B(1:iend,1,1) &
         &                   -A(1:iend,1,2)*B(1:iend,2,1))
    X(1:iend,2,1) = inv_det*( A(1:iend,1,1)*B(1:iend,2,1) &
         &                   -A(1:iend,2,1)*B(1:iend,1,1))
    X(1:iend,1,2) = inv_det*( A(1:iend,2,2)*B(1:iend,1,2) &
         &                   -A(1:iend,1,2)*B(1:iend,2,2))
    X(1:iend,2,2) = inv_det*( A(1:iend,1,1)*B(1:iend,2,2) &
         &                   -A(1:iend,2,1)*B(1:iend,1,2))

  end subroutine solve_mat_2


  !---------------------------------------------------------------------
  ! Solve Ax=b optimized for 3x3 matrices, using LU
  ! factorization and substitution without pivoting.
  pure subroutine solve_vec_3(n,iend,A,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    ! Some compilers unfortunately don't support assocate
    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3), &
    !         y1 => b(:,1), x1 => solve_vec3(:,1), &
    !         x2 => solve_vec3(:,2), x3 => solve_vec3(:,3))

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(1:iend,2,1) / A(1:iend,1,1)
    L31 = A(1:iend,3,1) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,1,2)
    U23 = A(1:iend,2,3) - L21*A(1:iend,1,3)
    L32 =(A(1:iend,3,2) - L31*A(1:iend,1,2)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(1:iend,2) - L21*b(1:iend,1)
    y3 = b(1:iend,3) - L31*b(1:iend,1) - L32*y2

    ! Solve Ux = y by back substitution
    x(1:iend,3) = y3/U33
    x(1:iend,2) = (y2 - U23*x(1:iend,3)) / U22
    x(1:iend,1) = (b(1:iend,1) - A(1:iend,1,2)*x(1:iend,2) &
         &         - A(1:iend,1,3)*x(1:iend,3)) / A(1:iend,1,1)
    !    end associate

  end subroutine solve_vec_3


  !---------------------------------------------------------------------
  ! Solve AX=B optimized for 3x3 matrices, using LU factorization and
  ! substitution with no pivoting.
  pure subroutine solve_mat_3(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)
    real(jprb), intent(out) :: X(:,:,:)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    integer :: j

    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(1:iend,2,1) / A(1:iend,1,1)
    L31 = A(1:iend,3,1) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,1,2)
    U23 = A(1:iend,2,3) - L21*A(1:iend,1,3)
    L32 =(A(1:iend,3,2) - L31*A(1:iend,1,2)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,1,3) - L32*U23

    do j = 1,3
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = B(:,1,j)
      y2 = B(1:iend,2,j) - L21*B(1:iend,1,j)
      y3 = B(1:iend,3,j) - L31*B(1:iend,1,j) - L32*y2
      ! Solve UX(:,:,j) = y by back substitution
      X(1:iend,3,j) = y3 / U33
      X(1:iend,2,j) = (y2 - U23*X(1:iend,3,j)) / U22
      X(1:iend,1,j) = (B(1:iend,1,j) - A(1:iend,1,2)*X(1:iend,2,j) &
           &          - A(1:iend,1,3)*X(1:iend,3,j)) / A(1:iend,1,1)
    end do

  end subroutine solve_mat_3


  !---------------------------------------------------------------------
  ! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
  ! is a diagonal matrix, using LU factorization and substitution with
  ! no pivoting.
  pure subroutine diag_mat_right_divide_3(n,A,B,X)
    integer,    intent(in)  :: n
    real(jprb), intent(in)  :: A(n,3,3)
    real(jprb), intent(in)  :: B(n,3)
    real(jprb), intent(out) :: X(n,3,3)

    real(jprb), dimension(n) :: L21, L31, L32
    real(jprb), dimension(n) :: U22, U23, U33
    real(jprb), dimension(n) :: y2, y3

    integer :: j

    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
    ! LU decomposition of the *transpose* of A:
    !       ( 1        )   (U11 U12 U13)
    ! A^T = (L21  1    ) * (    U22 U23)
    !       (L31 L32  1)   (        U33)
    L21 = A(:,1,2) / A(:,1,1)
    L31 = A(:,1,3) / A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,2,1)
    U23 = A(:,3,2) - L21*A(:,3,1)
    L32 =(A(:,2,3) - L31*A(:,2,1)) / U22
    U33 = A(:,3,3) - L31*A(:,3,1) - L32*U23

    ! Solve X(1,:) = A^-T ( B(1) )
    !                     (  0   )
    !                     (  0   )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = B(:,1)
    y2 = - L21*B(:,1)
    y3 = - L31*B(:,1) - L32*y2
    ! Solve UX(:,:,j) = y by back substitution
    X(:,1,3) = y3 / U33
    X(:,1,2) = (y2 - U23*X(:,1,3)) / U22
    X(:,1,1) = (B(:,1) - A(:,2,1)*X(:,1,2) &
         &          - A(:,3,1)*X(:,1,3)) / A(:,1,1)

    ! Solve X(2,:) = A^-T (  0   )
    !                     ( B(2) )
    !                     (  0   )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = B(:,2)
    y3 = - L32*B(:,2)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,2,3) = y3 / U33
    X(:,2,2) = (B(:,2) - U23*X(:,2,3)) / U22
    X(:,2,1) = (-A(:,2,1)*X(:,2,2) &
         &           -A(:,3,1)*X(:,2,3)) / A(:,1,1)

    ! Solve X(3,:) = A^-T (  0   )
    !                     (  0   )
    !                     ( B(3) )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = 0
    ! y3 = B(:,3)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,3,3) = B(:,3) / U33
    X(:,3,2) = -U23*X(:,3,3) / U22
    X(:,3,1) = (-A(:,2,1)*X(:,3,2) &
         &          - A(:,3,1)*X(:,3,3)) / A(:,1,1)

  end subroutine diag_mat_right_divide_3


  !---------------------------------------------------------------------
  ! Treat A as n m-by-m matrices and return the LU factorization of A
  ! compressed into a single matrice (with L below the diagonal and U
  ! on and above the diagonal; the diagonal elements of L are 1). No
  ! pivoting is performed.
  pure subroutine lu_factorization(n, iend, m, A, LU)
    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(out) :: LU(iend,m,m)

    real(jprb) :: s(iend)
    integer    :: j1, j2, j3

    ! This routine is adapted from an in-place one, so we first copy
    ! the input into the output.
    LU(1:iend,1:m,1:m) = A(1:iend,1:m,1:m)

    do j2 = 1, m
      do j1 = 1, j2-1
        s = LU(1:iend,j1,j2)
        do j3 = 1, j1-1
          s = s - LU(1:iend,j1,j3) * LU(1:iend,j3,j2)
        end do
        LU(1:iend,j1,j2) = s
      end do
      do j1 = j2, m
        s = LU(1:iend,j1,j2)
        do j3 = 1, j2-1
          s = s - LU(1:iend,j1,j3) * LU(1:iend,j3,j2)
        end do
        LU(1:iend,j1,j2) = s
      end do
      if (j2 /= m) then
        s = 1.0_jprb / LU(1:iend,j2,j2)
        do j1 = j2+1, m
          LU(1:iend,j1,j2) = LU(1:iend,j1,j2) * s
        end do
      end if
    end do

  end subroutine lu_factorization


  !---------------------------------------------------------------------
  ! Treat LU as an LU-factorization of an original matrix A, and
  ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
  ! m-element vectors.
  pure subroutine lu_substitution(n,iend,m,LU,b,x)
    ! CHECK: dimensions should be ":"?
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: LU(iend,m,m)
    real(jprb), intent(in) :: b(:,:)
    real(jprb), intent(out):: x(iend,m)

    integer :: j1, j2

    x(1:iend,1:m) = b(1:iend,1:m)

    ! First solve Ly=b
    do j2 = 2, m
      do j1 = 1, j2-1
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*LU(1:iend,j2,j1)
      end do
    end do
    ! Now solve Ux=y
    do j2 = m, 1, -1
      do j1 = j2+1, m
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*LU(1:iend,j2,j1)
      end do
      x(1:iend,j2) = x(1:iend,j2) / LU(1:iend,j2,j2)
    end do

  end subroutine lu_substitution


  !---------------------------------------------------------------------
  ! Return matrix X where AX=B. LU, A, X, B all consist of n m-by-m
  ! matrices.
  pure subroutine solve_mat_n(n,iend,m,A,B,X)
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: A(:,:,:)
    real(jprb), intent(in) :: B(:,:,:)
    real(jprb), intent(out):: X(iend,m,m)

    real(jprb) :: LU(iend,m,m)

    integer :: j

    call lu_factorization(n,iend,m,A,LU)

    do j = 1, m
      call lu_substitution(n,iend,m,LU,B(1:,1:m,j),X(1:iend,1:m,j))
!      call lu_substitution(n,iend,m,LU,B(1:n,1:m,j),X(1:iend,1:m,j))
    end do

  end subroutine solve_mat_n


  !---------------------------------------------------------------------
  ! Solve Ax=b, where A consists of n m-by-m matrices and x and b
  ! consist of n m-element vectors. For m=2 or m=3, this function
  ! calls optimized versions, otherwise it uses general LU
  ! decomposition without pivoting.
  function solve_vec(n,iend,m,A,b)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: A(:,:,:)
    real(jprb), intent(in) :: b(:,:)

    real(jprb)             :: solve_vec(iend,m)
    real(jprb)             :: LU(iend,m,m)
    real(jprb)             :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_vec',0,hook_handle)

    if (m == 2) then
      call solve_vec_2(n,iend,A,b,solve_vec)
    elseif (m == 3) then
      call solve_vec_3(n,iend,A,b,solve_vec)
    else
      call lu_factorization(n,iend,m,A,LU)
      call lu_substitution(n,iend,m,LU,b,solve_vec)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_vec',1,hook_handle)

  end function solve_vec


  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of n m-by-m matrices. For m=2
  ! or m=3, this function calls optimized versions, otherwise it uses
  ! general LU decomposition without pivoting.
  function solve_mat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)

    real(jprb)              :: solve_mat(iend,m,m)
    real(jprb)              :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_mat',0,hook_handle)

    if (m == 2) then
      call solve_mat_2(n,iend,A,B,solve_mat)
    elseif (m == 3) then
      call solve_mat_3(n,iend,A,B,solve_mat)
    else
      call solve_mat_n(n,iend,m,A,B,solve_mat)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_mat',1,hook_handle)

  end function solve_mat

  pure subroutine mat_square_sw(m,iend,A,C)

    integer,    intent(in)                        :: m,iend
    real(jprb), intent(in),  dimension(iend,m,m)  :: A
    real(jprb), intent(out), dimension(iend,m,m)  :: C
    integer    :: j1, j2, j3
    integer    :: mblock, m2block

    ! Array-wise assignment
    C = 0.0_jprb

    ! Matrix has a sparsity pattern
    !     (C D E)
    ! A = (F G H)
    !     (0 0 I)

    mblock = m/3       ! 
    m2block = 2*mblock ! 

    ! Do the top-left (C, D, F, G)
    do j2 = 1,m2block  !    1,6 
      do j1 = 1,m2block !   1,6
        do j3 = 1,m2block ! 1,6
          C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
        ! using sum was faster on GCC+AMD Zen platform but not on Intel
        ! C(:,j1,j2) = sum(A(:,j1,1:m2block)*A(:,1:m2block,j2),2)
      end do
    end do

    do j2 = m2block+1,m  ! 7,9
      ! Do the top-right (E & H)
      do j1 = 1,m2block  ! 1,6
        do j3 = 1,m      ! 1,9
          C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
      end do

      ! Do the bottom-right (I)
      do j1 = m2block+1,m   ! 7,9
        do j3 = m2block+1,m ! 7,9
          C(:,j1,j2) =  C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
      end do
    end do

  end subroutine mat_square_sw

  pure subroutine mat_square_sw_9(iend,A,C)

    integer,    intent(in)                        :: iend
    real(jprb), intent(in),  dimension(iend,9,9)  :: A
    real(jprb), intent(out), dimension(iend,9,9)  :: C
    integer    :: j1, j2, j3, jg

    ! First input matrix has pattern

    !     (C    D     E)
    ! A = (F=-D G=-C  H)
    !     (0    0     I)

    ! As a result, output and subsequent input matrices have pattern
    
    !     (C    D    E)
    ! A = (F=D  G=C  H)
    !     (0    0    I)

    do j2 = 1,3  !    1,3 
      do j1 = 1,6 !   1,6   C, F
        ! do j3 = 1,6 ! 1,6
        !   C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        ! end do
        ! Further speedup: flatten last loop, only one write SIMD instruction
        C(:,j1,j2) = A(:,j1,1)*A(:,1,j2) + A(:,j1,2)*A(:,2,j2) &
        &          + A(:,j1,3)*A(:,3,j2) + A(:,j1,4)*A(:,4,j2) &
        &          + A(:,j1,5)*A(:,5,j2) + A(:,j1,6)*A(:,6,j2)
      end do
    end do
    ! D,G
    C(:,1:3,4:6) = C(:,4:6,1:3) ! D = F
    C(:,4:6,4:6) = C(:,1:3,1:3) ! G = C

    ! Lower left corner with zeros
    C(:,7:9,1:6) = 0.0_jprb

    do j2 = 7,9  ! 7,9
      ! Do the top-right (E & H)
      do j1 = 1,6  ! 1,6
        C(:,j1,j2) = A(:,j1,1)*A(:,1,j2) + A(:,j1,2)*A(:,2,j2) &
        &          + A(:,j1,3)*A(:,3,j2) + A(:,j1,4)*A(:,4,j2) &
        &          + A(:,j1,5)*A(:,5,j2) + A(:,j1,6)*A(:,6,j2) &
        &          + A(:,j1,7)*A(:,7,j2) + A(:,j1,8)*A(:,8,j2) &
        &          + A(:,j1,9)*A(:,9,j2)
      end do
      ! Do the bottom-right (I)
      do j1 = 7,9   ! 7,9
        C(:,j1,j2) = A(:,j1,7)*A(:,7,j2) + A(:,j1,8)*A(:,8,j2) + A(:,j1,9)*A(:,9,j2)
      end do
    end do

  end subroutine mat_square_sw_9

  pure subroutine mat_x_mat_sw_9(iend,A,B,C)

    integer,    intent(in)                      :: iend
    real(jprb), intent(in), dimension(iend,9,9) :: A, B
    real(jprb), intent(out),dimension(iend,9,9) :: C
    integer    :: j1, j2, j3

    do j2 = 1,3  !    1,3 
      do j1 = 1,6 !   1,6   C, F
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) &
        &          + A(:,j1,3)*B(:,3,j2) + A(:,j1,4)*B(:,4,j2) &
        &          + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2)
      end do
    end do
    ! D,G
    C(:,1:3,4:6) = C(:,4:6,1:3) ! D = F
    C(:,4:6,4:6) = C(:,1:3,1:3) ! G = C

    ! Lower left corner with zeros
    C(:,7:9,1:6) = 0.0_jprb

    do j2 = 7,9  ! 7,9
      ! Do the top-right (E & H)
      do j1 = 1,6  ! 1,6
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) &
        &          + A(:,j1,3)*B(:,3,j2) + A(:,j1,4)*B(:,4,j2) &
        &          + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2) &
        &          + A(:,j1,7)*B(:,7,j2) + A(:,j1,8)*B(:,8,j2) &
        &          + A(:,j1,9)*B(:,9,j2)
      end do

      ! Do the bottom-right (I)
      do j1 = 7,9   ! 7,9
        C(:,j1,j2) = A(:,j1,7)*B(:,7,j2) + A(:,j1,8)*B(:,8,j2) + A(:,j1,9)*B(:,9,j2)
      end do
    end do

  end subroutine mat_x_mat_sw_9

  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of iend m-by-m matrices
  ! Overwrite B with X. A is corrupted
  pure subroutine solve_mat_sw_9(iend,A,B)
    integer,    intent(in) :: iend
    real(jprb), intent(inout) :: A(iend,9,9) ! A=LU is corrupted
    real(jprb), intent(inout) :: B(iend,9,9) ! X = B, both input and output
    ! real(jprb), intent(out):: X(iend,m,m)

    integer :: j,j1, j2, j3, mblock, m2block,m
    real(jprb) :: s(iend)

    !     (C   D  E)
    ! A = (-D -C  H)
    !     (0   0  I)

    m = 9
    mblock = 3 !m/3       
    m2block = 6 !2*mblock 

    ! factorization of A into LU

    ! First do columns 1-6, for which only rows 1-6 have non-negative entries
    do j2 = 1, m2block
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          A(:,j1,j2) = A(:,j1,j2)- A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      do j1 = j2, m2block
        do j3 = 1, j2-1
          A(:,j1,j2) = A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      s = 1.0_jprb / A(:,j2,j2)
      do j1 = j2+1, m2block
        A(:,j1,j2) = A(:,j1,j2) * s
      end do
    end do

    ! Remaining columns
    do j2 = m2block+1, m
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          A(:,j1,j2) = A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      do j1 = j2, m
        do j3 = 1, j2-1
          A(:,j1,j2)= A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      if (j2 /= m) then
        s = 1.0_jprb / A(:,j2,j2)
        do j1 = j2+1, m
          A(:,j1,j2) = A(:,j1,j2) * s
        end do
      end if
    end do

    !---------------------------------------------------------------------
    ! Treat LU as an LU-factorization of an original matrix A, and
    ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
    ! m-element vectors.
    ! Here B is both input b and output x, and A has been LU factorized, combining L and U
    ! into one matrix where the diagonal is the diagonal of U (L has ones in the diagonal)

    ! A and B both have following structure:
    !     (C   D  E)
    !     (F   G  H)
    !     (0   0  I)

    ! Separate j3 (columns) into two regions to avoid redundant operations with zero
    do j3 = 1,m2block ! in this region B(:,7:9),A(:,7:9) are 0
      ! First solve Ly=b
      do j2 = 2, m2block
        do j1 = 1, j2-1
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        ! No division because diagonal of L is unity
      end do
      ! Now solve Ux=y
      do j2 = m2block, 1, -1
        do j1 = j2+1, m2block
          B(:,j2,j3) = ( B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1) )
        end do
        B(:,j2,j3) = B(:,j2,j3) / A(:,j2,j2) ! Divide by diagonal of A=U
      end do
    end do

    do j3 = m2block+1,m ! columns 7-9: here B has nonzero values for all rows, but A doesn't 
      ! First solve Ly=b
      ! do j2 = 2, m
      do j2 = 2, m2block 
        do j1 = 1, j2-1
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        ! No division because diagonal of L is unity
      end do
      ! When j2 = 7, the A terms are all 0, because A(7,1:6)=0
      ! when j2 = 8, only the last j1 has nonzero A
      ! When j2 = 9, two last j1 are nonzero
      B(:,8,j3) = B(:,8,j3) - B(:,7,j3)*A(:,8,7)   ! j2 = 8
      B(:,9,j3) = B(:,9,j3) - B(:,8,j3)*A(:,9,8) - B(:,7,j3)*A(:,9,7) ! j2 = 9

      ! Now solve Ux=y
      do j2 = m, 1, -1
        do j1 = j2+1, m
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        B(:,j2,j3) = B(:,j2,j3) / A(:,j2,j2) ! Divide by diagonal of A=U
      end do
    end do

  end subroutine solve_mat_sw_9


  ! --- MATRIX EXPONENTIATION ---
  !---------------------------------------------------------------------
  ! Perform matrix exponential of n m-by-m matrices stored in A (where
  ! index n varies fastest) using the Higham scaling and squaring
  ! method. The result is placed in A. This routine is intended for
  ! speed so is accurate only to single precision.  For simplicity and
  ! to aid vectorization, the Pade approximant of order 7 is used for
  ! all input matrices, perhaps leading to a few too many
  ! multiplications for matrices with a small norm.
  subroutine expm(n,iend,m,A,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)      :: n, m, iend
    real(jprb), intent(inout)   :: A(n,m,m)
    integer,    intent(in)      :: i_matrix_pattern

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(iend,m,m) :: A2, A4, A6
    real(jprb), dimension(iend,m,m) :: U, V

    real(jprb) :: normA(iend), sum_column(iend)

    integer    :: j1, j2, j3
    real(jprb) :: frac(iend)
    integer    :: expo(iend)
    real(jprb) :: scaling(iend)

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:expm',0,hook_handle)

    normA = 0.0_jprb

    ! Compute the 1-norms of A
    do j3 = 1,m
      sum_column(:) = 0.0_jprb
      do j2 = 1,m
        do j1 = 1,iend
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      do j1 = 1,iend
        if (sum_column(j1) > normA(j1)) then
          normA(j1) = sum_column(j1)
        end if
      end do
    end do

    frac = fraction(normA/theta(3))
    expo = exponent(normA/theta(3))
    where (frac == 0.5_jprb)
      expo = expo - 1
    end where

    where (expo < 0)
      expo = 0
    end where

    ! Scale the input matrices by a power of 2
    scaling = 2.0_jprb**(-expo)
    do j3 = 1,m
      do j2 = 1,m
        A(1:iend,j2,j3) = A(1:iend,j2,j3) * scaling
      end do
    end do
    ! Pade approximant of degree 7
    A2 = mat_x_mat(n,iend,m,A, A, i_matrix_pattern)
    A4 = mat_x_mat(n,iend,m,A2,A2,i_matrix_pattern)
    A6 = mat_x_mat(n,iend,m,A2,A4,i_matrix_pattern)

    V = c(8)*A6 + c(6)*A4 + c(4)*A2
    do j3 = 1,m
      V(:,j3,j3) = V(:,j3,j3) + c(2)
    end do
    U = mat_x_mat(n,iend,m,A,V,i_matrix_pattern)
    V = c(7)*A6 + c(5)*A4 + c(3)*A2
    ! Add a multiple of the identity matrix
    do j3 = 1,m
      V(:,j3,j3) = V(:,j3,j3) + c(1)
    end do

    V = V-U
    U = 2.0_jprb*U
    A(1:iend,1:m,1:m) = solve_mat(n,iend,m,V,U)

    ! Add the identity matrix
    do j3 = 1,m
      A(1:iend,j3,j3) = A(1:iend,j3,j3) + 1.0_jprb
    end do

    ! Loop through the matrices
    do j1 = 1,iend
      if (expo(j1) > 0) then
        ! Square matrix j1 expo(j1) times          
        A(j1,:,:) = repeated_square(m,A(j1,:,:),expo(j1),i_matrix_pattern)
      end if
    end do

    if (lhook) call dr_hook('radiation_matrix:expm',1,hook_handle)

  end subroutine expm

  !---------------------------------------------------------------------
  ! Like expm, but optimized for the shortwave, which has
  ! a special matrix structure with zeros and repeated elements. 
  ! Further assumes nreg = 3  =>  m = 9
  subroutine expm_opt(iend,A)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)      :: iend
    real(jprb), intent(inout)   :: A(iend,9,9)

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(iend,9,9) :: A2, A4, A6
    real(jprb), dimension(iend,9,9) :: U, V
    real(jprb), dimension(9,9)      :: temp_in, temp_out
    real(jprb) :: normA(iend), sum_column(iend)

    integer    :: j1, j2, j3, j4,minexpo, nrepeat
    real(jprb) :: frac(iend)
    integer    :: expo(iend)
    real(jprb) :: scaling(iend)

    real(jprb) :: hook_handle
#ifdef USE_TIMING
    ret =  gptlstart('expm_opt')
#endif 
    if (lhook) call dr_hook('radiation_matrix:expm_opt',0,hook_handle)

    normA = 0.0_jprb

    ! Compute the 1-norms of A
    do j3 = 1,9
      sum_column(:) = 0.0_jprb
      do j2 = 1,9
        do j1 = 1,iend
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      do j1 = 1,iend
        if (sum_column(j1) > normA(j1)) then
          normA(j1) = sum_column(j1)
        end if
      end do
    end do

    frac = fraction(normA/theta(3))
    expo = exponent(normA/theta(3))
    where (frac == 0.5_jprb)
      expo = expo - 1
    end where

    where (expo < 0)
      expo = 0
    end where
    
    minexpo = minval(expo)
    ! Scale the input matrices by a power of 2
    scaling = 2.0_jprb**(-expo)
    do j3 = 1,9
      do j2 = 1,9
        A(:,j2,j3) = A(:,j2,j3) * scaling
      end do
    end do
    ! Pade approximant of degree 7
#ifdef USE_TIMING
    ret =  gptlstart('expm_Pade_mat_x_mat')
#endif 
#ifdef USE_TIMING
    ret =  gptlstart('expm_mat_squares')
#endif 
    call mat_square_sw_9(iend,A,A2)  ! These matrices have zeroes in the lower left corner AND repeated elements
    call mat_square_sw_9(iend,A2,A4)

#ifdef USE_TIMING
    ret =  gptlstop('expm_mat_squares')
#endif 
    call mat_x_mat_sw_9(iend,A2,A4,A6)     ! These matrices have zeroes in the lower left corner AND repeated elements

    V = c(8)*A6 + c(6)*A4 + c(4)*A2
    do j3 = 1,9
      V(:,j3,j3) = V(:,j3,j3) + c(2)
    end do

    U = mat_x_mat(iend,iend,9,A,V,IMatrixPatternShortwave)

#ifdef USE_TIMING
    ret =  gptlstop('expm_Pade_mat_x_mat')
#endif 
    V = c(7)*A6 + c(5)*A4 + c(3)*A2
    ! Add a multiple of the identity matrix
    do j3 = 1,9
      V(:,j3,j3) = V(:,j3,j3) + c(1)
    end do

    V = V-U
    A = 2.0_jprb*U
#ifdef USE_TIMING
    ret =  gptlstart('expm_solve_mat')
#endif 
    ! A = solve_mat(n,iend,m,V,U)
    call solve_mat_sw_9(iend,V,A)
#ifdef USE_TIMING
    ret =  gptlstop('expm_solve_mat')
#endif 
    ! Add the identity matrix
    do j3 = 1,9
      A(:,j3,j3) = A(:,j3,j3) + 1.0_jprb
    end do

    ! Loop through the matrices    
#ifdef USE_TIMING
    ret =  gptlstart('expm_square_matrix_repeated')
#endif 

#ifdef USE_TIMING
    ret =  gptlstart('expm_square_matrix_minexpo')
#endif 

    ! To improve efficiency, square all matrices with the minimum expo first, and then square individual matrices as needed
    do j1 = 1,minexpo
      call mat_square_sw(9,iend,A,A2)
      A = A2
    end do

#ifdef USE_TIMING
    ret =  gptlstop('expm_square_matrix_minexpo')
#endif 

    do j1 = 1,iend
      if (expo(j1) > minexpo) then
        nrepeat = expo(j1)-minexpo
        !Square matrix nrepeat times 
        temp_in =  A(j1,:,:)
        call repeated_square_sw_9(nrepeat,temp_in,temp_out)
        A(j1,:,:) = temp_out
      end if
    end do

#ifdef USE_TIMING
    ret =  gptlstop('expm_square_matrix_repeated')
#endif 
    if (lhook) call dr_hook('radiation_matrix:expm_opt',1,hook_handle)
#ifdef USE_TIMING
    ret =  gptlstop('expm_opt')
#endif 
  end subroutine expm_opt


  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 2x2 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b)
  !   ( a  -b)
  ! and a and b are assumed to be positive or zero.  The solution uses
  ! Putzer's algorithm - see the appendix of Hogan et al. (GMD 2018)
  subroutine fast_expm_exchange_2(n,a,b,R)

    use yomhook, only : lhook, dr_hook

    integer,                      intent(in)  :: n
    real(jprb), dimension(n),     intent(in)  :: a, b
    real(jprb), dimension(n,2,2), intent(out) :: R

    real(jprb), dimension(n) :: factor

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',0,hook_handle)

    ! Security to ensure that if a==b==0 then the identity matrix is returned
    factor = (1.0_jprb - exp(-(a(:)+b(:))))/max(1.0e-12_jprb,a(:)+b(:))

    R(:,1,1) = 1.0_jprb - factor*a(:)
    R(:,2,1) = factor*a(:)
    R(:,1,2) = factor*b(:)
    R(:,2,2) = 1.0_jprb - factor*b(:)

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',1,hook_handle)

  end subroutine fast_expm_exchange_2


  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 3x3 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b   0)
  !   ( a -b-c  d)
  !   ( 0   c  -d)
  ! and a-d are assumed to be positive or zero.  The solution uses the
  ! diagonalization method and is a slight generalization of the
  ! solution provided in the appendix of Hogan et al. (GMD 2018),
  ! which assumed c==d.
  subroutine fast_expm_exchange_3(n,a,b,c,d,R)

    use yomhook, only : lhook, dr_hook

    real(jprb), parameter :: my_epsilon = 1.0e-12_jprb

    integer,                      intent(in)  :: n
    real(jprb), dimension(n),     intent(in)  :: a, b, c, d
    real(jprb), dimension(n,3,3), intent(out) :: R

    ! Eigenvectors
    real(jprb), dimension(n,3,3) :: V

    ! Non-zero Eigenvalues
    real(jprb), dimension(n) :: lambda1, lambda2

    ! Diagonal matrix of the exponential of the eigenvalues
    real(jprb), dimension(n,3) :: diag

    ! Result of diag right-divided by V
    real(jprb), dimension(n,3,3) :: X

    ! Intermediate arrays
    real(jprb), dimension(n) :: y2, y3

    real(jprb), dimension(n) :: L21, L31, L32
    real(jprb), dimension(n) :: U22, U23, U33

    integer :: j1, j2

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',0,hook_handle)

    ! Eigenvalues
    y2 = 0.5_jprb * (a(:)+b(:)+c(:)+d(:))
    y3 = sqrt(y2*y2 - (a(:)*c(:) + a(:)*d(:) + b(:)*d(:)))
    lambda1 = -y2 + y3
    lambda2 = -y2 - y3

    ! Eigenvectors, with securities such taht if a--d are all zero
    ! then V is non-singular and the identity matrix is returned in R;
    ! note that lambdaX is typically negative so we need a
    ! sign-preserving security
    V(:,1,1) = max(my_epsilon, b(:)) &
         &  / sign(max(my_epsilon, abs(a(:) + lambda1)), a(:) + lambda1)
    V(:,1,2) = b(:) &
         &  / sign(max(my_epsilon, abs(a(:) + lambda2)), a(:) + lambda2)
    V(:,1,3) = b(:) / max(my_epsilon, a(:))
    V(:,2,:) = 1.0_jprb
    V(:,3,1) = c(:) &
         &  / sign(max(my_epsilon, abs(d(:) + lambda1)), d(:) + lambda1)
    V(:,3,2) = c(:) &
         &  / sign(max(my_epsilon, abs(d(:) + lambda2)), d(:) + lambda2)
    V(:,3,3) = max(my_epsilon, c(:)) / max(my_epsilon, d(:))
    
    diag(:,1) = exp(lambda1)
    diag(:,2) = exp(lambda2)
    diag(:,3) = 1.0_jprb

    ! ------ Compute X = diag * V^-1 ---------
    !  call diag_mat_right_divide_3(n,V,diag,X)

      !    associate (U11 => V(:,1,1), U12 => V(:,1,2), U13 => V(1,3))
    ! LU decomposition of the *transpose* of V:
    !       ( 1        )   (U11 U12 U13)
    ! V^T = (L21  1    ) * (    U22 U23)
    !       (L31 L32  1)   (        U33)
    L21 = V(:,1,2) / V(:,1,1)
    L31 = V(:,1,3) / V(:,1,1)
    U22 = V(:,2,2) - L21*V(:,2,1)
    U23 = V(:,3,2) - L21*V(:,3,1)
    L32 =(V(:,2,3) - L31*V(:,2,1)) / U22
    U33 = V(:,3,3) - L31*V(:,3,1) - L32*U23

    ! Solve X(1,:) = V^-T ( diag(1) )
    !                     (  0   )
    !                     (  0   )
    ! Solve Ly = diag(:,:,j) by forward substitution
    ! y1 = diag(:,1)
    y2 = - L21*diag(:,1)
    y3 = - L31*diag(:,1) - L32*y2
    ! Solve UX(:,:,j) = y by back substitution
    X(:,1,3) = y3 / U33
    X(:,1,2) = (y2 - U23*X(:,1,3)) / U22
    X(:,1,1) = (diag(:,1) - V(:,2,1)*X(:,1,2) &
         &          - V(:,3,1)*X(:,1,3)) / V(:,1,1)

    ! Solve X(2,:) = V^-T (  0   )
    !                     ( diag(2) )
    !                     (  0   )
    ! Solve Ly = diag(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = diag(:,2)
    y3 = - L32*diag(:,2)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,2,3) = y3 / U33
    X(:,2,2) = (diag(:,2) - U23*X(:,2,3)) / U22
    X(:,2,1) = (-V(:,2,1)*X(:,2,2) &
         &           -V(:,3,1)*X(:,2,3)) / V(:,1,1)

    ! Solve X(3,:) = V^-T (  0   )
    !                     (  0   )
    !                     ( diag(3) )
    ! Solve Ly = diag(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = 0
    ! y3 = diag(:,3)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,3,3) = diag(:,3) / U33
    X(:,3,2) = -U23*X(:,3,3) / U22
    X(:,3,1) = (-V(:,2,1)*X(:,3,2) &
         &          - V(:,3,1)*X(:,3,3)) / V(:,1,1)

    ! Compute V * X
    do j1 = 1,3
      do j2 = 1,3
        R(:,j2,j1) = V(:,j2,1)*X(:,1,j1) &
             &          + V(:,j2,2)*X(:,2,j1) &
             &          + V(:,j2,3)*X(:,3,j1)
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',1,hook_handle)

  end subroutine fast_expm_exchange_3



!   ! ------------------------------------------------------------------------------
!   ! ------ PROCEDURES FOR USING EIGENVALUE DECOMPOSITION IN SPARTACUS, -----------
!   ! ------ RATHER THAN MATRIX EXPONENTIALS. THIS WAS MUCH SLOWER WHEN TESTED -----
!   ! ------ THE MATRIX EXPONENTIAL METHOD SEEMS COMPUTATIONALLY MORE ATTRACTIVE ---
!   ! ------------------------------------------------------------------------------
!   !--------------------------------------------------------------------------------
!   ! Treat A as n k-by-l and B as n l-by-m matrices (with the n
!   ! dimension varying fastest) and perform matrix multiplications on
!   ! all n matrix pairs
!   function rect_mat_x_mat(n,k,l,m,A,B)

!     use yomhook, only : lhook, dr_hook

!     integer,    intent(in)                   :: n, k, l, m
!     real(jprb), intent(in), dimension(n,k,l) :: A
!     real(jprb), intent(in), dimension(n,l,m) :: B

!     real(jprb),             dimension(n,k,m) :: rect_mat_x_mat
!     integer    :: j1, j2, j3
!     real(jprb) :: hook_handle

!     if (lhook) call dr_hook('radtool_matrix:rect_mat_x_mat',0,hook_handle)

!     ! Array-wise assignment
!     rect_mat_x_mat = 0.0_jprb

!     ! Ordinary dense matrix
!     do j2 = 1,m
!        do j1 = 1,k
!           do j3 = 1,l
!              rect_mat_x_mat(:,j1,j2) = rect_mat_x_mat(:,j1,j2) &
!                   &                  + A(:,j1,j3)*B(:,j3,j2)
!           end do
!        end do
!     end do

!     if (lhook) call dr_hook('radtool_matrix:rect_mat_x_mat',1,hook_handle)

!   end function rect_mat_x_mat

!    !---------------------------------------------------------------------
!   ! Treat LU as an LU-factorization of an original matrix A, and
!   ! return the inverse of A. This is done using the same method as
!   ! lu_substitution but with the identity matrix on the right-hand
!   ! side. LU consists of n m-by-m matrices and b as n m-element
!   ! vectors.
!   pure subroutine lu_invert(n,iend,m,LU,X)
!     ! CHECK: dimensions should be ":"?
!     integer,    intent(in) :: n, m, iend
!     real(jprb), intent(in) :: LU(iend,m,m)
!     real(jprb), intent(out):: X(iend,m,m)

!     integer :: j1, j2, j3

!     do j3 = 1,m
!       ! Initialize with the identity matrix
!       do j2 = 1,j3-1
!         X(1:iend,j2,j3) = 0.0_jprb
!       end do
!       X(1:iend,j3,j3) = 1.0_jprb
!       do j2 = j3+1,m
!         X(1:iend,j2,j3) = 0.0_jprb
!       end do

!       ! First solve Ly=b
!       do j2 = 2, m
!         do j1 = 1, j2-1
!           X(1:iend,j2,j3) = X(1:iend,j2,j3) - X(1:iend,j1,j3)*LU(1:iend,j2,j1)
!         end do
!       end do
!       ! Now solve Ux=y
!       do j2 = m, 1, -1
!         do j1 = j2+1, m
!           X(1:iend,j2,j3) = X(1:iend,j2,j3) - X(1:iend,j1,j3)*LU(1:iend,j2,j1)
!         end do
!         X(1:iend,j2,j3) = X(1:iend,j2,j3) / LU(1:iend,j2,j2)
!       end do
!     end do

!   end subroutine lu_invert

!     !---------------------------------------------------------------------
!   ! Invert A, which consists of n m-by-m matrices.
!   function invert(n,iend,m,A)

!     use yomhook, only : lhook, dr_hook

!     integer,    intent(in)  :: n, m, iend
!     real(jprb), intent(in)  :: A(:,:,:)

!     real(jprb)              :: invert(iend,m,m)
!     real(jprb)              :: LU(iend,m,m)

!     integer :: j

!     real(jprb)              :: hook_handle

!     if (lhook) call dr_hook('radtool_matrix:invert',0,hook_handle)

!     ! if (m == 2) then
!     !   ! Reciprocal of determinant; use LU as scratch space
!     !   LU(1:iend,1,1) = 1.0_jprb / (A(1:iend,1,1)*A(1:iend,2,2) &
!     !        &                      -A(1:iend,2,1)*A(1:iend,1,2))
!     !   invert(1:iend,1,1) =  LU(1:iend,1,1) * A(1:iend,2,2)
!     !   invert(1:iend,2,1) = -LU(1:iend,1,1) * A(1:iend,2,1)
!     !   invert(1:iend,1,2) = -LU(1:iend,1,1) * A(1:iend,1,2)
!     !   invert(1:iend,2,2) =  LU(1:iend,1,1) * A(1:iend,1,1)
!     !   call invert_2(n,iend,A,invert)
!     ! else
!       call lu_factorization(n,iend,m,A,LU)
!       call lu_invert(n,iend,m,LU,invert)
!     ! end if

!     if (lhook) call dr_hook('radtool_matrix:invert',1,hook_handle)

!   end function invert

!    ! Brute-force method, without any clever optimizations using the
!   ! Schur complement
!   subroutine direct_diffuse_part(nmat, ndiff, ndir, mu0, &
!        &  exp_lambda_dz, exp_eigenval_dir_dz, &
!        &  g0, g1, g2, g3, g4, &
!        &  s_up, s_dn)

!     ! Number of input matrices
!     integer, intent(in) :: nmat

!     ! Size of diffuse gamma matrices
!     integer, intent(in) :: ndiff

!     ! Size of direct gamma matrices
!     integer, intent(in) :: ndir

!     ! Cosine of solar zenith angle
!     real(kind=jprb), intent(in) :: mu0

!     ! exp(-lambda*dz)
!     real(kind=jprb), intent(in) :: exp_lambda_dz(nmat,ndiff)

!     ! exp([Eigenvalues of gamma0]*dz)
!     real(kind=jprb), intent(in) :: exp_eigenval_dir_dz(nmat,ndir)

!     ! The unique submatrices of the eigenvector matrix of the full
!     ! Gamma matrix, "G_0", "G_1" etc. in Eq. 45
!     real(kind=jprb), intent(in), dimension(nmat,ndir,ndir)   :: g0
!     real(kind=jprb), intent(in), dimension(nmat,ndiff,ndiff) :: g1, g2
!     real(kind=jprb), intent(in), dimension(nmat,ndiff,ndir)  :: g3, g4

!     ! Radiation emerging from top and bottom of layer due to
!     ! scattering by the direct beam
!     real(kind=jprb), intent(out), dimension(nmat,ndiff, ndir) &
!          &  :: s_up, s_dn

!     real(kind=jprb) :: g_d(nmat,2*ndiff+ndir,2*ndiff+ndir)
!     real(kind=jprb) :: g_d2(nmat,ndiff,2*ndiff+ndir)
!     real(kind=jprb) :: cprime_dir(nmat,2*ndiff+ndir,ndir)
!     real(kind=jprb) :: rhs(nmat,2*ndiff+ndir,ndir)

!     integer :: jj

!     g_d = 0.0_jprb
!     ! Fill diffuse part of eigenvector matrix
!     g_d(:,1:ndiff,1:ndiff) = g1
!     g_d(:,ndiff+1:2*ndiff,1:ndiff) = g2 * spread(exp_lambda_dz,2,ndiff)
!     g_d(:,1:ndiff,ndiff+1:2*ndiff) = g_d(:,ndiff+1:2*ndiff,1:ndiff)
!     g_d(:,ndiff+1:2*ndiff,ndiff+1:2*ndiff) = g1
!     ! Fill direct part
!     g_d(:,2*ndiff+1:2*ndiff+ndir,2*ndiff+1:2*ndiff+ndir) = g0
!     ! Mixed part

! !    print *, g_d(:,1:ndiff,2*ndiff+1:2*ndiff+ndir)
! !    print *, g3
! !    print *, spread(exp_eigenval_dir_dz,2,ndiff)

!     g_d(:,1:ndiff,2*ndiff+1:2*ndiff+ndir) = g3 * spread(exp_eigenval_dir_dz,2,ndiff)
!     g_d(:,ndiff+1:2*ndiff,2*ndiff+1:2*ndiff+ndir) = g4

!     rhs = 0.0_jprb
!     !rhs = 1.0e-6_jprb
!     do jj = 1,ndir
!       rhs(:,2*ndiff+jj,jj) = 1.0_jprb
!     end do

!     cprime_dir = solve_rect_mat(nmat, 2*ndiff+ndir, ndir, g_d, rhs)

!     g_d2(:,1:ndiff,1:ndiff) = g1 * spread(exp_lambda_dz,2,ndiff)
!     g_d2(:,1:ndiff,ndiff+1:2*ndiff) = g2
!     g_d2(:,1:ndiff,2*ndiff+1:2*ndiff+ndir) = g3

!     s_up = rect_mat_x_mat(nmat, ndiff, 2*ndiff+ndir, ndir, g_d2, cprime_dir)
!     !s_up = s_up * mu0

!     g_d2(:,1:ndiff,1:ndiff) = g2
!     g_d2(:,1:ndiff,ndiff+1:2*ndiff) = g1 * spread(exp_lambda_dz,2,ndiff)
!     g_d2(:,1:ndiff,2*ndiff+1:2*ndiff+ndir) = g4 * spread(exp_eigenval_dir_dz,2,ndiff)

!     s_dn = rect_mat_x_mat(nmat, ndiff, 2*ndiff+ndir, ndir, g_d2, cprime_dir)
!     !s_dn = s_dn * mu0

!   end subroutine direct_diffuse_part


!     !---------------------------------------------------------------------
!   ! Return matrix X where AX=B. LU and A are n m-by-m matrices, and 
!   ! X and B are n m-by-l matrices.
!   function solve_rect_mat(n,m,l,A,B)
!     integer,    intent(in) :: n, m, l
!     real(jprb), intent(in) :: A(n,m,m)
!     real(jprb), intent(in) :: B(n,m,l)
!     real(jprb)             :: solve_rect_mat(n,m,l)

!     real(jprb) :: LU(n,m,m)

!     integer :: j

!     call lu_factorization(n,n,m,A,LU)

!     do j = 1, l
!       call lu_substitution(n,n,m,LU,B(1:n,1:m,j),solve_rect_mat(1:n,1:m,j))
!     end do

!   end function solve_rect_mat

!   subroutine schur_invert_sw(nmat, n0, n1, g0, g1, g2, g3, g0i, g1i, g2i, g3i)
!     integer(kind=jpim), intent(in)  :: nmat, n0, n1
!     real(kind=jprb),    intent(in)  :: g0(nmat,n0,n0), g1(nmat,n1,n1)
!     real(kind=jprb),    intent(in)  :: g2(nmat,n1,n1), g3(nmat,n1,n0)
!     real(kind=jprb),    intent(out) :: g0i(nmat,n0,n0), g1i(nmat,n1,n1)
!     real(kind=jprb),    intent(out) :: g2i(nmat,n1,n1), g3i(nmat,n1,n0)

!     g0i = invert(nmat,nmat,n0,g0)
!     g1i = invert(nmat,nmat,n1,g1 - mat_x_mat(nmat,nmat,n1,g2, &
!          &                           solve_mat(nmat,nmat,n1,g1,g2)))
!     g2i = mat_x_mat(nmat,nmat,n1,g1i,mat_x_mat(nmat,nmat,n1,g2, &
!          &                                  invert(nmat,nmat,n1,g1)))
!     g3i = rect_mat_x_mat(nmat,n1,n1,n0,g1i-g2i, &
!          &               rect_mat_x_mat(nmat,n1,n0,n0,g3,g0i))

!   end subroutine schur_invert_sw


!  subroutine calc_matrices_sw_eig(nmat, ndiff, ndir, dz, mu0, &
!        &  gamma0, gamma1, gamma2, gamma3, &
!        &  reflectance, transmittance, s_up, s_dn, trans_dir, &
!        &  int_dir, int_diff, int_dir_diff)
!     ! Inputs
!     ! Number of input matrices
!     integer, intent(in) :: nmat
!     ! Size of diffuse gamma matrices
!     integer, intent(in) :: ndiff
!     ! Size of direct gamma matrices
!     integer, intent(in) :: ndir
!     ! Layer thickness (m)
!     real(kind=jprb), intent(in) :: dz
!     ! Cosine of solar zenith angle
!     real(kind=jprb), intent(in) :: mu0
!     ! Exchange matrices; the unique parts of the matrix "Gamma" that
!     ! is not explicitly constructed:
!     !   Gamma = [-gamma1 -gamma2 -gamma3]
!     !           [+gamma2 +gamma1 +gamma3]
!     !           [                 gamma0]
!     ! Gamma in spartacus : (nmat, 9, 9) = (nmat, (3*nreg), (3*nreg))
!     real(kind=jprb), intent(in) :: gamma0(nmat,ndir,ndir)
!     real(kind=jprb), intent(in) :: gamma1(nmat,ndiff,ndiff)
!     real(kind=jprb), intent(in) :: gamma2(nmat,ndiff,ndiff)
!     real(kind=jprb), intent(in) :: gamma3(nmat,ndiff,ndir)
!     ! Outputs
!     ! Diffuse reflectance and transmittance matrices, "R" and "T" in
!     ! Eq. 49
!     real(kind=jprb), intent(out), dimension(nmat,ndiff,ndiff) &
!          &  :: reflectance, transmittance
  
!     ! Radiation emerging from top and bottom of layer due to
!     ! scattering by the direct beam, and the transmittance of the
!     ! layer to direct radiation with no scattering on the way
!     real(kind=jprb), intent(out), dimension(nmat,ndiff, ndir) &
!          &  :: s_up, s_dn
!     real(kind=jprb), intent(out), dimension(nmat,ndir,ndir) :: trans_dir
!     ! The following matrices are defined such that the integrated
!     ! diffuse flux (u_hat+v_hat in Eq. 29) is
!     !   int_diff * (u_conv+v_conv) + int_dir_diff * s_conv
!     ! where u_conv+v_conv is the convergence of diffuse fluxes into
!     ! the layer, i.e. the sum of the fluxes entering the layer from
!     ! the top or base, minus the fluxes leaving the layer from top or
!     ! base, and s_conv is the convergence of direct fluxes into the
!     ! layer, i.e. the direct flux at layer top minus the direct flux
!     ! at layer base. Likewise, the integrated direct flux is 
!     !   int_diff * s_conv
!     real(kind=jprb), intent(out), optional &
!          &  :: int_diff(nmat,ndiff,ndiff), &
!          &     int_dir(nmat,ndir,ndir), &
!          &     int_dir_diff(nmat,ndiff,ndir)
    
!     ! Local variables

!     ! Difference between gamma1 and gamma2
!     real(kind=jprb) :: gamma_diff(nmat,ndiff,ndiff)

!     ! The matrix product (gamma1-gamma2)*(gamma1+gamma2)
!     real(kind=jprb) :: gamma_product(nmat,ndiff,ndiff)

!     ! Eigenvalues/eigenvectors of gamma_product
!     real(kind=jprb) :: eigenval_prod(nmat,ndiff), &
!          &             eigenvec_prod(nmat,ndiff,ndiff)

!     ! Eigenvalues of gamma0
!     real(kind=jprb) :: eigenval_dir(nmat,ndir)
    
!     ! The unique submatrices of the eigenvector matrix of the full
!     ! Gamma matrix, "G_0", "G_1" etc. in Eq. 45
!     real(kind=jprb), dimension(nmat,ndir,ndir)   :: g0
!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: g1, g2
!     real(kind=jprb), dimension(nmat,ndiff,ndir)  :: g3, g4
    
!     ! g1 and g2 multiplied by the diagonal matrix "D" defined after
!     ! Eq. 46
!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: g1_d, g2_d

!     ! The square upper and lower parts of the rectangular matrix "C'"
!     ! in Eq. 47
!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: cprime_upper, cprime_lower

! #ifdef FAST_BUT_INCORRECT
!     ! Similar but for computing s_up and s_dn
!     real(kind=jprb), dimension(nmat,ndiff,ndir) :: cdir_upper, cdir_lower
! #endif

!     ! Positive eigenvalues of the full Gamma matrix
!     real(kind=jprb) :: lambda(nmat,ndiff)

!     ! exp(-lambda*dz)
!     real(kind=jprb) :: exp_lambda_dz(nmat,ndiff)

!     ! gamma2*g0
!     real(kind=jprb), dimension(nmat,ndiff,ndir) :: gamma3_g0

!     ! Inverse of g0
!     real(kind=jprb), dimension(nmat,ndir,ndir) :: g0_inv

!     ! DEFINE
!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: gamma1_d

!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: gamma2_inv_gamma1_d

! #ifdef FAST_BUT_INCORRECT
!     real(kind=jprb), dimension(nmat,ndiff,ndir) :: g3_d_inv_g0, g4_inv_g0
! #endif

!     ! If the full gamma matrix is
!     !        gamma = [ -gamma1  -gamma2  -gamma3  ]
!     !                [  gamma2   gamma1   gamma3  ]
!     !                [                    gamma0  ]
!     ! then its inverse is
!     !   inv(gamma) = [ -gamma1i -gamma2i -gamma3i ]
!     !                [  gamma2i  gamma1i -gamma3i ]
!     !                [                    gamma0i ]
!     ! and may be computed efficiently via the Schur-complement method
!     real(kind=jprb), dimension(nmat,ndir,ndir)   :: gamma0i
!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: gamma1i, gamma2i
!     real(kind=jprb), dimension(nmat,ndiff,ndir)  :: gamma3i

!     ! Temporary vector and matrix
!     real(kind=jprb), dimension(nmat,ndiff) :: tmp_vec
!     real(kind=jprb), dimension(nmat,ndiff,ndiff) :: tmp_mat

!     ! Loop index over matrix elements
!     integer :: jd, jo

    
!     ! Section 1: Compute the eigenvectors and eigenvalues of the
!     ! diffuse part of the Gamma matrix using the clever DISORT method
!     ! exploiting its regular structure

!     ! Compute the eigenvectors and eigenvalues of "gamma_product"
!     gamma_diff = gamma1-gamma2;
!     gamma_product = mat_x_mat(nmat,nmat,ndiff,gamma_diff,gamma1+gamma2)
! #ifdef USE_TIMING
!     ret =  gptlstart('eigen_decomposition_real')
! #endif 
!     call eigen_decomposition_real(ndiff,nmat,gamma_product, &
!          &                        eigenval_prod, eigenvec_prod)
! #ifdef USE_TIMING
!     ret =  gptlstop('eigen_decomposition_real')
! #endif 
!     ! EDIT PETER UKKONEN: problem with zero eigenvectors
!     ! print *, "EIGENVEC1", eigenvec_prod(:,1,1)
!     if (any(eigenvec_prod == 0.0_jprb)) then
!       call random_number(eigenvec_prod)
!     end if

!     ! Eigenvalues of the diffuse part of the full Gamma matrix are
!     ! (+/-)sqrt(eigenval_prod)
!     ! lambda = sqrt(max(0.0_jprb,eigenval_prod))
!     lambda = sqrt(max(0.0_jprb,eigenval_prod))
!     exp_lambda_dz = exp(-lambda*dz)

!     ! Compute the submatrices g1 and g2 of the eigenvector matrix of
!     ! the full Gamma matrix
!     tmp_mat = -solve_mat(nmat,nmat,ndiff,gamma_diff,eigenvec_prod)
!     tmp_mat = tmp_mat * spread(lambda,2,ndiff)
!     g1 = eigenvec_prod + tmp_mat
!     g2 = eigenvec_prod - tmp_mat
!     !g1 = -0.5*(eigenvec_prod + tmp_mat)
!     !g2 = -0.5*(eigenvec_prod - tmp_mat)
!     ! print *, "g1", g1(:,1,1)

!     ! Section 2: Compute diffuse reflectance and transmittance
!     ! matrices

!     ! Various equations require product of g1,g2 and exp(-lambda*dz)
!     g1_d = g1 * spread(exp_lambda_dz,2,ndiff)
!     g2_d = g2 * spread(exp_lambda_dz,2,ndiff)

!     ! Solve Eq. 48 to compute the upper and lower part of "C'", using
!     ! the Schur complement to exploit the regular structure of the
!     ! matrix on the left hand side of Eq. 48.
!     cprime_lower = invert(nmat,nmat,ndiff,(g1 &
!          &  -mat_x_mat(nmat,nmat,ndiff,g2_d, &
!          &             solve_mat(nmat,nmat,ndiff,g1,g2_d))))
!     cprime_upper = -solve_mat(nmat,nmat,ndiff,g1, &
!          &   mat_x_mat(nmat,nmat,ndiff,g2_d, cprime_lower))

!     ! Apply Eq. 49
!     reflectance   = mat_x_mat(nmat,nmat,ndiff,g1_d,cprime_upper) &
!          &        + mat_x_mat(nmat,nmat,ndiff,g2,  cprime_lower)
!     transmittance = mat_x_mat(nmat,nmat,ndiff,g2  ,cprime_upper) &
!          &        + mat_x_mat(nmat,nmat,ndiff,g1_d,cprime_lower)
! #ifdef USE_TIMING
!     ret =  gptlstart('eigen_decomposition_real2')
! #endif 
!     ! Section 3: Direct transmittance matrix: the matrix exponential
!     ! of gamma0, which we compute using eigen-decomposition
!     call eigen_decomposition_real(ndir,nmat,gamma0, &
!          &                        eigenval_dir, g0)
! #ifdef USE_TIMING
!     ret =  gptlstop('eigen_decomposition_real2')
! #endif 
!     !! EDIT PETER UKKONEN
!     if (any(g0 == 0.0_jprb)) then
!       call random_number(g0)
!     end if
!     g0_inv = invert(nmat,nmat,ndir,g0)
!     trans_dir = mat_x_mat(nmat,nmat,ndir, &
!          &  g0*spread(exp(eigenval_dir*dz),2,ndir), g0_inv)
   
!     ! Section 4: Mixed direct-diffuse part
!     gamma3_g0 = rect_mat_x_mat(nmat,ndiff,ndir,ndir,gamma3,g0)
!     gamma1_d = gamma1
!     do jd = 1,ndir
!        do jo = 1,ndiff
!           gamma1_d(:,jo,jo) = gamma1(:,jo,jo) + eigenval_dir(:,jd)
!        end do
!        gamma2_inv_gamma1_d = mat_x_mat(nmat,nmat,ndiff,gamma2,&
!             &                        invert(nmat,nmat,ndiff,gamma1_d))
!        tmp_mat = gamma1 - mat_x_mat(nmat,nmat,ndiff,&
!             &                       gamma2_inv_gamma1_d, gamma2)
!        do jo = 1,ndiff
!           tmp_mat(:,jo,jo) = tmp_mat(:,jo,jo) - eigenval_dir(:,jd)
!        end do
!        ! Subtract identity matrix
!        do jo = 1,ndiff
!           gamma2_inv_gamma1_d(:,jo,jo) = gamma2_inv_gamma1_d(:,jo,jo) - 1.0_jprb
!        end do
!        g4(:,:,jd) = solve_vec(nmat,nmat,ndiff,tmp_mat, &
!             &  mat_x_vec(nmat,nmat,ndiff,gamma2_inv_gamma1_d,gamma3_g0(:,:,jd)))
!        g3(:,:,jd) = -solve_vec(nmat,nmat,ndiff,gamma1_d, &
!             &  gamma3_g0(:,:,jd)+mat_x_vec(nmat,nmat,ndiff,gamma2,g4(:,:,jd)))
!     end do

!     call direct_diffuse_part(nmat, ndiff, ndir, mu0, &
!          &  exp_lambda_dz, exp(eigenval_dir*dz), &
!          &  g0, g1, g2, g3, g4, &
!          &  s_up, s_dn)

! #ifdef FAST_BUT_INCORRECT

!     g3_d_inv_g0 = rect_mat_x_mat(nmat,ndiff,ndir,ndir,g3*spread(exp(eigenval_dir*dz),2,ndir), &
!          &                  g0_inv)
!     g4_inv_g0 = rect_mat_x_mat(nmat,ndiff,ndir,ndir,g4, g0_inv)

!     ! Schur complement
!     tmp_mat = g1 - mat_x_mat(nmat,nmat,ndiff,g2_d, &
!          &    solve_mat(nmat,nmat,ndiff,g1,g2_d))
!     cdir_upper = -solve_rect_mat(nmat,ndiff,ndir,g1,g3_d_inv_g0 &
!          &  + rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2_d, &
!          &  solve_rect_mat(nmat,ndiff,ndir,tmp_mat, &
!          &  rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2_d, &
!          &  rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g1,g3_d_inv_g0) &
!          &  - g4_inv_g0))))
!     cdir_lower = solve_rect_mat(nmat,ndiff,ndir,tmp_mat, &
!          &  rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2_d, &
!          &  solve_rect_mat(nmat,ndiff,ndir,g1,g3_d_inv_g0) &
!          &  + g4_inv_g0))
!     s_up = mu0 * (rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g1_d,cdir_upper) &
!          &      + rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2,  cdir_lower) &
!          &      + rect_mat_x_mat(nmat,ndiff,ndir, ndir,g3,  g0_inv))
!     s_dn = mu0 * (rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2,  cdir_upper) &
!          &      + rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g1_d,cdir_lower) &
!          &      + rect_mat_x_mat(nmat,ndiff,ndir, ndir, &
!          &          g4*spread(exp(eigenval_dir*dz),2,ndir),  g0_inv))

! #endif
    
!     if (present(int_dir) .or. present(int_diff) &
!          &                    .or. present(int_dir_diff)) then
!       call schur_invert_sw(nmat,ndir,ndiff,gamma0, gamma1, gamma2, gamma3, &
!            &                               gamma0i,gamma1i,gamma2i,gamma3i)
!       int_dir      = -gamma0i
!       int_diff     = gamma2i - gamma1i
!       int_dir_diff = 2.0_jprb * gamma3i
!     end if

!   end subroutine calc_matrices_sw_eig

!  ! This routine performs an Eigen decomposition of "nmat" matrices
!   ! dimensioned "norder-x-norder", stored in the input array "amat".
!   ! The resulting Eigenvalues are stored in "eigenvalue" and
!   ! Eigenvectors as the column vectors comprising "eigenvector". This
!   ! routine assumes that the eigenvalues must be real, appropriate for
!   ! radiative transfer problems.  If EIGEN_OUTER_STACKING is defined
!   ! then the matrix number is the outer dimension (slowest varying in
!   ! memory) of the input and output arrays, otherwise it is the inner
!   ! dimension (fastest varying in memory).
!   subroutine eigen_decomposition_real(norder, nmat, amat, &
!        &                              eigenvalue, eigenvector, &
!        &                              nerror, ierror)

!     ! "jprb" is the precision in the input and output data and may be
!     ! single or double. "jprd" is double precision for local
!     ! variables.

!     implicit none

!     ! Constants

!     ! Tolerance
!     real(kind=jprd) :: Tol = epsilon(1.0_jprd)

!     ! Inputs

!     ! Order of matrix
!     integer,         intent(in)  :: norder
!     ! Number of matrices to decompose
!     integer,         intent(in)  :: nmat

!     ! Matrices to be decomposed
! #ifdef EIGEN_OUTER_STACKING
!     real(kind=jprb), intent(in)  :: amat(norder,norder,nmat)
! #else
!     real(kind=jprb), intent(in)  :: amat(nmat,norder,norder)
! #endif
    
!     ! Outputs

!     ! Returned eigenvalues and eigenvectors
! #ifdef EIGEN_OUTER_STACKING
!     real(kind=jprb), intent(out) :: eigenvalue(norder,nmat)
!     real(kind=jprb), intent(out) :: eigenvector(norder,norder,nmat)
! #else
!     real(kind=jprb), intent(out) :: eigenvalue(nmat,norder)
!     real(kind=jprb), intent(out) :: eigenvector(nmat,norder,norder)
! #endif
    
!     ! Return the number of matrices that could not be eigen-decomposed
!     integer, optional, intent(out) :: nerror

!     ! For each of the nmat matrices that could not be fully
!     ! eigen-decomposed, a positive value is added to ierror indicating
!     ! that eigenvalues ierror+1:norder were found but eigenvalues
!     ! 1:ierror were not found.  On success, a zero entry is added to
!     ! ierror.
!     integer, optional, intent(out) :: ierror(nmat)

!     ! Parameters
!     real(kind=jprd), parameter :: &
!          &  C1 = 0.4375_jprd, &
!          &  C2 = 0.5_jprd, &
!          &  C3 = 0.75_jprd, &
!          &  C4 = 0.95_jprd, &
!          &  C5 = 16.0_jprd, &
!          &  C6 = 256.0_jprd
    
!     ! Local variables

!     ! Balanced version of amat
!     real(kind=jprd) :: abal(norder,norder)

!     ! Double-precision eigenvalues and eigenvectors
!     real(kind=jprd) :: eigenval(norder)
!     real(kind=jprd) :: eigenvec(norder,norder)

!     ! Discriminant of polynomial, and a function of it
!     real(kind=jprd) :: discriminant, half_sqrt_disc

!     ! Loop index for matrices
!     integer :: jmat

!     ! Other loop indices
!     integer :: ji, jj, jn

!     ! Other indices
!     integer :: ka, kk, ll, lb, lll, n1, n2, in, kkk

!     ! Local scalars
!     real(kind=jprd) :: column, ff, gg, hh
!     real(kind=jprd) :: pp, qq, rr, tmp, rnorm, row, ss, scale, tt
!     real(kind=jprd) :: uu, vv, ww, xx, yy, zz

!     ! Work space
!     real(kind=jprd) :: wkd(2*norder)

!     logical :: not_finished, not_found, not_last, is_error

!     if (present(nerror)) then
!       nerror = 0;
!     end if
!     if (present(ierror)) then
!       ierror(:) = 0
!     end if

! #ifdef USE_TIMING
!      ret =  gptlstart('eigen_decomposition_real')
! #endif 
    
!     if (norder > 2) then
!       ! Nominal case: matrices of order larger than 2

!       ! Loop over input matrices
!       do jmat = 1,nmat

!         is_error = .false.

!         ! Initialize outputs
!         eigenval(:)    = 0.0_jprd
!         eigenvec(:,:)  = 0.0_jprd
!         do ji = 1,norder
!           eigenvec(ji,ji) = 1.0_jprd
!         end do

!         ! Balance the input matrix and reduce its norm by diagonal
!         ! similarity transformation stored in WK; then search for rows
!         ! isolating an eigenvalue and push them down

!         ! Initialize balanced matrix
! #ifdef EIGEN_OUTER_STACKING        
!         abal = amat(:,:,jmat)
! #else
!         abal = amat(jmat,:,:)
! #endif

!         rnorm = 0.0_jprd

!         ll = 1
!         kk = norder

!         ! The flow of the following is a bit complicated, but it is
!         ! better than the original which involved GOTO statements
!         not_finished = .true.

!         do while(not_finished) 
!           not_finished = .false.

!           kkk = kk

!           do jj = kkk,1,-1
!             row = 0.0_jprd

!             do ji = 1,kk
!               if (ji /= jj) then
!                 row = row + abs(abal(jj,ji))
!               end if
!             end do

!             if (row == 0.0_jprd) then

!               wkd(kk) = jj
!               if (ji /= kk) then
!                 do ji = 1,kk
!                   tmp         = abal(ji,jj)
!                   abal(ji,jj) = abal(ji,kk)
!                   abal(ji,kk) = tmp
!                 end do
!                 do ji = ll,norder
!                   tmp         = abal(jj,ji)
!                   abal(jj,ji) = abal(kk,ji)
!                   abal(kk,ji) = tmp
!                 end do
!               end if

!               kk = kk - 1
!               not_finished = .true.
!               exit ! Quit the jj loop and repeat
!             else
!               not_finished = .false.
!             end if

!           end do ! jj
!         end do ! while not_finished

!         ! Search for columns isolating an eigenvalue and push them
!         ! left

!         not_finished = .true.
!         do while(not_finished)
!           not_finished = .false.
!           lll = ll
!           do jj = lll,kk
!             column = 0.0_jprd

!             do ji = ll, kk
!               if (ji /= jj) then
!                 column = column + abs(abal(ji,jj))
!               end if
!             end do

!             if (column == 0.0_jprd) then
!               wkd(ll) = jj
!               if (jj /= ll) then
!                 do ji = 1,kk
!                   tmp         = abal(ji,jj)
!                   abal(ji,jj) = abal(ji,ll)
!                   abal(ji,ll) = tmp
!                 end do

!                 do ji = ll, norder
!                   tmp         = abal(jj,ji)
!                   abal(jj,ji) = abal(ll,ji)
!                   abal(ll,ji) = tmp               
!                 end do
!               end if
!               ll = ll + 1
!               not_finished = .true.
!               exit ! Quit the jj loop and repeat
!             else
!               not_finished = .false.
!             end if
!           end do
!         end do

!         ! Balance the submatrix in rows L through K
!         do ji = ll,kk
!           wkd(ji) = 1.0_jprd
!         end do

!         not_finished = .true.
!         do while (not_finished)
!           not_finished = .false.
!           do ji = ll,kk
!             column = 0.0_jprd
!             row    = 0.0_jprd
!             do jj = ll,kk
!               if (jj /= ji) then
!                 column = column + abs(abal(jj,ji))
!                 row    = row    + abs(abal(ji,jj))
!               end if
!             end do

!             ff = 1.0_jprd
!             gg = row / C5
!             hh = column + row

!             do while (column < gg)
!               ff = ff*C5
!               column = column*C6
!             end do

!             gg = row*C5

!             do while (column > gg)
!               ff = ff / C5
!               column = column/C6
!             end do

!             if ((column + row) / ff < C4*hh) then
!               wkd(ji) = wkd(ji)*ff
!               not_finished = .true.

!               do jj = ll,norder
!                 abal(ji,jj) = abal(ji,jj) / ff
!               end do
!               do jj = 1,kk
!                 abal(jj,ji) = abal(jj,ji) * ff
!               end do
!             end if
!           end do
!         end do ! not_finished

!         ! Is A already in Hessenberg form?
!         if (kk - 1 >= ll + 1) then
!           ! No: convert A to Hessenberg form

!           do jn = ll+1, kk-1
!             hh = 0.0_jprd
!             wkd(jn+norder) = 0.0_jprd
!             scale = 0.0_jprd

!             ! Scale column
!             do ji = jn,kk
!               scale = scale + abs(abal(ji,jn-1))
!             end do

!             if (scale /= 0.0_jprd) then

!               do ji = kk,jn,-1
!                 wkd(ji+norder) = abal(ji,jn-1) / scale
!                 hh = hh + wkd(ji+norder)*wkd(ji+norder)
!               end do

!               gg = -sign(sqrt(hh), wkd(jn+norder))
!               hh = hh - wkd(jn+norder)*gg
!               wkd(jn+norder) = wkd(jn+norder) - gg

!               hh = 1.0_jprd / hh
!               ! Form (I-(U*UT)*H)*A
!               do jj = jn,norder
!                 ff = 0.0_jprd
!                 do ji = kk,jn,-1
!                   ff = ff + wkd(ji+norder)*abal(ji,jj)
!                 end do
!                 do ji = jn,kk
!                   abal(ji,jj) = abal(ji,jj) - wkd(ji+norder)*ff*hh
!                 end do

!               end do

!               do ji = 1,kk
!                 ff = 0.0_jprd
!                 do jj = kk,jn,-1
!                   ff = ff + wkd(jj+norder)*abal(ji,jj)
!                 end do
!                 do jj = jn,kk
!                   abal(ji,jj) = abal(ji,jj) - wkd(jj+norder)*ff*hh
!                 end do
!               end do

!               wkd(jn+norder)  = scale*wkd(jn+norder)
!               abal(jn,jn-1) = scale*gg
!             end if
!           end do

!           do jn = kk-2,ll,-1 ! 360
!             n1 = jn+1
!             n2 = jn+2
!             ff = abal(n1,jn)

!             if (ff /= 0.0_jprd) then
!               ff = ff*wkd(jn+1+norder)
!               do ji = n2,kk
!                 wkd(ji+norder) = abal(ji,jn)
!               end do

!               if (n1 < kk) then
!                 do jj = 1,norder
!                   gg = 0.0_jprd
!                   do ji = n1,kk
!                     gg = gg + wkd(ji+norder)*eigenvec(ji,jj)
!                   end do

!                   gg = gg / ff

!                   do ji = n1,kk
!                     eigenvec(ji,jj) = eigenvec(ji,jj) + gg*wkd(ji+norder)
!                   end do
!                 end do

!               end if
!             end if

!           end do

!         end if

!         jn = 1
!         do ji = 1,norder
!           do jj = jn,norder
!             rnorm = rnorm + abs(abal(ji,jj))
!           end do
!           jn = ji
!           if (ji < ll .or. ji > kk) then
!             eigenval(ji) = abal(ji,ji)
!           end if
!         end do
!         jn = kk
!         tt = 0.0_jprd

!         ! Search for next eigenvalue
!         not_finished = .true.
!         do while (not_finished .and. jn >= ll)
!           not_finished = .false.
!           in = 0
!           n1 = jn-1
!           n2 = jn-2

!           ! Look for single small sub-diagonal element
!           not_found = .true.
!           do while (not_found) ! 410
!             !not_found = .false.

!             do ji = ll,jn
!               lb = jn + ll - ji
!               if (lb == ll) then
!                 exit
!               end if
!               ss = abs(abal(lb-1,lb-1)) + abs(abal(lb,lb))
!               if (ss == 0.0_jprd) then
!                 ss = rnorm
!               end if
!               if (abs(abal(lb,lb-1)) < tol*ss) then
!                 exit
!               end if
!             end do

!             xx = abal(jn,jn)
!             if (lb == jn) then
!               abal(jn,jn) = xx + tt
!               eigenval(jn) = abal(jn,jn)
!               jn = n1
!               not_finished = .true.
!               exit
!             end if

!             yy = abal(n1,n1)
!             ww = abal(jn,n1) * abal(n1,jn)

!             if (lb == n1) then
!               ! Two eigenvalues found
!               pp = (yy - xx) * C2
!               qq = pp*pp + ww
!               zz = sqrt(abs(qq))
!               abal(jn,jn) = xx + tt
!               xx = abal(jn,jn)
!               abal(n1,n1) = yy + tt
!               ! Real pair
!               zz = pp + sign(zz,pp)
!               eigenval(n1) = xx + zz
!               eigenval(jn) = eigenval(n1)

!               if (zz /= 0.0_jprd) then
!                 eigenval(jn) = xx - ww/zz
!               end if

!               xx = abal(jn,n1)
!               ! Employ scale factor in case xx and zz are very
!               ! small
!               rr = 1.0_jprd / sqrt(xx*xx + zz*zz)
!               pp = xx * rr
!               qq = zz * rr

!               ! Row modification
!               do jj = n1,norder
!                 zz = abal(n1,jj)
!                 abal(n1,jj) = qq*zz + pp*abal(jn,jj)
!                 abal(jn,jj) = qq*abal(jn,jj) - pp*zz
!               end do

!               ! Column modification
!               do ji = 1,jn
!                 zz = abal(ji,n1)
!                 abal(ji,n1) = qq*zz + pp*abal(ji,jn)
!                 abal(ji,jn) = qq*abal(ji,jn) - pp*zz
!               end do

!               ! Accumulate transforms
!               do ji = ll,kk
!                 zz = eigenvec(ji,n1)
!                 eigenvec(ji,n1) = qq*zz + pp*eigenvec(ji,jn)
!                 eigenvec(ji,jn) = qq*eigenvec(ji,jn) - pp*zz
!               end do
!               jn = n2
!               not_finished = .true.
!               exit
!             end if

!             if (in == 30) then
!               ! no convergence after 30 iterations; set error
!               ! indicator to the index of the current eigenvalue
!               if (present(ierror)) then
!                 ierror(jmat) = jn
!               end if
!               if (present(nerror)) then
!                 nerror = nerror + 1
!               end if
!               is_error = .true.
!               not_finished = .false.
!               exit
!             end if

!             ! Form shift
!             if (in == 10 .or. in == 20) then
!               tt = tt+xx
!               do ji = ll,jn
!                 abal(ji,ji) = abal(ji,ji) - xx
!               end do

!               ss = abs(abal(jn,n1)) + abs(abal(n1,n2))
!               xx = C3*ss
!               yy = xx
!               ww = -C1 * ss*ss
!             end if

!             in = in + 1
!             ! Look for two consecutive small sub-diagonal
!             ! elements
!             do jj = lb,n2
!               ji = n2+lb-jj
!               zz = abal(ji,ji)
!               rr = xx-zz
!               ss = yy-zz
!               pp = (rr*ss - ww) / abal(ji+1,ji) + abal(ji,ji+1)
!               qq = abal(ji+1,ji+1) -zz-rr-ss
!               rr = abal(ji+2,ji+1)
!               ss = 1.0_jprd / (abs(pp) + abs(qq) + abs(rr))
!               pp = pp * ss
!               qq = qq * ss
!               rr = rr * ss
!               if (ji == lb) then
!                 exit
!               end if
!               uu = abs(abal(ji,ji-1)) * (abs(qq) + abs(rr))
!               vv = abs(pp) * (abs(abal(ji-1,ji-1)) + abs(zz) &
!                    &              +abs(abal(ji+1,ji+1)))
!               if (uu <= tol*vv) then
!                 exit
!               end if
!             end do

!             abal(ji+2,ji) = 0.0_jprd
!             do jj = ji+3,jn
!               abal(jj,jj-2) = 0.0_jprd
!               abal(jj,jj-3) = 0.0_jprd
!             end do

!             ! Double QR step involving rows K to N and columns M
!             ! to N
!             do ka = ji,n1
!               not_last = ka /= n1
!               if (ka == ji) then
!                 ss = sign(sqrt(pp*pp + qq*qq + rr*rr), pp)
!                 if (lb /= ji) then
!                   abal(ka,ka-1) = -abal(ka,ka-1)
!                 end if
!               else
!                 pp = abal(ka,ka-1)
!                 qq = abal(ka+1,ka-1)
!                 rr = 0.0_jprd

!                 if (not_last) then
!                   rr = abal(ka+2,ka-1)
!                 end if

!                 xx = abs(pp) + abs(qq) + abs(rr)
!                 if (xx == 0.0_jprd) then
!                   cycle
!                 end if
!                 pp = pp / xx
!                 qq = qq / xx
!                 rr = rr / xx
!                 ss = sign(sqrt(pp*pp + qq*qq + rr*rr), pp)
!                 abal(ka,ka-1) = -ss*xx
!               end if

!               pp = pp + ss
!               ss = 1.0_jprd / ss
!               xx = pp*ss
!               yy = qq*ss
!               zz = rr*ss
!               pp = 1.0_jprd / pp
!               qq = qq*pp
!               rr = rr*pp

!               ! Row modification
!               do jj = ka,norder
!                 pp = abal(ka,jj) + qq*abal(ka+1,jj)

!                 if (not_last) then
!                   pp = pp + rr*abal(ka+2,jj)
!                   abal(ka+2,jj) = abal(ka+2,jj) - pp*zz
!                 end if

!                 abal(ka+1,jj) = abal(ka+1,jj) - pp*yy
!                 abal(ka,jj)   = abal(ka,jj)   - pp*xx
!               end do

!               ! Column modification
!               do ji = 1,min(jn,ka+3)
!                 pp = xx*abal(ji,ka) + yy*abal(ji,ka+1)

!                 if (not_last) then
!                   pp = pp + zz*abal(ji,ka+2)
!                   abal(ji,ka+2) = abal(ji,ka+2) - pp*rr
!                 end if

!                 abal(ji,ka+1) = abal(ji,ka+1) - pp*qq
!                 abal(ji,ka)   = abal(ji,ka)   - pp
!               end do

!               ! Accumulate transformations
!               do ji = ll,kk
!                 pp = xx*eigenvec(ji,ka) + yy*eigenvec(ji,ka+1)

!                 if (not_last) then
!                   pp = pp + zz*eigenvec(ji,ka+2)
!                   eigenvec(ji,ka+2) = eigenvec(ji,ka+2) - pp*rr
!                 end if
!                 eigenvec(ji,ka+1) = eigenvec(ji,ka+1) - pp*qq
!                 eigenvec(ji,ka)   = eigenvec(ji,ka)   - pp

!               end do
!             end do
!           end do ! do while not_found
!         end do ! do while not_finished
!         ! All eigenvalues found, now backsubstitute real vector
        
!         ! If an error occurred with this matrix, save what we have and
!         ! move to next one
!         if (.not. is_error) then

!           if (rnorm /= 0.0_jprd) then

!             do jn = norder,1,-1
!               n2 = jn
!               abal(jn,jn) = 1.0_jprd

!               do ji = jn-1,1,-1
!                 ww = abal(ji,ji) - eigenval(jn)
!                 if (abs(ww) < abs(Tol*rnorm)) then
!                   ww = sign(Tol*rnorm, ww)
!                 end if
!                 rr = abal(ji,jn)
!                 do jj = n2,jn-1
!                   rr = rr + abal(ji,jj) * abal(jj,jn)
!                 end do
!                 abal(ji,jn) = -rr / ww
!                 n2 = ji
!               end do
!             end do

!             ! End backsubstitution vectors of isolated evals
!             do ji = 1,norder
!               if (ji < ll .or. ji > kk) then
!                 do jj = ji,norder
!                   eigenvec(ji,jj) = abal(ji,jj)
!                 end do
!               end if
!             end do

!             ! Multiply by transformation matrix
!             if (kk /= 0.0_jprd) then
!               do jj = norder,ll,-1
!                 do ji = ll,kk
!                   zz = 0.0_jprd
!                   do jn = ll,min(jj,kk)
!                     zz = zz + eigenvec(ji,jn) * abal(jn,jj)
!                   end do
!                   eigenvec(ji,jj) = zz
!                 end do
!               end do
!             end if
!           end if

!           do ji = ll,kk
!             do jj = 1,norder
!               eigenvec(ji,jj) = eigenvec(ji,jj) * wkd(ji)
!             end do
!           end do

!           ! Interchange rows if permutations occurred

!           do ji = ll-1,1,-1
!             jj = nint(wkd(ji))

!             if (ji < jj) then
!               do jn = 1,norder
!                 tmp = eigenvec(ji,jn)
!                 eigenvec(ji,jn) = eigenvec(jj,jn)
!                 eigenvec(jj,jn) = tmp
!               end do
!             end if
!           end do

!           do ji = kk+1,norder
!             jj = nint(wkd(ji))

!             if (ji /= jj) then
!               do jn = 1,norder
!                 tmp = eigenvec(ji,jn)
!                 eigenvec(ji,jn) = eigenvec(jj,jn)
!                 eigenvec(jj,jn) = tmp
!               end do
!             end if

!           end do
!         end if

!         ! Put results into output arrays
! #ifdef EIGEN_OUTER_STACKING
!         eigenvalue(:,jmat) = eigenval(:)
!         eigenvector(:,:,jmat) = eigenvec(:,:)
! #else
!         eigenvalue(jmat,:) = eigenval(:)
!         eigenvector(jmat,:,:) = eigenvec(:,:)
! #endif
        
!       end do ! Loop over matrices

!     else if (norder == 2) then
!       ! Special case: 2x2 matrices
!       do jmat = 1,nmat
        
! #ifdef EIGEN_OUTER_STACKING
!         discriminant = (amat(1,1,jmat)-amat(2,2,jmat)) ** 2 &
!              & + 4.0_jprd*amat(1,2,jmat)*amat(2,1,jmat)
!         if (discriminant < 0.0_jprd) then
!           if (present(nerror)) then
!             nerror = nerror + 1
!           end if
!           if (present(ierror)) then
!             ierror(jmat) = 2
!           end if
!         end if

!         eigenvalue(1,jmat) = 0.5_jprd*(amat(1,1,jmat)+amat(2,2,jmat))
!         eigenvalue(2,jmat) = eigenvalue(1,jmat)
!         half_sqrt_disc = 0.5_jprd*sqrt(discriminant)
!         if (amat(1,1,jmat) >= amat(2,2,jmat)) then
!           eigenvalue(1,jmat) = eigenvalue(1,jmat) + half_sqrt_disc
!           eigenvalue(2,jmat) = eigenvalue(2,jmat) - half_sqrt_disc
!         else
!           eigenvalue(1,jmat) = eigenvalue(1,jmat) - half_sqrt_disc
!           eigenvalue(2,jmat) = eigenvalue(2,jmat) + half_sqrt_disc
!         end if
!         eigenvector(1,1,jmat) = 1.0_jprd
!         eigenvector(2,2,jmat) = 1.0_jprd
!         if (amat(1,1,jmat) == amat(2,2,jmat) .and. &
!              &  (amat(2,1,jmat) == 0.0_jprd .or. amat(1,2,jmat) == 0.0_jprd)) then
!           rnorm = 1.0_jprd / (Tol * &
!                &  abs(amat(1,1,jmat)) + abs(amat(2,1,jmat)) &
!                & +abs(amat(1,2,jmat)) + abs(amat(2,2,jmat)))
!           eigenvector(2,1,jmat) = amat(2,1,jmat) * rnorm
!           eigenvector(1,2,jmat) = amat(1,2,jmat) * rnorm
!         else
!           eigenvector(2,1,jmat) = amat(2,1,jmat) / (eigenvalue(1,jmat)-amat(2,2,jmat))
!           eigenvector(1,2,jmat) = amat(1,2,jmat) / (eigenvalue(2,jmat)-amat(1,1,jmat))
!         end if
! #else
!         discriminant = (amat(jmat,1,1)-amat(jmat,2,2)) ** 2 &
!              & + 4.0_jprd*amat(jmat,1,2)*amat(jmat,2,1)
!         if (discriminant < 0.0_jprd) then
!           if (present(nerror)) then
!             nerror = nerror + 1
!           end if
!           if (present(ierror)) then
!             ierror(jmat) = 2
!           end if
!         end if

!         eigenvalue(jmat,1) = 0.5_jprd*(amat(jmat,1,1)+amat(jmat,2,2))
!         eigenvalue(jmat,2) = eigenvalue(jmat,1)
!         half_sqrt_disc = 0.5_jprd*sqrt(discriminant)
!         if (amat(jmat,1,1) >= amat(jmat,2,2)) then
!           eigenvalue(jmat,1) = eigenvalue(jmat,1) + half_sqrt_disc
!           eigenvalue(jmat,2) = eigenvalue(jmat,2) - half_sqrt_disc
!         else
!           eigenvalue(jmat,1) = eigenvalue(jmat,1) - half_sqrt_disc
!           eigenvalue(jmat,2) = eigenvalue(jmat,2) + half_sqrt_disc
!         end if
!         eigenvector(jmat,1,1) = 1.0_jprd
!         eigenvector(jmat,2,2) = 1.0_jprd
!         if (amat(jmat,1,1) == amat(jmat,2,2) .and. &
!              &  (amat(jmat,2,1) == 0.0_jprd .or. amat(jmat,1,2) == 0.0_jprd)) then
!           rnorm = 1.0_jprd / (Tol * &
!                &  abs(amat(jmat,1,1)) + abs(amat(jmat,2,1)) &
!                & +abs(amat(jmat,1,2)) + abs(amat(jmat,2,2)))
!           eigenvector(jmat,2,1) = amat(jmat,2,1) * rnorm
!           eigenvector(jmat,1,2) = amat(jmat,1,2) * rnorm
!         else
!           eigenvector(jmat,2,1) = amat(jmat,2,1) / (eigenvalue(jmat,1)-amat(jmat,2,2))
!           eigenvector(jmat,1,2) = amat(jmat,1,2) / (eigenvalue(jmat,2)-amat(jmat,1,1))
!         end if
! #endif
        
!       end do ! Loop over matrices

!     else if (norder == 1) then
!       ! Special case: 1x1 matrices
! #ifdef EIGEN_OUTER_STACKING
!       eigenvalue(1,1:nmat) = amat(1,1,1:nmat)
!       eigenvector(1,1,1:nmat) = 1.0_jprd
! #else
!       eigenvalue(1:nmat,1) = amat(1:nmat,1,1)
!       eigenvector(1:nmat,1,1) = 1.0_jprd
! #endif     
!     else
!       ! norder not a positive number
!       if (present(nerror)) then
!         nerror = nerror + 1
!       end if

!       if (present(ierror)) then
!         ierror(:) = 1
!       end if
!     end if

! #ifdef USE_TIMING
!      ret =  gptlstop('eigen_decomposition_real')
! #endif 

!   end subroutine eigen_decomposition_real

! ------------------------------------------------------------------------------------------------------------------

! !!! EXPM WITH G-POINTS AS THE LAST DIMENSION; SO DOING MANY SMALL MATRIX-MATRIX MULTIPLICATIONS
! !!! dimensions are of order (9,9,112) rather than (112,9,9) in reference, where 112 is the independent g-point dim. 
! !!! and the matrix exponential is computed for the (9,9) size matrices in both cases.
! !!! this solution has so far been slower than the reference.
! !!! expm3 allows to use external libraries for the matrix-matrix multiplications, in particular,
! !!! LIBXSMM optimized for small matrices. However, this is not really faster than the handwritten
! !!! matrix kernels, because 1) they utilize the special matrix structure to avoid redundant computations
! !!! 2) the reference solution allows to vectorize over g-points which is the larger dimension
!   subroutine expm3(n,iend,m,A,i_matrix_pattern)

!     use yomhook, only : lhook, dr_hook

!     integer,    intent(in)      :: n, m, iend
!     real(jprb), intent(inout)   :: A(m,m,iend)
!     integer,    intent(in)      :: i_matrix_pattern

!     real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
!          &                                1.880152677804762e+00_jprb, &
!          &                                3.925724783138660e+00_jprb/) 
!     real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
!          &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
!          &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

!     real(jprb), dimension(m,m,iend) :: A2, A4, A6
!     real(jprb), dimension(m,m,iend) :: U, V

!     real(jprb) :: normA, sum_column(m)

!     integer    :: j1, j2, j3
!     real(jprb) :: frac
!     integer    :: expo(iend)
!     real(jprb) :: scaling

!     real(jprb) :: hook_handle

!     if (lhook) call dr_hook('radiation_matrix:expm',0,hook_handle)

!     ! normA = 0.0_jprb

!     ! Compute the 1-norms of A

!     do j3 = 1,iend
!         sum_column(:) = 0.0_jprb
!         do j2 = 1,m
!             sum_column(:) = sum_column(:) + abs(A(:,j2,j3))
!         end do
!         normA = 0.0_jprb
!         do j1 = 1,m
!             normA= max(normA,sum_column(j1))
!         end do

!         frac = fraction(normA/theta(3))
!         expo(j3) = exponent(normA/theta(3))
!         if (frac== 0.5_jprb) then
!             expo(j3) = expo(j3) - 1
!         end if

!         if (expo(j3) < 0) then
!             expo(j3) = 0
!         end if

!         ! Scale the input matrices by a power of 2
!         scaling = 2.0_jprb**(-expo(j3))
!         do j2 = 1,m
!             A(:,j2,j3) = A(:,j2,j3) * scaling
!         end do
!     end do
!       ! Pade approximant of degree 7

!       ! A2 = mat_x_mat3(n,iend,m,A, A, i_matrix_pattern)
!       ! A4 = mat_x_mat3(n,iend,m,A2,A2,i_matrix_pattern)
! #ifdef USE_TIMING
!     ret =  gptlstart('expm_Pade_mat_x_mat3')
! #endif 
!       ! A2 = mat_x_mat3(iend,m,A,A)
!       ! A4 = mat_x_mat3(iend,m,A2,A2)
! #ifdef USE_TIMING
!     ret =  gptlstart('mat_squares')
! #endif 
!       call mat_x_mat3_square(iend,m,A,A2)
!       call mat_x_mat3_square(iend,m,A2,A4)
! #ifdef USE_TIMING
!     ret =  gptlstop('mat_squares')
! #endif 
!       call mat_x_mat3(iend,m,A2,A4,A6)
!       ! A6 = mat_x_mat3(iend,m,A2,A4)
      
!       ! mat_x_mat3(nmat,nmat,ndiff,g2_d,
!       V = c(8)*A6 + c(6)*A4 + c(4)*A2
!       do j3 = 1,iend
!         do j1 = 1,m
!             V(j1,j1,j3) = V(j1,j1,j3) + c(2)
!         end do
!       end do
!       ! U = mat_x_mat3(iend,m,A,V)
!       call mat_x_mat3(iend,m,A,V,U)
! #ifdef USE_TIMING
!     ret =  gptlstop('expm_Pade_mat_x_mat3')
! #endif 
!       V = c(7)*A6 + c(5)*A4 + c(3)*A2
!       ! Add a multiple of the identity matrix
!       do j3 = 1,iend
!         do j1 = 1,m
!             V(j1,j1,j3) = V(j1,j1,j3) + c(1)
!         end do
!       end do

!       V = V-U
!       U = 2.0_jprb*U
! #ifdef USE_TIMING
!     ret =  gptlstart('expm_solve_mat3')
! #endif 
!       call solve_mat3(iend,m,V,U,A)
! #ifdef USE_TIMING
!     ret =  gptlstop('expm_solve_mat3')
! #endif 
!       ! Add the identity matrix
!       do j3 = 1,iend
!         do j1 = 1,m
!             A(j1,j1,j3) = A(j1,j1,j3) + 1.0_jprb
!         end do
!       end do

!       ! Loop through the matrices
!       ! print *, " MEAN EXPO:", sum(expo)/iend
! #ifdef USE_TIMING
!     ret =  gptlstart('expm_square_matrix_repeated2')
! #endif 
!       do j3 = 1,iend
!         do j1 = 1,expo(j3)
!             ! Square matrix j1 expo(j1) times  
!             A(:,:,j3) = mat_x_mat2(m,A(:,:,j3),A(:,:,j3))      
!         end do
!       end do
! #ifdef USE_TIMING
!     ret =  gptlstop('expm_square_matrix_repeated2')
! #endif 
!       if (lhook) call dr_hook('radiation_matrix:expm',1,hook_handle)

!     ! end do

    
!   end subroutine expm3

! pure subroutine mat_x_mat3(iend,h,A,B,C)
!     integer,    intent(in)                      :: iend,h
!     real(jprb), intent(in),   dimension(9,9,iend) :: A,B
!     real(jprb), intent(out),  dimension(9,9,iend) :: C
!     integer    :: j1, j2, j3, jg,m
!     integer    :: mblock, m2block

!     real(jprb) :: hook_handle
    
  
!     ! Array-wise assignment
!     ! mat_x_mat3 = 0.0_jprb

!     ! Matrix has a sparsity pattern
!     !     (C D E)
!     ! A = (F G H)
!     !     (0 0 I)
!     ! mblock = m/3
!     ! m2block = 2*mblock 
!     m = 9
!     mblock = 3
!     m2block = 6
!     ! Do the top-left (C, D, F, G)
!     do jg = 1,iend
!       C(:,:,jg) = 0.0_jprb
!       do j2 = 1,m2block
!           do j3 = 1,m2block
!               C(1:m2block,j2,jg) = C(1:m2block,j2,jg) + A(1:m2block,j3,jg)*B(j3,j2,jg)
!           end do
!       end do

!       do j2 = m2block+1,m
!       ! Do the top-right (E & H)
!           do j3 = 1,m
!             C(1:m2block,j2,jg) = C(1:m2block,j2,jg) + A(1:m2block,j3,jg)*B(j3,j2,jg)
!           end do
!       ! Do the bottom-right (I)
!           do j3 = m2block+1,m
!             C(m2block+1:m,j2,jg) = C(m2block+1:m,j2,jg) + A(m2block+1:m,j3,jg)*B(j3,j2,jg)
!          end do
!       end do

!     end do
!   end subroutine mat_x_mat3

!   subroutine mat_x_mat3_square(iend,h,A,C)
!     integer,    intent(in)                      :: iend,h
!     real(jprb), intent(in),   dimension(9,9,iend) :: A
!     real(jprb), intent(out),  dimension(9,9,iend) :: C
!     integer    :: j1, j2, j3, jg,m
!     integer    :: mblock, m2block
! #ifdef USE_LIBXSMM
!     TYPE(LIBXSMM_MMFUNCTION) :: xmm
!     INTEGER(LIBXSMM_BLASINT_KIND) :: mm = 9, n = 9, k = 9
!     integer :: LIBXSMM_PREFETCH_AUTO
! #endif  
!     ! Array-wise assignment
!     C = 0.0_jprb

!     ! Matrix has a sparsity pattern
!     !     (C D E)
!     ! A = (F G H)
!     !     (0 0 I)
!     m = 9
!     mblock = 3
!     m2block = 6

!    ! !GCC$ unroll 3
! #ifndef USE_LIBXSMM
!     ! Do the top-left (C, D, F, G)
!     do j1 = 1, iend
!       C(:,:,j1) = 0.0_jprb
!       do j2 = 1,m2block
!           do j3 = 1,m2block
!               C(1:m2block,j2,j1) = C(1:m2block,j2,j1) + A(1:m2block,j3,j1)*A(j3,j2,j1)
!           end do
!       end do

!       do j2 = m2block+1,m
!       ! Do the top-right (E & H)
!           do j3 = 1,m
!             C(1:m2block,j2,j1) = C(1:m2block,j2,j1) + A(1:m2block,j3,j1)*A(j3,j2,j1)
!           end do
!       ! Do the bottom-right (I)
!           do j3 = m2block+1,m
!             C(m2block+1:m,j2,j1) = C(m2block+1:m,j2,j1) + A(m2block+1:m,j3,j1)*A(j3,j2,j1)
!          end do
!       end do
!     end do
! #else

! ! generates and dispatches a matrix multiplication kernel
!       call libxsmm_mmdispatch(xmm, mm, n, k,                           &
!   &    alpha=1.0_jprb, beta=1.0_jprb, prefetch=LIBXSMM_PREFETCH_AUTO)
!      ! kernel multiplies and accumulates matrices: C += Ai * Bi
!      DO j1 = 1, iend
!        CALL libxsmm_mmcall(xmm, A(:,:,j1), A(:,:,j1), C(:,:,j1))
!      END DO

! #endif 

!   end subroutine mat_x_mat3_square

!     !---------------------------------------------------------------------
!   ! Return matrix X where AX=B. LU, A, X, B all consist of n m-by-m
!   ! matrices.
!   pure subroutine solve_mat3(iend,m,A,B,X)
!     integer,    intent(in) :: m,iend
!     real(jprb), intent(in) :: A(m,m,iend)
!     real(jprb), intent(in) :: B(m,m,iend)
!     real(jprb), intent(out):: X(m,m,iend)

!     real(jprb) :: LU(m,m)

!     integer :: j,i

!     do i = 1,iend
!       call lu_factorization3(m,A(:,:,i),LU)

!       do j = 1, m
!         call lu_substitution3(m,LU,B(1:m,j,i),X(1:m,j,i))
!   !      call lu_substitution(n,iend,m,LU,B(1:n,1:m,j),X(1:m,j))
!       end do
!     end do

!     ! X(:,7:9,1:6) always zero

!   end subroutine solve_mat3


  ! !---------------------------------------------------------------------
  ! ! Treat LU as an LU-factorization of an original matrix A, and
  ! ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
  ! ! m-element vectors.
  ! pure subroutine lu_substitution3(m,LU,b,x)
  !   ! CHECK: dimensions should be ":"?
  !   integer,    intent(in) :: m
  !   real(jprb), intent(in) :: LU(m,m)
  !   real(jprb), intent(in) :: b(m)
  !   real(jprb), intent(out):: x(m)

  !   integer :: j1, j2

  !   x(1:m) = b(1:m)

  !   ! First solve Ly=b
  !   do j2 = 2, m
  !     do j1 = 1, j2-1
  !       x(j2) = x(j2) - x(j1)*LU(j2,j1)
  !     end do
  !   end do
  !   ! Now solve Ux=y
  !   do j2 = m, 1, -1
  !     do j1 = j2+1, m
  !       x(j2) = x(j2) - x(j1)*LU(j2,j1)
  !     end do
  !     x(j2) = x(j2) / LU(j2,j2)
  !   end do

  ! end subroutine lu_substitution3

!   !---------------------------------------------------------------------
  ! ! Treat A as n m-by-m matrices and return the LU factorization of A
  ! ! compressed into a single matrice (with L below the diagonal and U
  ! ! on and above the diagonal; the diagonal elements of L are 1). No
  ! ! pivoting is performed.
  ! pure subroutine lu_factorization3(m, A, LU)
  !   integer,    intent(in)  :: m
  !   real(jprb), intent(in)  :: A(m,m)
  !   real(jprb), intent(out) :: LU(m,m)

  !   real(jprb) :: s
  !   integer    :: j1, j2, j3

  !   ! This routine is adapted from an in-place one, so we first copy
  !   ! the input into the output.
  !   LU(1:m,1:m) = A(1:m,1:m)
  !   !     (C   D  E)
  !   ! A = (-D -C  H)
  !   !     (0   0  I)
  !   do j2 = 1, m
  !     do j1 = 1, j2-1
  !       s = LU(j1,j2)
  !       do j3 = 1, j1-1
  !         s = s - LU(j1,j3) * LU(j3,j2)
  !       end do
  !       LU(j1,j2) = s
  !     end do
  !     do j1 = j2, m
  !       s = LU(j1,j2)
  !       do j3 = 1, j2-1
  !         s = s - LU(j1,j3) * LU(j3,j2)
  !       end do
  !       LU(j1,j2) = s
  !     end do
  !     if (j2 /= m) then
  !       s = 1.0_jprb / LU(j2,j2)
  !       do j1 = j2+1, m
  !         LU(j1,j2) = LU(j1,j2) * s
  !       end do
  !     end if
  !   end do

  ! end subroutine lu_factorization3

   ! subroutine expm_taylor(iend,m,A)

  !   use yomhook, only : lhook, dr_hook

  !   integer,    intent(in)      :: m, iend
  !   real(jprb), intent(inout)   :: A(m,m,iend)
  !   ! local variables
  !   real(jprb)    :: normA, sum_column(m)
  !   integer       :: m_vals(6)
  !   real(jprb)    :: theta(6)
  !   integer       :: j1, j2, j3
  !   integer       :: s(iend)

  !   call expmchk(m_vals, theta)

  !   do j3 = 1,iend ! g-point loop
  !     !  Compute the 1-norms of A
  !     sum_column(:) = 0.0_jprb
  !     do j2 = 1,m
  !         sum_column(:) = sum_column(:) + abs(A(:,j2,j3))
  !     end do
  !     normA = 0.0_jprb
  !     do j1 = 1,m
  !         normA = max(normA,sum_column(j1))
  !     end do

  !     if (normA <= theta(6)) then
  !       ! no scaling and squaring required
  !       s(j3) = 0

  !       do j1 = 1:6
  !         if (normA <= theta(j1)) then
  !             call taylor_approximant_of_degree(A(:,:,j3),m_vals(j1))
  !           exit
  !         end if
  !       end do

  !     else 

  !     end if

  !   end do


  ! end subroutine expm_taylor

  ! subroutine taylor_approximant_of_degree(order, A, E)
  !   ! !  Improved Paterson Stockmeyer scheme, with Pade at end point
  !           I = eye(size(A,1),classA);
  !           if (order>=2) then
  !               A2 = matmul(A,A)
  !           end if
  !           if (order>8) then
  !               A3 = matmul(A,A2)
  !           end if
  !           select case(order)
  !               case 1
  !                   E = I + A;
  !               case 2  
  !                   E = I + A + A2/2;
  !               case 4  
  !                   E = I + A + A2*(I/2 + A/6 + A2/24);
  !               case 8
  !                   !  Minimizes ||coefficients||_1
  !                   x3=2/3;    
  !                   a1=1/88*(1+sqrt(177))*x3;
  !                   a2=1/352*(1+sqrt(177))*x3;
  !                   u0=1; u1=1;
  !                   u2=1/630*(857-58*sqrt(177));
  !                   c0=(-271+29*sqrt(177))/(315*x3);
  !                   c1=(11*(-1+sqrt(177)))/(1260*x3);
  !                   c2=(11*(-9+sqrt(177)))/(5040*x3);
  !                   c4=-((-89+sqrt(177))/(5040*x3^2));
  !                   !  Matrix products
  !                   A4 = A2*(a1*A + a2*A2);
  !                   A8 = (x3*A2 + 1*A4)*(c0*I + c1*A + c2*A2 + c4*A4) ;
  !                   E = u0*I + u1*A + u2*A2 + + A8;
  !               case 12
  !                   !  lower values of ||_||_1
  !                   a01 = -0.0186023205146205532243437300433;
  !                   a02 = 4.60000000000000000000000000000;
  !                   a03 = 0.211693118299809442949323323336; a04 = 0;
  !                   a11 = -0.00500702322573317730979741843919;
  !                   a12 = 0.992875103538486836140479571505;
  !                   a13 = 0.158224384715726725371768893252; 
  !                   a14 = -0.131810610138301840156819349464; 
  !                   a21 = -0.573420122960522263905952420789; 
  !                   a22 = -0.132445561052799638845074997454; 
  !                   a23 = 0.165635169436727415011171668419;
  !                   a24 = -0.0202785554058925907933568229945; 
  !                   a31 = -0.133399693943892059700768926983;
  !                   a32 = 0.00172990000000000000000000000000; 
  !                   a33 = 0.0107862779315792425026320640108; 
  !                   a34 = -0.00675951846863086359778560766482;
            
  !                   q31 = a01*I+a11*A+a21*A2+a31*A3;
  !                   q32 = a02*I+a12*A+a22*A2+a32*A3;
  !                   q33 = a03*I+a13*A+a23*A2+a33*A3;
  !                   q34 = a04*I+a14*A+a24*A2+a34*A3;
  !                   !  Matrix products
  !                   q61 = q33 + q34*q34;
  !                   E = (q31 + (q32 + q61)*q61);
  !               case 18
  !                   !  Minimizes ||coefficients||_1
  !                   a01 = 0; 
  !                   a11 = -0.100365581030144620014614939405;
  !                   a21 =-0.00802924648241156960116919515244;
  !                   a31 =-0.000892138498045729955685466128049;
  !                   b01 = 0;
  !                   b11 =0.3978497494996450761451961277102845756965081084076856223845951607640145373149032030404660339703426170;
  !                   b21 =1.367837784604117199225237068782228242106453540654915795267462332707000476284638745738812082761860458;
  !                   b31 =0.4982896225253826775568588172622733561849181397319696269923450773179744578906675427707618377504305561;
  !                   b61 =-0.0006378981945947233092415500564919285518773827581013332055148653241353120789646323186965398317523194760;
  !                   b02 =-10.96763960529620625935174626753684124863041876254774214673058385106461743913502064396554589372626845;
  !                   b12 =1.680158138789061971827854017248105095278579593010566858091585875627364747895724070033586802947436157;
  !                   b22 =0.05717798464788655127028717132252481317274975093707600394506105236273081373356534970690710643085727120;
  !                   b32 =-0.006982101224880520842904664838015531664389838152077753258353447125605205335288966277257918925881337834;
  !                   b62 =0.00003349750170860705383133673406684398020225996444991565389728295589367037178816169580298011831485225359;
  !                   b03 =-0.09043168323908105619714688770183495499470866281162649735086602288456671216199491949073419844120202066;
  !                   b13 =-0.06764045190713819075600799003797074854562826013273060070581796689746210078264908587143648477465431597;
  !                   b23 =0.06759613017704596460827991950789370565071000768138331455747136830953224978586305376886338214283464385;
  !                   b33 =0.02955525704293155274260691822422312437473046472743510159951445245685893720550532525991666174704105350;
  !                   b63 =-0.00001391802575160607011247399793437424028074102305234870853220069004278819595344896368546425681168813708;
  !                   b04 = 0;
  !                   b14 = 0;
  !                   b24 =-0.0923364619367118592764570775143;
  !                   b34 =-0.0169364939002081717191385115723;
  !                   b64 =-0.0000140086798182036159794363205726;
    
  !                   !  Matrix products
  !                   A6 = A3*A3;
  !                   q31 = a01*I + a11*A + a21*A2 + a31*A3;
  !                   q61 = (b01 *I + b11*A + b21*A2 + b31*A3 + b61*A6); 
  !                   q62 = (b02 *I + b12*A + b22*A2 + b32*A3 + b62*A6);
  !                   q63 = (b03 *I + b13*A + b23*A2 + b33*A3 + b63*A6);
  !                   q64 = (b04 *I+  b14*A + b24*A2 + b34*A3 + b64*A6);
  !                   q91 = q31*q64 + q63;
  !                   q18 = q61 + (q62 + q91)*q91;
  !                   E = q18;
  !           end select
  !       end

  ! subroutine expmchk_singleprec(m_vals, theta)
  !   !EXPMCHK Check the class of input A and
  !   !    initialize M_VALS and THETA accordingly.
  !   integer,    dimension(6), intent(out) :: m_vals
  !   real(jprd), dimension(6), intent(out) :: theta

  !   m_vals = (/ 1,2,4,8,12,18/)
  !   theta  = ( / 1.192092800768788e-7_jprb, &    ! m_vals = 1
  !               5.978858893805233e-04_jprb, &   ! m_vals = 2 
  !               5.116619363445086e-02_jprb, &   ! m_vals = 4
  !               5.800524627688768e-01_jprb, &   ! m_vals = 8
  !               1.461661507209034e+00_jprb, &   ! m_vals = 12
  !               3.010066362817634e+00_jprb /)   ! m_vals = 18
  ! end subroutine expmchk_singleprec

  ! subroutine expmchk_doubleprec(m_vals, theta)
  !   !EXPMCHK Check the class of input A and
  !   !    initialize M_VALS and THETA accordingly.
  !   integer,    dimension(6), intent(out) :: m_vals
  !   real(jprd), dimension(6), intent(out) :: theta

  !   m_vals = (/ 1,2,4,8,12,18/)
  !   theta  = ( /  2.220446049250313e-16_jprb, &    ! m_vals = 1
  !               2.580956802971767e-08_jprb, &   ! m_vals = 2 
  !               3.397168839976962e-04_jprb, &   ! m_vals = 4
  !               4.991228871115323e-02_jprb, &   ! m_vals = 8
  !               2.996158913811580e-01_jprb, &   ! m_vals = 12
  !               1.090863719290036e+00_jprb /)   ! m_vals = 18
  ! end subroutine expmchk_doubleprec

end module radiation_matrix
