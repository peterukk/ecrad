! radiation_ifs_rrtmgp.F90 - Interface to IFS implementation of RRTMGP-NN
!
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Peter U kkonen
! Email:   peter.ukkonen@nbi.ku.dk
!
! Modifications
!   2020-01-30  P. Ukkonen

module radiation_ifs_rrtmgp

    implicit none
  
    ! public  :: setup_gas_optics_ifs_rrtmgp, gas_optics_ifs_rrtmgp
  
  contains
  
    !---------------------------------------------------------------------
    ! Setup the IFS implementation of RRTMGP gas absorption model
    subroutine setup_gas_optics_ifs_rrtmgp(config, directory)
  
      use radiation_config
      use yomhook,   only : lhook, dr_hook
      use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
      use mo_gas_concentrations, only: ty_gas_concs

      type(config_type), intent(inout), target :: config
      character(len=*), intent(in)     :: directory
  
      integer :: irep ! For implied do

      real(jprb) :: hook_handle

      type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array

  
      if (lhook) call dr_hook('radiation_ifs_rrtmgp:setup_gas_optics',0,hook_handle)
  
      ! The IFS implementation of RRTMG uses many global variables.  In
      ! the IFS these will have been set up already; otherwise set them
      ! up now.


    end subroutine setup_gas_optics_ifs_rrtmgp
  
  
    !---------------------------------------------------------------------
    ! Scale gas mixing ratios according to required units
    subroutine set_gas_units(gas)
  
      use radiation_gas,           only : gas_type, IMassMixingRatio
      type(gas_type),    intent(inout) :: gas
  
      call gas%set_units(IMassMixingRatio)
  
    end subroutine set_gas_units
  
  
    !---------------------------------------------------------------------
    ! Compute gas optical depths, shortwave scattering, Planck function
    ! and incoming shortwave radiation at top-of-atmosphere
    subroutine gas_optics_ifs_rrtmgp(ncol,nlev, &
         &  config, single_level, thermodynamics, gas, & 
         &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
         &  incoming_sw)
  
      use parkind1,                 only : jprb, jpim
      use yomhook,   only : lhook, dr_hook

      use radiation_config,         only : config_type, ISolverSpartacus
      use radiation_thermodynamics, only : thermodynamics_type
      use radiation_single_level,   only : single_level_type
      use radiation_gas
      use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp

  
      integer, intent(in) :: ncol               ! number of columns
      integer, intent(in) :: nlev               ! number of model levels
      type(config_type), intent(in) :: config
      type(single_level_type),  intent(in) :: single_level
      type(thermodynamics_type),intent(in) :: thermodynamics
      type(gas_type),           intent(in) :: gas
  
      ! Longwave albedo of the surface
      real(jprb), dimension(config%n_g_lw,ncol), &
           &  intent(in), optional :: lw_albedo
  
      ! Gaseous layer optical depth in longwave and shortwave, and
      ! shortwave single scattering albedo (i.e. fraction of extinction
      ! due to Rayleigh scattering) at each g-point
      real(jprb), dimension(config%n_g_lw,nlev,ncol), intent(out) :: &
           &   od_lw
      real(jprb), dimension(config%n_g_sw,nlev,ncol), intent(out) :: &
           &   od_sw, ssa_sw
  
      ! The Planck function (emitted flux from a black body) at half
      ! levels at each longwave g-point
      real(jprb), dimension(config%n_g_lw,nlev+1,ncol), &
           &   intent(out), optional :: planck_hl
      ! Planck function for the surface (W m-2)
      real(jprb), dimension(config%n_g_lw,ncol), &
           &   intent(out), optional :: lw_emission
  
      ! The incoming shortwave flux into a plane perpendicular to the
      ! incoming radiation at top-of-atmosphere in each of the shortwave
      ! g-points
      real(jprb), dimension(config%n_g_sw,ncol), &
           &   intent(out), optional :: incoming_sw
  
      real(jprb) :: incoming_sw_scale(ncol)
  

      integer :: jlev, jgreorder, jg, ig, iband, jcol
  
      real(jprb) :: hook_handle
  
      
      if (lhook) call dr_hook('radiation_ifs_rrtm:gas_optics',1,hook_handle)
      
    end subroutine gas_optics_ifs_rrtmgp
    
  
    !---------------------------------------------------------------------
    ! Compute Planck function of the atmosphere
    ! subroutine planck_function_atmos(nlev,istartcol,iendcol, &
    !      config, thermodynamics, PFRAC, &
    !      planck_hl)
  
    !   use parkind1,                 only : jprb, jpim
  
    !   use radiation_config,         only : config_type, ISolverSpartacus
    !   use radiation_thermodynamics, only : thermodynamics_type
    !   !use radiation_gas
  
    !   integer, intent(in) :: nlev               ! number of model levels
    !   integer, intent(in) :: istartcol, iendcol ! range of columns to process
    !   type(config_type), intent(in) :: config
    !   type(thermodynamics_type),intent(in) :: thermodynamics
    !   real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW,nlev)
  
    !   ! The Planck function (emitted flux from a black body) at half
    !   ! levels at each longwave g-point
    !   real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
    !        &   planck_hl
  
    !   ! Planck function values per band
    !   real(jprb), dimension(istartcol:iendcol, config%n_bands_lw) :: planck_store
  
    !   ! Look-up table variables for Planck function
    !   real(jprb), dimension(istartcol:iendcol) :: frac
    !   integer,    dimension(istartcol:iendcol) :: ind
  
    !   ! Temperature (K) of a half-level
    !   real(jprb) :: temperature
  
    !   real(jprb) :: factor
    !   real(jprb) :: ZFLUXFAC
  
    !   integer :: jlev, jgreorder, jg, ig, iband, jband, jcol, ilevoffset
  
    !   real(jprb) :: hook_handle
  
    !   if (lhook) call dr_hook('radiation_ifs_rrtmgp:planck_function_atmos',0,hook_handle)
  
    !   ZFLUXFAC = 2.0_jprb*ASIN(1.0_jprb) * 1.0e4_jprb
      
    !   ! nlev may be less than the number of original levels, in which
    !   ! case we assume that the user wants the lower part of the
    !   ! atmosphere
    !   ilevoffset = ubound(thermodynamics%temperature_hl,2)-nlev-1
  
    !   ! Work out interpolations: for each half level, the index of the
    !   ! lowest interpolation bound, and the fraction into interpolation
    !   ! interval
    !   do jlev = 1,nlev+1
    !     do jcol = istartcol,iendcol
    !       temperature = thermodynamics%temperature_hl(jcol,jlev+ilevoffset)
    !       if (temperature < 339.0_jprb .and. temperature >= 160.0_jprb) then
    !         ! Linear interpolation between -113 and 66 degC
    !         ind(jcol)  = int(temperature - 159.0_jprb)
    !         frac(jcol) = temperature - int(temperature)
    !       else if(temperature >= 339.0_jprb) then
    !         ! Extrapolation above 66 degC
    !         ind(jcol)  = 180
    !         frac(jcol) = temperature - 339.0_jprb
    !       else
    !         ! Cap below -113 degC (to avoid possible negative Planck
    !         ! function values)
    !         ind(jcol)  = 1
    !         frac(jcol) = 0.0_jprb
    !       end if
    !     end do
  
    !     ! Calculate Planck functions per band
    !     do jband = 1,config%n_bands_lw
    !       factor = zfluxfac * delwave(jband)
    !       do jcol = istartcol,iendcol
    !         planck_store(jcol,jband) = factor &
    !              &  * (totplnk(ind(jcol),jband) &
    !              &  + frac(jcol)*(totplnk(ind(jcol)+1,jband)-totplnk(ind(jcol),jband)))
    !       end do
    !     end do
  
    !     if (config%i_solver_lw == ISolverSpartacus) then
    !       ! We need to rearrange the gas optics info in memory:
    !       ! reordering the g points in order of approximately increasing
    !       ! optical depth (for efficient 3D processing on only the
    !       ! regions of the spectrum that are optically thin for gases)
    !       ! and reorder in pressure since the the functions above treat
    !       ! pressure decreasing with increasing index.
    !       if (jlev == 1) then
    !         ! Top-of-atmosphere half level - note that PFRAC is on model
    !         ! levels not half levels
    !         do jgreorder = 1,config%n_g_lw
    !           iband = config%i_band_from_reordered_g_lw(jgreorder)
    !           ig = config%i_g_from_reordered_g_lw(jgreorder)
    !           planck_hl(jgreorder,1,:) = planck_store(:,iband) &
    !                &   * PFRAC(:,ig,nlev)
    !         end do
    !       else
    !         do jgreorder = 1,config%n_g_lw
    !           iband = config%i_band_from_reordered_g_lw(jgreorder)
    !           ig = config%i_g_from_reordered_g_lw(jgreorder)
    !           planck_hl(jgreorder,jlev,:) &
    !                  &   = planck_store(:,iband) &
    !                  &   * PFRAC(:,ig,nlev+2-jlev)
    !         end do
    !       end if
    !     else
    !       ! G points have not been reordered 
    !       if (jlev == 1) then
    !         ! Top-of-atmosphere half level - note that PFRAC is on model
    !         ! levels not half levels
    !         do jg = 1,config%n_g_lw
    !           iband = config%i_band_from_g_lw(jg)
    !           planck_hl(jg,1,:) = planck_store(:,iband) * PFRAC(:,jg,nlev)
    !         end do
    !       else
    !         do jg = 1,config%n_g_lw
    !           iband = config%i_band_from_g_lw(jg)
    !           planck_hl(jg,jlev,:) = planck_store(:,iband) * PFRAC(:,jg,nlev+2-jlev)
    !         end do
    !       end if
    !     end if
    !   end do
  
    !   if (lhook) call dr_hook('radiation_ifs_rrtmgp:planck_function_atmos',1,hook_handle)
  
    ! end subroutine planck_function_atmos
  
  
    !---------------------------------------------------------------------
    ! Compute Planck function of the surface
    ! subroutine planck_function_surf(istartcol, iendcol, config, temperature, PFRAC, &
    !      &  planck_surf)
  
    !   use parkind1,                 only : jprb, jpim
  
    !   USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
    !   use yoerrtwn, only : totplnk, delwave
  
    !   use yomhook, only : lhook, dr_hook
  
    !   use radiation_config,         only : config_type, ISolverSpartacus
    !   !    use radiation_gas
  
    !   integer, intent(in) :: istartcol, iendcol ! range of columns to process
    !   type(config_type), intent(in) :: config
    !   real(jprb), intent(in) :: temperature(:)
  
    !   real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW)
  
    !   ! Planck function of the surface (W m-2)
    !   real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
    !        &  intent(out) :: planck_surf
  
    !   ! Planck function values per band
    !   real(jprb), dimension(istartcol:iendcol, config%n_bands_lw) :: planck_store
  
    !   ! Look-up table variables for Planck function
    !   real(jprb), dimension(istartcol:iendcol) :: frac
    !   integer,    dimension(istartcol:iendcol) :: ind
  
    !   ! Temperature (K)
    !   real(jprb) :: Tsurf
  
    !   real(jprb) :: factor
    !   real(jprb) :: ZFLUXFAC
  
    !   integer :: jgreorder, jg, ig, iband, jband, jcol
  
    !   real(jprb) :: hook_handle
  
    !   if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_surf',0,hook_handle)
  
    !   ZFLUXFAC = 2.0_jprb*ASIN(1.0_jprb) * 1.0e4_jprb
  
    !   ! Work out surface interpolations
    !   do jcol = istartcol,iendcol
    !     Tsurf = temperature(jcol)
    !     if (Tsurf < 339.0_jprb .and. Tsurf >= 160.0_jprb) then
    !       ! Linear interpolation between -113 and 66 degC
    !       ind(jcol)  = int(Tsurf - 159.0_jprb)
    !       frac(jcol) = Tsurf - int(Tsurf)
    !     else if(Tsurf >= 339.0_jprb) then
    !       ! Extrapolation above 66 degC
    !       ind(jcol)  = 180
    !       frac(jcol) = Tsurf - 339.0_jprb
    !     else
    !       ! Cap below -113 degC (to avoid possible negative Planck
    !       ! function values)
    !       ind(jcol)  = 1
    !       frac(jcol) = 0.0_jprb
    !     end if
    !   end do
  
    !   ! Calculate Planck functions per band
    !   do jband = 1,config%n_bands_lw
    !     factor = zfluxfac * delwave(jband)
    !     do jcol = istartcol,iendcol
    !       planck_store(jcol,jband) = factor &
    !            &  * (totplnk(ind(jcol),jband) &
    !            &  + frac(jcol)*(totplnk(ind(jcol)+1,jband)-totplnk(ind(jcol),jband)))
    !     end do
    !   end do
  
    !   if (config%i_solver_lw == ISolverSpartacus) then
    !     ! We need to rearrange the gas optics info in memory: reordering
    !     ! the g points in order of approximately increasing optical
    !     ! depth (for efficient 3D processing on only the regions of
    !     ! the spectrum that are optically thin for gases) and reorder
    !     ! in pressure since the the functions above treat pressure
    !     ! decreasing with increasing index.
    !     do jgreorder = 1,config%n_g_lw
    !       iband = config%i_band_from_reordered_g_lw(jgreorder)
    !       ig = config%i_g_from_reordered_g_lw(jgreorder)
    !       planck_surf(jgreorder,:) = planck_store(:,iband) * PFRAC(:,ig)
    !     end do
    !   else
    !     ! G points have not been reordered 
    !     do jg = 1,config%n_g_lw
    !       iband = config%i_band_from_g_lw(jg)
    !       planck_surf(jg,:) = planck_store(:,iband) * PFRAC(:,jg)
    !     end do
    !   end if
  
    !   if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_surf',1,hook_handle)
      
    ! end subroutine planck_function_surf
  
  
    ! !---------------------------------------------------------------------
    ! ! Externally facing function for computing the Planck function
    ! ! without reference to any gas profile; typically this would be used
    ! ! for computing the emission by facets of a complex surface.  Note
    ! ! that this uses fixed "PFRAC" values, obtained by averaging over
    ! ! those derived from RRTM-G for near-surface conditions over a line
    ! ! of meridian from the ECMWF model.
    ! subroutine planck_function(config, temperature, planck_surf)
  
    !   use parkind1,                 only : jprb, jpim
  
    !   use radiation_config,         only : config_type
  
    !   type(config_type), intent(in) :: config
    !   real(jprb), intent(in) :: temperature
  
    !   ! Planck function of the surface (W m-2)
    !   real(jprb), dimension(config%n_g_lw), &
    !        &  intent(out) :: planck_surf
  
    !   ! Fraction of each band contributed by each g-point within
    !   ! it. Since there are 16 bands, this array sums to 16
    !   real(jprb), parameter, dimension(1,140) :: frac &
    !        = reshape( (/ 0.21227E+00, 0.18897E+00, 0.25491E+00, 0.17864E+00, 0.11735E+00, 0.38298E-01, 0.57871E-02, &
    !        &    0.31753E-02, 0.53169E-03, 0.76476E-04, 0.16388E+00, 0.15241E+00, 0.14290E+00, 0.12864E+00, &
    !        &    0.11615E+00, 0.10047E+00, 0.80013E-01, 0.60445E-01, 0.44918E-01, 0.63395E-02, 0.32942E-02, &
    !        &    0.54541E-03, 0.15380E+00, 0.15194E+00, 0.14339E+00, 0.13138E+00, 0.11701E+00, 0.10081E+00, &
    !        &    0.82296E-01, 0.61735E-01, 0.41918E-01, 0.45918E-02, 0.37743E-02, 0.30121E-02, 0.22500E-02, &
    !        &    0.14490E-02, 0.55410E-03, 0.78364E-04, 0.15938E+00, 0.15146E+00, 0.14213E+00, 0.13079E+00, &
    !        &    0.11672E+00, 0.10053E+00, 0.81566E-01, 0.61126E-01, 0.41150E-01, 0.44488E-02, 0.36950E-02, &
    !        &    0.29101E-02, 0.21357E-02, 0.19609E-02, 0.14134E+00, 0.14390E+00, 0.13913E+00, 0.13246E+00, &
    !        &    0.12185E+00, 0.10596E+00, 0.87518E-01, 0.66164E-01, 0.44862E-01, 0.49402E-02, 0.40857E-02, &
    !        &    0.32288E-02, 0.23613E-02, 0.15406E-02, 0.58258E-03, 0.82171E-04, 0.29127E+00, 0.28252E+00, &
    !        &    0.22590E+00, 0.14314E+00, 0.45494E-01, 0.71792E-02, 0.38483E-02, 0.65712E-03, 0.29810E+00, &
    !        &    0.27559E+00, 0.11997E+00, 0.10351E+00, 0.84515E-01, 0.62253E-01, 0.41050E-01, 0.44217E-02, &
    !        &    0.36946E-02, 0.29113E-02, 0.34290E-02, 0.55993E-03, 0.31441E+00, 0.27586E+00, 0.21297E+00, &
    !        &    0.14064E+00, 0.45588E-01, 0.65665E-02, 0.34232E-02, 0.53199E-03, 0.19811E+00, 0.16833E+00, &
    !        &    0.13536E+00, 0.11549E+00, 0.10649E+00, 0.93264E-01, 0.75720E-01, 0.56405E-01, 0.41865E-01, &
    !        &    0.59331E-02, 0.26510E-02, 0.40040E-03, 0.32328E+00, 0.26636E+00, 0.21397E+00, 0.14038E+00, &
    !        &    0.52142E-01, 0.38852E-02, 0.14601E+00, 0.13824E+00, 0.27703E+00, 0.22388E+00, 0.15446E+00, &
    !        &    0.48687E-01, 0.98054E-02, 0.18870E-02, 0.11961E+00, 0.12106E+00, 0.13215E+00, 0.13516E+00, &
    !        &    0.25249E+00, 0.16542E+00, 0.68157E-01, 0.59725E-02, 0.49258E+00, 0.33651E+00, 0.16182E+00, &
    !        &    0.90984E-02, 0.95202E+00, 0.47978E-01, 0.91716E+00, 0.82857E-01, 0.77464E+00, 0.22536E+00 /), (/ 1,140 /) )
  
    !   call planck_function_surf(1, 1, config, spread(temperature,1,1), &
    !        &                    frac, planck_surf)
  
    ! end subroutine planck_function
  
  end module radiation_ifs_rrtmgp
  