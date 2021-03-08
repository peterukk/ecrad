! ecrad_driver.F90 - Driver for offline ECRAD radiation scheme
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
! ECRAD is the radiation scheme used in the ECMWF Integrated
! Forecasting System in cycle 43R3 and later. Several solvers are
! available, including McICA, Tripleclouds and SPARTACUS (the Speedy
! Algorithm for Radiative Transfer through Cloud Sides, a modification
! of the two-stream formulation of shortwave and longwave radiative
! transfer to account for 3D radiative effects). Gas optical
! properties are provided by the RRTM-G gas optics scheme.

! This program takes three arguments:
! 1) Namelist file to configure the radiation calculation
! 2) Name of a NetCDF file containing one or more atmospheric profiles
! 3) Name of output NetCDF file

program ecrad_driver

    ! --------------------------------------------------------
    ! Section 1: Declarations
    ! --------------------------------------------------------
    use parkind1,                 only : jprb, jprd ! Working/double precision
  
    use radiation_io,             only : nulout
    use radiation_interface,      only : setup_radiation, radiation, set_gas_units
    use radiation_config,         only : config_type, IGasModelRRTMGP
    use radiation_single_level,   only : single_level_type
    use radsurf_properties,       only : surface_type, print_surface_representation
    use radsurf_intermediate,     only : surface_intermediate_type
    use radsurf_flux,             only : surface_flux_type
    use radsurf_save,             only : save_surface_fluxes
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type, &
         &   IVolumeMixingRatio, IMassMixingRatio, &
         &   IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, &
         &   IHCFC22, ICCl4, GasName, GasLowerCaseName, NMaxGases
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_flux,           only : flux_type
    use radiation_save,           only : save_fluxes, save_inputs
    use ecrad_driver_config,      only : driver_config_type
  #ifdef BLOCK_DERIVED_TYPES
    use ecrad_driver_read_input,  only : read_input, read_input_blocked, unblock_fluxes
  #else
    use ecrad_driver_read_input,  only : read_input
  #endif
    use easy_netcdf
    use mo_gas_concentrations,    only : ty_gas_concs
    use mod_network
  #ifdef USE_TIMING
    ! Timing library
    use gptl,                  only: gptlstart, gptlstop, gptlinitialize, gptlpr, gptlfinalize, gptlsetoption, &
                                     gptlpercent, gptloverhead
  #endif
    implicit none
  #ifdef USE_PAPI  
  #include "f90papi.h"
  #endif  
  
    ! The NetCDF file containing the input profiles
    type(netcdf_file)         :: file
  
    ! Derived types for the inputs to the radiation scheme
    type(config_type)         :: config
  #ifdef BLOCK_DERIVED_TYPES
    type(single_level_type),    dimension(:), allocatable :: single_level_b
    type(thermodynamics_type),  dimension(:), allocatable :: thermodynamics_b
    type(gas_type),             dimension(:), allocatable :: gas_b
    type(cloud_type),           dimension(:), allocatable :: cloud_b
    type(aerosol_type),         dimension(:), allocatable :: aerosol_b
    ! Derived type for the surface inputs
    type(surface_type),         dimension(:), allocatable :: surface_b
  #else
    type(single_level_type)   :: single_level
    ! type(thermodynamics_type) :: thermodynamics
    type(gas_type)            :: gas
    type(cloud_type)          :: cloud
    type(aerosol_type)        :: aerosol
    ! Derived type for the surface inputs
    type(surface_type)        :: surface
  #endif
    type(thermodynamics_type) :: thermodynamics
  
    ! Derived types for the surface fluxes
    type(surface_intermediate_type) :: surface_intermediate
    type(surface_flux_type)         :: surface_flux
  
    ! Configuration specific to this driver
    type(driver_config_type)  :: driver_config
  
    ! Derived type containing outputs from the radiation scheme
    type(flux_type)                             :: flux
  #ifdef BLOCK_DERIVED_TYPES
    type(flux_type), dimension(:), allocatable  :: flux_b
  #endif
  
    integer :: ncol, nlev         ! Number of columns and levels
    integer :: istartcol, iendcol ! Range of columns to process
  
    ! Name of file names specified on command line
    character(len=512) :: file_name
    integer            :: istatus ! Result of command_argument_count
  
    ! For parallel processing of multiple blocks
    integer :: jblock, nblock, blocksize ! Block loop index and number
    integer, external :: omp_get_thread_num
    double precision, external :: omp_get_wtime
  
    ! Loop index for repeats (for benchmarking)
    integer :: jrepeat
  
    ! Are any variables out of bounds?
    logical :: is_out_of_bounds
  
    ! Are we using a complex surface representation stored in the
    ! "surface" structure?
    logical :: is_complex_surface
  
  !  integer    :: iband(20), nweights
  !  real(jprb) :: weight(20)
  
    ! Start/stop time in seconds
    real(kind=jprd) :: tstart, tstop
  
    character(len=10) :: file_id
  
    character(len=512) :: file_name_inp, file_name_ref ! optional, compare fluxes to some reference
    logical :: compare_to_ref = .false.
  
  #ifdef USE_TIMING
    integer :: ret
    ! print *, "using GPTL timing library"
    !
    ! Initialize timers
    !
    ret = gptlsetoption (gptlpercent, 1)        ! Turn on "% of" print
    ret = gptlsetoption (gptloverhead, 0)       ! Turn off overhead estimate
  
  #ifdef USE_PAPI  
  #ifdef SINGLE_PRECISION
    ret = GPTLsetoption (PAPI_SP_OPS, 1);
  #else
    ret = GPTLsetoption (PAPI_DP_OPS, 1);
  #endif
  #endif  
    ret = gptlinitialize()
    ret = gptlstart('ecrad_total')
  #endif
    ! --------------------------------------------------------
    ! Section 2: Configure
    ! --------------------------------------------------------
  
    ! Check program called with correct number of arguments
    if (command_argument_count() < 3) then
      stop 'Usage: ecrad config.nam input_file.nc output_file.nc [output_surface_file.nc]'
    end if
  
    ! Use namelist to configure the radiation calculation
    call get_command_argument(1, file_name, status=istatus)
    if (istatus /= 0) then
      stop 'Failed to read name of namelist file as string of length < 512'
    end if
  
    ! Read "radiation" namelist into radiation configuration type
    call config%read(file_name=file_name)
  
    ! Read "radiation_driver" namelist into radiation driver config type
    call driver_config%read(file_name)
  
    if (driver_config%iverbose >= 2) then
      write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
      write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
      write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
  #ifdef SINGLE_PRECISION
      write(nulout,'(a)') 'Floating-point precision: single'
  #else
      write(nulout,'(a)') 'Floating-point precision: double'
  #endif
      call config%print(driver_config%iverbose)
    end if
  
    ! Albedo/emissivity intervals may be specified like this
    !call config%define_sw_albedo_intervals(6, &
    !     &  [0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
    !     &     1.19_jprb, 2.38e-6_jprb], [1,2,3,4,5,6], &
    !     &   do_nearest=.false.)
    !call config%define_lw_emiss_intervals(3, &
    !     &  [8.0e-6_jprb, 13.0e-6_jprb], [1,2,1], &
    !     &   do_nearest=.false.)
  
    ! --------------------------------------------------------
    ! Section 3: Read input data file
    ! --------------------------------------------------------
  
    ! Get NetCDF input file name
    call get_command_argument(2, file_name, status=istatus)
    if (istatus /= 0) then
      stop 'Failed to read name of input NetCDF file as string of length < 512'
    end if
  
    ! Open the file and configure the way it is read
    call file%open(trim(file_name), iverbose=driver_config%iverbose)
  
    ! Get NetCDF output file name
    call get_command_argument(3, file_name, status=istatus)
    if (istatus /= 0) then
      stop 'Failed to read name of output NetCDF file as string of length < 512'
    end if
  
    ! 2D arrays are assumed to be stored in the file with height varying
    ! more rapidly than column index. Specifying "true" here transposes
    ! all 2D arrays so that the column index varies fastest within the
    ! program.
    call file%transpose_matrices(.true.)
  
    ! Read input variables from NetCDF file
  #ifdef USE_TIMING
    ret =  gptlstart('read_inputs')
  #endif
  
  #ifdef BLOCK_DERIVED_TYPES
    call read_input_blocked(file, config, driver_config, nblock, blocksize, ncol, nlev, &
         &          is_complex_surface, surface_b, single_level_b, thermodynamics_b, &
         &          gas_b, cloud_b, aerosol_b)
  #else
    call read_input(file, config, driver_config, ncol, nlev, &
          &          is_complex_surface, surface, single_level, thermodynamics, &
          &          gas, cloud, aerosol)
  #endif
  
  #ifdef USE_TIMING
    ret =  gptlstop('read_inputs')
    ret =  gptlstart('setup_radiation')
  #endif
  
    ! Setup the radiation scheme: load the coefficients for gas and
    ! cloud optics, currently from RRTMG
    if (config%i_gas_model == IGasModelRRTMGP) then
     ! call setup_radiation(config, gas_b(1)%gas_rrtmgp, single_level_b(1)%solar_irradiance)
    else
      call setup_radiation(config)
    end if
  
  #ifdef USE_TIMING
    ret =  gptlstop('setup_radiation')
  #endif
  
    if (is_complex_surface) then
      config%do_canopy_fluxes_sw = .true.
      config%do_canopy_fluxes_lw = .true.
    end if
  
  
    ! Close input file
    call file%close()
  
    ! if (is_complex_surface .and. driver_config%iverbose >= 2) then
    !   call print_surface_representation(surface%i_representation)
    ! end if
  
    ! Compute seed from skin temperature residual
  ! #ifdef BLOCK_DERIVED_TYPES
  !   do jblock = 1, nblock
  !     single_level_b(jblock)%iseed = int(1.0e9*(single_level_b(jblock)%skin_temperature &
  !           &                            -int(single_level_b(jblock)%skin_temperature)))
  !   end do
  ! #else
  !    single_level%iseed = int(1.0e9*(single_level%skin_temperature &
  !         &                            -int(single_level%skin_temperature)))
  ! #endif
  
    ! Set first and last columns to process
    if (driver_config%iendcol < 1 .or. driver_config%iendcol > ncol) then
      driver_config%iendcol = ncol
    end if
  
    if (driver_config%istartcol > driver_config%iendcol) then
      write(nulout,'(a,i0,a,i0,a,i0,a)') '*** Error: requested column range (', &
           &  driver_config%istartcol, &
           &  ' to ', driver_config%iendcol, ') is out of the range in the data (1 to ', &
           &  ncol, ')'
      stop 1
    end if
  
    ! --------------------------------------------------------
    ! Section 4: Call radiation scheme
    ! --------------------------------------------------------
  
  #ifdef BLOCK_DERIVED_TYPES
  
    allocate(flux_b(nblock))
  
    do jblock = 1, nblock
      if (driver_config%do_save_inputs) then
        write(file_id, '(i0)') jblock
        file_name_inp =  "inputs_" // trim(adjustl(file_id)) // ".nc"
      
        call save_inputs(file_name_inp, config, single_level_b(jblock), thermodynamics_b(jblock), &
            &                gas_b(jblock), cloud_b(jblock), aerosol_b(jblock), &
            &                lat=spread(0.0_jprb,1,blocksize), &
            &                lon=spread(0.0_jprb,1,blocksize), &
            &                iverbose=driver_config%iverbose)
      end if
    
        ! Ensure the units of the gas mixing ratios are what is required
      ! by the gas absorption model
      call set_gas_units(config, gas_b(jblock))
  
      ! Compute saturation with respect to liquid (needed for aerosol
      ! hydration) call
      call thermodynamics_b(jblock)%calc_saturation_wrt_liquid(1, blocksize)
  
      ! if (is_complex_surface) then
      !   call surface_intermediate%allocate(driver_config%istartcol, driver_config%iendcol, &
      !       &                             config, surface)
      !   call surface_flux%allocate(config, driver_config%istartcol, driver_config%iendcol, &
      !       &                     surface%i_representation)
      ! end if
  
      ! Check inputs are within physical bounds, printing message if not
      is_out_of_bounds =     gas_b(jblock)%out_of_physical_bounds() &
          & .or.   single_level_b(jblock)%out_of_physical_bounds() &
          & .or. thermodynamics_b(jblock)%out_of_physical_bounds() &
          & .or.          cloud_b(jblock)%out_of_physical_bounds() &
          & .or.        aerosol_b(jblock)%out_of_physical_bounds() 
  
      ! Allocate memory for the flux profiles, which may include arrays
      ! of dimension n_bands_sw/n_bands_lw, so must be called after
      ! setup_radiation
      call flux_b(jblock)%allocate(config, 1, blocksize, nlev)
  
    end do  
  
    call thermodynamics%allocate(ncol, nlev, allocated(thermodynamics_b(1)%h2o_sat_liq))
  
  #else
    
    ! Store inputs
    if (driver_config%do_save_inputs) then
      call save_inputs('inputs.nc', config, single_level, thermodynamics, &
           &                gas, cloud, aerosol, &
           &                lat=spread(0.0_jprb,1,ncol), &
           &                lon=spread(0.0_jprb,1,ncol), &
           &                iverbose=driver_config%iverbose)
    end if
  
  
    call set_gas_units(config, gas)
  
    ! Compute saturation with respect to liquid (needed for aerosol
    ! hydration) call
    call thermodynamics%calc_saturation_wrt_liquid(driver_config%istartcol,driver_config%iendcol)
  
    if (is_complex_surface) then
      call surface_intermediate%allocate(driver_config%istartcol, driver_config%iendcol, &
           &                             config, surface)
      call surface_flux%allocate(config, driver_config%istartcol, driver_config%iendcol, &
           &                     surface%i_representation)
    end if
  
    ! Check inputs are within physical bounds, printing message if not
    is_out_of_bounds =     gas%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
         &                                            driver_config%do_correct_unphysical_inputs) &
         & .or.   single_level%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
         &                                            driver_config%do_correct_unphysical_inputs) &
         & .or. thermodynamics%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
         &                                            driver_config%do_correct_unphysical_inputs) &
         & .or.          cloud%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
         &                                            driver_config%do_correct_unphysical_inputs) &
         & .or.        aerosol%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
         &                                            driver_config%do_correct_unphysical_inputs) 
    
  #endif
  
    ! Allocate memory for the flux profiles, which may include arrays
    ! of dimension n_bands_sw/n_bands_lw, so must be called after
    ! setup_radiation
    call flux%allocate(config, 1, ncol, nlev)
  
    if (driver_config%iverbose >= 2) then
      write(nulout,'(a)')  'Performing radiative transfer calculations'
    end if
  
    ! Option of repeating calculation multiple time for more accurate
    ! profiling
  #ifdef USE_TIMING
      ret =  gptlstart('radiation')
  #endif   
    tstart = omp_get_wtime() 
    do jrepeat = 1,driver_config%nrepeat
      
      if (driver_config%do_parallel) then
  
  #ifdef BLOCK_DERIVED_TYPES
        ! Run radiation scheme over blocks of columns in parallel
  
        ! Compute number of blocks to process
        ! nblock already calculated, in this case all columns are processed
        ! and nblock is simply ncol_tot / blocksize
       
        ! !$OMP PARALLEL DO shared(k_dist_lw, k_dist_sw) firstprivate(optical_props_lw, optical_props_sw, source)
        !$OMP PARALLEL DO 
        do jblock = 1, nblock
          ! Specify the range of columns to process.
          istartcol = 1
          iendcol = blocksize
          
          if (driver_config%iverbose >= 3) then
            write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
                 &  ' processing columns ', istartcol, '-', iendcol
          end if
          
          ! if (is_complex_surface) then
          !   call surface_intermediate%calc_boundary_conditions(driver_config%istartcol, &
          !        &  driver_config%iendcol, config, surface, thermodynamics, gas, single_level)
          ! end if
    
          ! Call the ECRAD radiation scheme
          ! if (config%i_gas_model == IGasModelRRTMGP) then
            
  
          call radiation(blocksize, nlev, istartcol, iendcol, config, &
                    &  single_level_b(jblock), thermodynamics_b(jblock), gas_b(jblock), cloud_b(jblock), &
                    aerosol_b(jblock), flux_b(jblock))
  
          call unblock_fluxes(flux, flux_b(jblock), thermodynamics, thermodynamics_b(jblock), jblock, blocksize)
  
          ! if (is_complex_surface) then
          !   call surface_intermediate%partition_fluxes(driver_config%istartcol, &
          !        &  driver_config%iendcol, config, surface, flux, surface_flux)
          ! end if
  
          ! write(file_id, '(i0)') jblock
          ! file_name =  "outputs_" // trim(adjustl(file_id)) // ".nc"
          ! call save_fluxes(file_name, config, thermodynamics(jblock), flux(jblock), &
          ! &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
          ! &   experiment_name=driver_config%experiment_name)
        end do
        !$OMP END PARALLEL DO
  
  #else
        ! No blocking
        
        ! Compute number of blocks to process
        nblock = (driver_config%iendcol - driver_config%istartcol &
             &  + driver_config%nblocksize) / driver_config%nblocksize
       
        !$OMP PARALLEL DO PRIVATE(istartcol, iendcol, ret) SCHEDULE(RUNTIME)
        do jblock = 1, nblock
          ! Specify the range of columns to process.
          istartcol = (jblock-1) * driver_config%nblocksize &
               &    + driver_config%istartcol
          iendcol = min(istartcol + driver_config%nblocksize - 1, &
               &        driver_config%iendcol)
            
          if (driver_config%iverbose >= 3) then
            write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
                 &  ' processing columns ', istartcol, '-', iendcol
          end if
          
          if (is_complex_surface) then
            call surface_intermediate%calc_boundary_conditions(driver_config%istartcol, &
                 &  driver_config%iendcol, config, surface, thermodynamics, gas, single_level)
          end if
          
          ! Call the ECRAD radiation scheme
          call radiation(ncol, nlev, istartcol, iendcol, config, &
               &  single_level, thermodynamics, gas, cloud, aerosol, flux)
          
          if (is_complex_surface) then
            call surface_intermediate%partition_fluxes(driver_config%istartcol, &
                 &  driver_config%iendcol, config, surface, flux, surface_flux)
          end if
          
        end do
        !$OMP END PARALLEL DO
  #endif
      else
        ! Run radiation scheme serially
  
  #ifdef BLOCK_DERIVED_TYPES
        write(nulout,'(a)')  'do_parallel must be turned on for BLOCK_DERIVED_TYPES=1'
        stop
  #else
  
        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a)')  'Processing ', ncol, ' columns'
        end if
        
        if (is_complex_surface) then
          call surface_intermediate%calc_boundary_conditions(driver_config%istartcol, &
               &  driver_config%iendcol, config, surface, thermodynamics, gas, single_level)
        end if
        
        ! Call the ECRAD radiation scheme
        call radiation(ncol, nlev, driver_config%istartcol, driver_config%iendcol, &
             &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)
        
        if (is_complex_surface) then
          call surface_intermediate%partition_fluxes(driver_config%istartcol, &
               &  driver_config%iendcol, config, surface, flux, surface_flux)
        end if
  #endif
      end if
      
    end do
  #ifdef USE_TIMING
    ret =  gptlstop('radiation')
  #endif
    tstop = omp_get_wtime()
    write(nulout, '(a,g11.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'
  
    ! do jrepeat = 1,config%n_g_sw
    !   print *, "SW band id for ig ", jrepeat, ":", config%i_band_from_g_sw(jrepeat)
    ! end do
    ! do jrepeat = 1,config%n_g_lw
    !   print *, "LW band id for ig ", jrepeat, ":", config%i_band_from_g_lw(jrepeat), config%i_band_from_reordered_g_lw(jrepeat)
    ! end do
    ! --------------------------------------------------------
    ! Section 5: Check and save output
    ! --------------------------------------------------------
  
    ! is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)
  
    ! Store the fluxes in the output file
  #ifdef USE_TIMING
    ret =  gptlstart('write_outputs')
  #endif
  
    ! compare_to_ref = .true.
    ! file_name_ref = file_name(1:index(file_name,'.nc')-1) // "_REFERENCE.nc"
  
    ! Store the fluxes in the output file
    if (compare_to_ref) then
      call save_fluxes(file_name, config, thermodynamics, flux, &
          &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
          &   experiment_name=driver_config%experiment_name, file_name_ref=file_name_ref)
    else
      call save_fluxes(file_name, config, thermodynamics, flux, &
      &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
      &   experiment_name=driver_config%experiment_name)
    end if
  
  #ifdef USE_TIMING
    ret =  gptlstop('write_outputs')
  #endif
      
    if (is_complex_surface) then
      ! Get NetCDF output file name for surface
      call get_command_argument(4, file_name, status=istatus)
      if (istatus /= 0) then
        write(nulout,'(a)') 'Warning: file name for surface-flux outputs not provided'
      else
        call save_surface_fluxes(file_name, config, surface_flux, iverbose=driver_config%iverbose)
      end if
    end if
  
    if (driver_config%iverbose >= 2) then
      write(nulout,'(a)') '------------------------------------------------------------------------------------'
    end if
    
  #ifdef USE_TIMING
    ! End timers
    !
    ret =  gptlstop('ecrad_total')
    ret = gptlpr(driver_config%nblocksize)
    ret = gptlfinalize()
  #endif
  
  end program ecrad_driver
  