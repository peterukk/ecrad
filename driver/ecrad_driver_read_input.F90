! ecrad_driver_read_input.F90 - Read input structures from NetCDF file
!
! (C) Copyright 2018- ECMWF.
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

module ecrad_driver_read_input

  public
  
contains

  subroutine read_input(file, config, driver_config, ncol, nlev, &
       &          is_complex_surface, surface, single_level, thermodynamics, &
       &          gas, cloud, aerosol)

    use parkind1,                 only : jprb, jpim
    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, ISolverSPARTACUS
    use ecrad_driver_config,      only : driver_config_type
    use radsurf_properties,       only : surface_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type, &
       &   IVolumeMixingRatio, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4, INO2, GasName, GasLowerCaseName, NMaxGases
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use easy_netcdf,              only : netcdf_file
    
    implicit none

    type(netcdf_file),         intent(in)    :: file
    type(config_type),         intent(in)    :: config
    type(driver_config_type),  intent(in)    :: driver_config
    type(surface_type),        intent(inout) :: surface
    type(single_level_type),   intent(inout) :: single_level
    type(thermodynamics_type), intent(inout) :: thermodynamics
    type(gas_type),            intent(inout) :: gas
    type(cloud_type),  target, intent(inout) :: cloud
    type(aerosol_type),        intent(inout) :: aerosol

    ! Number of columns and levels of input data
    integer, intent(out) :: ncol, nlev

    ! Are we using a complex surface representation stored in the
    ! "surface" structure?
    logical, intent(out) :: is_complex_surface

    integer :: ngases             ! Num of gases with concs described in 2D
    integer :: nwellmixedgases    ! Num of globally well-mixed gases

    ! Mixing ratio of gases described in 2D (ncol,nlev); this is
    ! volume mixing ratio (m3/m3) except for water vapour and ozone
    ! for which it is mass mixing ratio (kg/kg)
    real(jprb), allocatable, dimension(:,:) :: gas_mr

    ! Volume mixing ratio (m3/m3) of globally well-mixed gases
    real(jprb)                              :: well_mixed_gas_vmr

    ! Name of gas concentration variable in the file
    character(40)               :: gas_var_name

    ! Cloud overlap decorrelation length (m)
    real(jprb), parameter :: decorr_length_default = 2000.0_jprb

    ! General surface prop to be read and then modified before used in
    ! an ecRad structure
    real(jprb), allocatable, dimension(:)   :: prop_1d
    real(jprb), allocatable, dimension(:,:) :: prop_2d

    integer :: jgas               ! Loop index for reading gases
    integer :: irank              ! Dimensions of gas data

    ! Can we scale cloud size using namelist parameters?  No if the
    ! cloud size came from namelist parameters in the first place, yes
    ! if it came from the NetCDF file in the first place
    logical :: is_cloud_size_scalable

    ! The following calls read in the data, allocating memory for 1D and
    ! 2D arrays.  The program will stop if any variables are not found.
    
    ! Pressure and temperature (SI units) are on half-levels, i.e. of
    ! length (ncol,nlev+1)
    call file%get('pressure_hl',   thermodynamics%pressure_hl)
    call file%get('temperature_hl',thermodynamics%temperature_hl)

    ! Extract array dimensions
    ncol = size(thermodynamics%pressure_hl,1)
    nlev = size(thermodynamics%pressure_hl,2)-1

    if (driver_config%solar_irradiance_override > 0.0_jprb) then
      ! Optional override of solar irradiance
      single_level%solar_irradiance = driver_config%solar_irradiance_override
      if (driver_config%iverbose >= 2) then
        write(nulout,'(a,f10.1)')  '  Overriding solar irradiance with ', &
             &  driver_config%solar_irradiance_override
      end if
    else if (file%exists('solar_irradiance')) then
      call file%get('solar_irradiance', single_level%solar_irradiance)
    else
      single_level%solar_irradiance = 1366.0_jprb
      if (driver_config%iverbose >= 1 .and. config%do_sw) then
        write(nulout,'(a,g10.3,a)') 'Warning: solar irradiance set to ', &
             &  single_level%solar_irradiance, ' W m-2'
        end if
    end if

    if (driver_config%cos_sza_override >= 0.0_jprb) then
      ! Optional override of cosine of solar zenith angle
      allocate(single_level%cos_sza(ncol))
      single_level%cos_sza = driver_config%cos_sza_override
      if (driver_config%iverbose >= 2) then
        write(nulout,'(a,g10.3)') '  Overriding cosine of the solar zenith angle with ', &
             &  driver_config%cos_sza_override
      end if
    else if (file%exists('cos_solar_zenith_angle')) then
      ! Single-level variables, all with dimensions (ncol)
      call file%get('cos_solar_zenith_angle',single_level%cos_sza)
    else if (.not. config%do_sw) then
      ! If cos_solar_zenith_angle not present and shortwave radiation
      ! not to be performed, we create an array of zeros as some gas
      ! optics schemes still need to be run in the shortwave
      allocate(single_level%cos_sza(ncol))
      single_level%cos_sza = 0.0_jprb
    else
      write(nulout,'(a,a)') '*** Error: cos_solar_zenith_angle not provided'
      stop
    end if

    if (config%do_clouds) then

      ! --------------------------------------------------------
      ! Read cloud properties needed by most solvers
      ! --------------------------------------------------------

      ! Read cloud descriptors with dimensions (ncol, nlev)
      call file%get('cloud_fraction',cloud%fraction)

      ! Fractional standard deviation of in-cloud water content
      if (file%exists('fractional_std')) then
        call file%get('fractional_std', cloud%fractional_std)
      end if
      
      ! Cloud water content and effective radius may be provided
      ! generically, in which case they have dimensions (ncol, nlev,
      ! ntype)
      if (file%exists('q_hydrometeor')) then
        call file%get('q_hydrometeor',  cloud%mixing_ratio, ipermute=[2,1,3])     ! kg/kg
        call file%get('re_hydrometeor', cloud%effective_radius, ipermute=[2,1,3]) ! m
      else
        ! Ice and liquid properties provided in separate arrays
        allocate(cloud%mixing_ratio(ncol,nlev,2))
        allocate(cloud%effective_radius(ncol,nlev,2))
        call file%get('q_liquid', prop_2d)   ! kg/kg
        cloud%mixing_ratio(:,:,1) = prop_2d
        call file%get('q_ice', prop_2d)   ! kg/kg
        cloud%mixing_ratio(:,:,2) = prop_2d
        call file%get('re_liquid', prop_2d)   ! m
        cloud%effective_radius(:,:,1) = prop_2d
        call file%get('re_ice', prop_2d)   ! m
        cloud%effective_radius(:,:,2) = prop_2d
      end if
      ! For backwards compatibility, associate pointers for liquid and
      ! ice to the first and second slices of cloud%mixing_ratio and
      ! cloud%effective_radius
      cloud%q_liq  => cloud%mixing_ratio(:,:,1)
      cloud%q_ice  => cloud%mixing_ratio(:,:,2)
      cloud%re_liq => cloud%effective_radius(:,:,1)
      cloud%re_ice => cloud%effective_radius(:,:,2)

      ! Simple initialization of the seeds for the Monte Carlo scheme
      call single_level%init_seed_simple(1,ncol)
      ! Overwrite with user-specified values if available
      if (file%exists('iseed')) then
        call file%get('iseed', single_level%iseed)
      end if

      ! Cloud overlap parameter
      if (file%exists('overlap_param')) then
        call file%get('overlap_param', cloud%overlap_param)
      end if

      ! Optional scaling of liquid water mixing ratio
      if (driver_config%q_liq_scaling >= 0.0_jprb &
           &  .and. driver_config%q_liq_scaling /= 1.0_jprb) then
        cloud%q_liq = cloud%q_liq * driver_config%q_liq_scaling
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling liquid water mixing ratio by a factor of ', &
               &  driver_config%q_liq_scaling
        end if
      end if

      ! Optional scaling of ice water mixing ratio
      if (driver_config%q_ice_scaling >= 0.0_jprb .and. driver_config%q_ice_scaling /= 1.0_jprb) then
        cloud%q_ice = cloud%q_ice * driver_config%q_ice_scaling
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling ice water mixing ratio by a factor of ', &
               &  driver_config%q_ice_scaling
        end if
      end if

      ! Optional scaling of cloud fraction
      if (driver_config%cloud_fraction_scaling >= 0.0_jprb &
           &  .and. driver_config%cloud_fraction_scaling /= 1.0_jprb) then
        cloud%fraction = cloud%fraction * driver_config%cloud_fraction_scaling
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling cloud_fraction by a factor of ', &
               &  driver_config%cloud_fraction_scaling
        end if
      end if

      ! Cloud overlap is currently treated by an overlap decorrelation
      ! length (m) that is constant everywhere, and specified in one
      ! of the namelists
      if (driver_config%overlap_decorr_length_override > 0.0_jprb) then
        ! Convert overlap decorrelation length to overlap parameter between
        ! adjacent layers, stored in cloud%overlap_param
        call cloud%set_overlap_param(thermodynamics, &
             &    driver_config%overlap_decorr_length_override)
      else if (.not. allocated(cloud%overlap_param)) then 
        if (driver_config%iverbose >= 1) then
          write(nulout,'(a,g10.3,a)') 'Warning: overlap decorrelation length set to ', &
               &  decorr_length_default, ' m'
        end if
        call cloud%set_overlap_param(thermodynamics, decorr_length_default)
      else if (driver_config%overlap_decorr_length_scaling > 0.0_jprb) then
        ! Scale the overlap decorrelation length by taking the overlap
        ! parameter to a power
        !    where (cloud%overlap_param > 0.99_jprb) cloud%overlap_param = 0.99_jprb
        
        where (cloud%overlap_param > 0.0_jprb) 
          cloud%overlap_param = cloud%overlap_param**(1.0_jprb &
               &                             / driver_config%overlap_decorr_length_scaling)
        end where
        
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling overlap decorrelation length by a factor of ', &
               &  driver_config%overlap_decorr_length_scaling
        end if
      else if (driver_config%overlap_decorr_length_scaling == 0.0_jprb) then
        cloud%overlap_param = 0.0_jprb
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a)')  '  Setting overlap decorrelation length to zero (random overlap)'
        end if
      end if
      
      ! Cloud inhomogeneity is specified by the fractional standard
      ! deviation of cloud water content, that is currently constant
      ! everywhere (and the same for water and ice). The following copies
      ! this constant into the cloud%fractional_std array.
      if (driver_config%fractional_std_override >= 0.0_jprb) then
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3,a)') '  Overriding cloud fractional standard deviation with ', &
               &  driver_config%fractional_std_override
        end if
        call cloud%create_fractional_std(ncol, nlev, &
             &  driver_config%fractional_std_override)
      else if (.not. allocated(cloud%fractional_std)) then
        call cloud%create_fractional_std(ncol, nlev, 0.0_jprb)
        if (driver_config%iverbose >= 1) then
          write(nulout,'(a)') 'Warning: cloud optical depth fractional standard deviation set to zero'
        end if
      end if

      ! --------------------------------------------------------
      ! Read cloud properties needed by SPARTACUS
      ! --------------------------------------------------------

      if (config%i_solver_sw == ISolverSPARTACUS &
           &  .or.   config%i_solver_lw == ISolverSPARTACUS) then

        ! 3D radiative effects are governed by the length of cloud
        ! edge per area of gridbox, which is characterized by the
        ! inverse of the cloud effective size (m-1). Order of
        ! precedence: (1) effective size namelist overrides, (2)
        ! separation namelist overrides, (3) inv_cloud_effective_size
        ! present in NetCDF, (4) inv_cloud_effective_separation
        ! present in NetCDF. Only in the latter two cases may the
        ! effective size be scaled by the namelist variable
        ! "effective_size_scaling".

        is_cloud_size_scalable = .false. ! Default for cases (1) and (2)

        if (driver_config%low_inv_effective_size_override >= 0.0_jprb &
             &  .or. driver_config%middle_inv_effective_size_override >= 0.0_jprb &
             &  .or. driver_config%high_inv_effective_size_override >= 0.0_jprb) then
          ! (1) Cloud effective size specified in namelist

          ! First check all three ranges provided
          if (driver_config%low_inv_effective_size_override < 0.0_jprb &
             &  .or. driver_config%middle_inv_effective_size_override < 0.0_jprb &
             &  .or. driver_config%high_inv_effective_size_override < 0.0_jprb) then
            write(nulout,'(a,a)') '*** Error: if one of [low|middle|high]_inv_effective_size_override', &
                 & ' is provided then all must be'
            stop
          end if
          if (driver_config%iverbose >= 2) then
            write(nulout,'(a,g10.3,a)') '  Overriding inverse cloud effective size with:'
            write(nulout,'(a,g10.3,a)') '    ', driver_config%low_inv_effective_size_override, &
                 &       ' m-1 (low clouds)'
            write(nulout,'(a,g10.3,a)') '    ', driver_config%middle_inv_effective_size_override, &
                 &       ' m-1 (mid-level clouds)'
            write(nulout,'(a,g10.3,a)') '    ', driver_config%high_inv_effective_size_override, &
                 &       ' m-1 (high clouds)'
          end if
          call cloud%create_inv_cloud_effective_size_eta(ncol, nlev, &
               &  thermodynamics%pressure_hl, &
               &  driver_config%low_inv_effective_size_override, &
               &  driver_config%middle_inv_effective_size_override, &
               &  driver_config%high_inv_effective_size_override, 0.8_jprb, 0.45_jprb)

        else if (driver_config%cloud_separation_scale_surface > 0.0_jprb &
             &  .and. driver_config%cloud_separation_scale_toa > 0.0_jprb) then
          ! (2) Cloud separation scale provided in namelist

          if (driver_config%iverbose >= 2) then
            write(nulout,'(a)') '  Effective cloud separation parameterized versus eta:'
            write(nulout,'(a,f8.1,a)') '    ', &
                 &  driver_config%cloud_separation_scale_surface, ' m at the surface'
            write(nulout,'(a,f8.1,a)') '    ', &
                 &  driver_config%cloud_separation_scale_toa, ' m at top-of-atmosphere'
            write(nulout,'(a,f6.2)') '     Eta power is', &
                 &  driver_config%cloud_separation_scale_power
            write(nulout,'(a,f6.2)') '     Inhomogeneity separation scaling is', &
                 &  driver_config%cloud_inhom_separation_factor
          end if
          call cloud%param_cloud_effective_separation_eta(ncol, nlev, &
               &  thermodynamics%pressure_hl, &
               &  driver_config%cloud_separation_scale_surface, &
               &  driver_config%cloud_separation_scale_toa, &
               &  driver_config%cloud_separation_scale_power, &
               &  driver_config%cloud_inhom_separation_factor)
          
        else if (file%exists('inv_cloud_effective_size')) then
          ! (3) NetCDF file contains cloud effective size

          is_cloud_size_scalable = .true.

          call file%get('inv_cloud_effective_size', cloud%inv_cloud_effective_size)
          ! For finer control we can specify the effective size for
          ! in-cloud inhomogeneities as well
          if (file%exists('inv_inhom_effective_size')) then
            if (.not. driver_config%do_ignore_inhom_effective_size) then
              call file%get('inv_inhom_effective_size', cloud%inv_inhom_effective_size)
            else
              if (driver_config%iverbose >= 1) then
                write(nulout,'(a)') 'Ignoring inv_inhom_effective_size so treated as equal to inv_cloud_effective_size'
                write(nulout,'(a)') 'Warning: ...this is unlikely to be accurate for cloud fraction near one'
              end if
            end if
          else
            if (driver_config%iverbose >= 1) then
              write(nulout,'(a)') 'Warning: inv_inhom_effective_size not set so treated as equal to inv_cloud_effective_size'
              write(nulout,'(a)') 'Warning: ...this is unlikely to be accurate for cloud fraction near one'
            end if
          end if
          
        else if (file%exists('inv_cloud_effective_separation')) then
          ! (4) Alternative way to specify cloud scale

          is_cloud_size_scalable = .true.
          
          call file%get('inv_cloud_effective_separation', prop_2d)
          allocate(cloud%inv_cloud_effective_size(ncol,nlev))
          allocate(cloud%inv_inhom_effective_size(ncol,nlev))
          where (cloud%fraction > config%cloud_fraction_threshold &
               &  .and. cloud%fraction < 1.0_jprb - config%cloud_fraction_threshold)
            ! Convert effective cloud separation to effective cloud
            ! size, noting divisions rather than multiplications
            ! because we're working in terms of inverse sizes
            cloud%inv_cloud_effective_size = prop_2d / sqrt(cloud%fraction*(1.0_jprb-cloud%fraction))
          elsewhere
            cloud%inv_cloud_effective_size = 0.0_jprb
          end where
          if (file%exists('inv_inhom_effective_separation')) then
            if (driver_config%iverbose >= 2) then
              write(nulout,'(a)') '  Effective size of clouds and their inhomogeneities being computed from input'
              write(nulout,'(a)') '  ...variables inv_cloud_effective_separation and inv_inhom_effective_separation'
            end if
            call file%get('inv_inhom_effective_separation', prop_2d)
            where (cloud%fraction > config%cloud_fraction_threshold)
              ! Convert effective separation of cloud inhomogeneities
              ! to effective size of cloud inhomogeneities, assuming
              ! here that the Tripleclouds treatment of cloud
              ! inhomogeneity will divide the cloudy part of the area
              ! into regions of equal area
              cloud%inv_inhom_effective_size = prop_2d &
                   &  / sqrt(0.5_jprb*cloud%fraction * (1.0_jprb-0.5_jprb*cloud%fraction))
            elsewhere
              cloud%inv_inhom_effective_size = 0.0_jprb
            end where
          else
            ! Assume that the effective separation of cloud
            ! inhomogeneities is equal to that of clouds but
            ! multiplied by a constant provided by the user; note that
            ! prop_2d at this point contains
            ! inv_cloud_effective_separation
            if (driver_config%iverbose >= 2) then
              write(nulout,'(a)') '  Effective size of clouds being computed from inv_cloud_effective_separation'
              write(nulout,'(a,f6.2,a)') '  ...and multiplied by ', driver_config%cloud_inhom_separation_factor, &
                   &  ' to get effective size of inhomogeneities'
            end if
            where (cloud%fraction > config%cloud_fraction_threshold)
              ! Note divisions rather than multiplications because
              ! we're working in terms of inverse sizes
              cloud%inv_inhom_effective_size = (1.0_jprb / driver_config%cloud_inhom_separation_factor) * prop_2d &
                   &  / sqrt(0.5_jprb*cloud%fraction * (1.0_jprb-0.5_jprb*cloud%fraction))
            elsewhere
              cloud%inv_inhom_effective_size = 0.0_jprb
            end where
          end if ! exists inv_inhom_effective_separation
          deallocate(prop_2d)
          
        else

          write(nulout,'(a)') '*** Error: SPARTACUS solver specified but cloud size not, either in namelist or input file'
          stop

        end if ! Select method of specifying cloud effective size
        
        ! In cases (3) and (4) above the effective size obtained from
        ! the NetCDF may be scaled by a namelist variable
        if (is_cloud_size_scalable .and. driver_config%effective_size_scaling > 0.0_jprb) then
          ! Scale cloud effective size
          cloud%inv_cloud_effective_size = cloud%inv_cloud_effective_size &
               &                         / driver_config%effective_size_scaling
          if (allocated(cloud%inv_inhom_effective_size)) then
            if (driver_config%iverbose >= 2) then
              write(nulout, '(a,g10.3)') '  Scaling effective size of clouds and their inhomogeneities with ', &
                   &                           driver_config%effective_size_scaling
            end if
            cloud%inv_inhom_effective_size = cloud%inv_inhom_effective_size &
                 &                         / driver_config%effective_size_scaling
          else
            if (driver_config%iverbose >= 2) then
              write(nulout, '(a,g10.3)') '  Scaling cloud effective size with ', &
                   &                           driver_config%effective_size_scaling
            end if
          end if
        end if

      end if ! Using SPARTACUS solver

    end if ! do_cloud

    ! --------------------------------------------------------
    ! Read surface properties
    ! --------------------------------------------------------

    ! Surface properties
    if (file%exists('tile_representation')) then
      ! We have a complex representation
      single_level%is_simple_surface = .false.
      is_complex_surface = .true.
      call surface%read(file)
      if (config%use_canopy_full_spectrum_sw) then
        allocate(single_level%sw_albedo(ncol,config%n_g_sw))
        allocate(single_level%sw_albedo_direct(ncol,config%n_g_sw))
      else
        allocate(single_level%sw_albedo(ncol,surface%nalbedobands))
        allocate(single_level%sw_albedo_direct(ncol,surface%nalbedobands))
      end if

      ! Optional override of shortwave albedo
      if (driver_config%sw_albedo_override >= 0.0_jprb) then
        surface%sw_albedo = driver_config%sw_albedo_override
        if (allocated(surface%sw_albedo_direct)) then
          surface%sw_albedo_direct = driver_config%sw_albedo_override
        end if
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding shortwave albedo of each facet with ', &
               &  driver_config%sw_albedo_override
        end if
      end if

      if (config%use_canopy_full_spectrum_lw) then
        allocate(single_level%lw_emission(ncol,config%n_g_lw))
        allocate(single_level%lw_emissivity(ncol,config%n_g_lw))
      else
        allocate(single_level%lw_emission(ncol,surface%nemissbands))
        allocate(single_level%lw_emissivity(ncol,surface%nemissbands))
      end if

      ! Optional override of longwave emissivity
      if (driver_config%lw_emissivity_override >= 0.0_jprb) then
        surface%lw_emissivity = driver_config%lw_emissivity_override
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Overriding longwave emissivity of each facet with ', &
               &  driver_config%lw_emissivity_override
        end if
      end if

    else
      ! We have a "simple" representation with a single flat tile, so
      ! the "surface" structure is not used
      single_level%is_simple_surface = .true.
      is_complex_surface = .false.

      ! Single-level variable with dimensions (ncol)
      if (file%exists('skin_temperature')) then
        call file%get('skin_temperature',single_level%skin_temperature) ! K
      else
        allocate(single_level%skin_temperature(ncol))
        single_level%skin_temperature(1:ncol) = thermodynamics%temperature_hl(1:ncol,nlev+1)
        if (driver_config%iverbose >= 1 .and. config%do_lw &
             &  .and. driver_config%skin_temperature_override < 0.0_jprb) then 
          write(nulout,'(a)') 'Warning: skin temperature set equal to lowest air temperature'
        end if
      end if

      if (driver_config%sw_albedo_override >= 0.0_jprb) then
        ! Optional override of shortwave albedo
        allocate(single_level%sw_albedo(ncol,1))
        single_level%sw_albedo = driver_config%sw_albedo_override
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding shortwave albedo with ', &
               &  driver_config%sw_albedo_override
        end if
        !if (allocated(single_level%sw_albedo_direct)) then
        !  single_level%sw_albedo_direct = driver_config%sw_albedo_override
        !end if
      else
        ! Shortwave albedo is stored with dimensions (ncol,nalbedobands)
        if (file%get_rank('sw_albedo') == 1) then
          ! ...but if in the NetCDF file it has only dimension (ncol), in
          ! order that nalbedobands is correctly set to 1, we need to turn
          ! off transposition
          call file%get('sw_albedo',    single_level%sw_albedo, do_transp=.false.)
          if (file%exists('sw_albedo_direct')) then
            call file%get('sw_albedo_direct', single_level%sw_albedo_direct, do_transp=.false.)
          end if
        else
          call file%get('sw_albedo',    single_level%sw_albedo, do_transp=.true.)
          if (file%exists('sw_albedo_direct')) then
            call file%get('sw_albedo_direct', single_level%sw_albedo_direct, do_transp=.true.)
          end if
        end if
      end if

      ! Longwave emissivity
      if (driver_config%lw_emissivity_override >= 0.0_jprb) then
        ! Optional override of longwave emissivity
        allocate(single_level%lw_emissivity(ncol,1))
        single_level%lw_emissivity = driver_config%lw_emissivity_override
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Overriding longwave emissivity with ', &
               &  driver_config%lw_emissivity_override
        end if
      else
        if (file%get_rank('lw_emissivity') == 1) then
          call file%get('lw_emissivity',single_level%lw_emissivity, do_transp=.false.)
        else
          call file%get('lw_emissivity',single_level%lw_emissivity, do_transp=.true.)
        end if
      end if
    end if

    ! Optional override of skin temperature
    if (driver_config%skin_temperature_override >= 0.0_jprb) then
      if (is_complex_surface) then
        surface%skin_temperature = driver_config%skin_temperature_override
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding skin_temperature of each facet with ', &
               &  driver_config%skin_temperature_override
        end if
      else
        single_level%skin_temperature = driver_config%skin_temperature_override
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding skin_temperature with ', &
               &  driver_config%skin_temperature_override
        end if
      end if
    end if

    ! --------------------------------------------------------
    ! Read aerosol and gas concentrations
    ! --------------------------------------------------------

    if (config%use_aerosols) then
      ! Load aerosol data
      call file%get('aerosol_mmr', aerosol%mixing_ratio, ipermute=[2,3,1]);
      ! Store aerosol level bounds
      aerosol%istartlev = lbound(aerosol%mixing_ratio, 2)
      aerosol%iendlev   = ubound(aerosol%mixing_ratio, 2)
    end if

    ! Load in gas volume mixing ratios, which can be either 2D arrays
    ! (varying with height and column) or 0D scalars (constant volume
    ! mixing ratio everywhere).
    ngases          = 0 ! Gases with varying mixing ratio
    nwellmixedgases = 0 ! Gases with constant mixing ratio

    ! Water vapour and ozone are always in terms of mass mixing ratio
    ! (kg/kg) and always 2D arrays with dimensions (ncol,nlev), unlike
    ! other gases (see below)

    call gas%allocate(ncol, nlev)

    ! Loop through all radiatively important gases
    do jgas = 1,NMaxGases
      if (jgas == IH2O) then
        if (file%exists('q')) then
          call file%get('q', gas_mr)
          call gas%put(IH2O, IMassMixingRatio, gas_mr)
        else if (file%exists('h2o_mmr')) then
          call file%get('h2o_mmr', gas_mr)
          call gas%put(IH2O, IMassMixingRatio, gas_mr)
        else
          call file%get('h2o' // trim(driver_config%vmr_suffix_str), gas_mr);
          call gas%put(IH2O, IVolumeMixingRatio, gas_mr)
        end if
      else if (jgas == IO3) then
        if (file%exists('o3_mmr')) then
          call file%get('o3_mmr', gas_mr)
          call gas%put(IO3, IMassMixingRatio, gas_mr)
        else
          call file%get('o3' // trim(driver_config%vmr_suffix_str), gas_mr)
          call gas%put(IO3, IVolumeMixingRatio, gas_mr)
        end if
      else
        ! Find number of dimensions of the variable holding gas "jgas" in
        ! the input file, where the following function returns -1 if the
        ! gas is not found
        gas_var_name = trim(GasLowerCaseName(jgas)) // trim(driver_config%vmr_suffix_str)
        irank = file%get_rank(trim(gas_var_name))
        ! Note that if the gas is not present then a warning will have
        ! been issued, and irank will be returned as -1
        if (irank == 0) then
          ! Store this as a well-mixed gas
          call file%get(trim(gas_var_name), well_mixed_gas_vmr)
          call gas%put_well_mixed(jgas, IVolumeMixingRatio, well_mixed_gas_vmr)
        else if (irank == 2) then
          call file%get(trim(gas_var_name), gas_mr)
          call gas%put(jgas, IVolumeMixingRatio, gas_mr)
        else if (irank > 0) then
          write(nulout,'(a,a,a)')  '***  Error: ', trim(gas_var_name), ' does not have 0 or 2 dimensions'
          stop
        end if
      end if
      if (allocated(gas_mr)) deallocate(gas_mr)
    end do

    ! Scale gas concentrations if needed
    call gas%scale(IH2O,    driver_config%h2o_scaling,    driver_config%iverbose >= 2)
    call gas%scale(ICO2,    driver_config%co2_scaling,    driver_config%iverbose >= 2)
    call gas%scale(IO3,     driver_config%o3_scaling,     driver_config%iverbose >= 2)
    call gas%scale(IN2O,    driver_config%n2o_scaling,    driver_config%iverbose >= 2)
    call gas%scale(ICO,     driver_config%co_scaling,     driver_config%iverbose >= 2)
    call gas%scale(ICH4,    driver_config%ch4_scaling,    driver_config%iverbose >= 2)
    call gas%scale(IO2,     driver_config%o2_scaling,     driver_config%iverbose >= 2)
    call gas%scale(ICFC11,  driver_config%cfc11_scaling,  driver_config%iverbose >= 2)
    call gas%scale(ICFC12,  driver_config%cfc12_scaling,  driver_config%iverbose >= 2)
    call gas%scale(IHCFC22, driver_config%hcfc22_scaling, driver_config%iverbose >= 2)
    call gas%scale(ICCL4,   driver_config%ccl4_scaling,   driver_config%iverbose >= 2)
    call gas%scale(INO2,    driver_config%no2_scaling,    driver_config%iverbose >= 2)

  end subroutine read_input


  subroutine read_input_blocked(file, config, driver_config, nblocks, ncol, nlev, &
       &          is_complex_surface, surface, single_level, thermodynamics, &
       &          gas, cloud, aerosol)!, gas_rrtmgp)

    use parkind1,                 only : jprb, jpim
    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, ISolverSPARTACUS, IGasModelRRTMGP
    use ecrad_driver_config,      only : driver_config_type
    use radsurf_properties,       only : surface_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type, &
       &   IVolumeMixingRatio, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4, INO2, GasName, GasLowerCaseName, NMaxGases
    use radiation_gas_constants,  only : GasMolarMass, AirMolarMass
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use easy_netcdf,              only : netcdf_file
    use mo_gas_concentrations,    only : ty_gas_concs
    
    implicit none

    type(netcdf_file),         intent(in)    :: file
    type(config_type),         intent(in)    :: config
    type(driver_config_type),  intent(in)    :: driver_config
    type(surface_type),         dimension(:), allocatable,  intent(inout) :: surface
    type(single_level_type),    dimension(:), allocatable,  intent(inout) :: single_level
    type(thermodynamics_type),  dimension(:), allocatable,  intent(inout) :: thermodynamics
    type(gas_type),             dimension(:), allocatable,  intent(inout) :: gas
    type(cloud_type), target,   dimension(:), allocatable,  intent(inout) :: cloud
    type(aerosol_type),         dimension(:), allocatable,  intent(inout) :: aerosol
    ! type(ty_gas_concs),         dimension(:), allocatable,  optional, &
    !        &                                                intent(inout):: gas_rrtmgp

    ! Number of columns and levels of input data
    integer, intent(out) :: ncol, nlev, nblocks

    ! Are we using a complex surface representation stored in the
    ! "surface" structure?
    logical, intent(out) :: is_complex_surface

    integer :: ngases             ! Num of gases with concs described in 2D
    integer :: nwellmixedgases    ! Num of globally well-mixed gases
    integer :: blocksize    

    ! Mixing ratio of gases described in 2D (ncol,nlev); this is
    ! volume mixing ratio (m3/m3) except for water vapour and ozone
    ! for which it is mass mixing ratio (kg/kg)
    real(jprb), allocatable, dimension(:,:) :: gas_mr

    ! Volume mixing ratio (m3/m3) of globally well-mixed gases
    real(jprb)                              :: well_mixed_gas_vmr

    ! Name of gas concentration variable in the file
    character(40)               :: gas_var_name

    ! Cloud overlap decorrelation length (m)
    real(jprb), parameter :: decorr_length_default = 2000.0_jprb

    ! General surface prop to be read and then modified before used in
    ! an ecRad structure
    real(jprb) :: prop_scalar
    real(jprb), allocatable, dimension(:)   :: prop_1d, prop_1d_2
    real(jprb), allocatable, dimension(:,:) :: prop_2d, prop_2d_2, prop_2d_3,prop_2d_4
    real(jprb), allocatable, dimension(:,:,:)   :: prop_3d, prop_3d_2
    integer, allocatable, dimension(:)   :: prop_1d_int


    integer, allocatable, dimension(:) ::  istartcols, iendcols
    integer :: jgas               ! Loop index for reading gases
    integer :: irank              ! Dimensions of gas data
    integer :: b, istartcol, iendcol

    ! Can we scale cloud size using namelist parameters?  No if the
    ! cloud size came from namelist parameters in the first place, yes
    ! if it came from the NetCDF file in the first place
    logical :: is_cloud_size_scalable
    ! Use RRTMGP gas optics?
    logical :: use_rrtmgp = .false.

     ! Scaling factor for RRTMGP gases
    real(jprb) :: sf

    ! Array of integers for keeping track of which gases are present in input file
    integer, dimension(:), allocatable :: gas_inds
    character*6, dimension(:), allocatable :: rrtmgp_gas_names

    if (config%i_gas_model == IGasModelRRTMGP) use_rrtmgp = .true.

    ! The following calls read in the data, allocating memory for 1D and
    ! 2D arrays.  The program will stop if any variables are not found.
    
    ! Pressure and temperature (SI units) are on half-levels, i.e. of
    ! length (ncol,nlev+1)
    call file%get('pressure_hl',   prop_2d, do_transp=.false.)
    call file%get('temperature_hl',prop_2d_2, do_transp=.false.)

    ! Extract array dimensions
    ncol = size(prop_2d,2)
    nlev = size(prop_2d,1)-1

      ! Compute number of blocks to process
    blocksize = driver_config%nblocksize
    nblocks = ncol / blocksize ! no remainder allowed
    ! nblocks = (driver_config%iendcol - driver_config%istartcol &
    ! &  + driver_config%nblocksize) / driver_config%nblocksize
    
    if ( .not. ((mod(ncol, blocksize) == 0 ) .and. (driver_config%istartcol == 1) &
     & .and. (driver_config%iendcol == 0 .or. driver_config%iendcol == ncol))) then
      stop 'to use full blocking, columns must fit neatly into blocks and all columns must be used'
    end if
    
    ! Allocate all the class structures
    allocate(surface(nblocks))
    allocate(single_level(nblocks))
    allocate(thermodynamics(nblocks))
    allocate(gas(nblocks))
    allocate(cloud(nblocks))
    allocate(aerosol(nblocks))
    !if (present(gas_rrtmgp)) allocate(gas_rrtmgp(nblocks))
    !if (config%i_gas_model == IGasModelRRTMGP) allocate(gas%gas_rrtmgp)
    allocate(istartcols(nblocks),iendcols(nblocks))

    ! Start blocking variables and get the block start and end indices
    do b = 1, nblocks
      istartcol = (b-1) * blocksize + 1
      iendcol = istartcol + blocksize - 1
      thermodynamics(b)%pressure_hl     = transpose(prop_2d(:,istartcol:iendcol))
      thermodynamics(b)%temperature_hl  = transpose(prop_2d_2(:,istartcol:iendcol))
      istartcols(b) = istartcol
      iendcols(b)   = iendcol

      if (config%i_gas_model == IGasModelRRTMGP) then
        thermodynamics(b)%pressure_hl_reverse     = prop_2d(:,istartcol:iendcol)
        thermodynamics(b)%temperature_hl_reverse  = prop_2d_2(:,istartcol:iendcol)

        thermodynamics(b)%pressure_fl_reverse &
        &  = 0.5_jprb * (thermodynamics(b)%pressure_hl_reverse(1:nlev, :) &
        &               +thermodynamics(b)%pressure_hl_reverse(2:nlev+1,:))
        thermodynamics(b)%temperature_fl_reverse &
        &  = 0.5_jprb * (thermodynamics(b)%temperature_hl_reverse(1:nlev,:) &
        &               +thermodynamics(b)%temperature_hl_reverse(2:nlev+1,:))
   
      end if
    end do

    deallocate(prop_2d, prop_2d_2)

    if (driver_config%solar_irradiance_override > 0.0_jprb) then
      ! Optional override of solar irradiance
      prop_scalar = driver_config%solar_irradiance_override
      if (driver_config%iverbose >= 2) then
        write(nulout,'(a,f10.1)')  '  Overriding solar irradiance with ', &
             &  driver_config%solar_irradiance_override
      end if
    else if (file%exists('solar_irradiance')) then
      call file%get('solar_irradiance', prop_scalar)
    else
      prop_scalar = 1366.0_jprb
      if (driver_config%iverbose >= 1 .and. config%do_sw) then
        write(nulout,'(a,g10.3,a)') 'Warning: solar irradiance set to ', &
             &  prop_scalar, ' W m-2'
        end if
    end if
    do b = 1, nblocks
      single_level(b)%solar_irradiance = prop_scalar
      allocate(single_level(b)%cos_sza(blocksize))
    end do 

    if (file%exists('cos_solar_zenith_angle')) then
      ! Single-level variables, all with dimensions (ncol)
      call file%get('cos_solar_zenith_angle',prop_1d)
    end if 
    

    do b = 1, nblocks
      istartcol = (b-1) * blocksize + driver_config%istartcol
      iendcol = istartcol + blocksize - 1

      if (driver_config%cos_sza_override >= 0.0_jprb) then
        ! Optional override of cosine of solar zenith angle
        single_level(b)%cos_sza = driver_config%cos_sza_override
        if (driver_config%iverbose >= 2 .and. b==1) then
          write(nulout,'(a,g10.3)') '  Overriding cosine of the solar zenith angle with ', &
               &  driver_config%cos_sza_override
        end if
      else if (allocated(prop_1d)) then
        ! Single-level variables, all with dimensions (ncol)
        single_level(b)%cos_sza = prop_1d(istartcol:iendcol)
      else if (.not. config%do_sw) then
        ! If cos_solar_zenith_angle not present and shortwave radiation
        ! not to be performed, we create an array of zeros as some gas
        ! optics schemes still need to be run in the shortwave
        single_level(b)%cos_sza = 0.0_jprb
      else
        write(nulout,'(a,a)') '*** Error: cos_solar_zenith_angle not provided'
        stop
      end if
    end do

    if (config%do_clouds) then

      ! --------------------------------------------------------
      ! Read cloud properties needed by most solvers
      ! --------------------------------------------------------

      ! Read cloud descriptors, all with dimensions (ncol, nlev)
      call file%get('cloud_fraction',prop_2d, do_transp=.false.)
      do b = 1, nblocks
        cloud(b)%fraction = transpose(prop_2d(:,istartcols(b):iendcols(b)))
      end do

      ! Fractional standard deviation of in-cloud water content
      if (file%exists('fractional_std')) then
        call file%get('fractional_std', prop_2d, do_transp=.false.)
        do b = 1, nblocks
          cloud(b)%fractional_std = transpose(prop_2d(:,istartcols(b):iendcols(b)))
        end do
      end if

      ! Cloud water content and effective radius may be provided
      ! generically, in which case they have dimensions (ncol, nlev,
      ! ntype)
      if (file%exists('q_hydrometeor')) then                ! originally nlev, ncol, ntype
        call file%get('q_hydrometeor',  prop_3d)     ! kg/kg
        call file%get('re_hydrometeor', prop_3d_2) ! m

        do b = 1, nblocks
          cloud(b)%mixing_ratio  = &
            & reshape(prop_3d(:,istartcols(b):iendcols(b),:),[blocksize,nlev,size(prop_3d,3)],order=[2,1,3])
          cloud(b)%effective_radius  = &
            & reshape(prop_3d_2(:,istartcols(b):iendcols(b),:),[blocksize,nlev,size(prop_3d,3)],order=[2,1,3])
        end do
        deallocate(prop_3d, prop_3d_2)
      else
        ! Ice and liquid properties provided in separate arrays

        call file%get('q_liquid', prop_2d, do_transp = .false.)   ! kg/kg
        call file%get('q_ice', prop_2d_2, do_transp = .false.)   ! kg/kg
        call file%get('re_liquid', prop_2d_3, do_transp = .false.)   ! m
        call file%get('re_ice',  prop_2d_4, do_transp = .false.)   ! m

        do b = 1, nblocks
          allocate(cloud(b)%mixing_ratio(blocksize,nlev,2))
          allocate(cloud(b)%effective_radius(blocksize,nlev,2))

          cloud(b)%mixing_ratio(:,:,1) = transpose(prop_2d(:,istartcols(b):iendcols(b)))
          cloud(b)%mixing_ratio(:,:,2) = transpose(prop_2d_2(:,istartcols(b):iendcols(b)))
          cloud(b)%effective_radius(:,:,1) = transpose(prop_2d_3(:,istartcols(b):iendcols(b)))
          cloud(b)%effective_radius(:,:,2) = transpose(prop_2d_4(:,istartcols(b):iendcols(b)))
        end do
        deallocate(prop_2d_2, prop_2d_3, prop_2d_4)
      end if

      deallocate(prop_2d)

      ! For backwards compatibility, associate pointers for liquid and
      ! ice to the first and second slices of cloud%mixing_ratio and
      ! cloud%effective_radius
      do b = 1, nblocks
        cloud(b)%q_liq  => cloud(b)%mixing_ratio(:,:,1)
        cloud(b)%q_ice  => cloud(b)%mixing_ratio(:,:,2)
        cloud(b)%re_liq => cloud(b)%effective_radius(:,:,1)
        cloud(b)%re_ice => cloud(b)%effective_radius(:,:,2)
      end do
      
      ! Simple initialization of the seeds for the Monte Carlo scheme
      do b = 1, nblocks
        call single_level(b)%init_seed_simple(1,blocksize)
      end do
      ! Overwrite with user-specified values if available
      if (file%exists('iseed')) then
        call file%get('iseed', prop_1d_int)              ! read iseed as real
        do b = 1, nblocks
          single_level(b)%iseed = prop_1d_int(istartcols(b):iendcols(b))
        end do
      end if

      ! Cloud overlap parameter
      if (file%exists('overlap_param')) then
        call file%get('overlap_param', prop_2d, do_transp =.false.)
        do b = 1, nblocks
          cloud(b)%overlap_param = transpose(prop_2d(:,istartcols(b):iendcols(b)))
        end do
        deallocate(prop_2d)
      end if

      ! Optional scaling of liquid water mixing ratio
      if (driver_config%q_liq_scaling >= 0.0_jprb &
           &  .and. driver_config%q_liq_scaling /= 1.0_jprb) then
        do b = 1,nblocks
          cloud(b)%q_liq = cloud(b)%q_liq * driver_config%q_liq_scaling
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling liquid water mixing ratio by a factor of ', &
               &  driver_config%q_liq_scaling
        end if
      end if

      ! Optional scaling of ice water mixing ratio
      if (driver_config%q_ice_scaling >= 0.0_jprb .and. driver_config%q_ice_scaling /= 1.0_jprb) then
        do b = 1,nblocks
          cloud(b)%q_ice = cloud(b)%q_ice * driver_config%q_ice_scaling
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling ice water mixing ratio by a factor of ', &
               &  driver_config%q_ice_scaling
        end if
      end if

      ! Optional scaling of cloud fraction
      if (driver_config%cloud_fraction_scaling >= 0.0_jprb &
           &  .and. driver_config%cloud_fraction_scaling /= 1.0_jprb) then
        do b = 1,nblocks
          cloud(b)%fraction = cloud(b)%fraction * driver_config%cloud_fraction_scaling
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling cloud_fraction by a factor of ', &
               &  driver_config%cloud_fraction_scaling
        end if
      end if

      ! Cloud overlap is currently treated by an overlap decorrelation
      ! length (m) that is constant everywhere, and specified in one
      ! of the namelists
      if (driver_config%overlap_decorr_length_override > 0.0_jprb) then
        ! Convert overlap decorrelation length to overlap parameter between
        ! adjacent layers, stored in cloud%overlap_param
        do b = 1,nblocks
          call cloud(b)%set_overlap_param(thermodynamics(b), &
              &    driver_config%overlap_decorr_length_override)
        end do
      else if (.not. allocated(cloud(1)%overlap_param)) then 
        if (driver_config%iverbose >= 1) then
          write(nulout,'(a,g10.3,a)') 'Warning: overlap decorrelation length set to ', &
               &  decorr_length_default, ' m'
        end if
        do b = 1, nblocks
          call cloud(b)%set_overlap_param(thermodynamics(b), decorr_length_default)
        end do
      else if (driver_config%overlap_decorr_length_scaling > 0.0_jprb) then
        ! Scale the overlap decorrelation length by taking the overlap
        ! parameter to a power
        !    where (cloud%overlap_param > 0.99_jprb) cloud%overlap_param = 0.99_jprb
        do b = 1, nblocks
          where (cloud(b)%overlap_param > 0.0_jprb) 
            cloud(b)%overlap_param = cloud(b)%overlap_param**(1.0_jprb / driver_config%overlap_decorr_length_scaling)
          end where
        end do
        
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Scaling overlap decorrelation length by a factor of ', &
               &  driver_config%overlap_decorr_length_scaling
        end if
      else if (driver_config%overlap_decorr_length_scaling == 0.0_jprb) then
        do b = 1, nblocks
          cloud(b)%overlap_param = 0.0_jprb
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a)')  '  Setting overlap decorrelation length to zero (random overlap)'
        end if
      end if
      
      ! Cloud inhomogeneity is specified by the fractional standard
      ! deviation of cloud water content, that is currently constant
      ! everywhere (and the same for water and ice). The following copies
      ! this constant into the cloud%fractional_std array.
      if (driver_config%fractional_std_override >= 0.0_jprb) then
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3,a)') '  Overriding cloud fractional standard deviation with ', &
               &  driver_config%fractional_std_override
        end if
        prop_scalar = driver_config%fractional_std_override
        do b = 1,nblocks
          call cloud(b)%create_fractional_std(blocksize, nlev, prop_scalar)
        end do
      else if (.not. allocated(cloud(1)%fractional_std)) then
        prop_scalar = 0.0_jprb
        if (driver_config%iverbose >= 1) then
          write(nulout,'(a)') 'Warning: cloud optical depth fractional standard deviation set to zero'
        end if
        do b = 1,nblocks
          call cloud(b)%create_fractional_std(blocksize, nlev, prop_scalar)
        end do
      end if


      ! --------------------------------------------------------
      ! Read cloud properties needed by SPARTACUS
      ! --------------------------------------------------------

      if (config%i_solver_sw == ISolverSPARTACUS &
           &  .or.   config%i_solver_lw == ISolverSPARTACUS) then

        ! 3D radiative effects are governed by the length of cloud
        ! edge per area of gridbox, which is characterized by the
        ! inverse of the cloud effective size (m-1). Order of
        ! precedence: (1) effective size namelist overrides, (2)
        ! separation namelist overrides, (3) inv_cloud_effective_size
        ! present in NetCDF, (4) inv_cloud_effective_separation
        ! present in NetCDF. Only in the latter two cases may the
        ! effective size be scaled by the namelist variable
        ! "effective_size_scaling".

        is_cloud_size_scalable = .false. ! Default for cases (1) and (2)

        if (driver_config%low_inv_effective_size_override >= 0.0_jprb &
             &  .or. driver_config%middle_inv_effective_size_override >= 0.0_jprb &
             &  .or. driver_config%high_inv_effective_size_override >= 0.0_jprb) then
          ! (1) Cloud effective size specified in namelist

          ! First check all three ranges provided
          if (driver_config%low_inv_effective_size_override < 0.0_jprb &
             &  .or. driver_config%middle_inv_effective_size_override < 0.0_jprb &
             &  .or. driver_config%high_inv_effective_size_override < 0.0_jprb) then
            write(nulout,'(a,a)') '*** Error: if one of [low|middle|high]_inv_effective_size_override', &
                 & ' is provided then all must be'
            stop
          end if
          if (driver_config%iverbose >= 2) then
            write(nulout,'(a,g10.3,a)') '  Overriding inverse cloud effective size with:'
            write(nulout,'(a,g10.3,a)') '    ', driver_config%low_inv_effective_size_override, &
                 &       ' m-1 (low clouds)'
            write(nulout,'(a,g10.3,a)') '    ', driver_config%middle_inv_effective_size_override, &
                 &       ' m-1 (mid-level clouds)'
            write(nulout,'(a,g10.3,a)') '    ', driver_config%high_inv_effective_size_override, &
                 &       ' m-1 (high clouds)'
          end if
          do b = 1, nblocks
            call cloud(b)%create_inv_cloud_effective_size_eta(blocksize, nlev, &
                &  thermodynamics(b)%pressure_hl, &
                &  driver_config%low_inv_effective_size_override, &
                &  driver_config%middle_inv_effective_size_override, &
                &  driver_config%high_inv_effective_size_override, 0.8_jprb, 0.45_jprb)
          end do

        else if (driver_config%cloud_separation_scale_surface > 0.0_jprb &
             &  .and. driver_config%cloud_separation_scale_toa > 0.0_jprb) then
          ! (2) Cloud separation scale provided in namelist

          if (driver_config%iverbose >= 2) then
            write(nulout,'(a)') '  Effective cloud separation parameterized versus eta:'
            write(nulout,'(a,f8.1,a)') '    ', &
                 &  driver_config%cloud_separation_scale_surface, ' m at the surface'
            write(nulout,'(a,f8.1,a)') '    ', &
                 &  driver_config%cloud_separation_scale_toa, ' m at top-of-atmosphere'
            write(nulout,'(a,f6.2)') '     Eta power is', &
                 &  driver_config%cloud_separation_scale_power
            write(nulout,'(a,f6.2)') '     Inhomogeneity separation scaling is', &
                 &  driver_config%cloud_inhom_separation_factor
          end if
          do b = 1, nblocks
            call cloud(b)%param_cloud_effective_separation_eta(blocksize, nlev, &
                &  thermodynamics(b)%pressure_hl, &
                &  driver_config%cloud_separation_scale_surface, &
                &  driver_config%cloud_separation_scale_toa, &
                &  driver_config%cloud_separation_scale_power, &
                &  driver_config%cloud_inhom_separation_factor)
          end do
        else if (file%exists('inv_cloud_effective_size')) then
          ! (3) NetCDF file contains cloud effective size

          is_cloud_size_scalable = .true.

          call file%get('inv_cloud_effective_size', prop_2d, do_transp = .false.)

          do b = 1, nblocks
            cloud(b)%inv_cloud_effective_size = transpose(prop_2d(:,istartcols(b):iendcols(b)))
          end do
          deallocate(prop_2d)
          ! For finer control we can specify the effective size for
          ! in-cloud inhomogeneities as well
          if (file%exists('inv_inhom_effective_size')) then
            if (.not. driver_config%do_ignore_inhom_effective_size) then
              call file%get('inv_inhom_effective_size', prop_2d, do_transp = .false.)  
              do b = 1, nblocks
                cloud(b)%inv_inhom_effective_size = transpose(prop_2d(:,istartcols(b):iendcols(b)))
              end do
            else
              if (driver_config%iverbose >= 1) then
                write(nulout,'(a)') 'Ignoring inv_inhom_effective_size so treated as equal to inv_cloud_effective_size'
                write(nulout,'(a)') 'Warning: ...this is unlikely to be accurate for cloud fraction near one'
              end if
            end if
          else
            if (driver_config%iverbose >= 1) then
              write(nulout,'(a)') 'Warning: inv_inhom_effective_size not set so treated as equal to inv_cloud_effective_size'
              write(nulout,'(a)') 'Warning: ...this is unlikely to be accurate for cloud fraction near one'
            end if
          end if
          
        else if (file%exists('inv_cloud_effective_separation')) then
          ! (4) Alternative way to specify cloud scale

          is_cloud_size_scalable = .true.
          
          call file%get('inv_cloud_effective_separation', prop_2d, do_transp = .false.)  
          do b = 1, nblocks
            allocate(cloud(b)%inv_cloud_effective_size(blocksize,nlev))
            allocate(cloud(b)%inv_inhom_effective_size(blocksize,nlev))
            where (cloud(b)%fraction > config%cloud_fraction_threshold &
            &  .and. cloud(b)%fraction < 1.0_jprb - config%cloud_fraction_threshold)
              ! Convert effective cloud separation to effective cloud
              ! size, noting divisions rather than multiplications
              ! because we're working in terms of inverse sizes
              cloud(b)%inv_cloud_effective_size = transpose(prop_2d(:,istartcols(b):iendcols(b))) &
                  &  / sqrt(cloud(b)%fraction*(1.0_jprb-cloud(b)%fraction))
            elsewhere
              cloud(b)%inv_cloud_effective_size = 0.0_jprb
            end where
          end do

          if (file%exists('inv_inhom_effective_separation')) then
            if (driver_config%iverbose >= 2) then
              write(nulout,'(a)') '  Effective size of clouds and their inhomogeneities being computed from input'
              write(nulout,'(a)') '  ...variables inv_cloud_effective_separation and inv_inhom_effective_separation'
            end if
            call file%get('inv_inhom_effective_separation', prop_2d, do_transp = .false.)  
            do b = 1, nblocks
              where (cloud(b)%fraction > config%cloud_fraction_threshold)
                ! Convert effective separation of cloud inhomogeneities
                ! to effective size of cloud inhomogeneities, assuming
                ! here that the Tripleclouds treatment of cloud
                ! inhomogeneity will divide the cloudy part of the area
                ! into regions of equal area
                cloud(b)%inv_inhom_effective_size =  transpose(prop_2d(:,istartcols(b):iendcols(b))) &
                    &  / sqrt(0.5_jprb*cloud(b)%fraction * (1.0_jprb-0.5_jprb*cloud(b)%fraction))
              elsewhere
                cloud(b)%inv_inhom_effective_size = 0.0_jprb
              end where
            end do
          else
            ! Assume that the effective separation of cloud
            ! inhomogeneities is equal to that of clouds but
            ! multiplied by a constant provided by the user; note that
            ! prop_2d at this point contains
            ! inv_cloud_effective_separation
            if (driver_config%iverbose >= 2) then
              write(nulout,'(a)') '  Effective size of clouds being computed from inv_cloud_effective_separation'
              write(nulout,'(a,f6.2,a)') '  ...and multiplied by ', driver_config%cloud_inhom_separation_factor, &
                   &  ' to get effective size of inhomogeneities'
            end if
            do b = 1, nblocks
              where (cloud(b)%fraction > config%cloud_fraction_threshold)
                ! Note divisions rather than multiplications because
                ! we're working in terms of inverse sizes
               cloud(b)%inv_inhom_effective_size = (1.0_jprb / driver_config%cloud_inhom_separation_factor) &
                    & * transpose(prop_2d(:,istartcols(b):iendcols(b))) &
                    &  / sqrt(0.5_jprb*cloud(b)%fraction * (1.0_jprb-0.5_jprb*cloud(b)%fraction))
              elsewhere
                cloud(b)%inv_inhom_effective_size = 0.0_jprb
              end where
            end do
          end if ! exists inv_inhom_effective_separation
          deallocate(prop_2d)
          
        else

          write(nulout,'(a)') '*** Error: SPARTACUS solver specified but cloud size not, either in namelist or input file'
          stop

        end if ! Select method of specifying cloud effective size
        
        ! In cases (3) and (4) above the effective size obtained from
        ! the NetCDF may be scaled by a namelist variable
        if (is_cloud_size_scalable .and. driver_config%effective_size_scaling > 0.0_jprb) then
          ! Scale cloud effective size
          do b = 1, nblocks
            cloud(b)%inv_cloud_effective_size = cloud(b)%inv_cloud_effective_size &
                &                         / driver_config%effective_size_scaling
            if (allocated(cloud(b)%inv_inhom_effective_size)) then
              if (driver_config%iverbose >= 2) then
                write(nulout, '(a,g10.3)') '  Scaling effective size of clouds and their inhomogeneities with ', &
                    &                           driver_config%effective_size_scaling
              end if
              cloud(b)%inv_inhom_effective_size = cloud(b)%inv_inhom_effective_size &
                  &                         / driver_config%effective_size_scaling
            else
              if (driver_config%iverbose >= 2) then
                write(nulout, '(a,g10.3)') '  Scaling cloud effective size with ', &
                    &                           driver_config%effective_size_scaling
              end if
            end if
         end do
        end if

      end if ! Using SPARTACUS solver

    end if ! do_cloud

    ! --------------------------------------------------------
    ! Read surface properties
    ! --------------------------------------------------------

    ! Surface properties
    if (file%exists('tile_representation')) then
      ! We have a complex representation
      write(nulout,'(a)') 'Error: blocked derived type not yet implemented for complex surface'
      stop
      ! single_level%is_simple_surface = .false.
      ! is_complex_surface = .true.
      ! call surface%read(file)
      ! if (config%use_canopy_full_spectrum_sw) then
      !   allocate(single_level%sw_albedo(ncol,config%n_g_sw))
      !   allocate(single_level%sw_albedo_direct(ncol,config%n_g_sw))
      ! else
      !   allocate(single_level%sw_albedo(ncol,surface%nalbedobands))
      !   allocate(single_level%sw_albedo_direct(ncol,surface%nalbedobands))
      ! end if

      ! ! Optional override of shortwave albedo
      ! if (driver_config%sw_albedo_override >= 0.0_jprb) then
      !   surface%sw_albedo = driver_config%sw_albedo_override
      !   if (allocated(surface%sw_albedo_direct)) then
      !     surface%sw_albedo_direct = driver_config%sw_albedo_override
      !   end if
      !   if (driver_config%iverbose >= 2) then
      !     write(nulout,'(a,g10.3)') '  Overriding shortwave albedo of each facet with ', &
      !          &  driver_config%sw_albedo_override
      !   end if
      ! end if

      ! if (config%use_canopy_full_spectrum_lw) then
      !   allocate(single_level%lw_emission(ncol,config%n_g_lw))
      !   allocate(single_level%lw_emissivity(ncol,config%n_g_lw))
      ! else
      !   allocate(single_level%lw_emission(ncol,surface%nemissbands))
      !   allocate(single_level%lw_emissivity(ncol,surface%nemissbands))
      ! end if

      ! ! Optional override of longwave emissivity
      ! if (driver_config%lw_emissivity_override >= 0.0_jprb) then
      !   surface%lw_emissivity = driver_config%lw_emissivity_override
      !   if (driver_config%iverbose >= 2) then
      !     write(nulout,'(a,g10.3)')  '  Overriding longwave emissivity of each facet with ', &
      !          &  driver_config%lw_emissivity_override
      !   end if
      ! end if

    else
      ! We have a "simple" representation with a single flat tile, so
      ! the "surface" structure is not used
      do b = 1, nblocks
        single_level(b)%is_simple_surface = .true.
      end do
      is_complex_surface = .false.

      ! Single-level variable with dimensions (ncol)
      if (file%exists('skin_temperature')) then
        call file%get('skin_temperature', prop_1d) ! K
        do b = 1,nblocks
          single_level(b)%skin_temperature = prop_1d(istartcols(b):iendcols(b))
        end do
      else
        do b = 1, nblocks
          allocate(single_level(b)%skin_temperature(blocksize))
          single_level(b)%skin_temperature = thermodynamics(b)%temperature_hl(istartcols(b):iendcols(b),nlev+1)
        end do
        if (driver_config%iverbose >= 1 .and. config%do_lw &
             &  .and. driver_config%skin_temperature_override < 0.0_jprb) then 
          write(nulout,'(a)') 'Warning: skin temperature set equal to lowest air temperature'
        end if
      end if

      if (driver_config%sw_albedo_override >= 0.0_jprb) then
        ! Optional override of shortwave albedo
        do b = 1,nblocks
          allocate(single_level(b)%sw_albedo(blocksize,1))
          single_level(b)%sw_albedo = driver_config%sw_albedo_override
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding shortwave albedo with ', &
               &  driver_config%sw_albedo_override
        end if
        !if (allocated(single_level%sw_albedo_direct)) then
        !  single_level%sw_albedo_direct = driver_config%sw_albedo_override
        !end if
      else
        ! Shortwave albedo is stored with dimensions (ncol,nalbedobands)
        if (file%get_rank('sw_albedo') == 1) then
          ! ...but if in the NetCDF file it has only dimension (ncol), in
          ! order that nalbedobands is correctly set to 1, we need to turn
          ! off transposition
          call file%get('sw_albedo',  prop_1d)!, do_transp=.false.)
          do b = 1, nblocks
            allocate(single_level(b)%sw_albedo(blocksize,1))
            single_level(b)%sw_albedo(:,1)         = prop_1d(istartcols(b):iendcols(b))
          end do

          if (file%exists('sw_albedo_direct')) then
            call file%get('sw_albedo_direct', prop_1d_2)!, do_transp=.false.)
            do b = 1, nblocks
              allocate(single_level(b)%sw_albedo_direct(blocksize,1))
              single_level(b)%sw_albedo_direct(:,1)  = prop_1d_2(istartcols(b):iendcols(b))
            end do
          end if

        else
          call file%get('sw_albedo',    prop_2d, do_transp=.false.)
          if (file%exists('sw_albedo_direct')) then
            call file%get('sw_albedo_direct', prop_2d_2, do_transp=.false.)
          end if
          do b = 1, nblocks
            single_level(b)%sw_albedo         = transpose(prop_2d(:,istartcols(b):iendcols(b)))
            if(allocated(prop_2d_2)) single_level(b)%sw_albedo_direct = transpose(prop_2d_2(:,istartcols(b):iendcols(b)))
          end do
          deallocate(prop_2d)
          if (allocated(prop_2d_2)) deallocate(prop_2d_2)
        end if
      end if

      ! Longwave emissivity
      if (driver_config%lw_emissivity_override >= 0.0_jprb) then
        ! Optional override of longwave emissivity
        do b = 1, nblocks
          allocate(single_level(b)%lw_emissivity(blocksize,1))
          single_level(b)%lw_emissivity = driver_config%lw_emissivity_override
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)')  '  Overriding longwave emissivity with ', &
               &  driver_config%lw_emissivity_override
        end if
      else
        if (file%get_rank('lw_emissivity') == 1) then
          call file%get('lw_emissivity',prop_1d)
          do b = 1,nblocks
            allocate(single_level(b)%lw_emissivity(blocksize,1))
            single_level(b)%lw_emissivity(:,1) = prop_1d(istartcols(b):iendcols(b))
          end do
        else
          call file%get('lw_emissivity',prop_2d, do_transp=.false.)
          do b = 1, nblocks
            single_level(b)%lw_emissivity = transpose(prop_2d(:,istartcols(b):iendcols(b)))
          end do
          deallocate(prop_2d)
        end if
      end if
    end if ! complex surface

    ! Optional override of skin temperature
    if (driver_config%skin_temperature_override >= 0.0_jprb) then
      if (is_complex_surface) then
        do b = 1, nblocks
          surface(b)%skin_temperature = driver_config%skin_temperature_override
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding skin_temperature of each facet with ', &
               &  driver_config%skin_temperature_override
        end if
      else
        do b = 1,nblocks
          single_level(b)%skin_temperature = driver_config%skin_temperature_override
        end do
        if (driver_config%iverbose >= 2) then
          write(nulout,'(a,g10.3)') '  Overriding skin_temperature with ', &
               &  driver_config%skin_temperature_override
        end if
      end if
    end if

    ! --------------------------------------------------------
    ! Read aerosol and gas concentrations
    ! --------------------------------------------------------

    if (config%use_aerosols) then
      ! Load aerosol data
      call file%get('aerosol_mmr', prop_3d, ipermute=[2,3,1]);
      ! call file%get('aerosol_mmr', prop_3d, ipermute=(/2,1,3/));  !  ntype, nlev, ncol

      do b = 1,nblocks
        aerosol(b)%mixing_ratio = prop_3d(istartcols(b):iendcols(b),:,:)
        ! Store aerosol level bounds
        aerosol(b)%istartlev = lbound(aerosol(b)%mixing_ratio, 2)
        aerosol(b)%iendlev   = ubound(aerosol(b)%mixing_ratio, 2)
      end do
      deallocate(prop_3d)
    end if

    ! Load in gas volume mixing ratios, which can be either 2D arrays
    ! (varying with height and column) or 0D scalars (constant volume
    ! mixing ratio everywhere).
    ngases          = 0 ! Gases with varying mixing ratio
    nwellmixedgases = 0 ! Gases with constant mixing ratio

    ! Water vapour and ozone are always in terms of mass mixing ratio
    ! (kg/kg) and always 2D arrays with dimensions (ncol,nlev), unlike
    ! other gases (see below)

    ! For RRTMGP: which gases in GasLowerCaseName are present in file?
    ! H2O and O3 are assumed to always be present
    !if (present(gas_rrtmgp)) then
    if (use_rrtmgp) then
      gas_inds = [1]
      do jgas = 2,NMaxGases
        if ((jgas == IO3)) then
          gas_inds = [gas_inds, jgas]
        else
          gas_var_name = trim(GasLowerCaseName(jgas)) // trim(driver_config%vmr_suffix_str)
          if (file%exists(trim(gas_var_name))) gas_inds = [gas_inds, jgas]
        end if
      end do
      allocate(rrtmgp_gas_names(size(gas_inds)+1))
      rrtmgp_gas_names(1) = "n2"
      rrtmgp_gas_names(2:size(rrtmgp_gas_names)) = GasLowerCaseName(gas_inds)
    end if

    ! if (present(gas_rrtmgp)) then
    !   allocate(rrtmgp_gas_names(size(gas_inds)+1))
    !   rrtmgp_gas_names(1) = "n2"
    !   rrtmgp_gas_names(2:size(rrtmgp_gas_names)) = GasLowerCaseName(gas_inds)
    ! end if

    ! Initialize derived types
    do b = 1,nblocks
      call gas(b)%allocate(blocksize, nlev)
      !if(present(gas_rrtmgp)) call stop_on_err(gas_rrtmgp(b)%init(rrtmgp_gas_names))
      if (use_rrtmgp) call stop_on_err(gas(b)%gas_rrtmgp%init(rrtmgp_gas_names))

      ! if(present(gas_rrtmgp)) call stop_on_err(gas_rrtmgp(b)%init(GasLowerCaseName(gas_inds)))
    end do
    
    if (use_rrtmgp) then
      ! Add nitrogen to RRTMGP gas type
      gas_var_name = "n2" // trim(driver_config%vmr_suffix_str)
      if(file%exists('gas_var_name')) then
        irank = file%get_rank(trim(gas_var_name))
        ! Note that if the gas is not present then a warning will have
        ! been issued, and irank will be returned as -1
        if (irank == 0) then
          ! Store this as a well-mixed gas
          call file%get(trim(gas_var_name), well_mixed_gas_vmr)
          do b = 1,nblocks
            call stop_on_err(gas(b)%gas_rrtmgp%set_vmr("n2",well_mixed_gas_vmr))
          end do
        else if (irank == 2) then
          call file%get(trim(gas_var_name), gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call stop_on_err(gas(b)%gas_rrtmgp%set_vmr("n2", gas_mr(:,istartcols(b):iendcols(b))))
          end do
        else if (irank > 0) then
          write(nulout,'(a,a,a)')  '***  Error: ', trim(gas_var_name), ' does not have 0 or 2 dimensions'
          stop
        end if
      else
        do b = 1,nblocks
          call stop_on_err(gas(b)%gas_rrtmgp%set_vmr("n2", 0.781000018_jprb))
        end do
      end if 
   end if

    ! Loop through all radiatively important gases
    do jgas = 1,NMaxGases
      if (jgas == IH2O) then
        if (file%exists('q')) then
          call file%get('q', gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call gas(b)%put(IH2O, IMassMixingRatio, transpose(gas_mr(:,istartcols(b):iendcols(b))))
            sf = AirMolarMass / GasMolarMass(IH2O)
          end do
        else if (file%exists('h2o_mmr')) then
          call file%get('h2o_mmr', gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call gas(b)%put(IH2O, IMassMixingRatio, transpose(gas_mr(:,istartcols(b):iendcols(b))))
            sf = AirMolarMass / GasMolarMass(IH2O)
          end do
        else
          call file%get('h2o' // trim(driver_config%vmr_suffix_str), gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call gas(b)%put(IH2O, IVolumeMixingRatio, transpose(gas_mr(:,istartcols(b):iendcols(b))))
            sf = 1.0_jprb
          end do
        end if
        if(use_rrtmgp) then
          do b = 1, nblocks
            call stop_on_err(gas(b)%gas_rrtmgp%set_vmr(GasLowerCaseName(IH2O), sf*gas_mr(:,istartcols(b):iendcols(b))))
          end do
        end if
      else if (jgas == IO3) then
        if (file%exists('o3_mmr')) then
          call file%get('o3_mmr', gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call gas(b)%put(IO3, IMassMixingRatio, transpose(gas_mr(:,istartcols(b):iendcols(b))))
          end do
          sf = AirMolarMass / GasMolarMass(IO3)
        else
          call file%get('o3' // trim(driver_config%vmr_suffix_str), gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call gas(b)%put(IO3, IVolumeMixingRatio, transpose(gas_mr(:,istartcols(b):iendcols(b))))
          end do
          sf = 1.0_jprb
        end if
        if(use_rrtmgp) then
          do b = 1, nblocks
            call stop_on_err(gas(b)%gas_rrtmgp%set_vmr(GasLowerCaseName(IO3), sf*gas_mr(:,istartcols(b):iendcols(b))))
          end do
        end if
      else
        ! Find number of dimensions of the variable holding gas "jgas" in
        ! the input file, where the following function returns -1 if the
        ! gas is not found
        gas_var_name = trim(GasLowerCaseName(jgas)) // trim(driver_config%vmr_suffix_str)
        irank = file%get_rank(trim(gas_var_name))
        ! Note that if the gas is not present then a warning will have
        ! been issued, and irank will be returned as -1
        if (irank == 0) then
          ! Store this as a well-mixed gas
          call file%get(trim(gas_var_name), well_mixed_gas_vmr)
          do b = 1,nblocks
            call gas(b)%put_well_mixed(jgas, IVolumeMixingRatio, well_mixed_gas_vmr)
            if(use_rrtmgp) call stop_on_err(gas(b)%gas_rrtmgp%set_vmr(GasLowerCaseName(jgas),well_mixed_gas_vmr))
          end do
        else if (irank == 2) then
          call file%get(trim(gas_var_name), gas_mr, do_transp=.false.)
          do b = 1,nblocks
            call gas(b)%put(jgas, IVolumeMixingRatio, transpose(gas_mr(:,istartcols(b):iendcols(b))))
            if(use_rrtmgp) then
              call stop_on_err(gas(b)%gas_rrtmgp%set_vmr(GasLowerCaseName(jgas), gas_mr(:,istartcols(b):iendcols(b))))
            end if
          end do
        else if (irank > 0) then
          write(nulout,'(a,a,a)')  '***  Error: ', trim(gas_var_name), ' does not have 0 or 2 dimensions'
          stop
        end if
      end if
      if (allocated(gas_mr)) deallocate(gas_mr)
    end do

    ! do b = 1, size(gas_rrtmgp(1)%gas_name)
    !   print *," max of gas ",gas_rrtmgp(1)%gas_name(b), maxval(gas_rrtmgp(1)%concs(b)%conc)
    ! end do

    ! Scale gas concentrations if needed
    do b = 1,nblocks
      call gas(b)%scale(IH2O,    driver_config%h2o_scaling,    driver_config%iverbose >= 2)
      call gas(b)%scale(ICO2,    driver_config%co2_scaling,    driver_config%iverbose >= 2)
      call gas(b)%scale(IO3,     driver_config%o3_scaling,     driver_config%iverbose >= 2)
      call gas(b)%scale(IN2O,    driver_config%n2o_scaling,    driver_config%iverbose >= 2)
      call gas(b)%scale(ICO,     driver_config%co_scaling,     driver_config%iverbose >= 2)
      call gas(b)%scale(ICH4,    driver_config%ch4_scaling,    driver_config%iverbose >= 2)
      call gas(b)%scale(IO2,     driver_config%o2_scaling,     driver_config%iverbose >= 2)
      call gas(b)%scale(ICFC11,  driver_config%cfc11_scaling,  driver_config%iverbose >= 2)
      call gas(b)%scale(ICFC12,  driver_config%cfc12_scaling,  driver_config%iverbose >= 2)
      call gas(b)%scale(IHCFC22, driver_config%hcfc22_scaling, driver_config%iverbose >= 2)
      call gas(b)%scale(ICCL4,   driver_config%ccl4_scaling,   driver_config%iverbose >= 2)
      call gas(b)%scale(INO2,    driver_config%no2_scaling,    driver_config%iverbose >= 2)
    end do
  end subroutine read_input_blocked

    !---------------------------------------------------------------------
  ! Deallocate flux arrays
  subroutine unblock_fluxes(ib, blocksize, thermodynamics_b, flux_b, thermodynamics, flux )
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_flux,   only : flux_type
    use parkind1,         only : jprb, jpim
    use yomhook,          only : lhook, dr_hook
    integer,                              intent(in)    :: ib, blocksize
    type(thermodynamics_type),            intent(in)    :: thermodynamics_b
    type(flux_type),                      intent(in)    :: flux_b
    type(thermodynamics_type),            intent(inout) :: thermodynamics
    type(flux_type),                      intent(inout) :: flux

    real(jprb)                      :: hook_handle
    integer                         :: i1, i2

    if (lhook) call dr_hook('radiation_flux:unblock_fluxes',0,hook_handle)

    ! do ib = 1, nblocks
      i1 = (ib-1) * blocksize + 1
      i2 = i1 + blocksize - 1

      thermodynamics%pressure_hl(i1:i2,:)     = thermodynamics_b%pressure_hl
      thermodynamics%temperature_hl(i1:i2,:)  = thermodynamics_b%temperature_hl

      if (allocated(thermodynamics%h2o_sat_liq)) thermodynamics%h2o_sat_liq(i1:i2,:)  = thermodynamics_b%h2o_sat_liq

      if (allocated(flux%lw_up))       flux%lw_up(i1:i2,:) = flux_b%lw_up
      if (allocated(flux%lw_dn))       flux%lw_dn(i1:i2,:) = flux_b%lw_dn
      if (allocated(flux%lw_up_clear)) flux%lw_up_clear(i1:i2,:) = flux_b%lw_up_clear
      if (allocated(flux%lw_dn_clear)) flux%lw_up_clear(i1:i2,:) = flux_b%lw_up_clear
    
      if (allocated(flux%sw_up))        flux%sw_up(i1:i2,:) = flux_b%sw_up
      if (allocated(flux%sw_dn))        flux%sw_dn(i1:i2,:) = flux_b%sw_dn
      if (allocated(flux%sw_up_clear))  flux%sw_up_clear(i1:i2,:) = flux_b%sw_up_clear
      if (allocated(flux%sw_dn_clear))  flux%sw_dn_clear(i1:i2,:) = flux_b%sw_dn_clear
      if (allocated(flux%sw_dn_direct)) flux%sw_dn_direct(i1:i2,:) = flux_b%sw_dn_direct
      if (allocated(flux%sw_dn_direct_clear)) flux%sw_dn_direct_clear(i1:i2,:) = flux_b%sw_dn_direct_clear
  
      if (allocated(flux%lw_up_band))       flux%lw_up_band(:,i1:i2,:) = flux_b%lw_up_band
      if (allocated(flux%lw_dn_band))       flux%lw_dn_band(:,i1:i2,:) = flux_b%lw_dn_band
      if (allocated(flux%lw_up_clear_band)) flux%lw_up_clear_band(:,i1:i2,:) = flux_b%lw_up_clear_band
      if (allocated(flux%lw_dn_clear_band)) flux%lw_dn_clear_band(:,i1:i2,:) = flux_b%lw_dn_clear_band
    
      if (allocated(flux%sw_up_band))       flux%sw_up_band(:,i1:i2,:) = flux_b%sw_up_band 
      if (allocated(flux%sw_dn_band))       flux%sw_dn_band(:,i1:i2,:) = flux_b%sw_dn_band
      if (allocated(flux%sw_up_clear_band)) flux%sw_up_clear_band(:,i1:i2,:) = flux_b%sw_up_clear_band
      if (allocated(flux%sw_dn_clear_band))  flux%sw_dn_clear_band(:,i1:i2,:) = flux_b%sw_dn_clear_band
      if (allocated(flux%sw_dn_direct_band)) flux%sw_dn_direct_band(:,i1:i2,:) = flux_b%sw_dn_direct_band
      if (allocated(flux%sw_dn_direct_clear_band)) flux%sw_dn_direct_clear_band(:,i1:i2,:) = flux_b%sw_dn_direct_clear_band
  
      if (allocated(flux%sw_dn_surf_band)) then
        flux%sw_dn_surf_band(:,i1:i2) = flux_b%sw_dn_surf_band 
        flux%sw_dn_direct_surf_band(:,i1:i2) = flux_b%sw_dn_direct_surf_band 
      end if

      if (allocated(flux%sw_dn_surf_clear_band)) then
        flux%sw_dn_surf_clear_band(:,i1:i2) = flux_b%sw_dn_surf_clear_band 
        flux%sw_dn_direct_surf_clear_band(:,i1:i2) = flux_b%sw_dn_direct_surf_clear_band 
      end if
  
      if (allocated(flux%lw_dn_surf_canopy)) flux%lw_dn_surf_canopy(:,i1:i2) = flux_b%lw_dn_surf_canopy
      if (allocated(flux%sw_dn_diffuse_surf_canopy)) flux%sw_dn_diffuse_surf_canopy(:,i1:i2) = flux_b%sw_dn_diffuse_surf_canopy
      if (allocated(flux%sw_dn_direct_surf_canopy)) flux%sw_dn_direct_surf_canopy(:,i1:i2) = flux_b%sw_dn_direct_surf_canopy
  
      if (allocated(flux%cloud_cover_sw)) flux%cloud_cover_sw(i1:i2) = flux_b%cloud_cover_sw 
      if (allocated(flux%cloud_cover_lw)) flux%cloud_cover_lw(i1:i2) = flux_b%cloud_cover_lw 
  
      if (allocated(flux%lw_derivatives)) flux%lw_derivatives(i1:i2,:) = flux_b%lw_derivatives
  

    ! end do
  
    if (lhook) call dr_hook('radiation_flux:deallocate',1,hook_handle)

  end subroutine unblock_fluxes

  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    use iso_c_binding
    character(len=*), intent(in) :: error_msg
  
    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "read_input stopping"
      stop
    end if
  end subroutine stop_on_err

end module ecrad_driver_read_input
