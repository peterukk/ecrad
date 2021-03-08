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
  
    public  :: setup_gas_optics_ifs_rrtmgp, gas_optics_ifs_rrtmgp
  
  contains
  
    !---------------------------------------------------------------------
    ! Setup the IFS implementation of RRTMGP gas absorption model
    subroutine setup_gas_optics_ifs_rrtmgp(config, gas_rrtmgp, solar_irradiance)

     use parkind1,         only : jprb
     use yomhook,          only : lhook, dr_hook
     use radiation_config
     use radiation_monochromatic,  only : &
          &   setup_cloud_optics_mono   => setup_cloud_optics, &
          &   setup_aerosol_optics_mono => setup_aerosol_optics
     use radiation_cloud_optics,   only :  setup_cloud_optics
     use radiation_aerosol_optics, only :  setup_aerosol_optics
     use mo_load_coefficients,     only:   load_and_init
     use mo_gas_concentrations,    only :  ty_gas_concs
     use mo_gas_optics_rrtmgp,     only :  ty_gas_optics_rrtmgp

     type(config_type),      intent(inout), target    :: config
     type(ty_gas_concs),     intent(in)               :: gas_rrtmgp
     real(jprb),             intent(in)               :: solar_irradiance  

     integer :: irep ! For implied do

     integer, parameter :: RRTMGP_GPOINT_REORDERING_LW(256)  = (/ & 
          & 225, 241, 226, 242, 227, 193, 228, 113, 114, 115, 97, 116, 194, 98, 81, 117, 229,  & 
          & 243, 129, 82, 83, 118, 99, 84, 195, 85, 230, 86, 87, 100, 65, 119, 231, 177, 244, 88,  & 
          & 66, 196, 130, 101, 232, 67, 102, 245, 233, 120, 89, 178, 197, 68, 131, 234, 103, 235,  & 
          & 246, 198, 236, 90, 69, 179, 237, 33, 104, 91, 132, 238, 121, 34, 199, 92, 247, 239, 35,  & 
          & 70, 93, 180, 133, 36, 240, 94, 105, 200, 122, 17, 95, 209, 37, 71, 123, 181, 248, 96,  & 
          & 134, 124, 106, 38, 145, 201, 18, 107, 182, 72, 49, 135, 108, 125, 39, 146, 109, 210,  & 
          & 202, 249, 161, 19, 110, 203, 183, 147, 50, 126, 204, 111, 73, 136, 148, 112, 51, 162,  & 
          & 40, 205, 20, 250, 1, 127, 52, 128, 251, 163, 149, 211, 184, 74, 206, 252, 21, 53, 164,  & 
          & 75, 253, 41, 137, 150, 54, 254, 76, 22, 207, 165, 2, 255, 212, 77, 256, 208, 185, 42,  & 
          & 55, 78, 151, 166, 3, 23, 138, 43, 213, 139, 4, 79, 44, 56, 167, 186, 214, 140, 152, 45,  & 
          & 5, 24, 80, 187, 141, 215, 46, 168, 6, 188, 57, 47, 142, 48, 153, 189, 216, 7, 25, 190,  & 
          & 169, 143, 58, 154, 155, 191, 144, 59, 156, 8, 170, 157, 192, 26, 171, 217, 60, 158, & 
          & 172, 61, 27, 159, 173, 160, 62, 174, 9, 63, 28, 175, 218, 64, 176, 219, 29, 220, 10, & 
          & 11, 221, 30, 12, 222, 13, 223, 14, 224, 31, 15, 32, 16 /)
     integer, parameter :: RRTMGP_GPOINT_REORDERING_SW(224)  = (/ & 
          & 97, 81, 65, 113, 129, 49, 130, 82, 114, 131, 98, 115, 50, 66, 145, 116, 132, 146, & 
          & 17, 117, 147, 51, 148, 1, 83, 118, 149, 133, 67, 150, 119, 161, 52, 162, 18, 151, 163, & 
          & 120, 2, 164, 99, 165, 134, 68, 152, 166, 53, 167, 168, 169, 177, 170, 19, 171, 178, & 
          & 172, 173, 179, 121, 69, 3, 180, 54, 174, 181, 84, 135, 182, 183, 184, 20, 185, 186, & 
          & 187, 188, 189, 190, 191, 192, 175, 153, 100, 193, 70, 122, 55, 194, 176, 136, 21, 123, & 
          & 33, 71, 195, 124, 4, 154, 56, 101, 22, 125, 155, 85, 72, 137, 126, 196, 34, 23, 102, & 
          & 127, 156, 57, 128, 86, 5, 73, 138, 197, 24, 103, 157, 139, 209, 58, 87, 35, 210, 198, & 
          & 74, 140, 59, 211, 75, 104, 199, 60, 25, 212, 6, 141, 200, 76, 158, 88, 36, 201, 202, & 
          & 203, 204, 205, 206, 207, 208, 213, 61, 77, 142, 214, 215, 216, 217, 37, 218, 219, 220, & 
          & 221, 222, 224, 223, 26, 105, 7, 27, 62, 78, 143, 28, 159, 38, 89, 29, 144, 63, 106, 30, & 
          & 160, 79, 39, 8, 64, 31, 107, 32, 90, 108, 80, 91, 40, 109, 92, 93, 110, 9, 94, 111, 41, & 
          & 112, 95, 96, 10, 42, 11, 43, 12, 44, 45, 13, 46, 47, 48, 14, 15, 16 /)
     real(jprb) :: hook_handle

     if (lhook) call dr_hook('radiation_ifs_rrtmgp:setup_gas_optics_ifs_rrtmgp',0,hook_handle)


    ! Load neural nets (optional)
     if (config%do_rrtmgp_neural_nets) then
          call config%rrtmgp_neural_nets(1) % load(trim(config%rrtmgp_neural_net_sw_tau))
          call config%rrtmgp_neural_nets(2) % load(trim(config%rrtmgp_neural_net_sw_ray))
          call config%rrtmgp_neural_nets(3) % load(trim(config%rrtmgp_neural_net_lw_tau))
          call config%rrtmgp_neural_nets(4) % load(trim(config%rrtmgp_neural_net_lw_pfrac))
     end if
     ! Load k-distributions
     call load_and_init(config%k_dist_lw, trim(config%rrtmgp_gas_optics_file_name_lw), gas_rrtmgp)
     call load_and_init(config%k_dist_sw, trim(config%rrtmgp_gas_optics_file_name_sw), gas_rrtmgp)

     ! Scale the spectral solar source function with user-provided solar irradiance
     call stop_on_err(config%k_dist_sw%set_tsi(solar_irradiance ))

     ! Cloud and aerosol properties can only be defined per band
     config%do_cloud_aerosol_per_sw_g_point = .false.
     config%do_cloud_aerosol_per_lw_g_point = .false.

     config%n_g_sw = config%k_dist_sw%get_ngpt()
     config%n_g_lw = config%k_dist_lw%get_ngpt()

     config%wavenumber1_sw = config%k_dist_sw%band_lims_wvn(1,:)
     config%wavenumber2_sw = config%k_dist_sw%band_lims_wvn(2,:)

     config%wavenumber1_lw = config%k_dist_lw%band_lims_wvn(1,:)
     config%wavenumber2_lw = config%k_dist_lw%band_lims_wvn(2,:)

     config%i_band_from_g_sw =  config%k_dist_sw%get_gpoint_bands() 
     config%i_band_from_g_lw =  config%k_dist_lw%get_gpoint_bands() 

     config%n_bands_sw = config%k_dist_sw%get_nband()
     config%n_bands_lw = config%k_dist_lw%get_nband()

     allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
     allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
     allocate(config%i_g_from_reordered_g_sw(config%n_g_sw))
     allocate(config%i_g_from_reordered_g_lw(config%n_g_lw))

     ! Store band positions if using generalized cloud or aerosol
     call config%gas_optics_sw%spectral_def%allocate_bands_only(config%wavenumber1_sw, &
          &                                                     config%wavenumber2_sw)
     call config%gas_optics_lw%spectral_def%allocate_bands_only(config%wavenumber1_lw, &
          &                                                     config%wavenumber2_lw)

     if (config%i_solver_sw == ISolverSpartacus) then
          ! SPARTACUS requires g points ordered in approximately
          ! increasing order of optical depth
          config%i_g_from_reordered_g_sw = RRTMGP_GPOINT_REORDERING_SW
     else
          ! Implied-do for no reordering
          config%i_g_from_reordered_g_sw = (/ (irep, irep=1,config%n_g_sw) /)
     end if

     if (config%i_solver_lw == ISolverSpartacus) then
          ! SPARTACUS requires g points ordered in approximately
          ! increasing order of optical depth
          config%i_g_from_reordered_g_lw = RRTMGP_GPOINT_REORDERING_LW
     else
          ! Implied-do for no reordering
          config%i_g_from_reordered_g_lw = (/ (irep, irep=1,config%n_g_lw) /)
     end if

     config%i_band_from_reordered_g_sw &
          = config%i_band_from_g_sw(config%i_g_from_reordered_g_sw)

     config%i_band_from_reordered_g_lw &
          = config%i_band_from_g_lw(config%i_g_from_reordered_g_lw)

     ! The i_spec_* variables are used solely for storing spectral
     ! data, and this can either be by band or by g-point
     if (config%do_save_spectral_flux) then
          if (config%do_save_gpoint_flux) then
               config%n_spec_sw = config%n_g_sw
               config%n_spec_lw = config%n_g_lw
               config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
               config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
          else
               config%n_spec_sw = config%n_bands_sw
               config%n_spec_lw = config%n_bands_lw
               config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
               config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
          end if
     else
          config%n_spec_sw = 0
          config%n_spec_lw = 0
          nullify(config%i_spec_from_reordered_g_sw)
          nullify(config%i_spec_from_reordered_g_lw)
     end if

    end subroutine setup_gas_optics_ifs_rrtmgp
  
  
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
#ifdef USE_TIMING
    ! Timing library
     use gptl,                  only: gptlstart, gptlstop, gptlinitialize, gptlpr, &
          & gptlfinalize, gptlsetoption, gptlpercent, gptloverhead
#endif
  
     integer, intent(in) :: ncol               ! number of columns
     integer, intent(in) :: nlev               ! number of model levels
     type(config_type), intent(in) :: config
     type(single_level_type),  intent(in) :: single_level
     type(thermodynamics_type),intent(in) :: thermodynamics
     type(gas_type),           intent(in) :: gas

     ! Longwave albedo of the surface
     real(jprb), dimension(config%n_g_lw,ncol), intent(in)  :: lw_albedo

     ! Gaseous layer optical depth in longwave and shortwave, and
     ! shortwave single scattering albedo (i.e. fraction of extinction
     ! due to Rayleigh scattering) at each g-point
     real(jprb), dimension(config%n_g_lw,nlev,ncol), intent(out) :: od_lw
     real(jprb), dimension(config%n_g_sw,nlev,ncol), intent(out) :: od_sw, ssa_sw

     ! The Planck function (emitted flux from a black body) at half
     ! levels at each longwave g-point
     real(jprb), dimension(config%n_g_lw,nlev+1,ncol), intent(out)   :: planck_hl
     ! Planck function for the surface (W m-2)
     real(jprb), dimension(config%n_g_lw,ncol),      intent(out)     :: lw_emission

     ! The incoming shortwave flux into a plane perpendicular to the
     ! incoming radiation at top-of-atmosphere in each of the shortwave
     ! g-points
     real(jprb), dimension(config%n_g_sw,ncol),  intent(out)         :: incoming_sw

     ! RRTMGP sources are provided in W/m2-str; factor of pi converts to flux units
     real(jprb), parameter :: pi = acos(-1._jprb)

     integer :: ret, jcol, jlev

     real(jprb) :: hook_handle

     
     if (lhook) call dr_hook('radiation_ifs_rrtm:gas_optics',1,hook_handle)

#ifdef USE_TIMING
    ret =  gptlstart('gas_optics_sw')
#endif     

     if (config%do_sw) then

     ! Scale the spectral solar source function with user-provided solar irradiance
     ! call stop_on_err(config%k_dist_sw%set_tsi(single_level%solar_irradiance ))

          if (config%do_rrtmgp_neural_nets) then
          
               call stop_on_err(config%k_dist_sw%gas_optics_ext_ecrad( &
                    &   ncol, nlev, config%n_g_sw, &
                    &   thermodynamics%pressure_fl_reverse, &
                    &   thermodynamics%pressure_hl_reverse, &
                    &   thermodynamics%temperature_fl_reverse, &
                    &   gas%gas_rrtmgp, &
                    &   od_sw, ssa_sw, incoming_sw, &
                    &   neural_nets = config%rrtmgp_neural_nets(1:2))) 
          else 
               call stop_on_err(config%k_dist_sw%gas_optics_ext_ecrad( &
                    &   ncol, nlev, config%n_g_sw, &
                    &   thermodynamics%pressure_fl_reverse, &
                    &   thermodynamics%pressure_hl_reverse, &
                    &   thermodynamics%temperature_fl_reverse, &
                    &   gas%gas_rrtmgp, &
                    &   od_sw, ssa_sw, incoming_sw)) 
          end if
          
          if (config%i_solver_lw == ISolverSpartacus) then
               !    if (.true.) then
               ! We need to rearrange the gas optics info in memory: reordering
               ! the g points in order of approximately increasing optical
               ! depth (for efficient 3D processing on only the regions of the
               ! spectrum that are optically thin for gases) and reorder in
               ! pressure since the the functions above treat pressure
               ! decreasing with increasing index.  Note that the output gas
               ! arrays have dimensions in a different order to the inputs,
               ! so there is some inefficiency here.
               do jcol = 1,ncol
                    incoming_sw(:,jcol) = incoming_sw(config%i_g_from_reordered_g_sw,jcol)
                    do jlev = 1, nlev
                         od_sw(:,jlev,jcol) = od_sw(config%i_g_from_reordered_g_sw,jlev,jcol)
                         ssa_sw(:,jlev,jcol) = ssa_sw(config%i_g_from_reordered_g_sw,jlev,jcol)
                    end do
               end do
          end if

     end if
#ifdef USE_TIMING
    ret =  gptlstop('gas_optics_sw')
    ret =  gptlstart('gas_optics_lw')
#endif  

     if (config%do_lw) then
          if (config%do_rrtmgp_neural_nets) then

               call stop_on_err(config%k_dist_lw%gas_optics_int_ecrad( &
                    &   ncol, nlev, config%n_g_lw, &
                    &   thermodynamics%pressure_fl_reverse, &
                    &   thermodynamics%pressure_hl_reverse, &
                    &   thermodynamics%temperature_fl_reverse, &
                    &   thermodynamics%temperature_hl_reverse, &
                    &   single_level%skin_temperature, &
                    &   gas%gas_rrtmgp, &
                    &   od_lw, lw_emission, planck_hl, &
                    &   neural_nets = config%rrtmgp_neural_nets(3:4))) 
          else 

               call stop_on_err(config%k_dist_lw%gas_optics_int_ecrad( &
               &   ncol, nlev, config%n_g_lw, &
               &   thermodynamics%pressure_fl_reverse, &
               &   thermodynamics%pressure_hl_reverse, &
               &   thermodynamics%temperature_fl_reverse, &
               &   thermodynamics%temperature_hl_reverse, &
               &   single_level%skin_temperature, &
               &   gas%gas_rrtmgp, &
               &   od_lw, lw_emission, planck_hl))
          end if
          if (single_level%is_simple_surface) then
               ! lw_emission at this point is actually the planck function of
               ! the surface
               lw_emission = lw_emission * (1.0_jprb - lw_albedo)
          else
          ! Longwave emission has already been computed
               if (config%use_canopy_full_spectrum_lw) then
               lw_emission = transpose(single_level%lw_emission(1:ncol,:))
               else
               lw_emission = transpose(single_level%lw_emission(1:ncol, &
                    & config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw)))
               end if
          end if
          ! RRTMGP sources are provided in W/m2-str; factor of pi converts to flux units
          planck_hl = planck_hl * pi 
          ! source%lay_source = source%lay_source * pi 
          lw_emission = lw_emission * pi

          if (config%i_solver_lw == ISolverSpartacus) then
               !    if (.true.) then
               ! We need to rearrange the gas optics info in memory: reordering
               ! the g points in order of approximately increasing optical
               ! depth (for efficient 3D processing on only the regions of the
               ! spectrum that are optically thin for gases) and reorder in
               ! pressure since the the functions above treat pressure
               ! decreasing with increasing index.  Note that the output gas
               ! arrays have dimensions in a different order to the inputs,
               ! so there is some inefficiency here.
               do jcol = 1,ncol
                    lw_emission(:,jcol) = lw_emission(config%i_g_from_reordered_g_lw,jcol)
                    do jlev = 1, nlev
                         od_lw(:,jlev,jcol) = od_lw(config%i_g_from_reordered_g_lw,jlev,jcol)
                         planck_hl(:,jlev,jcol) = planck_hl(config%i_g_from_reordered_g_lw,jlev,jcol)
                    end do
                    planck_hl(:,nlev+1,jcol) = planck_hl(config%i_g_from_reordered_g_lw,nlev+1,jcol)
               end do
          end if 
      end if

#ifdef USE_TIMING
    ret =  gptlstop('gas_optics_lw')
#endif 
      
    end subroutine gas_optics_ifs_rrtmgp


  subroutine stop_on_err(error_msg)
     use iso_fortran_env, only : error_unit
     use iso_c_binding
     character(len=*), intent(in) :: error_msg
   
     if(error_msg /= "") then
       write (error_unit,*) trim(error_msg)
       write (error_unit,*) "radiation_interface stopping"
       stop
     end if
   end subroutine stop_on_err

  
  end module radiation_ifs_rrtmgp
  