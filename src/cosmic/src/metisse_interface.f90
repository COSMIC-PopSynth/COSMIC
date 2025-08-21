module metisse_interface
    implicit none
    integer, parameter :: f2py_strlen = 256

    ! -----------------------
    ! Python track info
    ! -----------------------
    ! Arrays
    integer :: ntracks_h_in, max_ncol_h_in, max_neep_h_in, max_points_h_in
    logical :: is_he_in
    character(len=f2py_strlen), allocatable :: filenames_h_in(:)
    real(8), allocatable :: initial_mass_h_in(:)
    real(8), allocatable :: initial_Y_h_in(:)
    real(8), allocatable :: initial_Z_h_in(:)
    real(8), allocatable :: Fe_div_H_h_in(:)
    real(8), allocatable :: alpha_div_Fe_h_in(:)
    real(8), allocatable :: v_div_vcrit_h_in(:)
    integer, allocatable :: ntrack_arr_h_in(:)
    integer, allocatable :: neep_arr_h_in(:)
    integer, allocatable :: ncol_arr_h_in(:)
    logical, allocatable :: is_he_arr_h_in(:)
    integer, allocatable :: eep_data_h_in(:,:)       ! (max_neep, ntracks)
    real(8), allocatable :: tr_data_h_in(:,:)       ! (max_ncol, max_points)
    character(len=f2py_strlen), allocatable :: col_names_h_in(:,:) ! (max_ncol, ntracks)

    ! Similarly for helium tracks
    integer :: ntracks_he_in, max_ncol_he_in, max_neep_he_in, max_points_he_in
    character(len=f2py_strlen), allocatable :: filenames_he_in(:)
    real(8), allocatable :: initial_mass_he_in(:)
    real(8), allocatable :: initial_Y_he_in(:)
    real(8), allocatable :: initial_Z_he_in(:)
    real(8), allocatable :: Fe_div_H_he_in(:)
    real(8), allocatable :: alpha_div_Fe_he_in(:)
    real(8), allocatable :: v_div_vcrit_he_in(:)
    integer, allocatable :: ntrack_arr_he_in(:)
    integer, allocatable :: neep_arr_he_in(:)
    integer, allocatable :: ncol_arr_he_in(:)
    logical, allocatable :: is_he_arr_he_in(:)
    integer, allocatable :: eep_data_he_in(:,:)       
    real(8), allocatable :: tr_data_he_in(:,:)       
    character(len=f2py_strlen), allocatable :: col_names_he_in(:,:) 
    ! -----------------------
    ! File lists
    ! -----------------------
    character(len=f2py_strlen), allocatable :: py_track_list(:)
    character(len=f2py_strlen), allocatable :: py_track_list_he(:)
    character(len=f2py_strlen), allocatable :: py_metallicity_file_list(:)
    character(len=f2py_strlen), allocatable :: py_metallicity_file_list_he(:)

    ! -----------------------
    ! Metallicity values
    ! -----------------------
    real(8), allocatable :: py_Z_list(:)
    real(8), allocatable :: py_Z_list_he(:)

    ! -----------------------
    ! Hydrogen format controls
    ! -----------------------
    logical :: py_read_eep_files_H
    integer :: py_PreMS_EEP_H, py_ZAMS_EEP_H, py_IAMS_EEP_H, py_TAMS_EEP_H, py_BGB_EEP_H
    integer :: py_cHeIgnition_EEP_H, py_cHeBurn_EEP_H, py_TA_cHeB_EEP_H, py_TPAGB_EEP_H
    integer :: py_cCBurn_EEP_H, py_post_AGB_EEP_H, py_Initial_EEP_H, py_Final_EEP_H
    integer :: py_low_mass_final_eep_H, py_high_mass_final_eep_H
    logical :: py_fix_track_H
    character(len=f2py_strlen) :: py_age_colname_H, py_mass_colname_H, py_log_L_colname_H
    character(len=f2py_strlen) :: py_log_T_colname_H, py_log_R_colname_H
    character(len=f2py_strlen) :: py_he_core_mass_H, py_co_core_mass_H
    character(len=f2py_strlen) :: py_he_core_radius_H, py_co_core_radius_H
    character(len=f2py_strlen) :: py_mass_conv_envelope_H, py_radius_conv_envelope_H
    character(len=f2py_strlen) :: py_log_Tc_H, py_He4_mass_frac_H, py_c12_mass_frac_H, py_o16_mass_frac_H

    ! -----------------------
    ! Helium format controls
    ! -----------------------
    logical :: py_read_eep_files_He
    integer :: py_PreMS_EEP_He, py_ZAMS_EEP_He, py_IAMS_EEP_He, py_TAMS_EEP_He, py_BGB_EEP_He
    integer :: py_cHeIgnition_EEP_He, py_cHeBurn_EEP_He, py_TA_cHeB_EEP_He, py_TPAGB_EEP_He
    integer :: py_cCBurn_EEP_He, py_post_AGB_EEP_He, py_Initial_EEP_He, py_Final_EEP_He
    integer :: py_low_mass_final_eep_He, py_high_mass_final_eep_He
    logical :: py_fix_track_He
    character(len=f2py_strlen) :: py_age_colname_He, py_mass_colname_He, py_log_L_colname_He
    character(len=f2py_strlen) :: py_log_T_colname_He, py_log_R_colname_He
    character(len=f2py_strlen) :: py_he_core_mass_He, py_co_core_mass_He
    character(len=f2py_strlen) :: py_he_core_radius_He, py_co_core_radius_He
    character(len=f2py_strlen) :: py_mass_conv_envelope_He, py_radius_conv_envelope_He
    character(len=f2py_strlen) :: py_log_Tc_He, py_He4_mass_frac_He, py_c12_mass_frac_He, py_o16_mass_frac_He

contains

    ! --------------------------
    ! Python interface: set format controls for hydrogen
    ! --------------------------
    subroutine set_format_controls_H(read_eep, prems, zams, iams, tams, bgb, &
                                     cHeIgn, cHeBurn, ta_cHeB, tpagb, cCBurn, postAGB, &
                                     initEEP, finalEEP, fixtrack, lowEEP, highEEP, &
                                     age_col, mass_col, logL_col, logT_col, logR_col, &
                                     he_mass_col, co_mass_col, he_radius_col, co_radius_col, &
                                     mass_env_col, radius_env_col, logTc_col, He4_col, c12_col, o16_col)
                                     
        logical, intent(in) :: read_eep, fixtrack
        integer, intent(in) :: prems, zams, iams, tams, bgb
        integer, intent(in) :: cHeIgn, cHeBurn, ta_cHeB, tpagb
        integer, intent(in) :: cCBurn, postAGB, initEEP, finalEEP
        integer, intent(in) :: lowEEP, highEEP
        character(len=*), intent(in) :: age_col, mass_col, logL_col, logT_col, logR_col
        character(len=*), intent(in) :: he_mass_col, co_mass_col, he_radius_col, co_radius_col
        character(len=*), intent(in) :: mass_env_col, radius_env_col
        character(len=*), intent(in) :: logTc_col, He4_col, c12_col, o16_col

        py_read_eep_files_H = read_eep
        py_PreMS_EEP_H = prems
        py_ZAMS_EEP_H = zams
        py_IAMS_EEP_H = iams
        py_TAMS_EEP_H = tams
        py_BGB_EEP_H = bgb
        py_cHeIgnition_EEP_H = cHeIgn
        py_cHeBurn_EEP_H = cHeBurn
        py_TA_cHeB_EEP_H = ta_cHeB
        py_TPAGB_EEP_H = tpagb
        py_cCBurn_EEP_H = cCBurn
        py_post_AGB_EEP_H = postAGB
        py_Initial_EEP_H = initEEP
        py_Final_EEP_H = finalEEP
        py_fix_track_H = fixtrack
        py_low_mass_final_eep_H = lowEEP
        py_high_mass_final_eep_H = highEEP

        py_age_colname_H = age_col
        py_mass_colname_H = mass_col
        py_log_L_colname_H = logL_col
        py_log_T_colname_H = logT_col
        py_log_R_colname_H = logR_col
        py_he_core_mass_H = he_mass_col
        py_co_core_mass_H = co_mass_col
        py_he_core_radius_H = he_radius_col
        py_co_core_radius_H = co_radius_col
        py_mass_conv_envelope_H = mass_env_col
        py_radius_conv_envelope_H = radius_env_col
        py_log_Tc_H = logTc_col
        py_He4_mass_frac_H = He4_col
        py_c12_mass_frac_H = c12_col
        py_o16_mass_frac_H = o16_col
    end subroutine set_format_controls_H

    ! --------------------------
    ! Python interface: set format controls for helium
    ! --------------------------
    subroutine set_format_controls_He(read_eep, bgb, cHeBurn, ta_cHeB, tpagb, cCBurn, postAGB, &
                                      initEEP, finalEEP, fixtrack, lowEEP, highEEP, &
                                      age_col, mass_col, logL_col, logT_col, logR_col, &
                                      he_mass_col, co_mass_col, he_radius_col, co_radius_col, &
                                      mass_env_col, radius_env_col, logTc_col, He4_col, c12_col, o16_col)
        logical, intent(in) :: read_eep, fixtrack
        integer, intent(in) :: bgb, cHeBurn, ta_cHeB, tpagb
        integer, intent(in) :: cCBurn, postAGB, initEEP, finalEEP
        integer, intent(in) :: lowEEP, highEEP
        character(len=*), intent(in) :: age_col, mass_col, logL_col, logT_col, logR_col
        character(len=*), intent(in) :: he_mass_col, co_mass_col, he_radius_col, co_radius_col
        character(len=*), intent(in) :: mass_env_col, radius_env_col
        character(len=*), intent(in) :: logTc_col, He4_col, c12_col, o16_col

        py_read_eep_files_He = read_eep
        py_BGB_EEP_He = bgb
        py_cHeBurn_EEP_He = cHeBurn
        py_TA_cHeB_EEP_He = ta_cHeB
        py_TPAGB_EEP_He = tpagb
        py_cCBurn_EEP_He = cCBurn
        py_post_AGB_EEP_He = postAGB
        py_Initial_EEP_He = initEEP
        py_Final_EEP_He = finalEEP
        py_fix_track_He = fixtrack
        py_low_mass_final_eep_He = lowEEP
        py_high_mass_final_eep_He = highEEP

        py_age_colname_He = age_col
        py_mass_colname_He = mass_col
        py_log_L_colname_He = logL_col
        py_log_T_colname_He = logT_col
        py_log_R_colname_He = logR_col
        py_he_core_mass_He = he_mass_col
        py_co_core_mass_He = co_mass_col
        py_he_core_radius_He = he_radius_col
        py_co_core_radius_He = co_radius_col
        py_mass_conv_envelope_He = mass_env_col
        py_radius_conv_envelope_He = radius_env_col
        py_log_Tc_He = logTc_col
        py_He4_mass_frac_He = He4_col
        py_c12_mass_frac_He = c12_col
        py_o16_mass_frac_He = o16_col
    end subroutine set_format_controls_He

    ! --------------------------
    ! File/metallicity routines 
    ! --------------------------
    subroutine set_file_lists(nmet, met_files, nmet_he, met_he_files, &
                              nh_track, h_tracks, nhe_track, he_tracks)
        integer, intent(in) :: nmet, nmet_he, nh_track, nhe_track
        character(len=*), intent(in) :: met_files(nmet)
        character(len=*), intent(in) :: met_he_files(nmet_he)
        character(len=*), intent(in) :: h_tracks(nh_track)
        character(len=*), intent(in) :: he_tracks(nhe_track)

        if (allocated(py_metallicity_file_list)) deallocate(py_metallicity_file_list)
        allocate(py_metallicity_file_list(nmet))
        py_metallicity_file_list = met_files

        if (allocated(py_metallicity_file_list_he)) deallocate(py_metallicity_file_list_he)
        allocate(py_metallicity_file_list_he(nmet_he))
        py_metallicity_file_list_he = met_he_files

        if (allocated(py_track_list)) deallocate(py_track_list)
        allocate(py_track_list(nh_track))
        py_track_list = h_tracks

        if (allocated(py_track_list_he)) deallocate(py_track_list_he)
        allocate(py_track_list_he(nhe_track))
        py_track_list_he = he_tracks
    end subroutine set_file_lists

    subroutine set_mets(nmet, Z_values, nmet_he, Z_values_he)
        integer, intent(in) :: nmet, nmet_he
        real(8), intent(in) :: Z_values(nmet), Z_values_he(nmet_he)

        if (allocated(py_Z_list)) deallocate(py_Z_list)
        allocate(py_Z_list(nmet))
        py_Z_list = Z_values

        if (allocated(py_Z_list_he)) deallocate(py_Z_list_he)
        allocate(py_Z_list_he(nmet_he))
        py_Z_list_he = Z_values_he
    end subroutine set_mets

    subroutine set_tracks_from_python(ntracks, max_ncol, max_neep, max_points, filenames, &
                                      initial_mass, initial_Y, initial_Z, &
                                      Fe_div_H, alpha_div_Fe, v_div_vcrit, ntrack_arr, &
                                      neep_arr, ncol_arr, is_he_arr, eep_data, tr_data, col_names, is_he)
        implicit none
        integer, intent(in) :: ntracks, max_ncol, max_neep, max_points
        character(len=*), intent(in) :: filenames(ntracks)
        real(8), intent(in) :: initial_mass(ntracks), initial_Y(ntracks), initial_Z(ntracks)
        real(8), intent(in) :: Fe_div_H(ntracks), alpha_div_Fe(ntracks), v_div_vcrit(ntracks)
        integer, intent(in) :: ntrack_arr(ntracks), neep_arr(ntracks), ncol_arr(ntracks)
        logical, intent(in) :: is_he_arr(ntracks)
        integer, intent(in) :: eep_data(max_neep, ntracks)
        real(8), intent(in) :: tr_data(max_ncol, max_points)
        character(len=*), intent(in) :: col_names(max_ncol, ntracks)
        logical, intent(in) :: is_he
    
        integer :: i
    
        if (is_he) then
            ! Store in He arrays
            ntracks_he_in = ntracks
            max_ncol_he_in = max_ncol
            max_neep_he_in = max_neep
            max_points_he_in = max_points
            is_he_in = is_he
    
            if (allocated(filenames_he_in)) deallocate(filenames_he_in)
            allocate(filenames_he_in(ntracks)); filenames_he_in = filenames
            if (allocated(initial_mass_he_in)) deallocate(initial_mass_he_in)
            allocate(initial_mass_he_in(ntracks)); initial_mass_he_in = initial_mass
            if (allocated(initial_Y_he_in)) deallocate(initial_Y_he_in)
            allocate(initial_Y_he_in(ntracks)); initial_Y_he_in = initial_Y
            if (allocated(initial_Z_he_in)) deallocate(initial_Z_he_in)
            allocate(initial_Z_he_in(ntracks)); initial_Z_he_in = initial_Z
            if (allocated(Fe_div_H_he_in)) deallocate(Fe_div_H_he_in)
            allocate(Fe_div_H_he_in(ntracks)); Fe_div_H_he_in = Fe_div_H
            if (allocated(alpha_div_Fe_he_in)) deallocate(alpha_div_Fe_he_in)
            allocate(alpha_div_Fe_he_in(ntracks)); alpha_div_Fe_he_in = alpha_div_Fe
            if (allocated(v_div_vcrit_he_in)) deallocate(v_div_vcrit_he_in)
            allocate(v_div_vcrit_he_in(ntracks)); v_div_vcrit_he_in = v_div_vcrit
            if (allocated(ntrack_arr_he_in)) deallocate(ntrack_arr_he_in)
            allocate(ntrack_arr_he_in(ntracks)); ntrack_arr_he_in = ntrack_arr
            if (allocated(neep_arr_he_in)) deallocate(neep_arr_he_in)
            allocate(neep_arr_he_in(ntracks)); neep_arr_he_in = neep_arr
            if (allocated(ncol_arr_he_in)) deallocate(ncol_arr_he_in)
            allocate(ncol_arr_he_in(ntracks)); ncol_arr_he_in = ncol_arr
            if (allocated(is_he_arr_he_in)) deallocate(is_he_arr_he_in)
            allocate(is_he_arr_he_in(ntracks)); is_he_arr_he_in = is_he_arr
            if (allocated(eep_data_he_in)) deallocate(eep_data_he_in)
            allocate(eep_data_he_in(max_neep, ntracks)); eep_data_he_in = eep_data
            if (allocated(tr_data_he_in)) deallocate(tr_data_he_in)
            allocate(tr_data_he_in(max_ncol, max_points)); tr_data_he_in = tr_data
            if (allocated(col_names_he_in)) deallocate(col_names_he_in)
            allocate(col_names_he_in(max_ncol, ntracks)); col_names_he_in = col_names
    
        else
            ! Store in H arrays
            ntracks_h_in = ntracks
            max_ncol_h_in = max_ncol
            max_neep_h_in = max_neep
            max_points_h_in = max_points
            is_he_in = is_he
    
            if (allocated(filenames_h_in)) deallocate(filenames_h_in)
            allocate(filenames_h_in(ntracks)); filenames_h_in = filenames
            if (allocated(initial_mass_h_in)) deallocate(initial_mass_h_in)
            allocate(initial_mass_h_in(ntracks)); initial_mass_h_in = initial_mass
            if (allocated(initial_Y_h_in)) deallocate(initial_Y_h_in)
            allocate(initial_Y_h_in(ntracks)); initial_Y_h_in = initial_Y
            if (allocated(initial_Z_h_in)) deallocate(initial_Z_h_in)
            allocate(initial_Z_h_in(ntracks)); initial_Z_h_in = initial_Z
            if (allocated(Fe_div_H_h_in)) deallocate(Fe_div_H_h_in)
            allocate(Fe_div_H_h_in(ntracks)); Fe_div_H_h_in = Fe_div_H
            if (allocated(alpha_div_Fe_h_in)) deallocate(alpha_div_Fe_h_in)
            allocate(alpha_div_Fe_h_in(ntracks)); alpha_div_Fe_h_in = alpha_div_Fe
            if (allocated(v_div_vcrit_h_in)) deallocate(v_div_vcrit_h_in)
            allocate(v_div_vcrit_h_in(ntracks)); v_div_vcrit_h_in = v_div_vcrit
            if (allocated(ntrack_arr_h_in)) deallocate(ntrack_arr_h_in)
            allocate(ntrack_arr_h_in(ntracks)); ntrack_arr_h_in = ntrack_arr
            if (allocated(neep_arr_h_in)) deallocate(neep_arr_h_in)
            allocate(neep_arr_h_in(ntracks)); neep_arr_h_in = neep_arr
            if (allocated(ncol_arr_h_in)) deallocate(ncol_arr_h_in)
            allocate(ncol_arr_h_in(ntracks)); ncol_arr_h_in = ncol_arr
            if (allocated(is_he_arr_h_in)) deallocate(is_he_arr_h_in)
            allocate(is_he_arr_h_in(ntracks)); is_he_arr_h_in = is_he_arr
            if (allocated(eep_data_h_in)) deallocate(eep_data_h_in)
            allocate(eep_data_h_in(max_neep, ntracks)); eep_data_h_in = eep_data
            if (allocated(tr_data_h_in)) deallocate(tr_data_h_in)
            allocate(tr_data_h_in(max_ncol, max_points)); tr_data_h_in = tr_data
            if (allocated(col_names_h_in)) deallocate(col_names_h_in)
            allocate(col_names_h_in(max_ncol, ntracks)); col_names_h_in = col_names
    
        end if

end subroutine set_tracks_from_python



end module metisse_interface
