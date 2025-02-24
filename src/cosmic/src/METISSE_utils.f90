    subroutine assign_commons()
        use track_support
        implicit none
        
        !to assign common variables when METISSE is used with COSMIC
          
        REAL(dp) :: ecsn,ecsn_mlow
        COMMON /SNVARS1/ ecsn,ecsn_mlow
         
        real(dp) :: d
        
        if(front_end == COSMIC) then
        ! use inputs from COSMIC
        
            if (Mec_core > 0.d0) ecsn = Mec_core
            d = (Mec_core-Mup_core)
            if (Mup_core > 0.d0 .and. d>tiny ) ecsn_mlow = Mup_core
            
        else
            print*,'Error: Front end mismtach in assign commons'
            print*,'expected 2 (COSMIC); got ', front_end
        endif

    end subroutine

    subroutine get_bhspin(bhspin,id)
        use track_support, only: tarr,dp
        implicit none
        integer, intent(in) :: id
        real(dp), intent(out) :: bhspin

        bhspin = tarr(id)% pars% bhspin
    end subroutine
    
    subroutine check_error(err)
        use track_support, only: code_error
        integer, intent(out) :: err
        err = 0
        if(code_error) err = 1
    end subroutine
    
    
    subroutine assign_error()
        use track_support, only: code_error
        code_error = .true.
    end subroutine
      
    subroutine initialize_metisse_front_cmc()
    ! passing strings with c/cmc is not very realiable
    ! we set front end like this avoid possible seg faults
        call initialize_front_end('cosmic')
    end subroutine

    subroutine get_COSMIC_input()
        use track_support
        use z_support, only: Z_accuracy_limit, get_csafe_string
        
        ! takes inputs from cosmic and assigns them
        ! to appropiate variables in METISSE
        
        character(len=strlen) :: path_to_tracks, path_to_he_tracks
        real(dp) :: z_match_limit
        LOGICAL METISSE_verbose
        COMMON/ METISSEVARS/ path_to_tracks,path_to_he_tracks,&
                     z_match_limit, METISSE_verbose
        
        ! remove the null charcater if any
        call get_csafe_string(path_to_tracks,METALLICITY_DIR)
        call get_csafe_string(path_to_he_tracks,METALLICITY_DIR_HE)
        Z_accuracy_limit = z_match_limit
        verbose = METISSE_verbose
    
    end subroutine

    logical function check_path_change() result (load_tracks)
        use track_support, only: strlen,METALLICITY_DIR, METALLICITY_DIR_HE
        use z_support, only: get_csafe_string

        character(len=strlen) :: path_to_tracks, path_to_he_tracks
        COMMON/ METISSEVARS/ path_to_tracks,path_to_he_tracks
            
        INTEGER :: using_cmc
        COMMON /CMCPASS/ using_cmc
        
        character(len=strlen) :: string1,string2
        load_tracks = .false.

        ! remove the null charcater if any
        call get_csafe_string(path_to_tracks,string1)
        call get_csafe_string(path_to_he_tracks, string2)

        if((trim(path_to_tracks)/=trim(METALLICITY_DIR)) .or. &
            (trim(path_to_he_tracks)/=trim(METALLICITY_DIR_HE))) load_tracks = .true.
    end function

    
