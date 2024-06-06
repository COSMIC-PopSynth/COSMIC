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
      

