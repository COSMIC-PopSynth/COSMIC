subroutine zcnsts_METISSE(z,zpars)
    use track_support
    use z_support
!    use remnant_support

    real(dp) ,intent(in) :: z
    real(dp) ,intent(out) :: zpars(20)

    integer :: i,ierr
    logical :: debug
    character(len=strlen):: path

    ierr = 0
    debug = .false.
    
    call read_metisse_input()
    initial_Z = z
    
     if (.not. defined(initial_Z))then
        print*,"Error: Define initial_Z "
        STOP
    endif
    
    if (trim(INPUT_FILES_DIR) == '' )then
        print*,"Error: Define INPUT_FILES_DIR "
        STOP
    endif
    
    if (read_files_from_Z) then         !get folder name if read_files_from_Z
        if (Z_folder_list == '') then
            print*,"Error: Z_folder_list not defined for read_files_from_Z"
            STOP
        else
            call get_folder_from_Z(INPUT_FILES_DIR,initial_Z,path)
        endif
    else
        path = INPUT_FILES_DIR
    endif
    
    !reading format file
    call read_format(format_file)

    !getting filenames
    call get_files_from_path(path)

    if (verbose) print*,"Number of input tracks", num_tracks

    !determine key columns 
    if (key_columns_file /= '') call process_columns(key_columns_file,key_cols,ierr)

    !get column numbers for core related quantities
!    call get_core_columns()

    if (.not. read_eep_files .and. header_location<=0)then
            call process_columns(column_name_file,temp_cols,ierr)
            if(ierr/=0) then
                print*,"Check if header location and column_name_file are correct "
                STOP
            endif

            if (size(temp_cols) /= total_cols) then
                print*,'Erorr reading number of columns'
                print*,'Check if column_name_file and total_cols are correct'
                STOP
            endif
        if (debug) print*,"read column names from file"
    end if

    call read_key_eeps()
    if (debug) print*, "key eeps", key_eeps

        if (read_eep_files) then
            if (debug) print*,"reading eep files"
            do i=1,num_tracks
                call read_eep(s(i))
                if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% ncol
            end do
        else
            do i=1,num_tracks
                call read_input_file(s(i))
                if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% ncol
            end do

        endif

        ! print*, 'key cols size',size(key_cols)
        if (size(key_cols)<=1) then
            allocate(key_cols(s(1)% ncol))
            key_cols% name = s(1)% cols% name
        end if

    if(debug) print*, s(1)% cols% name, s(1)% tr(:,1)
    if(debug) print*,"using key columns", key_cols(1) 

    !locates key columns of mass, age etc.
    call locate_column_numbers(s,key_cols)
    do i = 1,size(s)
        call set_star_type_from_history(s(i))
        s(i)% has_mass_loss = check_mass_loss(s(i))
    end do

    !TODO: check for monotonicity of initial masses
    if (debug) print*,s% initial_mass

    !sets z parameters and cutoff masses
    call set_zparameters(zpars)
    
    !TODO: make user definable?
    zpars(11) = 0.76d0 - 3.d0*z
    zpars(12) = 0.24d0 + 2.d0*z
    Z04 = initial_Z**0.4
    zpars(14) = initial_Z**0.4

    ! if (direct_call) then
    !     call set_remnant_scheme()
    ! else
    !     call assign_commons()
    ! endif

end subroutine zcnsts_METISSE

