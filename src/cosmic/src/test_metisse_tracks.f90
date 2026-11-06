! Debugging helper for test_bse.f (not part of the python build)
!
! When COSMIC uses METISSE, python reads the format/eep files and passes them to
! the fortran via c_m_interface (see read_tracks_for_METISSE in evolve.py).
! debug/create_binary_in.py records exactly what python would pass and writes it
! to a stream file, which this subroutine replays through the same c_m_interface
! calls so the fortran ends up in the same state as when run through python.
!
! File layout (native endian, see write_metisse_tracks in debug/create_binary_in.py):
!   char(8) 'METTRK01', int32 has_he
!   then a block for the hydrogen tracks and (if has_he) one for the helium tracks:
!     int32 n_int, int32 format_ints(n_int)    (logicals stored as 0/1)
!     char(256) format_strings(15)
!     int32 ntracks, max_ncol, max_neep, max_points
!     char(256) filenames(ntracks)
!     real64 initial_mass, initial_Y, initial_Z, Fe_div_H, alpha_div_Fe, v_div_vcrit (each ntracks)
!     int32 ntrack_arr, neep_arr, ncol_arr (each ntracks)
!     int32 eep_data(max_neep, ntracks)
!     real64 tr_data(max_ncol, max_points)
!     char(32) col_names(max_ncol, ntracks)

subroutine load_metisse_tracks(filename)
    use c_m_interface
    implicit none
    character(len=*), intent(in) :: filename

    integer, parameter :: n_str = 15, col_width = 32
    character(len=8) :: tag
    integer :: u, ios, has_he

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print*, 'load_metisse_tracks: cannot open ', trim(filename)
        print*, 'run debug/create_binary_in.py on a METISSE initC first'
        stop 1
    endif

    read(u) tag
    if (tag /= 'METTRK01') then
        print*, 'load_metisse_tracks: ', trim(filename), ' is not a track file (or is an old format)'
        stop 1
    endif
    read(u) has_he

    call read_block(.false.)
    if (has_he /= 0) call read_block(.true.)

    close(u)

contains

    subroutine read_block(is_he)
        logical, intent(in) :: is_he
        integer :: n_int, ntracks, max_ncol, max_neep, max_points
        integer, allocatable :: ints(:)
        character(len=f2py_strlen) :: strs(n_str)
        character(len=f2py_strlen), allocatable :: filenames(:)
        real(dp), allocatable :: initial_mass(:), initial_Y(:), initial_Z(:)
        real(dp), allocatable :: Fe_div_H(:), alpha_div_Fe(:), v_div_vcrit(:)
        integer, allocatable :: ntrack_arr(:), neep_arr(:), ncol_arr(:)
        integer, allocatable :: eep_data(:,:)
        real(dp), allocatable :: tr_data(:,:)
        character(len=col_width), allocatable :: col_names(:,:)

        ! format controls, in the same order as the set_format_controls_* arguments
        read(u) n_int
        allocate(ints(n_int))
        read(u) ints
        read(u) strs

        if (is_he) then
            if (n_int /= 12) stop 'load_metisse_tracks: unexpected number of He format controls'
            call set_format_controls_He(ints(1) /= 0, ints(2), ints(3), ints(4), ints(5), &
                                        ints(6), ints(7), ints(8), ints(9), ints(10) /= 0, &
                                        ints(11), ints(12), &
                                        strs(1), strs(2), strs(3), strs(4), strs(5), &
                                        strs(6), strs(7), strs(8), strs(9), strs(10), &
                                        strs(11), strs(12), strs(13), strs(14), strs(15))
        else
            if (n_int /= 17) stop 'load_metisse_tracks: unexpected number of H format controls'
            call set_format_controls_H(ints(1) /= 0, ints(2), ints(3), ints(4), ints(5), &
                                       ints(6), ints(7), ints(8), ints(9), ints(10), &
                                       ints(11), ints(12), ints(13), ints(14), ints(15) /= 0, &
                                       ints(16), ints(17), &
                                       strs(1), strs(2), strs(3), strs(4), strs(5), &
                                       strs(6), strs(7), strs(8), strs(9), strs(10), &
                                       strs(11), strs(12), strs(13), strs(14), strs(15))
        endif

        ! the tracks themselves
        read(u) ntracks, max_ncol, max_neep, max_points
        allocate(filenames(ntracks))
        allocate(initial_mass(ntracks), initial_Y(ntracks), initial_Z(ntracks))
        allocate(Fe_div_H(ntracks), alpha_div_Fe(ntracks), v_div_vcrit(ntracks))
        allocate(ntrack_arr(ntracks), neep_arr(ntracks), ncol_arr(ntracks))
        allocate(eep_data(max_neep, ntracks), tr_data(max_ncol, max_points))
        allocate(col_names(max_ncol, ntracks))

        read(u) filenames
        read(u) initial_mass, initial_Y, initial_Z, Fe_div_H, alpha_div_Fe, v_div_vcrit
        read(u) ntrack_arr, neep_arr, ncol_arr
        read(u) eep_data
        read(u) tr_data
        read(u) col_names

        call set_tracks_from_python(ntracks, max_ncol, max_neep, max_points, &
                                    filenames, initial_mass, initial_Y, initial_Z, &
                                    Fe_div_H, alpha_div_Fe, v_div_vcrit, ntrack_arr, &
                                    neep_arr, ncol_arr, eep_data, tr_data, col_names, is_he)

        if (is_he) then
            write(*,'(a,i5,a)') ' Loaded ', ntracks, ' helium METISSE tracks'
        else
            write(*,'(a,i5,a)') ' Loaded ', ntracks, ' hydrogen METISSE tracks'
        endif
    end subroutine read_block

end subroutine load_metisse_tracks
