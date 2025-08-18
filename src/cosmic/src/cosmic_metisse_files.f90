module cosmic_metisse_files
    implicit none
    ! Dummy storage for Python-provided lists
    integer, parameter :: f2py_strlen = 256
    character(len=f2py_strlen), allocatable :: py_track_list(:)
    character(len=f2py_strlen), allocatable :: py_track_list_he(:)
    character(len=f2py_strlen), allocatable :: py_metallicity_file_list(:)
    character(len=f2py_strlen), allocatable :: py_metallicity_file_list_he(:)
contains
    subroutine set_file_lists_from_python(nmet, met_files, nmet_he, met_he_files, &
                                          nh_track, h_tracks, nhe_track, he_tracks)
        integer, intent(in) :: nmet, nmet_he, nh_track, nhe_track
        integer, parameter :: f2py_strlen = 256
        character(len=f2py_strlen), intent(in) :: met_files(1)
        character(len=f2py_strlen), intent(in) :: met_he_files(1)
        character(len=f2py_strlen), intent(in) :: h_tracks(100)
        character(len=f2py_strlen), intent(in) :: he_tracks(55)
        
        ! Save Python-provided lists in module variables
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
    end subroutine set_file_lists_from_python
end module cosmic_metisse_files
