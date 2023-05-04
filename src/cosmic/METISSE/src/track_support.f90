module track_support

    !some subroutines of module have been adpated from iso_eep_support module of ISO package (Dotter 2016).

    implicit none

    integer, parameter :: min_io_unit = 29
    integer, parameter :: max_io_unit = 99
    logical :: assigned(max_io_unit) = .false.

    !----from mesa const_def.f90
    ! real number precision options: single, double
    integer, parameter :: sp = selected_real_kind(p=5)
    integer, parameter :: dp = selected_real_kind(p=15)

    integer, parameter :: strlen = 256 ! for character (len=strlen)

    logical :: verbose
    logical :: write_track_to_file, write_eep_file
    logical :: direct_call = .false.


    character(len=strlen) :: METISSE_DIR, eep_dir
    real(dp) :: pts_1,pts_2,pts_3
    integer :: low_mass_final_eep, high_mass_final_eep
    integer, allocatable :: key_eeps(:)

    logical, parameter :: old_core_mass_names=.false.
    integer, parameter :: col_width = 32

    real(dp), parameter :: ln10=log(1.0d1)
    real(sp), parameter :: ln10_sp = log(10.0)
    real(dp), parameter :: tiny = 1.0d-6
    real(dp), parameter :: undefined  =  -1.0
    integer, parameter :: undefined_i = -1


    ! for use when constructing EEP distance
    logical :: weight_center_rho_T_by_Xc
    real(dp) :: Teff_scale=2d0
    real(dp) :: logL_scale=0.125d0
    real(dp) :: age_scale=0.05d0
    real(dp) :: Rhoc_scale=1d0
    real(dp) :: Tc_scale=1d0

    !stellar types for handling primary eeps
    integer, parameter :: unknown           =  1 !for initialization only
    integer, parameter :: sub_stellar       =  2 !no fusion = brown dwarf
    integer, parameter :: star_low_mass     =  3 !ends as a WD
    integer, parameter :: star_high_mass    =  4 !does not end as a WD

    character(len=10) :: star_label(4) = ['   unknown', 'substellar', '  low-mass', ' high-mass']
    character(len=5) :: phase_label(16) = ['lm_MS','   MS','   HG','  FGB',' CHeB',' EAGB',&
    'TPAGB','He_MS', 'He_HG','He_GB','He_WD','CO_WD','ONeWD','   NS','   BH','   MR']

    ! default column format specs
    integer :: head !=29
    integer :: main !=28
    integer :: xtra !=0

    !sse phases

    integer, parameter :: low_mass_MS = 0
    integer, parameter :: MS = 1
    integer, parameter :: HG = 2
    integer, parameter :: RGB = 3
    integer, parameter :: HeBurn = 4
    integer, parameter :: EAGB = 5
    integer, parameter :: TPAGB = 6
    integer, parameter :: He_MS = 7
    integer, parameter :: He_HG = 8
    integer, parameter :: He_GB = 9
    integer, parameter :: HeWD = 10
    integer, parameter :: CO_WD = 11
    integer, parameter :: ONeWD = 12
    integer, parameter :: NS = 13
    integer, parameter :: BH = 14
    integer, parameter :: Massless_REM = 15

    !EEPs

    integer :: PreMS_EEP = -1
    integer :: ZAMS_EEP = -1
    integer :: IAMS_EEP = -1
    integer :: TAMS_EEP = -1
    integer :: BGB_EEP = -1
    integer :: cHeIgnition_EEP = -1
    integer :: cHeBurn_EEP = -1
    integer :: TA_cHeB_EEP = -1

    integer :: cCBurn_EEP = -1
    integer :: TPAGB_EEP = -1
    integer :: post_AGB_EEP  = -1
    integer :: WD_EEP  = -1

    integer :: Initial_EEP = -1       !files will be read from this line number
    integer :: Final_EEP = -1        !to this line
    integer :: Extra_EEP1 = -1
    integer :: Extra_EEP2= -1
    integer :: Extra_EEP3 = -1


    ! min quantities from history file that need to be identified

    character(len=col_width) :: age_colname, mass_colname, log_L_colname,log_T_colname, &
                                log_R_colname, log_mdot_colname,he_core_mass,c_core_mass, &
                                log_Tc,c12_mass_frac,o16_mass_frac, he4_mass_frac, &
                                Lum_colname,Teff_colname,Radius_colname, mdot_colname, &
                                he_core_radius

    integer :: i_age, i_mass, i_logLH, i_logLHe, i_logTe, i_logL, i_logR
    integer :: i_logg, i_Tc, i_Rhoc, i_Xc, i_Yc, i_he_core, i_co_core
    integer :: i_Cc, i_gamma, i_surfH, i_c12,i_o16,i_he4, i_lum, i_rad, i_mdot
    integer :: number_of_core_columns, i_age2, i_core_radius

    integer, allocatable :: core_cols(:)!, surface_cols(:)
    !for columns
    integer, parameter :: max_col = 180
    integer, parameter :: column_int=0
    integer, parameter :: column_dbl=1
    character(len=strlen) :: core_columns        !TODO: make it flexible

    type column
     character(len=col_width) :: name
     integer :: type, loc
    end type column

    !EEP arrays
    integer, parameter :: primary = 10 ! number of primary EEPs !TODO: --change this
    ! as set by primary_eep
    integer :: eep_interval(primary-1) ! number of secondary EEPs
    ! between the primaries

    real(dp), allocatable :: t_incomplete(:), t_notfound(:)

  !holds an evolutionary track for input, use an array of these for multiple tracks

    type eep_track
        character(len=strlen) :: filename, cmd_suffix
        character(len=8) :: version_string
        type(column), allocatable :: cols(:)

        logical :: has_phase = .false., ignore=.false.
        logical :: has_mass_loss
        integer :: ncol, ntrack, neep, MESA_revision_number
        integer :: star_type = unknown

        integer, allocatable :: eep(:), phase(:)
        real(dp) :: initial_mass, initial_Y, Fe_div_H, initial_Z, v_div_vcrit, alpha_div_Fe
        real(dp), allocatable :: tr(:,:)

    end type eep_track

    !holds current parameters of star-- used by track
    type star_parameters
        integer :: phase,extra
        real(dp) :: mass,core_mass,McHe, McCO,luminosity,Teff,radius
        real(dp) :: log_L,log_Teff,log_R                !log values
        real(dp) :: epoch, age, age_old
        real(dp) :: delta, dt, dms

    end type star_parameters
    

    !holds values of agb parameters for constructing AGB to WD track
    type agb_parameters
        real(dp) :: age,radius,lum
        real(dp) :: t1,t2,t_post_agb
        integer :: phase_wd
    end type agb_parameters

    !holds values of agb parameters for constructing AGB to WD track
    type sse_parameters
        real(dp) :: D,Mx,Lx,LtMS
        real(dp) :: Rzams, Lzams      !zams values
    end type

    !holds interpolated track
    type track
        type(column), allocatable :: cols(:)
        logical :: has_RGB=.false., complete=.true.
        logical :: has_mass_loss
        integer :: ncol, ntrack, neep
        integer :: star_type = unknown,irecord
        integer, allocatable :: eep(:), phase(:)
        real(dp) :: initial_mass, initial_Z , initial_Y, Fe_div_H,  v_div_vcrit, alpha_div_Fe
        real(dp), allocatable :: tr(:,:)

        real(dp), allocatable :: times(:), times_new(:)           !timescales
        logical :: lost_envelope = .false., post_agb = .false.
        real(dp) :: zams_mass!, zams_radius, zams_lum      !zams values
        real(dp) :: MS_time, nuc_time
        type(star_parameters) :: pars!, old_pars    ! parameters at any instant
        type(agb_parameters) :: agb
        type(sse_parameters) :: He_pars
    end type track
    
    !defining array for input tracks
    type(eep_track), allocatable,target :: s(:)
    type(track), allocatable,target :: tarr(:)
    integer :: num_tracks
    real(dp) :: initial_Z

!variable declaration-- for main
    integer :: number_of_tracks
    character(len=strlen) :: input_mass_file
    logical :: read_mass_from_file
    
    !for z_support
    real(dp) :: Mhook, Mhef,Mfgb, Mup, Mec, Mextra,Mup_core,Mec_core
    real(dp) :: Z04

    !for interp_support
    logical:: fix_track
    real(dp) :: lookup_index, accuracy_limit
    
    !for remnant support
    real(dp) :: max_NS_mass         !maximum NS mass
    logical:: construct_wd_track, allow_electron_capture, use_Initial_final_mass_relation
    character (len = col_width) :: BHNS_mass_scheme, WD_mass_scheme
!    real(dp) :: mc1, mc2 !mass cutoffs for Belczynski methods
    
            
    contains
    

    subroutine get_track_ptr(id,t,ierr)
    integer :: ierr
    type(track), pointer :: t
    integer, optional :: id


    t => NULL()
    if(present(id))then
        t => tarr(id)
    else
        t => tarr(1)
        id = 1
    endif
    end subroutine

    subroutine deallocate_arrays(t)
    type(track), pointer :: t
        deallocate(t% eep)
        deallocate(t% cols)
        deallocate(t% tr)
        deallocate(t% phase)
        deallocate(t% times)
    end subroutine deallocate_arrays
    
    !from COMPAS (Team Compas 2020)
    real(dp) function quadratic(a,b,c) result(x)
      real(dp),intent(in) :: a,b,c
      real(dp) :: D,sqrtD,x1,x2

        D = (B*B)-(4.d0*A*C)
        x = 0.0

        if (D< 0.0) then
            print*,"fatal error: Non-real roots"
            STOP
        else if (D > 0.0) then
            sqrtD = sqrt(D)
            x1 = (-B + sqrtD)/(2*A)
            x2 = (-B - sqrtD)/(2*A)
            x= max(x1, x2)
        else
            x = -B/(2*A)
        endif
    end function quadratic

    !from mesa/utils_lib.f
    integer function alloc_iounit(ierr)
        !use utils_def
        integer, intent(out) :: ierr
        integer :: i
        ierr = 0
        alloc_iounit = -1
        do i = min_io_unit, max_io_unit
            if (.not. assigned(i)) then
                assigned(i) = .true.
                alloc_iounit = i
                exit
            end if
        end do
        if (alloc_iounit == -1) then
            ierr = -1
        end if
    end function alloc_iounit

    subroutine free_iounit(iounit)
        !use utils_def
        integer, intent(in) :: iounit
        logical :: bad_iounit
        bad_iounit = .false.
        if (iounit >= min_io_unit .and. iounit <= max_io_unit) then
            assigned(iounit) = .false.
        else
            bad_iounit = .true.
        end if
        if (bad_iounit) then
            write(*,*) 'called free_iounit with invalid arg', iounit
            stop 'free_iounit'
        end if
    end subroutine free_iounit

    ! from ISO (Dotter et al. 2016), adapted from MESA 

    ! if vec contains decreasing values,
    ! returns 1 if val > vec(1); returns n if val <= vec(n)
    ! else returns k between 1 and n-1 such that
    !     vec(k) >= val > vec(k+1)

    ! if vec contains increasing values,
    ! returns 0 if val < vec(1); returns n if val >= vec(n)
    ! else returns k between 1 and n-1 such that
    !     vec(k) <= val < vec(k+1)

    integer function binary_search(n, vec, val) result(loc)
     integer, intent(in) :: n
     real(dp), intent(in) :: val
     real(dp), intent(in) :: vec(:) ! (n)
!     real(dp), parameter :: tiny = 1.0d-13
     integer :: first, last, mid

     if (n <= 1) then
        loc = n; return
     end if

     if (vec(n) < vec(1)) then ! decreasing values

        if (val > vec(1)) then
           loc = 0; return
        else if (abs(val - vec(1)) < tiny ) then
           loc = 1; return
        else if (val <= vec(n)) then
           loc = n; return
        end if


        first = 1
        last = n-1
        loc = -1
        do while (first <= last)
           mid = (first + last)/2
           if (vec(mid) >= val) then
              if (val > vec(mid+1)) then
                 loc = mid
                 exit
              end if
              first = mid + 1
           else
              last = mid - 1
           end if
        end do

     else ! increasing values

        if (val < vec(1)) then
           loc = 0; return
        else if (abs(val - vec(1)) < tiny) then
           loc = 1; return
        else if (val >= vec(n)) then
           loc = n; return
        end if

        first = 1
        last = n-1
        loc = -1
        do while (first <= last)
           mid = (first + last)/2
           if (vec(mid) <= val) then
              if (val < vec(mid+1)) then
                 loc = mid
                 exit
              end if
              first = mid + 1
           else
              last = mid - 1
           end if
        end do

     end if

    end function binary_search
    
    subroutine write_eep_track(eep_filename,x)
    type(track), intent(in) :: x
    character(len=strlen), intent(in) :: eep_filename

    integer :: io, ierr, j, ncol
    real(dp) :: real_phase
    character(len=8) :: have_phase
    io=alloc_iounit(ierr)

    open(io,file=trim(eep_filename),action='write',status='unknown')
    have_phase = 'YES'

    !write(io,'(a25,a8)') '# MIST version number  = ', x% version_string
    !write(io,'(a25,i8)') '# MESA revision number = ', x% MESA_revision_number
!                     123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
   ! write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a88)') '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit                                        '
    write(io,'(a2,f6.4,1p1e13.5,0p3f9.2)') '# ', x% initial_Y, x% initial_Z, x% Fe_div_H, x% alpha_div_Fe, x% v_div_vcrit
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'
    write(io,'(a1,1x,a16,4a8,2x,a10)') '#','initial_mass', 'N_pts', 'N_EEP', 'N_col', 'phase', 'type'
    write(io,'(a1,1x,1p1e16.10,3i8,a8,2x,a10)') '#', x% initial_mass, x% ntrack, x% neep, ncol, have_phase, &
         star_label(x% star_type)
    write(io,'(a8,20i8)') '# EEPs: ', x% eep
    write(io,'(a88)') '# --------------------------------------------------------------------------------------'

       write(io,'(a1,299(27x,i5))') '#', (j,j=1,x% ncol)
       write(io,'(a1,299a32)') '#', adjustr(x% cols(:)% name), 'phase'
       do j=x% eep(1),x% ntrack
          real_phase=real(x% phase(j))            !added by Poojan
          write(io,'(1x,299(1pes32.16e3))') x% tr(:,j), real_phase
       enddo
    close(io)
    call free_iounit(io)
    end subroutine write_eep_track

    subroutine write_dat_track(tphys, pars)
        real(dp), intent(in) :: tphys
        type(star_parameters), intent(in) :: pars
        character(LEN=*), PARAMETER  :: FMT= '(1p9e15.6,2i10)'
        write(120,FMT) tphys,pars% age,pars% mass,pars% core_mass,pars% McHe, pars% McCO &
                    ,pars% log_L,pars% log_Teff,pars% log_R,pars% phase ,pars% extra
    end subroutine write_dat_track

!    subroutine alloc_track(filename,x)
!        character(len=strlen), intent(in) :: filename
!        type(eep_track), pointer :: x
!        allocate(x)
!        x% neep = primary
!        x% filename = trim(filename)
!        allocate(x% eep(x% neep))
!      end subroutine alloc_track

    
    elemental function pow10_sg(x) result(y)
        real(sp), intent(in) :: x
        real(sp) :: y
        y = exp(ln10_sp*x)
    end function pow10_sg

    elemental function pow10(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = exp(ln10*x)
    end function pow10
    
    !linear search alogorithm

    subroutine index_search(size_list,list,value,min_index,debug)
        integer, intent(in) :: size_list
        real(dp), intent(in) :: list(:)
        real(dp), intent(in) :: value
        integer, intent(out) :: min_index
        logical, optional :: debug
        real(dp) :: difference(size_list)

        if (present(debug)) then
            if (debug)  then
!            print*, "in index search",value,size_list
            print*, "in index search; list(1),list(size_list-1),list(size_list)"
            print*, list(1),list(size_list-1),list(size_list)
            endif
        endif
        
        if (value < list(1)) then             !from num_binary_search.inc
            min_index = 1; return
        elseif (check_equal(value, list(size_list)))then
            min_index = size_list
        else if (value > list(size_list)) then
            min_index = size_list+1; return
        end if

        difference = 0.d0
        difference = abs(list-value)
        min_index = minloc(difference, dim=1)
     
   
    end subroutine index_search
    
    subroutine save_values(new_line,pars)
        type(star_parameters) :: pars
        real(dp),intent (in) :: new_line(:,:)

        pars% mass = new_line(i_mass,1)
        pars% McHe = new_line(i_he_core,1)
        pars% McCO = new_line(i_co_core,1)
        pars% log_L = new_line(i_logL,1)
        pars% luminosity = 10**pars% log_L
!        pars% luminosity =  new_line(i_lum,1)
        pars% log_Teff = new_line(i_logTe,1)
        pars% Teff = 10**(pars% log_Teff)
        pars% log_R = new_line(i_logR,1)
!        pars% log_R = 2*(3.762+(0.25*pars% log_L)-pars% log_Teff )
        pars% radius = 10**pars% log_R
        if (pars% phase <= EAGB) then
            pars% core_mass = pars% McHe
        else
            pars% core_mass = pars% McCO
        endif
    end subroutine
    
    logical function check_ge(x,y) result(z)
    !TODO: needs to be checked before use
    real(dp), intent(in) :: x,y
        if (x>y .or. abs(x-y)<tiny) then
            z = .true.
        else
            z = .false.
        endif
        return
    end function check_ge
    
    logical function check_le(x,y) result(z)
    !TODO: needs to be checked before use
    real(dp), intent(in) :: x,y
        if (x<y .or. abs(x-y)<tiny) then
            z = .true.
        else
            z = .false.
        endif
        return
    end function check_le
    
    logical function check_equal(x,y,limit) result(z)
    real(dp), intent(in) :: x,y
    real(dp), intent(in), optional :: limit
    real(dp) :: diff, threshold

    if (present(limit)) then
        threshold = limit
        else
    threshold = 1d-4
    endif
        diff = abs(x-y)
        if (diff<=threshold) then
            z = .true.
        else
            z = .false.
        endif
        return
    end function check_equal

    logical function defined(x) result(y)
        real(dp), intent(in) :: x
        if (abs(x-undefined)<=tiny) then
            y = .false.
        else
            y = .true.
        endif
    return
    end function defined

    !same as the function 'defined' above but for integers
    logical function identified(x) result(y)
        integer, intent(in) :: x
        if (abs(x-undefined_i)<=tiny) then
            y = .false.
        else
            y = .true.
        endif
        return
    end function identified
    
    subroutine uniform_distribution(n, minval, maxval, marray)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: minval, maxval
        real(dp), intent(out) :: marray(n)
        real(dp) :: h
        integer :: i

        marray = 0.d0
        !linearly spaced
        h = abs(maxval- minval)/(n-1)
        do i= 1,n
            marray(i) = minval + (i-1)*h
        end do
        
        !for log spaced
    end subroutine uniform_distribution


    subroutine distance_along_track(t)
      type(track), intent(inout) :: t
      real(dp) :: tmp_dist, weight, max_center_h1
      integer :: j

      if(weight_center_rho_T_by_Xc)then
         max_center_h1 = maxval(t% tr(i_Xc,:))
         if(max_center_h1 <= 0d0) max_center_h1 = 1d0
      else
         max_center_h1 = 1d0
         weight = 1d0
      endif

!      t% dist(1) = 0d0
      if(t% ntrack > 3)then
         do j = 2, t% ntrack
            
            if(weight_center_rho_T_by_Xc)then
               weight = max(0d0, t% tr(i_Xc,j)/max_center_h1)
            endif
            
            !build up the distance between EEPs piece by piece
            tmp_dist =            Teff_scale*sqdiff(t% tr(i_logTe,j) , t% tr(i_logTe,j-1))
            tmp_dist = tmp_dist + logL_scale*sqdiff(t% tr(i_logL, j) , t% tr(i_logL, j-1))
            tmp_dist = tmp_dist + weight * Rhoc_scale * sqdiff(t% tr(i_Rhoc, j) , t% tr(i_Rhoc, j-1))
            tmp_dist = tmp_dist + weight * Tc_scale*  sqdiff(t% tr(i_Tc,   j) , t% tr(i_Tc,   j-1))
            tmp_dist = tmp_dist + age_scale* sqdiff(log10(t% tr(i_age,j)) , log10(t% tr(i_age,j-1)))

!            t% dist(j) = t% dist(j-1) + sqrt(tmp_dist)
         enddo
      endif
    end subroutine distance_along_track

    elemental function sqdiff(x0,x1) result(y) !square of x, y=x*x
      real(dp), intent(in) :: x0, x1
      real(dp) :: y, dx
      dx=x0-x1
      y = dx*dx
    end function sqdiff


end module track_support
