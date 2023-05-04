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

!    real(dp), parameter :: ln10=log(1.0d1)
!    real(sp), parameter :: ln10_sp = log(10.0)
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

    character(len=10) :: star_label(4)
    character(len=5) :: phase_label(16)
    
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
