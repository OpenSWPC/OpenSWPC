!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! SAC-formatted seismograms
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_sac

  !! -- Dependency
  use m_std

  !! -- Declarations
  implicit none
  private

  !! -- Public Procedures
  public :: sac__hdr     ! sac data type
  public :: sac__write   ! write sac datafile
  public :: sac__read    ! read  sac datafile  
  public :: sac__init    ! initialize sac data type
  public :: sac__whdr    ! read header
  public :: csf__write   ! write concatenated sac format file

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Sac Header Type Definition
  !<
  !! --
  type sac__hdr

    !!               var name          description                   record#
    real(dp)      :: delta           ! sampling interval             (001)
    real(dp)      :: depmin          ! minimum value                 (002)
    real(dp)      :: depmax          ! maximum value                 (003)
    real(dp)      :: scale           ! multiplying scale factor      (004)
    real(dp)      :: odelta          ! Observed increment            (005)
    real(dp)      :: b               ! begenning independent value   (006)
    real(dp)      :: e               ! ending independent value      (007)
    real(dp)      :: o               ! event origin time             (008)
    real(dp)      :: a               ! first arrival time            (009)
    real(dp)      :: t0              ! time picks                    (011)
    real(dp)      :: t1              ! time picks                    (012)
    real(dp)      :: t2              ! time picks                    (013)
    real(dp)      :: t3              ! time picks                    (014)
    real(dp)      :: t4              ! time picks                    (015)
    real(dp)      :: t5              ! time picks                    (016)
    real(dp)      :: t6              ! time picks                    (017)
    real(dp)      :: t7              ! time picks                    (018)
    real(dp)      :: t8              ! time picks                    (019)
    real(dp)      :: t9              ! time picks                    (020)
    real(dp)      :: f               ! fini or end of event time     (021)
    real(dp)      :: resp0           ! instrument response param.    (022)
    real(dp)      :: resp1           ! instrument response param.    (023)
    real(dp)      :: resp2           ! instrument response param.    (024)
    real(dp)      :: resp3           ! instrument response param.    (025)
    real(dp)      :: resp4           ! instrument response param.    (026)
    real(dp)      :: resp5           ! instrument response param.    (027)
    real(dp)      :: resp6           ! instrument response param.    (028)
    real(dp)      :: resp7           ! instrument response param.    (029)
    real(dp)      :: resp8           ! instrument response param.    (030)
    real(dp)      :: resp9           ! instrument response param.    (031)
    real(dp)      :: stla            ! station latitude              (032)
    real(dp)      :: stlo            ! station longitude             (033)
    real(dp)      :: stel            ! station elevation (m)         (034)
    real(dp)      :: stdp            ! station depth (m)             (035)
    real(dp)      :: evla            ! event latitude                (036)
    real(dp)      :: evlo            ! event longitude               (037)
    real(dp)      :: evel            ! event elevation (m)           (038)
    real(dp)      :: evdp            ! event depth (m)               (039)
    real(dp)      :: mag             ! event magnitude               (040)
    real(dp)      :: user0           ! user header                   (041)
    real(dp)      :: user1           ! user header                   (042)
    real(dp)      :: user2           ! user header                   (043)
    real(dp)      :: user3           ! user header                   (044)
    real(dp)      :: user4           ! user header                   (045)
    real(dp)      :: user5           ! user header                   (046)
    real(dp)      :: user6           ! user header                   (047)
    real(dp)      :: user7           ! user header                   (048)
    real(dp)      :: user8           ! user header                   (049)
    real(dp)      :: user9           ! user header                   (050)
    real(dp)      :: dist            ! distance (km)                 (051)
    real(dp)      :: az              ! azimuth (deg)                 (052)
    real(dp)      :: baz             ! back azimuth (deg)            (053)
    real(dp)      :: gcarc           ! angular distance (deg)        (054)
    real(dp)      :: depmen          ! mean value                    (057)
    real(dp)      :: cmpaz           ! component azimuth             (058)
    real(dp)      :: cmpinc          ! component incident angle      (059) 
    real(dp)      :: xminimum        ! minimum value of x (spec)     (060)    
    real(dp)      :: xmaximum        ! maximum value of x (spec)     (061)    
    real(dp)      :: yminimum        ! minimum value of y (spec)     (062)    
    real(dp)      :: ymaximum        ! maximum value of y (spec)     (063)    
    integer       :: nzyear          ! reference time, year          (071)
    integer       :: nzjday          ! reference time, julian day    (072)
    integer       :: nzhour          ! reference time, hour          (073)
    integer       :: nzmin           ! reference time, minute        (074)
    integer       :: nzsec           ! reference time, second        (075)
    integer       :: nzmsec          ! reference time, millisecond   (076)
    integer       :: nvhdr           ! header version                (077)
    integer       :: norid           ! origin ID (CSS3.0)            (078)
    integer       :: nevid           ! event ID (CSS3.0)             (079)
    integer       :: npts            ! number of data points         (080)
    integer       :: nwfid           ! waveform ID (CSS3.0)          (082)
    integer       :: nxsize          ! spectral length               (083)
    integer       :: nysize          ! spectral width                (084)
    integer       :: iftype          ! type of file                  (086)
    integer       :: idep            ! type of dependent var.        (087)
    integer       :: iztype          ! reference time equivallence   (088)
    integer       :: iinst           ! instrument type               (090)
    integer       :: istreg          ! station region                (091)
    integer       :: ievreg          ! event region                  (092)
    integer       :: ievtyp          ! event type                    (093)
    integer       :: iqual           ! data quality                  (094)
    integer       :: isynth          ! synthetic data flag real=49   (095)
    integer       :: imagtyp         ! magnitude type                (096)
    integer       :: imagsrc         ! source of magnitude info.     (097)
    logical       :: leven           ! is evenly spaced file         (106)
    logical       :: lpspol          ! is positive polarity          (107)
    logical       :: lovrok          ! is overwrite ok?              (108)
    logical       :: lcalda          ! is calc distance azimuth      (109)
    character(8)  :: kstnm           ! station name                  (111)
    character(16) :: kevnm           ! event name                    (113)
    character(8)  :: khole           ! hole name                     (117)
    character(8)  :: ko              ! origin time identification    (119)
    character(8)  :: ka              ! time pick name                (121)
    character(8)  :: kt0             ! time pick name                (123)
    character(8)  :: kt1             ! time pick name                (125)
    character(8)  :: kt2             ! time pick name                (127)
    character(8)  :: kt3             ! time pick name                (129)
    character(8)  :: kt4             ! time pick name                (131)
    character(8)  :: kt5             ! time pick name                (133)
    character(8)  :: kt6             ! time pick name                (135)
    character(8)  :: kt7             ! time pick name                (137)
    character(8)  :: kt8             ! time pick name                (139)
    character(8)  :: kt9             ! time pick name                (141)
    character(8)  :: kf              ! fini identification           (143)
    character(8)  :: kuser0          ! user area                     (145)
    character(8)  :: kuser1          ! user area                     (147)
    character(8)  :: kuser2          ! user area                     (149)
    character(8)  :: kcmpnm          ! component name                (151)
    character(8)  :: knetwk          ! network name                  (153)
    character(8)  :: kdatrd          ! date data onto comp.          (155)
    character(8)  :: kinst           ! instrument                    (157)

    !! Unofficial headers at unused blocks
    real(dp) :: user10 ! user header (064)
    real(dp) :: user11 ! user header (065)
    real(dp) :: user12 ! user header (066)
    real(dp) :: user13 ! user header (067)
    real(dp) :: user14 ! user header (068)
    real(dp) :: user15 ! user header (069)
    real(dp) :: user16 ! user header (070)
    integer  :: iuser0 ! user header (098)
    integer  :: iuser1 ! user header (099)
    integer  :: iuser2 ! user header (100)
    integer  :: iuser3 ! user header (101)
    integer  :: iuser4 ! user header (102)
    integer  :: iuser5 ! user header (103)
    integer  :: iuser6 ! user header (104)
    integer  :: iuser7 ! user header (105)
    logical  :: luser0 ! user header (110)

    !! original header
    integer :: nzmonth
    integer :: nzday
    integer :: tim

  end type sac__hdr

  interface sac__write
    module procedure sac_d, sac_s
  end interface sac__write

  interface sac__read
     module procedure rsac_d, rsac_s
  end interface

  interface csf__write
    module procedure wcsf_d, wcsf_s
  end interface csf__write

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Read SAC file
  !<
  !! --
  subroutine rsac_d ( fn_sac, ss, dat )

    character(*),          intent(in)    :: fn_sac  !< sac filename
    type(sac__hdr),        intent(out)   :: ss      !< header info
    real(dp), allocatable, intent(inout) :: dat(:)  !< waveform data
    real(sp), allocatable                :: fdat(:)
    !----
    
    call rsac_s( fn_sac, ss, fdat )

    if( .not. allocated( dat ) )  allocate( dat(1:ss%npts) )
    dat = dble(fdat)
    deallocate(fdat)
    
  end subroutine rsac_d

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Read SAC file
  !<
  !! --  
  subroutine rsac_s( fn_sac, ss, dat )

    character(*),                    intent(in)    :: fn_sac  !< sac filename
    type(sac__hdr),                  intent(out)   :: ss      !< header info
    real(sp), allocatable, optional, intent(inout) :: dat(:)  !< waveform data
    !--
    integer :: io
    integer :: i
    integer :: nmax
    logical :: same_endian
    integer :: ierr
    !----
    call std__getio(io)
    open( io, file=fn_sac, &
          action='read', access='stream', form='unformatted', status='old', iostat=ierr)
    
    if( ierr /= 0 ) then
      write(stderr,*) '[sac__read]: file '// trim(fn_sac) // ' not opened. '
      return
    end if

    call sac__rhdr(io, ss, same_endian)
        
    if( present( dat ) ) then       
       
       if( .not. allocated( dat ) ) then
          allocate( dat(1:ss%npts))
          nmax = ss%npts
       else
          nmax = max( ss%npts, size(dat) )
       end if
       
       if( ss%npts > nmax ) then
          write(stderr,*) '[sac__read]: data array does not have enough size' 
       end if
       
       read(io) dat
       if( .not. same_endian ) then
          do i=1, ss%npts
             call change_endian_r( dat(i) )
          end do
       end if
       
    end if
    
    close( io )
    
  end subroutine rsac_s  

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC file
  !<
  !! --
  subroutine sac_d( fn_sac, ss, dat, overwrite )

    !! -- Arguments
    character(*),   intent(in)           :: fn_sac
    type(sac__hdr), intent(in)           :: ss
    real(DP),       intent(in)           :: dat(1:ss%npts)
    logical,        intent(in), optional :: overwrite

    !! ----

    if( present( overwrite) ) then
      call sac_s( fn_sac, ss, real(dat), overwrite)
    else
      call sac_s( fn_sac, ss, real(dat) )
    end if

  end subroutine sac_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC file
  !<
  subroutine sac_s( fn_sac, ss, dat, overwrite )

    !! -- Arguments
    character(*),   intent(in)           :: fn_sac
    type(sac__hdr), intent(in)           :: ss
    real(SP),       intent(in)           :: dat(1:ss%npts)
    logical,        intent(in), optional :: overwrite

    logical        :: isexist
    integer        :: io
    character(1)   :: yn
    !! ----

    !! overwrite check
    inquire( file = fn_sac, exist=isexist )
    if( isexist ) then
      if( present( overwrite) ) then
        if( .not. overwrite ) then
          write(STDERR,*) 'sac: file '//trim(fn_sac)//' exists.'
          write(STDERR,*) 'sac: could not overwrite the file.'
          write(STDERR,*) 'sac: return without success'
          write(STDERR,*)
          return
        end if
      else
        write(STDERR,*) 'sac: file '//trim(fn_sac)//' exists.'
        write(STDERR,*) 'sac: Overwrite ? (y/n)'
        read(STDIN,'(A)') yn
        if( yn /= 'y' .and. yn /='Y' ) then
          write(STDERR,*) 'sac: could not overwrite the file.'
          write(STDERR,*) 'sac: return without success'
          write(STDERR,*)
          return
        end if
      end if
    end if

    call std__getio( io )
    open( io, file=trim(fn_sac), action='write', access='stream', form='unformatted', status='replace' )

    call sac__whdr(io, ss)

    write( io ) dat(1:ss%npts)
    close( io )

  end subroutine sac_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC data header from pre-opened file io
  !! No endian conversion will be made. Always write in machine-endian.
  !<
  !! --
  subroutine sac__whdr(io, ss)

    !! -- Arguments
    integer,        intent(in) :: io
    type(sac__hdr), intent(in) :: ss

    real(SP)       :: fheader(70)
    integer        :: iheader(71:105)
    logical        :: lheader(106:110)
    character(4)   :: aheader(111:158)
    integer        :: i
    !! ----

    !! header initialize
    fheader(1:70) = -12345.0
    iheader(71:105) = -12345
    lheader(106:110) = .false.
    do i=111, 157, 2
      aheader( i ) = '-123'
      aheader( i+1 ) = '45'
    end do


    ! Copy header data to temprary arrays
    fheader(  1) = real( int( ss % delta * 1d7 ) ) / 1e7
    fheader(  2) = real( ss % depmin )
    fheader(  3) = real( ss % depmax )
    fheader(  6) = real( ss % b      )
    fheader(  7) = real( ss % e      )
    fheader(  8) = real( ss % o      )
    fheader(  9) = real( ss % a      )
    fheader( 11) = real( ss % t0     )
    fheader( 12) = real( ss % t1     )
    fheader( 13) = real( ss % t2     )
    fheader( 14) = real( ss % t3     )
    fheader( 15) = real( ss % t4     )
    fheader( 16) = real( ss % t5     )
    fheader( 17) = real( ss % t6     )
    fheader( 18) = real( ss % t7     )
    fheader( 19) = real( ss % t8     )
    fheader( 20) = real( ss % t9     )
    fheader( 32) = real( ss % stla   )
    fheader( 33) = real( ss % stlo   )
    fheader( 34) = real( ss % stel   )
    fheader( 35) = real( ss % stdp   )
    fheader( 36) = real( ss % evla   )
    fheader( 37) = real( ss % evlo   )
    fheader( 38) = real( ss % evel   )
    fheader( 39) = real( ss % evdp   )
    fheader( 40) = real( ss % mag    )
    fheader( 41) = real( ss % user0  )
    fheader( 42) = real( ss % user1  )
    fheader( 43) = real( ss % user2  )
    fheader( 44) = real( ss % user3  )
    fheader( 45) = real( ss % user4  )
    fheader( 46) = real( ss % user5  )
    fheader( 47) = real( ss % user6  )
    fheader( 48) = real( ss % user7  )
    fheader( 49) = real( ss % user8  )
    fheader( 50) = real( ss % user9  )
    fheader( 51) = real( ss % dist   )
    fheader( 52) = real( ss % az     )
    fheader( 53) = real( ss % baz    )
    fheader( 54) = real( ss % gcarc  )
    fheader( 57) = real( ss % depmen )
    fheader( 58) = real( ss % cmpaz  )
    fheader( 59) = real( ss % cmpinc )
    fheader( 64) = real( ss % user10 )
    fheader( 65) = real( ss % user11 )
    fheader( 66) = real( ss % user12 )
    fheader( 67) = real( ss % user13 )
    fheader( 68) = real( ss % user14 )
    fheader( 69) = real( ss % user15 )
    fheader( 70) = real( ss % user16 )

    iheader( 71) = ss % nzyear
    iheader( 72) = ss % nzjday
    iheader( 73) = ss % nzhour
    iheader( 74) = ss % nzmin
    iheader( 75) = ss % nzsec
    iheader( 76) = ss % nzmsec
    iheader( 77) = ss % nvhdr
    iheader( 80) = ss % npts
    iheader( 86) = ss % iftype
    iheader( 87) = ss % idep
    iheader( 93) = ss % ievtyp

    iheader( 98) = ss % iuser0
    iheader( 99) = ss % iuser1
    iheader(100) = ss % iuser2
    iheader(101) = ss % iuser3
    iheader(102) = ss % iuser4
    iheader(103) = ss % iuser5
    iheader(104) = ss % iuser6
    iheader(105) = ss % iuser7

    lheader(106) = ss % leven
    lheader(107) = ss % lpspol
    lheader(108) = ss % lovrok
    lheader(109) = ss % lcalda

    lheader(110) = ss % luser0

    aheader(111) = ss%kstnm(1:4);  aheader(112) = ss%kstnm(5:8)
    aheader(113) = ss%kevnm(1:4);  aheader(114) = ss%kevnm(5:8)
    aheader(115) = ss%kevnm(9:12); aheader(116) = ss%kevnm(13:16)
    aheader(117) = ss%khole(1:4);  aheader(118) = ss%khole(5:8)
    aheader(119) = ss%ko(1:4);     aheader(120) = ss%ko(5:8)
    aheader(121) = ss%ka(1:4);     aheader(122) = ss%ka(5:8)
    aheader(123) = ss%kt0(1:4);    aheader(124) = ss%kt0(5:8)
    aheader(125) = ss%kt1(1:4);    aheader(126) = ss%kt1(5:8)
    aheader(127) = ss%kt2(1:4);    aheader(128) = ss%kt2(5:8)
    aheader(129) = ss%kt3(1:4);    aheader(130) = ss%kt3(5:8)
    aheader(131) = ss%kt4(1:4);    aheader(132) = ss%kt4(5:8)
    aheader(133) = ss%kt5(1:4);    aheader(134) = ss%kt5(5:8)
    aheader(135) = ss%kt6(1:4);    aheader(136) = ss%kt6(5:8)
    aheader(137) = ss%kt7(1:4);    aheader(138) = ss%kt7(5:8)
    aheader(139) = ss%kt8(1:4);    aheader(140) = ss%kt8(5:8)
    aheader(141) = ss%kt9(1:4);    aheader(142) = ss%kt9(5:8)
    aheader(143) = ss%kf(1:4);     aheader(143) = ss%kf(5:8)
    aheader(145) = ss%kuser0(1:4); aheader(146) = ss%kuser0(5:8)
    aheader(147) = ss%kuser1(1:4); aheader(148) = ss%kuser1(5:8)
    aheader(149) = ss%kuser2(1:4); aheader(150) = ss%kuser2(5:8)
    aheader(151) = ss%kcmpnm(1:4); aheader(152) = ss%kcmpnm(5:8)
    aheader(153) = ss%knetwk(1:4); aheader(154) = ss%knetwk(5:8)
    aheader(155) = ss%kdatrd(1:4); aheader(156) = ss%kdatrd(5:8)
    aheader(157) = ss%kinst(1:4);  aheader(158) = ss%kinst(5:8)

    !! write
    write( io ) fheader(1:70)
    write( io ) iheader(71:105)
    write( io ) lheader(106:110)
    write( io ) aheader(111:158)

  end subroutine sac__whdr
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Initialize SAC header
  !<
  !! --
  subroutine sac__init( ss )

    !! -- Arguments
    type(sac__hdr), intent(inout) :: ss

    real(SP)     :: ferr = -12345.0_SP
    integer      :: ierr = -12345
    character(6) :: cerr = '-12345'
    !! ----

    ss%delta   = ferr
    ss%depmin  = ferr
    ss%depmax  = ferr
    ss%b       = ferr
    ss%e       = ferr
    ss%o       = ferr
    ss%a       = ferr
    ss%t0      = ferr
    ss%t1      = ferr
    ss%t2      = ferr
    ss%t3      = ferr
    ss%t4      = ferr
    ss%t5      = ferr
    ss%t6      = ferr
    ss%t7      = ferr
    ss%t8      = ferr
    ss%t9      = ferr
    ss%stla    = ferr
    ss%stlo    = ferr
    ss%stel    = ferr
    ss%stdp    = ferr
    ss%evla    = ferr
    ss%evlo    = ferr
    ss%evel    = ferr
    ss%evdp    = ferr
    ss%mag     = ferr
    ss%user0   = ferr
    ss%user1   = ferr
    ss%user2   = ferr
    ss%user3   = ferr
    ss%user4   = ferr
    ss%user5   = ferr
    ss%user6   = ferr
    ss%user7   = ferr
    ss%user8   = ferr
    ss%user9   = ferr
    ss%dist    = ferr
    ss%az      = ferr
    ss%baz     = ferr
    ss%gcarc   = ferr
    ss%depmen  = ferr
    ss%cmpaz   = ferr
    ss%cmpinc  = ferr
    ss%nzyear  = ierr
    ss%nzjday  = ierr
    ss%nzhour  = ierr
    ss%nzmin   = ierr
    ss%nzsec   = ierr
    ss%nzmsec  = ierr
    ss%nvhdr   = 6
    ss%npts    = ierr
    ss%iftype  = 01
    ss%ievtyp  = ierr
    ss%idep    = ierr
    ss%leven   = .true.
    ss%lpspol  = .false.
    ss%lovrok  = .true.
    ss%lcalda  = .true.
    ss%kstnm   = cerr
    ss%kcmpnm  = cerr
    ss%kevnm   = cerr
    ss%khole   = cerr
    ss%ko      = cerr
    ss%ka      = cerr
    ss%kt0     = cerr
    ss%kt1     = cerr
    ss%kt2     = cerr
    ss%kt3     = cerr
    ss%kt4     = cerr
    ss%kt5     = cerr
    ss%kt6     = cerr
    ss%kt7     = cerr
    ss%kt8     = cerr
    ss%kt9     = cerr
    ss%kf      = cerr
    ss%kuser0  = cerr
    ss%kuser1  = cerr
    ss%kuser2  = cerr
    ss%knetwk  = cerr
    ss%kdatrd  = cerr
    ss%kinst   = cerr

    !! inoficial headers
    ss%user10 = ferr
    ss%user11 = ferr
    ss%user12 = ferr
    ss%user13 = ferr
    ss%user14 = ferr
    ss%user15 = ferr
    ss%user16 = ferr
    ss%iuser0 = ierr
    ss%iuser1 = ierr
    ss%iuser2 = ierr
    ss%iuser3 = ierr
    ss%iuser4 = ierr
    ss%iuser5 = ierr
    ss%iuser6 = ierr
    ss%iuser7 = ierr
    ss%luser0 = .false.    

  end subroutine sac__init
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write csf format
  !<
  !! --
  subroutine wcsf_d( fn_csf, ntrace, npts, sh, dat, overwrite )

    !! Arguments
    character(*),      intent(in) :: fn_csf
    integer,           intent(in) :: ntrace
    integer,           intent(in) :: npts
    type(sac__hdr),    intent(in) :: sh(ntrace)
    real(DP),          intent(in) :: dat(npts,ntrace)
    logical, optional, intent(in) :: overwrite
    !! ----

    if( present(overwrite) ) then
      call wcsf_s( fn_csf, ntrace, npts, sh, real(dat), overwrite )
    else
      call wcsf_s( fn_csf, ntrace, npts, sh, real(dat) )
    end if

  end subroutine wcsf_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write csf format
  !<
  !! --
  subroutine wcsf_s( fn_csf, ntrace, npts, sh, dat, overwrite )

    !! Arguments
    character(*),      intent(in) :: fn_csf
    integer,           intent(in) :: ntrace
    integer,           intent(in) :: npts
    type(sac__hdr),    intent(in) :: sh(ntrace)
    real(SP),          intent(in) :: dat(npts,ntrace)
    logical, optional, intent(in) :: overwrite
    !! --
    logical :: isexist
    integer :: io
    integer :: i
    character(1) :: yn
    integer, parameter :: NVHDR = 1
    !! ----

    !! overwrite
    inquire( file = fn_csf, exist=isexist )
    if( isexist ) then
      if( present( overwrite) ) then
        if( .not. overwrite ) then
          write(STDERR,*) 'wcsf: file '//trim(fn_csf)//' exists.' 
          write(STDERR,*) 'wcsf: could not overwrite the file.'
          write(STDERR,*) 'wcsf: return without success'
          write(STDERR,*)
          return
        end if
      else
        write(STDERR,*) 'wcsf: file '//trim(fn_csf)//' exists.' 
        write(STDERR,*) 'wcsf: Overwrite ? (y/n)' 
        read(STDIN,'(A)') yn
        if( yn /= 'y' .and. yn /='Y' ) then
          write(STDERR,*) 'wcsf: could not overwrite the file.'
          write(STDERR,*) 'wcsf: return without success'
          write(STDERR,*)
          return
        end if
      end if
    end if

    !file check
    do i=1, ntrace
      if( sh(i)%npts /= npts ) then
        write(STDERR,*) "wcdf: npts mismatch"
        write(STDERR,*) "wcdf: return without success"
        return
      end if
    end do

    call std__getio(io)
    open(io, file=trim(fn_csf), action='write', access='stream', form='unformatted', status='unknown')

    write(io) 'CSFD'
    write(io) ntrace
    write(io) npts

    do i=1, ntrace
      call sac__whdr(io, sh(i))
      write(io) dat(1:npts, i)
    end do

    close(io)

  end subroutine wcsf_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !-----------------------------------------------------------------------------------------------!
  !> Read SAC data header from pre-opened file io
  !--
  subroutine sac__rhdr( io, ss, same_endian )

    integer,        intent(in) :: io
    type(sac__hdr), intent(out) :: ss
    logical,        intent(out), optional :: same_endian
    !--
    integer :: i
    real(sp)     :: fheader(1:70)
    integer      :: iheader(71:105)
    logical      :: lheader(106:110)
    integer      :: nvhdr
    character(4) :: aheader(111:158)
    logical      :: same_endian0
    !----
    
    read(io) fheader
    read(io) iheader
    read(io) lheader
    read(io) aheader
    nvhdr = iheader( 77) 
    same_endian0 = ( 1 <= nvhdr .and. nvhdr <= 9 )
    
    if( .not. same_endian0 ) then
      do i=  1, 70
        call change_endian_r(fheader(i))
      end do
      do i= 71,105
        call change_endian_i(iheader(i))
      end do
      do i=106,110
        call change_endian_l(lheader(i))
      end do
    end if
    
    ss % delta  = dble( int( fheader(  1) * 1d7 + 0.5 ) ) / 1d7
    ss % depmin = dble( fheader(  2)  )
    ss % depmax = dble( fheader(  3)  ) 
    ss % b      = dble( fheader(  6)  )
    ss % e      = dble( fheader(  7)  )
    ss % o      = dble( fheader(  8)  )
    ss % a      = dble( fheader(  9)  )
    ss % t0     = dble( fheader( 11)  )
    ss % t1     = dble( fheader( 12)  )
    ss % t2     = dble( fheader( 13)  )
    ss % t3     = dble( fheader( 14)  )
    ss % t4     = dble( fheader( 15)  )
    ss % t5     = dble( fheader( 16)  )
    ss % t6     = dble( fheader( 17)  )
    ss % t7     = dble( fheader( 18)  )
    ss % t8     = dble( fheader( 19)  )
    ss % t9     = dble( fheader( 20)  )
    ss % stla   = dble( fheader( 32)  )
    ss % stlo   = dble( fheader( 33)  )
    ss % stel   = dble( fheader( 34)  )
    ss % stdp   = dble( fheader( 35)  )
    ss % evla   = dble( fheader( 36)  )
    ss % evlo   = dble( fheader( 37)  )
    ss % evel   = dble( fheader( 38)  )
    ss % evdp   = dble( fheader( 39)  )
    ss % mag    = dble( fheader( 40)  )
    ss % user0  = dble( fheader( 41)  )
    ss % user1  = dble( fheader( 42)  )
    ss % user2  = dble( fheader( 43)  )
    ss % user3  = dble( fheader( 44)  )
    ss % user4  = dble( fheader( 45)  )
    ss % user5  = dble( fheader( 46)  )
    ss % user6  = dble( fheader( 47)  )
    ss % user7  = dble( fheader( 48)  )
    ss % user8  = dble( fheader( 49)  )
    ss % user9  = dble( fheader( 50)  )
    ss % dist   = dble( fheader( 51)  )
    ss % az     = dble( fheader( 52)  )
    ss % baz    = dble( fheader( 53)  )
    ss % gcarc  = dble( fheader( 54)  )
    ss % depmen = dble( fheader( 57)  )
    ss % cmpaz  = dble( fheader( 58)  )
    ss % cmpinc = dble( fheader( 59)  )
    
    ss % nzyear = iheader( 71) 
    ss % nzjday = iheader( 72) 
    ss % nzhour = iheader( 73) 
    ss % nzmin  = iheader( 74) 
    ss % nzsec  = iheader( 75) 
    ss % nzmsec = iheader( 76) 
    ss % nvhdr  = iheader( 77) 
    ss % npts   = iheader( 80) 
    ss % iftype = iheader( 86) 
    ss % idep   = iheader( 87) 
    ss % ievtyp = iheader( 93) 
    
    ss % leven  = lheader(106)
    ss % lpspol = lheader(107)
    ss % lovrok = lheader(108)
    ss % lcalda = lheader(109)
    
    ss % kstnm  = aheader(111) // aheader(112)
    ss % kevnm  = aheader(113) // aheader(114) // aheader(115) // aheader(116)
    ss % khole  = aheader(117) // aheader(118)
    ss % ko     = aheader(119) // aheader(120)
    ss % ka     = aheader(121) // aheader(122)
    ss % kt0    = aheader(123) // aheader(124)
    ss % kt1    = aheader(125) // aheader(126)
    ss % kt2    = aheader(127) // aheader(128)
    ss % kt3    = aheader(129) // aheader(130)
    ss % kt4    = aheader(131) // aheader(132)
    ss % kt5    = aheader(133) // aheader(134)
    ss % kt6    = aheader(135) // aheader(136)
    ss % kt7    = aheader(137) // aheader(138)
    ss % kt8    = aheader(139) // aheader(140)
    ss % kt9    = aheader(141) // aheader(142)
    ss % kf     = aheader(143) // aheader(144)
    ss % kuser0 = aheader(145) // aheader(146)
    ss % kuser1 = aheader(147) // aheader(148)
    ss % kuser2 = aheader(149) // aheader(150)
    ss % kcmpnm = aheader(151) // aheader(152)
    ss % knetwk = aheader(153) // aheader(154)
    ss % kdatrd = aheader(155) // aheader(156)
    ss % kinst  = aheader(157) // aheader(158)
    
    call char_zeropad( ss%kstnm  )
    call char_zeropad( ss%kevnm  )
    call char_zeropad( ss%kt0    )
    call char_zeropad( ss%kt1    )
    call char_zeropad( ss%kt2    )
    call char_zeropad( ss%kt3    )
    call char_zeropad( ss%kt4    )
    call char_zeropad( ss%kt5    )
    call char_zeropad( ss%kt6    )
    call char_zeropad( ss%kt7    )
    call char_zeropad( ss%kt8    )
    call char_zeropad( ss%kt9    )
    call char_zeropad( ss%kuser0 )
    call char_zeropad( ss%kuser1 )
    call char_zeropad( ss%kuser2 )
    call char_zeropad( ss%kcmpnm )
    call char_zeropad( ss%knetwk )
    
    !! unofficieal headers
    ss % user10 = dble( fheader( 64) )
    ss % user11 = dble( fheader( 65) )
    ss % user12 = dble( fheader( 66) )
    ss % user13 = dble( fheader( 67) )
    ss % user14 = dble( fheader( 68) )
    ss % user15 = dble( fheader( 69) )
    ss % user16 = dble( fheader( 70) )
    ss % iuser0 = iheader( 98)
    ss % iuser1 = iheader( 99)
    ss % iuser2 = iheader(100)
    ss % iuser3 = iheader(101)
    ss % iuser4 = iheader(102)
    ss % iuser5 = iheader(103)
    ss % iuser6 = iheader(104)
    ss % iuser7 = iheader(105)
    ss % luser0 = lheader(110)

    if( present(same_endian) ) same_endian = same_endian0

  end subroutine sac__rhdr

  !-----------------------------------------------------------------------------------------------!
  subroutine byteswap(nbyte, foo)

    integer,          intent(in)    :: nbyte ! must be even
    character(nbyte), intent(inout) :: foo
    !--
    integer  :: i, j
    character(1) :: wk
    !----

    do i=1, nbyte/2
      j = nbyte - i + 1
      wk       = foo(i:i)
      foo(i:i) = foo(j:j)
      foo(j:j) = wk
    end do    

  end subroutine byteswap

  !-----------------------------------------------------------------------------------------------!  
  subroutine change_endian_r( var )

    real, intent(inout) :: var
    character(4) :: c
    
    c = transfer( var, c )
    call byteswap( 4, c )
    var = transfer( c, var )
    
  end subroutine change_endian_r

  !-----------------------------------------------------------------------------------------------!  
  subroutine change_endian_i( var )

    integer, intent(inout) :: var
    character(4) :: c
    
    c = transfer( var, c )
    call byteswap( 4, c)
    var = transfer( c, var )
    
  end subroutine change_endian_i

  !-----------------------------------------------------------------------------------------------!  
  subroutine change_endian_l( var )

    logical, intent(inout) :: var
    character(4) :: c

    c = transfer( var, c )
    call byteswap( 4, c )
    var = transfer( c, var )
    
  end subroutine change_endian_l

  !-----------------------------------------------------------------------------------------------!
  !> Substitute blanks in character ch by padding '0'
  !--
  subroutine char_zeropad( ch )

    character(*), intent(inout) :: ch
    !--
    integer :: i
    !----
    
    do i=1, len(ch)
       if( ichar(ch(i:i))== 0 )then
          ch(i:i) = ' ' 
       end if
    end do
    
  end subroutine char_zeropad  

end module m_sac
!! ----------------------------------------------------------------------------------------------------------------------------- !!
