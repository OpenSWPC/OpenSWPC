!!
!! waveform output
!!
!! copyright
!!   Copyright 2024 Takuto Maeda. All rights reseaved. This project is released under the MIT license. 
!!
#include "../shared/m_debug.h"
module m_wav

    use m_std
    use m_debug
    use m_global
    use m_pwatch
    use m_sac
    use m_readini
    use m_system
    use m_geomap
    implicit none
    private
    save

    public :: wav__setup
    public :: wav__store
    public :: wav__write
    public :: wav__stquery

    integer :: ntdec_w
    character(3) :: wav_format !! sac or csf
    logical :: wav_calc_dist = .false.

    logical :: sw_wav_v      = .false.
    logical :: sw_wav_u      = .false.
    logical :: sw_wav_stress = .false.
    logical :: sw_wav_strain = .false.

      !! waveform
    integer :: ntw ! number of wave samples
    real(SP), allocatable :: wav_vel(:,:,:), wav_disp(:,:,:)
    real(SP), allocatable :: wav_stress(:,:,:), wav_strain(:,:,:) ! (6,nt,nst)
    type(sac__hdr), allocatable :: sh_vel   (:,:), sh_disp  (:,:)
    type(sac__hdr), allocatable :: sh_stress(:,:), sh_strain(:,:)
    real(SP), allocatable :: ux(:), uy(:), uz(:)
    real(SP), allocatable :: exx(:), eyy(:), ezz(:), eyz(:), exz(:), exy(:)
    
    integer :: nst = 0
    real(SP),     allocatable :: xst(:), yst(:), zst(:)
    integer,      allocatable :: ist(:), jst(:), kst(:)
    real(SP),     allocatable :: stlo(:), stla(:)
    character(8), allocatable :: stnm(:)

    real(MP) :: r40x, r40y, r40z, r41x, r41y, r41z

contains

    subroutine wav__setup(io_prm)

        integer, intent(in) :: io_prm
        !! --
        logical :: green_mode
        character(256) :: fn_stloc
        character(2) :: st_format
        !! ----

        call pwatch__on('wav__setup')

        call readini( io_prm, 'green_mode', green_mode, .false. )        
        if (.not. green_mode) then
            sw_wav_v      = .false.
            sw_wav_u      = .false.
            sw_wav_stress = .false.
            sw_wav_strain = .false.
            nst = 0
            return
        end if

        call readini( io_prm, 'ntdec_w',       ntdec_w,     10      )
        call readini( io_prm, 'sw_wav_v',      sw_wav_v,   .false. )
        call readini( io_prm, 'sw_wav_u',      sw_wav_u,   .false. )
        call readini( io_prm, 'sw_wav_stress', sw_wav_stress,   .false. )
        call readini( io_prm, 'sw_wav_strain', sw_wav_strain,   .false. )    
        call readini( io_prm, 'wav_format',    wav_format, 'sac' ) 
        call readini( io_prm, 'wav_calc_dist', wav_calc_dist, .false. )
        call readini( io_prm, 'st_format',     st_format,  'xy'    )
        call readini( io_prm, 'fn_stloc',      fn_stloc,   ''      )

        if(.not. ( sw_wav_v .or. sw_wav_u .or. sw_wav_stress .or. sw_wav_strain)) then
            nst = 0
            call pwatch__off('wav__setup')
            return
        end if

        ntw = floor( real(nt-1)/real(ntdec_w) + 1.0 )
      
        call set_stinfo(fn_stloc, st_format)

        call system__call('mkdir -p '//trim(odir)//'/wav > /dev/null 2>&1' )

        if (sw_wav_v) then
            allocate(wav_vel(ntw,3,nst))
            allocate(sh_vel(3,nst))
            wav_vel(:,:,:) = 0.0
        end if
          
        if (sw_wav_u ) then
            allocate(wav_disp(ntw,3,nst))
            allocate(sh_disp(3,nst))
            wav_disp(:,:,:) = 0.0
        end if
          
        if (sw_wav_stress) then
            allocate(wav_stress(ntw,6,nst))
            allocate(sh_stress(6,nst))
            wav_stress(:,:,:) = 0.0
        end if
      
        if (sw_wav_strain) then
            allocate(wav_strain(ntw,6,nst))
            allocate(sh_strain(6,nst))
            wav_strain(:,:,:) = 0.0
        end if        

        call set_sac_header()


        !! FDM coefficients
        r40x = 9.0_MP /  8.0_MP / dx
        r40y = 9.0_MP /  8.0_MP / dy
        r40z = 9.0_MP /  8.0_MP / dz
        r41x = 1.0_MP / 24.0_MP / dx
        r41y = 1.0_MP / 24.0_MP / dy
        r41z = 1.0_MP / 24.0_MP / dz
                
        call pwatch__off('wav__setup')

    end subroutine wav__setup

    subroutine set_stinfo(fn_stloc, st_format)

        character(*), intent(in) :: fn_stloc
        character(*), intent(in) :: st_format
        !! --
        integer :: io_stlst
        integer :: err
        integer :: nst_g
        character(256) :: abuf
        real(SP) :: xst_g, yst_g, zst_g, stlo_g, stla_g
        integer  :: ist_g, jst_g, kst_g
        character(3) :: zsw_g
        character(8) :: stnm_g
        !! ----

        open(newunit=io_stlst, file=trim(fn_stloc), action='read', status='old', iostat=err)

        if(err /= 0) then
            call info('no station location file found')
            nst = 0
            return
        end if

        xst  = [real(SP) ::]
        yst  = [real(SP) ::]
        zst  = [real(SP) ::]
        stlo = [real(SP) ::]
        stla = [real(SP) ::]
        ist  = [integer  ::]
        jst  = [integer  ::]
        kst  = [integer  ::]      

        do
            read(io_stlst, '(a256)', iostat=err) abuf
            if( err /= 0 ) exit

            abuf = trim(adjustl(abuf))
            if( abuf(1:1)           == "#" ) cycle
            if( trim(adjustl(abuf)) == ""  ) cycle

            select case(st_format)

            case('xy')

                read(abuf,*) xst_g, yst_g, zst_g, stnm_g, zsw_g
                call geomap__c2g(xst_g, yst_g, clon, clat, phi, stlo_g, stla_g)
    
            case('ll')
                read(abuf,*) stlo_g, stla_g, zst_g, stnm_g, zsw_g
                call assert( -360. <= stlo_g .and. stlo_g <= 360. )
                call assert(  -90. <= stla_g .and. stla_g <=  90. )
                call geomap__g2c(stlo_g, stla_g, clon, clat, phi, xst_g, yst_g)

            case default 
                
                call info('unknown st_format: ' // st_format)
                call assert( .false. )

            end select

            ist_g = x2i ( xst_g, xbeg, real(dx) )
            jst_g = y2j ( yst_g, ybeg, real(dy) )
      
            if(i2x(1, xbeg, real(dx)) < xst_g .and. xst_g < i2x(nx, xbeg, real(dx)) .and. &
               j2y(1, ybeg, real(dy)) < yst_g .and. yst_g < j2y(ny, ybeg, real(dy)) .and. &
                      kbeg            < kst_g .and. kst_g <         kend                  )  then

                nst_g = nst_g + 1

                if(ibeg <= ist_g .and. ist_g <= iend .and. jbeg <= jst_g .and. jst_g <= jend) then

                    !! station depth setting
                    select case ( zsw_g )
                    case( 'dep' );  kst_g = z2k (zst_g, zbeg, real(dz))
                    case( 'fsb' );  kst_g = kfs(ist_g, jst_g) + 1  !! free surface: one-grid below for avoiding vacuum
                    case( 'obb' );  kst_g = kob(ist_g, jst_g) + 1   !! ocean column: below seafloor
                    case( 'oba' );  kst_g = kob(ist_g, jst_g) - 1   !! ocean column: above seafloor
                    case( 'bd0' );  kst_g = z2k(bddep(ist_g, jst_g, 0), zbeg, real(dz)) 
                    case( 'bd1' );  kst_g = z2k(bddep(ist_g, jst_g, 1), zbeg, real(dz)) 
                    case( 'bd2' );  kst_g = z2k(bddep(ist_g, jst_g, 2), zbeg, real(dz)) 
                    case( 'bd3' );  kst_g = z2k(bddep(ist_g, jst_g, 3), zbeg, real(dz)) 
                    case( 'bd4' );  kst_g = z2k(bddep(ist_g, jst_g, 4), zbeg, real(dz)) 
                    case( 'bd5' );  kst_g = z2k(bddep(ist_g, jst_g, 5), zbeg, real(dz)) 
                    case( 'bd6' );  kst_g = z2k(bddep(ist_g, jst_g, 6), zbeg, real(dz)) 
                    case( 'bd7' );  kst_g = z2k(bddep(ist_g, jst_g, 7), zbeg, real(dz)) 
                    case( 'bd8' );  kst_g = z2k(bddep(ist_g, jst_g, 8), zbeg, real(dz)) 
                    case( 'bd9' );  kst_g = z2k(bddep(ist_g, jst_g, 9), zbeg, real(dz)) 
                    case default 
                        call info( "unknown zsw type in station file. Assume 'dep'")
                        kst_g = z2k (zst_g, zbeg, real(dz))
                    end select

                    if(kst_g > kend) then
                        call info("station depth exceeds kend at station " // trim(stnm_g))
                        kst_g = kend - 1
                      end if
                      if( kst_g < kbeg ) then
                        call info("station depth fall short of kbeg at station" // trim(stnm_g))
                        kst_g = kbeg + 1
                      end if                    

                    nst  = nst + 1
                    xst  = [xst,  xst_g ]
                    yst  = [yst,  yst_g ]
                    zst  = [zst,  zst_g ]
                    stnm = [stnm, stnm_g]
                    stlo = [stlo, stlo_g]
                    stla = [stla, stla_g]
                    ist  = [ist,  ist_g ]
                    jst  = [jst,  jst_g ]
                    kst  = [kst,  kst_g ]

                end if
            else
                if( myid == 0 )  call info("station " // trim(stnm_g) // " is out of the region")                
            end if
        end do
        
        close(io_stlst)        

        if( nst_g == 0 ) then
            nst = 0
            if( myid == 0 ) call info("no station is detected. waveform files will not be created")
            return
        end if
      
    end subroutine set_stinfo

    subroutine set_sac_header()

        integer :: i, j
        real :: mag
        !! ----

        mag = moment_magnitude(M0)

        do i=1, nst

            if (sw_wav_v) then
                do j=1, 3
                    call initialize_sac_header(sh_vel(j,i), stnm(i), stlo(i), stla(i), xst(i), yst(i), zst(i), mag)
                end do
                sh_vel(1,i)%kcmpnm = "Vx"; sh_vel(1,i)%cmpinc = 90.0; sh_vel(1,i)%cmpaz =  0.0 + phi
                sh_vel(2,i)%kcmpnm = "Vy"; sh_vel(2,i)%cmpinc = 90.0; sh_vel(2,i)%cmpaz = 90.0 + phi
                sh_vel(3,i)%kcmpnm = "Vz"; sh_vel(3,i)%cmpinc = 90.0; sh_vel(3,i)%cmpaz =  0.0 

                sh_vel(:,i)%idep = 7 
            end if

            if (sw_wav_u) then
                do j=1, 3
                    call initialize_sac_header(sh_disp(j,i), stnm(i), stlo(i), stla(i), xst(i), yst(i), zst(i), mag)
                end do
                sh_vel(1,i)%kcmpnm = "Ux"; sh_vel(1,i)%cmpinc = 90.0; sh_vel(1,i)%cmpaz =  0.0 + phi
                sh_vel(2,i)%kcmpnm = "Uy"; sh_vel(2,i)%cmpinc = 90.0; sh_vel(2,i)%cmpaz = 90.0 + phi
                sh_vel(3,i)%kcmpnm = "Uz"; sh_vel(3,i)%cmpinc = 90.0; sh_vel(3,i)%cmpaz =  0.0 

                sh_vel(:,i)%idep = 6 
            end if

            if (sw_wav_stress) then
                do j=1, 6
                    call initialize_sac_header(sh_stress(j,i), stnm(i), stlo(i), stla(i), xst(i), yst(i), zst(i), mag)
                end do

                sh_stress(1,i)%kcmpnm = "Sxx"
                sh_stress(2,i)%kcmpnm = "Syy"
                sh_stress(3,i)%kcmpnm = "Szz"
                sh_stress(4,i)%kcmpnm = "Syz"
                sh_stress(5,i)%kcmpnm = "Sxz"
                sh_stress(6,i)%kcmpnm = "Sxy"
        
                sh_stress(:,i)%idep = 5
            end if

            if (sw_wav_strain) then
                do j=1, 6
                    call initialize_sac_header(sh_strain(j,i), stnm(i), stlo(i), stla(i), xst(i), yst(i), zst(i), mag)
                end do

                sh_stress(1,i)%kcmpnm = "Exx"
                sh_stress(2,i)%kcmpnm = "Eyy"
                sh_stress(3,i)%kcmpnm = "Ezz"
                sh_stress(4,i)%kcmpnm = "Eyz"
                sh_stress(5,i)%kcmpnm = "Exz"
                sh_stress(6,i)%kcmpnm = "Exy"
        
                sh_stress(:,i)%idep = 5
            end if

        end do
    end subroutine set_sac_header

    subroutine initialize_sac_header(sh, stnm, stlo, stla, xst, yst, zst, mag)

        use m_daytim
        type(sac__hdr), intent(inout) :: sh
        character(*), intent(in) :: stnm
        real, intent(in) :: stlo, stla, xst, yst, zst, mag
        logical, save :: first_call = .true.
        type(sac__hdr), save :: sh0
        !! ----

        if (first_call) then
            call sac__init(sh0)
            sh0%evlo    = evlo
            sh0%evla    = evla
            sh0%evdp    = evdp  !! km-unit after SWPC5.0
            sh0%tim     = exedate
            sh0%b       = tbeg
            sh0%delta   = ntdec_w * dt
            sh0%npts    = ntw
            sh0%mag     = mag

            if( bf_mode ) then
                sh0%user0   = fx0
                sh0%user1   = fy0
                sh0%user2   = fz0
            else
                sh0%user0   = mxx0
                sh0%user1   = myy0
                sh0%user2   = mzz0
                sh0%user3   = myz0
                sh0%user4   = mxz0
                sh0%user5   = mxy0
            end if
                    
            sh0%user6   = clon  !< coordinate
            sh0%user7   = clat  !< coordinate
            sh0%user8   = phi
            sh0%o       = otim
    
            call daytim__localtime( exedate, sh0%nzyear, sh0%nzmonth, sh0%nzday, sh0%nzhour, sh0%nzmin, sh0%nzsec )
            call daytim__ymd2jul  ( sh0%nzyear, sh0%nzmonth, sh0%nzday, sh0%nzjday )
            sh%nzmsec  = 0

            first_call = .false.
        end if

        sh = sh0
        sh%kevnm = trim(adjustl( title(1:16) ))
        sh%kstnm = trim(stnm)
        sh%stlo  = stlo
        sh%stla  = stla
        sh%stdp  = zst * 1000 ! in meter unit

        if (wav_calc_dist) then
            sh%lcalda = .false. 
            sh%dist = sqrt((sx0 - xst)**2 + (sy0 - yst)**2)
            sh%az   = std__rad2deg(atan2(yst - sy0, xst - sx0))
            sh%baz  = std__rad2deg(atan2(sy0 - yst, sx0 - xst))
        end if
            
    end subroutine initialize_sac_header

    subroutine wav__store(it)

        integer, intent(in) :: it
        !! --
        integer :: n, itw
        real(MP) :: dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz
        integer :: i, j, k
        !! ----

        call pwatch__on("wav__store")
        if (nst == 0) then
            call pwatch__off("wav__store")
            return
        end if

        if( it == 1 ) then
            if (sw_wav_u) then
                allocate(ux(nst), source=0.0)
                allocate(uy(nst), source=0.0)
                allocate(uz(nst), source=0.0)
            end if

            if( sw_wav_strain ) then
                allocate(exx(nst), source=0.0)
                allocate(eyy(nst), source=0.0)
                allocate(ezz(nst), source=0.0)
                allocate(eyz(nst), source=0.0)
                allocate(exz(nst), source=0.0)
                allocate(exy(nst), source=0.0)
            end if
        end if

        if (sw_wav_u) then

           !$omp parallel do private(n,i,j,k)
            do n=1, nst
                i = ist(n); j = jst(n); k = kst(n)
                ux(n) = ux(n) + (Vx(k, i, j) + Vx(k,   i-1, j  )) * 0.5 * dt
                uy(n) = uy(n) + (Vy(k, i, j) + Vy(k,   i  , j-1)) * 0.5 * dt
                uz(n) = uz(n) - (Vz(k, i, j) + Vz(k-1, i  , j  )) * 0.5 * dt ! positive upward
            end do
            !$omp end parallel do

        end if

        if (sw_wav_strain) then

            !$omp parallel do private(n, i, j, k, dxVx, dyVy, dzVz, dyVz, dzVy, dxVz, dzVx, dxVy, dyVx)
            do n=1, nst
                i = ist(n); j = jst(n); k = kst(n)
      
                dxVx =   (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r40x - (Vx(k  ,i+1,j  ) - Vx(k  ,i-2,j  )) * r41x
                dyVy =   (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r40y - (Vy(k  ,i  ,j+1) - Vy(k  ,i  ,j-2)) * r41y
                dzVz =   (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r40z - (Vz(k+1,i  ,j  ) - Vz(k-2,i  ,j  )) * r41z
      
                dxVy = ( (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r40x - (Vy(k  ,i+2,j  ) - Vy(k  ,i-1,j  )) * r41x &
                       + (Vy(k  ,i+1,j-1) - Vy(k  ,i  ,j-1)) * r40x - (Vy(k  ,i+2,j-1) - Vy(k  ,i-1,j-1)) * r41x &
                       + (Vy(k  ,i  ,j  ) - Vy(k  ,i-1,j  )) * r40x - (Vy(k  ,i+1,j  ) - Vy(k  ,i-2,j  )) * r41x &
                       + (Vy(k  ,i  ,j-1) - Vy(k  ,i-1,j-1)) * r40x - (Vy(k  ,i+1,j-1) - Vy(k  ,i-2,j-1)) * r41x ) / 4.0
       
                dxVz = ( (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r40x - (Vz(k  ,i+2,j  ) - Vz(k  ,i-1,j  )) * r41x &
                       + (Vz(k-1,i+1,j  ) - Vz(k-1,i  ,j  )) * r40x - (Vz(k-1,i+2,j  ) - Vz(k-1,i-1,j  )) * r41x &
                       + (Vz(k  ,i  ,j  ) - Vz(k  ,i-1,j  )) * r40x - (Vz(k  ,i+1,j  ) - Vz(k  ,i-2,j  )) * r41x &
                       + (Vz(k-1,i  ,j  ) - Vz(k-1,i-1,j  )) * r40x - (Vz(k-1,i+1,j  ) - Vz(k-1,i-2,j  )) * r41x ) / 4.0
        
                dyVx = ( (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r40y - (Vx(k  ,i  ,j+2) - Vx(k  ,i  ,j-1)) * r41y &
                       + (Vx(k  ,i-1,j+1) - Vx(k  ,i-1,j  )) * r40y - (Vx(k  ,i-1,j+2) - Vx(k  ,i-1,j-1)) * r41y &
                       + (Vx(k  ,i  ,j  ) - Vx(k  ,i  ,j-1)) * r40y - (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j-2)) * r41y &
                       + (Vx(k  ,i-1,j  ) - Vx(k  ,i-1,j-1)) * r40y - (Vx(k  ,i-1,j+1) - Vx(k  ,i-1,j-2)) * r41y ) / 4.0
        
                dyVz = ( (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r40y - (Vz(k  ,i  ,j+2) - Vz(k  ,i  ,j-1)) * r41y &
                       + (Vz(k-1,i  ,j+1) - Vz(k-1,i  ,j  )) * r40y - (Vz(k-1,i  ,j+2) - Vz(k-1,i  ,j-1)) * r41y &
                       + (Vz(k  ,i  ,j  ) - Vz(k  ,i  ,j-1)) * r40y - (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j-2)) * r41y &
                       + (Vz(k-1,i  ,j  ) - Vz(k-1,i  ,j-1)) * r40y - (Vz(k-1,i  ,j+1) - Vz(k-1,i  ,j-2)) * r41y ) / 4.0
        
                dzVx = ( (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r40z - (Vx(k+2,i  ,j  ) - Vx(k-1,i  ,j  )) * r41z &
                       + (Vx(k+1,i-1,j  ) - Vx(k  ,i-1,j  )) * r40z - (Vx(k+2,i-1,j  ) - Vx(k-1,i-1,j  )) * r41z &
                       + (Vx(k  ,i  ,j  ) - Vx(k-1,i  ,j  )) * r40z - (Vx(k+1,i  ,j  ) - Vx(k-2,i  ,j  )) * r41z &
                       + (Vx(k  ,i-1,j  ) - Vx(k-1,i-1,j  )) * r40z - (Vx(k+1,i-1,j  ) - Vx(k-2,i-1,j  )) * r41z ) / 4.0
        
                dzVy = ( (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r40z - (Vy(k+2,i  ,j  ) - Vy(k-1,i  ,j  )) * r41z &
                       + (Vy(k+1,i  ,j-1) - Vy(k  ,i  ,j-1)) * r40z - (Vy(k+2,i  ,j-1) - Vy(k-1,i  ,j-1)) * r41z &
                       + (Vy(k  ,i  ,j  ) - Vy(k-1,i  ,j  )) * r40z - (Vy(k+1,i  ,j  ) - Vy(k-2,i  ,j  )) * r41z &
                       + (Vy(k  ,i  ,j-1) - Vy(k-1,i  ,j-1)) * r40z - (Vy(k+1,i  ,j-1) - Vy(k-2,i  ,j-1)) * r41z ) / 4.0
      
                exx(i) = exx(i) + dxVx * dt
                eyy(i) = eyy(i) + dyVy * dt
                ezz(i) = ezz(i) + dzVz * dt
                eyz(i) = eyz(i) + (dyVz + dzVy) / 2.0 * dt
                exz(i) = exz(i) + (dxVz + dzVx) / 2.0 * dt
                exy(i) = exy(i) + (dxVy + dyVx) / 2.0 * dt
      
            end do
            !$omp end parallel do
        end if

        if (mod(it-1, ntdec_w) == 0) then

            itw = (it-1)/ntdec_w + 1
            if( sw_wav_v ) then

                !$omp parallel do private(n,i,j,k)
                do n=1, nst
                    i = ist(n); j = jst(n); k = kst(n)
                    wav_vel(itw,1,i) =  ( Vx(k, i, j) + Vx(k  , i-1, j  )) / 2.0 * M0 * UC * 1e9
                    wav_vel(itw,2,i) =  ( Vy(k, i, j) + Vy(k  , i  , j-1)) / 2.0 * M0 * UC * 1e9
                    wav_vel(itw,3,i) = -( Vz(k, i, j) + Vz(k-1, i  , j  )) / 2.0 * M0 * UC * 1e9
                end do
                !$omp end parallel do

            end if
          
            if( sw_wav_u ) then

                !$omp parallel do private(n)
                do n=1, nst
                    wav_disp(itw,1,n) = ux(n) * M0 * UC * 1e9
                    wav_disp(itw,2,n) = uy(n) * M0 * UC * 1e9
                    wav_disp(itw,3,n) = uz(n) * M0 * UC * 1e9
                end do
                !$omp end parallel do

            end if
          
            if( sw_wav_stress ) then

                !$omp parallel do private(n, i, j, k)
                do n=1, nst
                    i = ist(n); j = jst(n); k = kst(n)
                    wav_stress(itw,1,n) = Sxx(k, i, j) * M0 * UC * 1e6
                    wav_stress(itw,2,n) = Syy(k, i, j) * M0 * UC * 1e6
                    wav_stress(itw,3,n) = Szz(k, i, j) * M0 * UC * 1e6
                    wav_stress(itw,4,n) = ( Syz(k,   i,   j) + Syz(k,   i,   j-1)  &
                                          + Syz(k-1, i,   j) + Syz(k-1, i,   j-1)) / 4.0 * M0 * UC * 1e6
                    wav_stress(itw,5,n) = ( Sxz(k,   i,   j) + Sxz(k,   i-1, j  )  &
                                          + Sxz(k-1, i,   j) + Sxz(k-1, i-1, j  )) / 4.0 * M0 * UC * 1e6
                    wav_stress(itw,6,n) = ( Sxy(k,   i,   j) + Sxy(k,   i,   j-1)  &
                                          + Sxy(k,   i-1, j) + Sxy(k,   i-1, j-1)) / 4.0 * M0 * UC * 1e6
                end do
                !$omp end parallel do
      
            end if
          
            if( sw_wav_strain ) then

                !$omp parallel do private(n)
                do n=1, nst
                    wav_strain(itw,1,n) = exx(n) * M0 * UC * 1e-3
                    wav_strain(itw,2,n) = eyy(n) * M0 * UC * 1e-3
                    wav_strain(itw,3,n) = ezz(n) * M0 * UC * 1e-3
                    wav_strain(itw,4,n) = eyz(n) * M0 * UC * 1e-3
                    wav_strain(itw,5,n) = exz(n) * M0 * UC * 1e-3
                    wav_strain(itw,6,n) = exy(n) * M0 * UC * 1e-3
                end do
                !$omp end parallel do
        
            end if

        end if
        
        call pwatch__off("wav__store")

    end subroutine wav__store

    subroutine wav__stquery(stnm1, is_exist, ist1, jst1, kst1, xst1, yst1, zst1, stlo1, stla1)
    !! Return location and indices of the requested station

        character(*), intent(in)    :: stnm1
        logical,      intent(out)   :: is_exist
        integer,      intent(out)   :: ist1, jst1, kst1 !< station location indices
        real(SP),     intent(out)   :: xst1, yst1, zst1 !< station location coordinate
        real(SP),     intent(inout) :: stlo1, stla1
        !! --
        integer :: i
        !! ----
    
        is_exist = .false.
        do i=1, nst
            if (trim(stnm1) == trim(stnm(i))) then
    
                ist1     = ist(i)
                jst1     = jst(i)
                kst1     = kst(i)
                xst1     = xst(i)
                yst1     = yst(i)
                zst1     = zst(i)
                stlo1    = stlo(i)
                stla1    = stla(i)
                is_exist = .true.
                exit
            end if
        end do
    
    end subroutine wav__stquery

    subroutine wav__write()

        integer :: i, j
        character(6) :: cid
        !! ----

        call pwatch__on("wav__write")

        if (nst == 0) then
            call pwatch__off("wav__write")
            return
        end if

        if (wav_format == 'sac') then

            do i=1, nst
      
                if (sw_wav_v ) then
                    do j=1, 3
                        call export_wav__sac(sh_vel(j,i), wav_vel(:,j,i))
                    end do
                end if
        
                if( sw_wav_u ) then
                    do j=1, 3
                        call export_wav__sac(sh_disp(j,i), wav_disp(:,j,i))
                    end do
                end if
      
                if( sw_wav_stress ) then
                    do j=1, 6
                        call export_wav__sac(sh_stress(j,i), wav_stress(:,j,i))
                    end do
                end if
      
                if( sw_wav_strain ) then
                    do j=1, 6
                        call export_wav__sac(sh_strain(j,i), wav_strain(:,j,i))
                    end do
                end if
              
            end do
      
        else if ( wav_format == 'csf' ) then

            write(cid,'(I6.6)') myid
      
            if( sw_wav_v )      call export_wav__csf(nst, 3, sh_vel,    wav_vel)
            if( sw_wav_u )      call export_wav__csf(nst, 3, sh_disp,   wav_disp)
            if( sw_wav_stress ) call export_wav__csf(nst, 6, sh_stress, wav_stress)
            if( sw_wav_strain ) call export_wav__csf(nst, 6, sh_strain, wav_strain)
      
        end if 

        call pwatch__off("wav__write")

    end subroutine wav__write
    
    subroutine export_wav__sac( sh, dat )

        type(sac__hdr), intent(in) :: sh
        !! --
        real(SP), intent(in) :: dat(:)
        character(256) :: fn
        !! ----
        
        fn = trim(odir) // '/wav/' // trim(title) // '.' // trim(sh%kstnm) // '.' // trim(sh%kcmpnm) // '.sac'
        call sac__write(fn, sh, dat, .true.)
        
    end subroutine export_wav__sac
  
    subroutine export_wav__csf(nst, ncmp, sh, dat)
        
        integer, intent(in) :: nst, ncmp
        type(sac__hdr), intent(in) :: sh(ncmp, nst)
        real(SP), intent(in) :: dat(ntw, ncmp, nst)
        !! --
        character(5) :: cid
        character(256) :: fn
        !! ----

        write(cid,'(I5.5)') myid
        fn = trim(odir) // '/wav/' // trim(title) // '__' // cid // '__.csf'
        call csf__write(fn, nst*ncmp, ntw, reshape(sh,(/ncmp*nst/)), reshape(dat, (/ntw, ncmp*nst/)))
  
    end subroutine export_wav__csf
          
end module m_wav
