#include "../shared/m_debug.h"
program swpc_3d

    !! SWPC: Seismic Wave Propagation Code
    !!
    !! This software simulate seismic wave propagation, by solving equations of motion with constitutive equations of
    !! elastic/visco-elastic medium by finite difference method (FDM).
    !!
    !! Copyright 2013-2023 Takuto Maeda. All rights reserved. This project is released under the MIT license.
    !!
    !! #### References
    !!   - Noguchi et al.     (2016) GJI     doi:10.1093/gji/ggw074
    !!   - Maeda et al.       (2014) EPSL    doi:10.1016/j.epsl.2014.04.037
    !!   - Maeda et al.       (2013) BSSA    doi:10.1785/0120120118
    !!   - Maeda and Furumura (2013) PAGEOPH doi:10.1007/s00024-011-0430-z
    !!   - Noguchi et al.     (2013) PAGEOPH doi:10.1007/s00024-011-0412-1
    !!   - Furumura et al.    (2008) PAGEOPH doi:10.1007/s00024-008-0318-8
    !!   - Furumura ad Chen   (2005) PARCO   doi:10.1016/j.parco.2005.02.003

    use m_std
    use m_debug
    use m_global
    use m_kernel
    use m_getopt
    use m_source
    use m_medium
    use m_report
    use m_pwatch
    use m_snap
    use m_absorb
    use m_green
    use m_readini
    use m_version
    use m_wav
    use mpi
 
    implicit none

    character(256) :: fn_prm
    integer :: it
    integer :: ierr
    logical :: is_opt
    logical :: stopwatch_mode
    integer :: io_prm, io_watch
    logical :: strict_mode

    !! Version
    call getopt('v', is_opt); if (is_opt) call version__display('swpc_3d')
    call getopt('-version', is_opt); if (is_opt) call version__display('swpc_3d')

    !! Launch MPI process
    call mpi_init(ierr)

    !! option processing
    call getopt('i', is_opt, fn_prm, './in/input.inf')
    open (newunit=io_prm, file=trim(fn_prm), action='read', status='old')

    call readini(io_prm, 'stopwatch_mode', stopwatch_mode, .true.)

    call readini(io_prm, 'strict_mode', strict_mode, .false.)
    call readini__strict_mode(strict_mode)

    !! Read control parameters
    call global__setup(io_prm)

    !! stopwatch start
    call pwatch__setup(stopwatch_mode)

    !! set-up each module
    call global__setup2()
    call medium__setup(io_prm)
    call mpi_barrier(mpi_comm_world, ierr) !! wait until deallocation by medium__setup
    call kernel__setup()
    call source__setup(io_prm)
    call absorb__setup(io_prm)
    call snap__setup(io_prm)
    call wav__setup(io_prm)
    call green__setup(io_prm)

    !$acc enter data copyin(&
    !$acc       Vx (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Vy (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Vz (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Sxx(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Syy(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Szz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Syz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Sxz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Sxy(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       Rxx(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), &
    !$acc       Ryy(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), &
    !$acc       Rzz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), &
    !$acc       Ryz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), &
    !$acc       Rxz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), &
    !$acc       Rxy(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k)  )
    
    !$acc enter data copyin(&    
    !$acc       rho (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       bx  (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), & 
    !$acc       by  (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), & 
    !$acc       bz  (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), & 
    !$acc       lam (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       mu  (kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), & 
    !$acc       muyz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       muxz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       muxy(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       taup(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       taus(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), &   
    !$acc       ts(nm) )

    !$acc enter data copyin(&
    !$acc       kfs(ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       kob(ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       kfs_top(ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       kfs_bot(ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       kob_top(ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       kob_bot(ibeg_m:iend_m, jbeg_m:jend_m), &
    !$acc       bddep(ibeg_m:iend_m, jbeg_m:jend_m, 0:NBD), &
    !$acc       kbeg_a(ibeg_m:iend_m, jbeg_m:jend_m))

    call report__setup(io_prm)
    close (io_prm)

    !! mainloop
    do it = 1, nt
        call report__progress(it)

        call green__store(it)
        call wav__store(it)
        call snap__write(it)

        call kernel__update_stress()
        call absorb__update_stress()

        call source__stressglut(it)
        call global__comm_stress()

        call kernel__update_vel()
        call absorb__update_vel()

        call source__bodyforce(it)

        call green__source(it)
        call global__comm_vel()
    end do

    call green__export()
    call wav__write()
    call snap__closefiles()

    !! ending message
    call report__terminate()

    !! stopwatch report from 0-th node
    if (stopwatch_mode) then
        if (myid == 0) then
            open (newunit=io_watch, file=trim(odir)//'/'//trim(title)//'.tim', action='write', status='unknown')
        end if
        call pwatch__report(io_watch, 0)
    end if

    !! Program termination
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)

end program swpc_3d
