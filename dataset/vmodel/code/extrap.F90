!! smoothly extrapolate JIVSM data to rectangular region

program extrap

  use m_std
  use m_system
  implicit none

  integer :: N_NS
  integer :: N_EW

  real(DP), parameter :: NOT_A_NUMBER = -99999.0_DP
  real(DP), allocatable :: lon(:,:), lat(:,:), dep(:,:), dep2(:,:)
  character(99) :: fn_in
  real(DP) :: dtmp
  integer :: i, j
  integer :: rh, bh
  logical, allocatable :: isdat(:,:)
  real :: dlon, dlat, lon_beg, lon_end, lat_beg, lat_end
  integer :: A_E, A_S, A_N, A_W, B_E, B_S, B_N, B_W
  real(DP) :: w1, w2, l1, l2
  character(99) :: fn_out
  real(DP) :: dep_ave1, dep_ave2
  integer :: ii, jj, k
  integer :: di, dj
  integer :: imode !< 1: interpolate both, 2: se-only

  call system__getarg(1, fn_in)
  call system__getarg(2, fn_out)
  call system__getarg(3, dlon)
  call system__getarg(4, dlat)
  call system__getarg(5, lon_beg)
  call system__getarg(6, lon_end)
  call system__getarg(7, lat_beg)
  call system__getarg(8, lat_end)
  call system__getarg(9, imode)

  n_ns = floor( ( lat_end - lat_beg ) / dlat + 0.5 ) + 1
  n_ew = floor( ( lon_end - lon_beg ) / dlon + 0.5 ) + 1
  write(STDERR,*) n_ns, n_ew
  allocate(lon(n_ns,n_ew), lat(n_ns,n_ew), dep(n_ns,n_ew), dep2(n_ns,n_ew))
  allocate(isdat(n_ns,n_ew))
  open(10,file=fn_in,access='stream')

  ! grd2xyzのデータ出力方向に合わせた順番で読み込み
  do i=N_NS, 1, -1
     do j=1, N_EW, 1
        read(10) lon(i,j), lat(i,j), dtmp

        ! Not a Numberの検知
        !        if( dtmp == dtmp ) then
        if( dtmp > -100000. ) then
           dep(i,j) = dtmp
        else
           dep(i,j) = NOT_A_NUMBER
        end if
     end do
  end do
  write(STDERR,*) dep(N_NS/2, N_EW/2 )
  write(STDERR,*) maxval(dep), minval(dep)
  close(10)

  !! detection of edges of NaN-region A
  isdat(:,:) = .true.
  do i=1, N_EW
     if( abs( dep(N_NS-1,i) - NOT_A_NUMBER ) > epsilon(1.0) ) then
        exit
     end if
     A_E = i
  end do
  do i=N_NS,1,-1
     if( abs( dep(i,2) - NOT_A_NUMBER ) > epsilon(1.0) ) then
        exit
     end if
     A_S = i
  end do
  A_N = N_NS-1
  A_W = 2
  isdat( A_S:A_N,A_W:A_E ) = .false.

  write(STDERR,*) A_S, A_N, A_E, A_W
  !! same for NaN-region B
  do i=N_EW, 1, -1
     if( abs( dep(2,i) -  NOT_A_NUMBER ) > epsilon(1.0) ) then
        exit
     end if
     B_W = i
  end do
  do i=1, N_NS, 1
     if( abs( dep(i,N_EW-1) - NOT_A_NUMBER ) > epsilon(1.0) ) then
        exit
     end if
     B_N = i
  end do
  B_E = N_EW-1
  B_S = 2
  isdat( B_S:B_N,B_W:B_E ) = .false.

  write(STDERR,*) B_S, B_N, B_E, B_W



  !! interpolation: weighted sum
  ! region A
  if( imode == 1 ) then
    do j=A_W, A_E
      do i=A_S, A_N

        !! distance from edge
        l1 = dble( i - A_S + 1 ) / dble( A_N-A_S+1 )
        l2 = dble( A_E - j + 1 ) / dble( A_E-A_W+1 )
        !! weight
        w1 = l2 / ( l1 + l2 )
        w2 = l1 / ( l1 + l2 )

        !! averaged depth along the edges
        di = i-A_S + 1
        dj = A_E - j + 1

        dep_ave1 = 0
        k = 0
        do jj=max(j-di,A_W), min(j+di,A_E)
           dep_ave1 = dep_ave1 + dep(A_S-1,jj)
           k = k + 1
        end do
        dep_ave1 = dep_ave1 / k

        dep_ave2 = 0
        k = 0
        do ii=max(i-dj,A_S), min(i+dj,A_N)
           dep_ave2 = dep_ave2 + dep(ii,A_E+1)
           k = k + 1
        end do
        dep_ave2 = dep_ave2 / k

        dep(i,j) = w1 * dep_ave1 + w2 * dep_ave2
      end do
    end do
  else

    do j=A_W, A_E
      do i=A_S, A_N
        dep(i,j) = dep(i,A_E+1)
      end do
    end do

  end if
  ! region B
  do j=B_W, B_E
     do i=B_S, B_N
        l1 = dble( B_N - i + 1 ) / dble( B_N - B_S + 1 )
        l2 = dble( j - B_W + 1 ) / dble( B_E - B_W + 1 )
        w1 = l2 / ( l1 + l2 )
        w2 = l1 / ( l1 + l2 )

        di = B_N - i + 1
        dj = j - B_W + 1


        dep_ave1 = 0
        k = 0
        do jj=max(j-di,B_W), min(j+di,B_E)
           dep_ave1 = dep_ave1 + dep(B_N+1,jj)
           k = k + 1
        end do
        dep_ave1 = dep_ave1 / k

        dep_ave2 = 0
        k = 0
        do ii=max(i-dj,B_S), min(i+dj,B_N)
           dep_ave2 = dep_ave2 + dep(ii,B_W-1)
           k = k + 1
        end do
        dep_ave2 = dep_ave2 / k
        dep(i,j) = w1 * dep_ave1 + w2 * dep_ave2
     end do
  end do

  !! edge interpolation
  dep( :   ,1    ) = dep(:,2)
  dep( :   ,N_EW ) = dep(:,N_EW-1)
  dep(    1,:    ) = dep(2,:)
  dep( N_NS,:    ) = dep(N_NS-1,:)



  open(11,file=trim(fn_out),access='stream')
  do j=1, N_EW
     do i=1, N_NS
        write(11) lon(i,j), lat(i,j), dep(i,j)
     end do
  end do
  close(11)

end program extrap
