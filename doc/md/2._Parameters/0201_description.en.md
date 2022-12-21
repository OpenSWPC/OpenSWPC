# Notation of the Parameter File

In the parameter files, one parameter is defined on each line in the
following format.

```Fortran
variable_name = value
```

The description of the values should follow Fortran notation. For
example, logical (Boolean) values are denoted as `.true.` or `.false.`.

Lines that do not contain an equal sign (`=`) will be neglected; in
addition, lines starting with `!` or `#` are regarded as comment lines
and will be skipped. Comments can be followed by variable definitions,
that is, comments can be written on the same line as the parameter
definition. For example, the following parameter line will work without
errors.

```Fortran
nx = 1024 ! number of grids
```

The order of the parameter definition can be changed freely. If a parameter is not specified, `OpenSWPC` may use a pre-defined default variable if the missing parameter is not critical for executing computation. 
Note that there are parameters that must be defined explicitly. Multiple definitions of the same parameters in a single parameter file are not recommended; however, in this case, the first definition may be adopted. 
It is acceptable to leave blanks before and after the equal sign; however, **it is not permitted to have a blank space between the minus character and succeeding numbers** (e.g., `’- 35.0’` is not allowed). 
It is recommended to **use quotation marks around string parameters**. 
Without them, the directory path character ('/') may be unexpectedly interpreted as a termination of the string parameter.


!!! Caution 
    If a non-essential parameter is not given in the parameter file, an information message will be displayed in the standard error output as an example below: 
    ```txt
    [info] key XXXX is not found. 
    [info]     Use default value YYYY instead.
    ```

## An Example Parameter File

The following is a full set of example parameters. In the following
sections, detailed descriptions of each parameter will be given.

```Fortran

  !! ----------------------------------------------------------------------- !!
  !!
  !!  SWPC input file
  !!
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  !! Control
  !!

  title            = 'swpc'           !! exe title: used for output filenames
  odir             = './out'          !! output directory
  ntdec_r          = 50               !! screen report timing (1/cycle)
  strict_mode      = .false.          !! all parameters to be explicitly definied

  !! ----------------------------------------------------------------------- !!
  !! Model/Grid Size and Area
  !!

  nproc_x          = 2                !! parallelization in x-dir
  nproc_y          = 2                !! parallelization in x-dir
  nx               = 384              !! total grid number in x-dir
  ny               = 384              !! total grid number in y-dir
  nz               = 384              !! total grid number in z-dir
  nt               = 1000             !! time step number

  dx               = 0.5              !! grid width in x-dir
  dy               = 0.5              !! grid width in y-dir
  dz               = 0.5              !! grid width in z-dir
  dt               = 0.02             !! time step width

  vcut             = 1.5              !! minimum velocity
                                      !- smaller velocities will be increased

  xbeg             = -96.0            !! minimum in x-dir
  ybeg             = -96.0            !! minimum in y-dir
  zbeg             = -10.0            !! minimum in z-dir
  tbeg             = 0.0              !! start time

  clon             = 139.7604         !! center longitude
  clat             = 35.7182          !! center latitude
  phi              = 0.0              !! horizontal coordinate rotation
                                      !- measured clockwise from the north

  fq_min           = 0.02             !! minimum freq. for Q-constant model
  fq_max           = 2.00             !! maximum freq. for Q-constant model
  fq_ref           = 1.0              !! ref. freq. for physical dispersion

  !! ----------------------------------------------------------------------- !!
  !! Snapshot Output
  !!

  snp_format       = 'netcdf'         !! snapshot format (netcdf)

  xy_ps%sw         = .false.          !! P&S amp. for xy section
  xz_ps%sw         = .true.           !! P&S amp. for xz section
  yz_ps%sw         = .false.          !! P&S amp. for yz section
  fs_ps%sw         = .false.          !! P&S amp. for free surface
  ob_ps%sw         = .true.           !! P&S amp. for ocean bottom

  xy_v%sw          = .false.          !! 3-comp. velocity for xy section
  xz_v%sw          = .true.           !! 3-comp. velocity for xz section
  yz_v%sw          = .false.          !! 3-comp. velocity for yz section
  fs_v%sw          = .false.          !! 3-comp. velocity for free surface
  ob_v%sw          = .true.           !! 3-comp. velocity for ocean bottom

  xy_u%sw          = .false.          !! 3-comp. disp. for xy section
  xz_u%sw          = .true.           !! 3-comp. disp. for xz section
  yz_u%sw          = .false.          !! 3-comp. disp. for yz section
  fs_u%sw          = .false.          !! 3-comp. disp. for free surface
  ob_u%sw          = .true.           !! 3-comp. disp. for ocean bottom


  z0_xy            =  7.0             !! depth for xy cross section
  x0_yz            =  0.0             !! x-value for yz cross section
  y0_xz            =  0.0             !! y-value for xz cross section

  ntdec_s          = 5                !! time decimation of snapshot
                                      !- (specify 1 for no decimation)

  idec             = 2                !! x-decimation for snapshot
  jdec             = 2                !! y-decimation for snapshot
  kdec             = 2                !! z-decimation for snapshot

  !! ----------------------------------------------------------------------- !!
  !! Waveform Output
  !!

  sw_wav_v         = .true.           !! velocity trace output at stations
  sw_wav_u         = .false.          !! displacement trace output at stations
  sw_wav_stress    = .false.           !! stress tensor trace
  sw_wav_strain    = .false.           !! strain tansor trace
  ntdec_w          = 5                !! time decimation of waveform output
  st_format        = 'xy'             !! station format: 'xy' or 'll'
  fn_stloc         = './example/stloc.xy'  !! station location file
  wav_format       = 'sac'            !! 'sac' or 'csf' ('sac' recommended)
  wav_calc_dist    = .false.          !! Calculate epicentral distance
  
  !! ----------------------------------------------------------------------- !!
  !! Earthquake Source
  !!

  !! Moment tensor source format:
  !!   xymoij / xym0dc / llm0ij / llm0dc / xymwij / xymwdc / llmwij / llmwdc
  !! Body force source format:
  !!   xy or ll
  stf_format       = 'xym0ij'

  !! Basis source time function
  !! 'boxcar' / 'triangle' / 'herrmann' / 'kupper' / 'cosine' / 'texp'
  stftype          = 'kupper'

  fn_stf           = "./example/source.dat"   !! Source grid file name

  !! source depth correction
  !! 'asis':use z value, 'bd{i}': i-th boundary (i=0...9)
  sdep_fit         = 'asis'

    !! --------------------------------------------------------------------- !!
    !! Body force source mode
    !!
    bf_mode          = .false.


    !! --------------------------------------------------------------------- !!
    !! Plane wave source mode
    !!
    pw_mode          = .false.   !! plane wave input. neglects fn_stf
    pw_ztop          = 100.      !! top z-coordinate of the initial plane wave
    pw_zlen          = 30.       !! wavelength of the initial plane wave
    pw_ps            = 'p'       !! 'p' P-wave 's' S-wave
    pw_strike        = 0.0       !! strike direction of plane wave (deg.)
    pw_dip           = 0.0       !! dip of plane wave (deg.)
    pw_rake          = 0.0       !! rake of plane S-wave polarization (deg.)

  !! ----------------------------------------------------------------------- !!
  !! Absorbing Boundary Condition
  !!

  abc_type         = 'pml'            !! 'pml' or 'cerjan'
  na               = 20               !! absorbing layer thickness
  stabilize_pml    = .false.           !! avoid low-v layer in PML region

  !! ----------------------------------------------------------------------- !!
  !! Velocity model
  !!

  vmodel_type      = 'lhm'            !! velocity model type 'uni'/'grd'/'lhm'
  is_ocean         = .true.           !! topography z<0 is covered by ocean
  topo_flatten     = .false.          !! Force topography variation to zero (formerly is_flatten)
  munk_profile     = .true.           !! velocity gradient inside the seawater column
  earth_flattening = .false.          !! Earth-flattening tranformation
  
    !! --------------------------------------------------------------------- !!
    !! For uniform velocity model 'uni'
    !!
    vp0              = 5.0              !! P-wave velocity [km/s]
    vs0              = 3.0              !! S-wave velocity [km/s]
    rho0             = 2.7              !! mass density    [g/cm^3]
    qp0              = 200              !! Qp
    qs0              = 200              !! Qs
    topo0            = 0                !! topography location

    !! --------------------------------------------------------------------- !!
    !! For GMT grid file input 'grd' ( requires netcdf library )
    !!
    dir_grd          = '${DATASET}/vmodel/ejivsm/'    !! directory for grd file
    fn_grdlst        = './example/grd.lst'            !! grd file list
    node_grd         = 0                              !! input MPI node

    !! --------------------------------------------------------------------- !!
    !! For layered homogeneous medium model ('lhm')
    !!
    fn_lhm           = 'example/lhm.dat'    !! 1D velocity structure

    !! --------------------------------------------------------------------- !!
    !! For random medium models
    !!
    dir_rmed         = './in/'             !! location of random medium file
    fn_grdlst_rmed   = './example/grd.lst' !! grd file list
    rhomin           = 1.0                 !! minimum density threshold
    fn_rmed0         = 'dummy.ns'          !! vel. purturb. on a uniform media

  !! ----------------------------------------------------------------------- !!
  !! Checkpoint/Restart
  !!
  is_ckp           = .false.          !! perform checkpoint/restart
  ckpdir           = './out/ckp'      !! output directory
  ckp_interval     = 1000000          !! interval for checkpoint check（1/cycle）
  ckp_time         = 1000000.         !! checkpoint time
  ckp_seq          = .true.           !! sequential output mode

  !! ----------------------------------------------------------------------- !!
  !! Reciprocity Green's Function Mode
  !!
  green_mode       = .false.          !! reciprocity Green's function mode
  green_stnm       = 'st01'           !! virtual station name from fn_stlst
  green_cmp        = 'z'              !! virtual source direction 'x', 'y', 'z'
  green_trise      = 1.0              !! rise time
  green_bforce     = .false.          !! also calc. body force Green's function
  green_maxdist    = 550.             !! horizontal limit of source grid
  green_fmt        = 'llz'            !! list file format: 'xyz' or 'llz'
  fn_glst          = 'example/green.lst'   !! Green's function grid point list

  !! ----------------------------------------------------------------------- !!
  !! MISC
  !!

  stopwatch_mode   = .false.          !! measure computation time at routines
  benchmark_mode   = .false.          !! benchmark mode

  ipad             = 0                !! memory padding size for tuning
  jpad             = 0                !! memory padding size for tuning
  kpad             = 0                !! memory padding size for tuning
```