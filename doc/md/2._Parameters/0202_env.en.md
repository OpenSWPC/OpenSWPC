# Controlling Parameters

!!! Info "Parameters"
    

    **`title`**
    : Title of the computation to be used for the output filename.
    
    **`odir`**
    : Name of output directory. This is a relative directory path from the
    location of the program execution. If this directory does not exist
    at the time of run, `OpenSWPC` will automatically create it.
      
    **`ntdec_r`**
    : **N**umber of **T**ime-step **DEC**imation factors for screen
    `R`eporting. The maximum amplitudes of the velocity components are
    reported in the standard error output every `ntdec_r` steps. This
    screen output is generally used to confirm that the model is working
    correctly. A cycle that is too short (this parameter is too small)
    may slow down the computation.

    **`stopwatch_mode`**
    : Measure the computation times at major subroutines and export the accumulated times to `(odir)/(title).tim`. This function is used for benchmarking and perforance tuning. The output `*.tim` file can be visualized with GMT (versions 4 or 5) via `tools/timvis.gmt?`. 

    **`benchmark_mode`**
    : If this flag is `.true.`, the fixed homogeneous medium and
    single-point moment tensor source will be selected irrespective of
    the parameter specification. This is used for validation and
    performance measurements.

    **`strict_mode`** **(New in v5.1)**
    : If this flag is `.true.`, the all necessary parameter must be explicitly written in the input file. 
    By default this parameter is set to `.false.`, so that `OpenSWPC` will use parameter-specific default values if the parameter is not found in the input file. 

    **`ipad, jpad, kpad`**
    : Expand the Fortran array sizes along the x-, y-, and z-directions.
    In some computer architectures, the computation speed is very
    sensitive to the array size. In such cases, slightly changing the
    array size using these parameters may improve the performance. The
    expanded array will not be used for the simulation. Therefore, the
    simulation result is not affected by changing this option.


!!! Caution "Benchmark Mode"

    Until version 5.0.2, an example input file `example/input.inf` contained `benchmark_mode = .true.` which neglects many changes of parameters related to velocity structure and earthquake source. After version 5.1, this parameter is set to `.false.` in the example file. 