# Checkpointing and Restarting

Some large-scale computers limit the computational time of a single job.
To achieve long-duration computation, `OpenSWPC` can export all memory
contents to files at specific times (checkpointing), and then continue
the simulation as another job (restarting).

If this function is turned on, `OpenSWPC` will terminate the computation
after an elapsed time of `ckp_time` (in seconds) and will export all
memory images.

For the next job, `OpenSWPC` first tries to find the directory cdir to
locate the checkpointing file. If there are checkpointing files,
`OpenSWPC` reads them to continue the simulation. Otherwise, `OpenSWPC`
starts the simulation from scratch.

After finishing the computation of all time steps, `OpenSWPC` removes
most of the contents of the checkpointing files. However, it does not
delete the checkpointing files. This is to avoid unexpectedly starting
the computation from the beginning and overwriting the output files.

This function is only available for the three-dimensional simulation
code (`swpc_3d.x`).

!!! Info "Parameters"
    **`is_ckp`**
    : The flag to use checkpointing/restarting.
    
    **`cdir`**
    : Output directory name of the checkpointing file. At restart, the
    checkpointing files are assumed to be in this directory.

    **`ckp_time`**
    : Checkpointing time in seconds.

    **`ckp_interval`**
    : Investigate if the computation time exceeds `ckp_time` periodically
    at this interval. Setting this interval step as too small may affect
    the performance of the computation.
