# Execution

To run the program, the MPI program is necessary, such that

```bash
$ mpirun -np ${NP} ./bin/swpc_3d.x -i ${input}
```

where `${NP}` is the number of MPI processes and `${input}` is the name
of the input file. Note that the `mpirun` command may be different for
different computational systems.

A file `example/input.inf` describes 2x2=4 parallel 3D example simulation in a homogeneous medium. To run the 3D simulation with this, try

```bash
$ mpirun -np 4 ./bin/swpc_3d.x -i example/input.inf
```

Then, one will have the standard error output as follows. This result is contained in `example/example.out` file.  

```text

 ------------------------------------------------------------------------------
  SWPC_3D                                                                     
 ------------------------------------------------------------------------------

  Grid Size               :      384 x    384 x    384
  MPI Partitioning        :        2 x    2
  Total Memory Size       :          13.292  [GiB]
  Node Memory Size        :           3.323  [GiB]
  Stability  Condition c  :           0.645  (c<1)
  Wavelength Condition r  :          12.488  (r>5-10)
  Minimum velocity        :           3.122  [km/s]
  Maximum velocity        :           7.977  [km/s]
  Maximum frequency       :           0.500  [Hz]

 ------------------------------------------------------------------------------

  it=0000050, 1.024 s/loop, eta 000:16:12, ( 2.94E-05  2.94E-05  1.59E-04 )
  it=0000100, 1.028 s/loop, eta 000:15:25, ( 1.95E-04  1.95E-04  8.87E-04 )
  it=0000150, 1.023 s/loop, eta 000:14:29, ( 1.04E-04  1.04E-04  4.87E-04 )
  it=0000200, 1.014 s/loop, eta 000:13:31, ( 2.62E-05  2.62E-05  4.89E-05 )
  it=0000250, 1.010 s/loop, eta 000:12:37, ( 1.32E-05  1.32E-05  3.59E-05 )
  it=0000300, 1.001 s/loop, eta 000:11:40, ( 1.13E-05  1.13E-05  3.16E-05 )
  it=0000350, 0.996 s/loop, eta 000:10:47, ( 1.26E-05  1.26E-05  2.67E-05 )
  it=0000400, 0.992 s/loop, eta 000:09:55, ( 1.13E-05  1.13E-05  2.33E-05 )
  it=0000450, 0.991 s/loop, eta 000:09:05, ( 9.68E-06  9.68E-06  2.07E-05 )
  it=0000500, 0.987 s/loop, eta 000:08:13, ( 9.22E-06  9.22E-06  1.98E-05 )
  it=0000550, 0.984 s/loop, eta 000:07:22, ( 8.90E-06  8.90E-06  1.90E-05 )
  it=0000600, 0.981 s/loop, eta 000:06:32, ( 8.37E-06  8.37E-06  1.78E-05 )
  it=0000650, 0.978 s/loop, eta 000:05:42, ( 7.78E-06  7.78E-06  1.68E-05 )
  it=0000700, 0.975 s/loop, eta 000:04:52, ( 7.63E-06  7.63E-06  1.59E-05 )
  it=0000750, 0.973 s/loop, eta 000:04:03, ( 7.47E-06  7.47E-06  1.53E-05 )
  it=0000800, 0.972 s/loop, eta 000:03:14, ( 7.01E-06  7.01E-06  1.45E-05 )
  it=0000850, 0.974 s/loop, eta 000:02:26, ( 6.79E-06  6.79E-06  1.35E-05 )
  it=0000900, 0.974 s/loop, eta 000:01:37, ( 7.07E-06  7.07E-06  1.28E-05 )
  it=0000950, 0.973 s/loop, eta 000:00:48, ( 7.38E-06  7.38E-06  1.20E-05 )
  it=0001000, 0.973 s/loop, eta 000:00:00, ( 7.22E-06  7.22E-06  1.14E-05 )

 ------------------------------------------------------------------------------
 
  Total time             :         972.827 s

 ------------------------------------------------------------------------------```
```

The first part of the message contains information such as the estimated
memory usage, stability condition, and wavelength condition. As shown in
the above example, a stability condition of $c<1$ is mandatory to
execute; if the specified parameter violates this condition, the program
aborts immediately. In addition, the wavelength condition (the ratio of
spatial grid size and minimum wavelength) is recommended to $r>5-10$.
During the computation, the computational speed, remaining time (eta;
estimated time of arrival), and maximum velocity amplitude of the
components are shown. Among them, numbers related to the computational speed may vary depending on the environment. If the computation is successful, the maximum amplitude (the three numbers in the last `()`) should be the same.

The above simulation requires around ~14 GB of computer memory. If this size of the simulation is difficult, try the 2D P-SV simulation with the same input parameter file along the $xz$ cross section: 

```bash
$ mpirun -np 2 ./bin/swpc_psv.x -i example/input.inf
```

Note that the number of parallelization (`-np` option) is different between 2D and 3D even for the same input parameter file. 

The output in this case is as follows. The result is almost the same as the 3D code, only the maximum amplitude is changed from 3 components to 2 components.


```text
------------------------------------------------------------------------------
  SWPC_PSV                                                                    
 ------------------------------------------------------------------------------


  Grid Size               :      384 x    384
  MPI Partitioning        :               2
  Total Memory Size       :           0.020  [GiB]
  Node Memory Size        :           0.010  [GiB]
  Stability  Condition c  :           0.526  (c<1)
  Wavelength Condition r  :          12.488  (r>5-10)
  Minimum velocity        :           3.122  [km/s]
  Maximum velocity        :           7.977  [km/s]
  Maximum frequency       :           0.500  [Hz]

 ------------------------------------------------------------------------------

  it=0000050, 0.003 s/loop, eta 000:00:02, ( 7.40E-02  2.73E-01 )
  it=0000100, 0.003 s/loop, eta 000:00:02, ( 8.13E-01  2.01E+00 )
  it=0000150, 0.003 s/loop, eta 000:00:02, ( 7.71E-01  1.46E+00 )
  it=0000200, 0.003 s/loop, eta 000:00:02, ( 2.90E-01  2.56E-01 )
  it=0000250, 0.003 s/loop, eta 000:00:01, ( 2.21E-01  3.03E-01 )
  it=0000300, 0.003 s/loop, eta 000:00:01, ( 1.90E-01  3.22E-01 )
  it=0000350, 0.003 s/loop, eta 000:00:01, ( 1.76E-01  3.27E-01 )
  it=0000400, 0.003 s/loop, eta 000:00:01, ( 1.93E-01  3.23E-01 )
  it=0000450, 0.003 s/loop, eta 000:00:01, ( 1.92E-01  3.08E-01 )
  it=0000500, 0.003 s/loop, eta 000:00:01, ( 1.89E-01  3.19E-01 )
  it=0000550, 0.003 s/loop, eta 000:00:01, ( 1.83E-01  3.33E-01 )
  it=0000600, 0.003 s/loop, eta 000:00:01, ( 1.75E-01  3.34E-01 )
  it=0000650, 0.003 s/loop, eta 000:00:00, ( 1.70E-01  3.36E-01 )
  it=0000700, 0.003 s/loop, eta 000:00:00, ( 1.69E-01  3.37E-01 )
  it=0000750, 0.003 s/loop, eta 000:00:00, ( 1.73E-01  3.48E-01 )
  it=0000800, 0.003 s/loop, eta 000:00:00, ( 1.76E-01  3.53E-01 )
  it=0000850, 0.003 s/loop, eta 000:00:00, ( 1.77E-01  3.45E-01 )
  it=0000900, 0.003 s/loop, eta 000:00:00, ( 1.81E-01  3.41E-01 )
  it=0000950, 0.003 s/loop, eta 000:00:00, ( 1.77E-01  3.39E-01 )
  it=0001000, 0.003 s/loop, eta 000:00:00, ( 1.63E-01  3.36E-01 )

 ------------------------------------------------------------------------------
 
  Total time             :           2.652 s

 ------------------------------------------------------------------------------
```

The computation speed depends on the environment, but in general, the computation time is much shorter for 2D calculations than for 3D calculations.

Note that the message 
```text
Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG
```
may appear after the calculation is completed, but this does not affect the calculation results.

