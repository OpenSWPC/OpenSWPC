# Execution

To run the program, the MPI program is necessary, such that

```bash
$ mpirun -np ${NP} ./bin/swpc_3d.x -i ${input}
```

where `${NP}` is the number of MPI processes and `${input}` is the name
of the input file. Note that the `mpirun` command may be different for
different computational systems.

If the program runs properly, the following message will appear in the
standard error output. The result may be slightly different for
different programs (3D/P-SV/SH) or execution modes.

```text
 ------------------------------------------------------------------------------
  SWPC_3D (benchmark mode)
 ------------------------------------------------------------------------------

  Grid Size               :      384 x  384 x  384
  MPI Partitioning        :        4 x    6
  Total Memory Size       :          12.705  [GiB]
  Node Memory Size        :           0.529  [GiB]
  Stability  Condition c  :           0.980  (c<1)
  Wavelength Condition r  :           7.000  (r>5-10)
  Minimum velocity        :           3.500  [km/s]
  Maximum velocity        :           6.062  [km/s]
  Maximum frequency       :           1.000  [Hz]

 ------------------------------------------------------------------------------

  it=0000050, 1.877 s/loop, eta 000:29:43, ( 5.00E-05  5.00E-05  4.96E-05 )
  it=0000100, 1.887 s/loop, eta 000:28:18, ( 1.75E-05  1.75E-05  1.05E-05 )
  it=0000150, 1.932 s/loop, eta 000:27:22, ( 1.02E-05  1.02E-05  5.41E-06 )
  it=0000200, 1.943 s/loop, eta 000:25:54, ( 6.59E-06  6.59E-06  4.35E-06 )
          
           ・
           ・
           ・

  it=0000950, 1.986 s/loop, eta 000:01:39, ( 4.89E-07  4.89E-07  1.81E-06 )
  it=0001000, 1.982 s/loop, eta 000:00:00, ( 1.65E-07  1.65E-07  1.54E-07 )

 ------------------------------------------------------------------------------

  Total time             :        1982.348 s

 ------------------------------------------------------------------------------  
```

The first part of the message contains information such as the estimated
memory usage, stability condition, and wavelength condition. As shown in
the above example, a stability condition of $c<1$ is mandatory to
execute; if the specified parameter violates this condition, the program
aborts immediately. In addition, the wavelength condition (the ratio of
spatial grid size and minimum wavelength) is recommended to $r>5-10$.
During the computation, the computation speed, remaining time (eta;
estimated time of arrival), and maximum velocity amplitude of the
components are shown.

