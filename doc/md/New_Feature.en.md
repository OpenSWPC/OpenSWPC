# New features

## Version 5.3.0 (2023-02-02)

### Better parallel partitioning

In previous versions, automatic MPI area allocation sometimes failed when the number of grids in the X or Y direction was not divisible by the number of MPI partitions and the number of MPI partitions was very large. 

This problem has been fixed in version 5.3.0. 

In new version, the partitioning algorithm is as follows: 

Let $N$ be the number of grid, and $P$ the number of MPI partitions. If $N$ is divisible by $P$, i.e., $\mod(N,P)= 0$, the number of grids assigned to a node is $N_P = N/P$. If this is not the case, the following rule is applied: 

| Node ID | Number of Grids $N_P$ |
| ------- | -------------------- |
| 0 | $N_P = (N-M)/P$ |
| 1 to $M$ | $N_P = (N-M)/P + 1$ |
| $M+1$ to $P-1$ | $N_P = (N-M)/P$ |

(where $M = \mod(N,P)$)

#### Example

![](./fig/5.3.0_partition.png)


### Python integration

An example of processing OpenSWPC input/output in Python is included in [this manual](./3._Tools/0305_python.en.md).

### Updated documentation
#### Try OpenSWPC on cloud!

See [this example](./1._SetUp/0100_trial.en.md). This is also would be a nice guide to compile the OpenSWPC in Ubuntu Linux. 

#### Better switching between EN/JP documentation

- One can switch between Japanese and English documentation using the button to the left of the search box.

- Untranslated documents (for example, this page) are displayed in English even in Japanese mode. 

![](./fig/demo-en-jp-switch.gif)


### Others

- Tune-up in some supercomputers
- Updated version of Japanese community model [JIVSM](./1._SetUp/0104_dataset.md).


## Previous Revision Histories

2021-08-27 (v5.2.0)

:   Earth-flattening transformation, psmeca (a new source representation), seawater velocity structure w/ SOFAR channel, a new tool that converts 2D simulation output into pseudo-3D, supporting new Japanese supercomputers, and some minor bugfixes. 

2020-08-13 (v5.1.0)

:   Strict-mode

2019-08-27 (v5.0.0)

:   New web-based documentation, stress/strain waveform output, slip-based source specification, and many bugfixes. 

2017-09-21 (v4.0)

:   Minor bugfixes, new binary output for waveform, updated references.

2016-08-20 (v3.0)

:   Hybrid parallel simulation for 2D codes.

2016-06-19 (v2.0)

:   Hybrid parallel simulation for 2D codes.

2016-05-05 (v1.0)

:   Official first release as an open-source software

2016-02-03

:   Output in \texttt{NetCDF} format.

2016-01-14

:   Body force and reciprocity modes.

2015-07-14

:   MPI/OpenMP hybrid parallel simulation mode.

2015-06-29

:   Added random media.

2015-06-10

:   Revision for the new Earth Simulator.

2015-06-04

:   First closed version for the ERI/UT joint usage program.

