# For Parameter Settings

The 3D simulation is bounded by the total memory size. The code requires

\begin{align}
\begin{split}
    &m_\text{MP} = 116 + 24 \mathtt{NM} = 188 \quad (\mathtt{NM}=3) \quad \text{byte} \\
\end{split}\end{align}

of memory for the case of mixed precision
(`MP=DP`) with a GNZ viscoelastic body of `NM=3`. Note that this is a
coarse estimate excluding the effect of an absorbing boundary.

The computation time can be roughly estimated by the parameter
$\mathtt{n_G}$, which is defined as the number of spatial and/or
temporal grids that one CPU can process in a second. This value depends
on the CPU, as shown in the following table. 

| Architecture Name              | CPU                           | \#Core | $\mathtt{n_G}$      |
| ------------------------------ | ----------------------------- | ------ | ------------------- |
| Mac Pro 2010                   | Intel Xeon X5670 2.93GHz      | 6      | $6.7 \times 10^{6}$ |
| EIC2015 (ERI,UTokyo)           | Intel Xeon E5-2680 v3 2.5 GHz | 12     | $7.0 \times 10^{6}$ |
| The earth simulator (3rd gen.) | NEC SX-ACE                    | 4      | $57 \times 10^{6}$  |

The total computation time can be estimated by 

\begin{align}
    t_\text{comp} = \frac{ {\tt nx}\times {\tt ny}\times {\tt nz}}{ \mathtt{n_G} \times \tt nproc} \times {\tt nt}
    \quad \text{[s]}\end{align}

where `ncore` is the number of
CPU cores used in the computation. If the estimated time exceeds that
provided by the computer system, it is recommended to make the model
size smaller and/or to use checkpointing/restarting.


