# Reciprocity Mode

This mode excites the seismic wave at a specified station location and
exports the velocity and/or strain velocity of multiple virtual source
locations. Based on the reciprocity theorem, this result corresponds to
the body force and/or moment tensor response from virtual source
locations observed at specified stations. If the time duration of the
source time function is sufficiently short, they can be treated as
Green's functions.

If we denote the Green's tensor, from the virtual source
$\mathbf{\xi}$ to the receiver $\mathbf{r}$, as
$G_{ij}\left( \mathbf{r}, t; \mathbf{\xi}\right)$, this mode
simulates the convolution of the spatial derivatives of Green's tensor
with the source time function $s(t)$ as

\begin{align}
\begin{split}
  G^{M1}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_x} \ast s(t) =   \frac{\partial G_{ix} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x} \ast s(t) 
  \\
  G^{M2}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \frac{\partial G_{iy} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} \ast s(t) 
   =     \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_y} \ast s(t)    
  \\
  G^{M3}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \frac{\partial G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_z} \ast s(t) 
  =      \frac{\partial G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_z} \ast s(t) 
  \\
  G^{M4}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \left( 
    \frac{\partial G_{iy} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} \right) \ast s(t) 
    \\
  &= \left( 
    \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_y} \right) \ast s(t) 
  \\
  G^{M5}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \left(
    \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_x}\right) \ast s(t)  
    \\& = \left(
    \frac{\partial G_{ix} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x}\right) \ast s(t) 
  \\
  G^{M6}_{i} \left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \left(
    \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} + 
    \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x} \right) \ast s(t) 
    \\
  &= \left(
    \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} + 
    \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x} \right) \ast s(t)
 \end{split} 
\end{align}

which corresponds to the moment tensor
response. Optionally, the body-force response

\begin{align}
\begin{split}
G^{B1}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv  G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) \ast s(t)
  =       G_{ix} \left( \mathbf{\xi}, t; \mathbf{r}\right) \ast s(t)  
  \\
G^{B2}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv G_{iy} \left( \mathbf{r}, t; \mathbf{\xi}\right) \ast s(t) 
  =      G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) \ast s(t) 
\\
  G^{B3}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right)\ast s(t)
  =      G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right)\ast s(t)
\end{split}
\end{align}

can be calculated.


To use this mode, the users should specify the station name `green_stnm`
of the receiver. This station name should be contained in the station
list file. `OpenSWPC` radiates the seismic wave by an excitation force
with a direction specified by the `green_cmp` parameter and a source
time function of the rise time, `green_trise`. To obtain the full
response of all components, three independent simulations with
`green_cmp=’x’`, `’y’`, and `’z’` are necessary.

The virtual source location should be given in the Cartesian or
geographical coordinates and depth (the format is described in Table
[\[table:greenf\]](#table:greenf){reference-type="ref"
reference="table:greenf"}) with unique integer ID numbers (`gid`).
Multiple virtual source locations can be specified in the simulation.
The `gid`s do not need to be sequential.

The output file is stored in the directory `(odir)/green/(gid)` in the
`SAC` format with the name convention `(title)__(green_cmp)__mij__.sac`
(for the moment tensor response) or `(title)__(green_cmp)__fi__.sac`
(for the body force response).

The amplitudes of the output files are multiplied by $10^9$ to compare
the `SAC`-formatted files in nm or nm/s units. The vertical component of
the output file is changed to be positive upward. However, the
derivative with respect to depth is performed according to the original
definition of positive downward.


!!! Info "Parameters"
    `green_mode`
    :Flags to turn the reciprocity mode on. If this is `.true.`, the
    other earthquake source parameters will be ignored.

    `green_stnm`
    : Name of the virtual station. This name must be included in the
    station list.

    `green_cmp`
    : Component at the virtual receiver. Choose from `’x’`, `’y’`, or
    `’z’`.

    `green_trise`
    : Rise time of the source time function convolved with the simulated
    Green's function.

    `green_bforce`
    : If `.true.`, calculate the body force response as well as the moment
    tensor response. The default setting is `.false.`.

    `green_fmt`
    : Format specification of the virtual source location. Choose from
    `’xyz’` (Cartesian coordinate; default) or `’llz’` (longitude,
    latitude, and depth).

    `green_maxdist`
    : The reciprocity wave will only be calculated if the horizontal
    distance is shorter than this parameter. Specify in units of km.

    `fn_glst`
    : Name of the virtual source location file.

    `stftype`
    : Source time function type. Same as the case with the moment tensor
    source.

    `ntdec_w`
    : Temporal decimation factor of the output waveforms. Same as the case
    with the normal waveform output.


The virtual source location file takes the following format according to the settings of `green_fmt`: 

 | `green_fmt` | format |     |     |       |
 | ----------- | ------ | --- | --- | ----- |
 | `’xyz’`     | x      | y   | z   | `gid` |
 | `’llz’`     | lon    | lat | z   | `gid` |
 
