# Velocity Structure

## `qmodel_tau.x`

Calculate the frequency dependence of $Q^{-1}$ and the body wave
dispersion from the input parameter file.

``` bash
$ qmodel_tau.x -nm [nm] -i [prm_file] -f0 [min_freq] -f1 [max_freq] -nf [ngrid]
```

This discretizes the frequency range from `min_freq` to `max_freq` into
`ngrid` and exports $Q^{-1} ( f )$ and physical dispersion. The latter
is normalized to 1 at the reference frequency. The parameters related to
the viscoelastic body are read from the input parameter file; however,
the number of bodies `nm` should be specified separately because it is
hard-coded into the program.

## `grdsnp.x`

From the input parameter file, calculate and print the discontinuity of
the input `NetCDF` file in Cartesian coordinates for the simulation
($x$, $y$, depth) in the standard output. This program is used to
confirm the coordinate transformation and the detailed digital model,
and to visualize the model in the computational domain.

``` 
$ grdsnp.x -i [prm_file] -g [grd_file]
```
