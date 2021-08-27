# Programs for Supporting Parameter Settings


## `fdmcond.x` 

The grid width in space and time in the finite difference method is
controlled by the stability condition. The wavelength condition will
affect the allowed maximum frequency radiated from the source.

The tool fdmcond.x can help determine these parameters to satisfy the
conditions. After the user specifies several parameters, such as the
grid width, maximum frequency (fmax), rise time (Tr), and minimum and
maximum velocities in the medium (vmin, vmax), the program can suggest
the other parameters.

``` text
$  ./fdmcond.x 
 
----------------------------------------------------------------------
                           FDM CONDITION                           
----------------------------------------------------------------------
 
 
  Model Dimension ? --> 3
   1) 2D
   2) 3D
 
 
 Source Type ? --> 3
   1) Triangle
   2) Herrmann
   3) Kupper
 
 
 Parameter Combination ? --> 5
   1) dh   (space grid),  fmax (max freq.),  vmax (max vel.)
   2) dh   (space grid),  Tr   (rise time),  vmax (max vel.)
   3) dh   (space grid),  fmax (max freq.),  dt (time grid)
   4) dh   (space grid),  Tr   (rise time),  dt (time grid)
   5) dh   (space grid),  vmin (min vel.),   vmax (max vel.)
   6) dh   (space grid),  vmin (min vel.),   dt (time grid) 
   7) fmax (max freq.) ,  vmax (max vel.),   dt (time grid)
   8) Tr   (rise time) ,  vmax (max vel.),   dt (time grid)
   9) vmin (min vel.)  ,  vmax (max vel.),   dt (time grid)
 
 
 Assumed Parameters: 
   dx     =   0.25
   dy     =   0.25
   dz     =   0.25
   vmin   =   0.3
   vmax   =   8.0
 
 Derivaed Parameters: 
   dt    <=   0.01546
   fmax  <=   0.17143
   Tr    >=  13.41667
 
 
```

## `mapregion.x`

The geographical region of the simulation will be automatically
determined by the parameters `clon`, `clat`, `phi`, `xbeg`, `ybeg`, `nx`, `ny`, `dx`,
and `dy`. The `mapregion.x` program reads the parameter file and exports the
outer edge of the region in longitude and latitude.


``` bash
$ mapregion.x -i input.inf -o region.dat
```
If the option `-o` is omitted, the result will be printed to the
standard output on the screen. This program will also estimate the total
memory usage in the standard error output.


## `mapregion.gmt`

These scripts use `mapregion.x` to visualize the region by using GMT. GMT version 5 or later is necessary. By default, these scripts plot only the region around the Japanese Islands.