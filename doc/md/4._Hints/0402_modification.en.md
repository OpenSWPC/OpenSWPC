# For Modifying the Code

## Defining Your own Velocity Model

The velocity structure is defined by the subroutine `vmodel_*`,
called by the module `m_medium.F90`. These subroutines commonly have the
input/output parameters defined in the following table. 

  | 変数名                     | in/out | type    | 説明                                       |
  | -------------------------- | ------ | ------- | ------------------------------------------ |
  | `io_prm`                   | in     | integer | I/O number of the input parameter file     |
  | `i0,i1`                    | in     | integer | Start/end indices of arrays in x-direction |
  | `j0,j1`                    | in     | integer | Start/end indices of arrays in y-direction |
  | `k0,k1`                    | in     | integer | Start/end indices of arrays in z-direction |
  | `xc(i0:i1)`                | in     | real    | x grid locations                           |
  | `yc(i0:i1)`                | in     | real    | y grid locations                           |
  | `zc(i0:i1)`                | in     | real    | z grid locations                           |
  | `vcut`                     | in     | real    | Cut-off velocity                           |
  | `rho(k0:k1,i0:i1,j0:j1)`   | out    | real    | Mass density \[g/cm${}^3$\]                |
  | `lam(k0:k1,i0:i1,j0:j1)`   | out    | real    | Lame coefficient $\lambda$ \[g/cm${}^3$\]  |
  | `mu(k0:k1,i0:i1,j0:j1)`    | out    | real    | Lame coefficient $\mu$ \[g/cm${}^3$\]      |
  | `qp(k0:k1,i0:i1,j0:j1)`    | out    | real    | $Q_P$                                      |
  | `qs(k0:k1,i0:i1,j0:j1)`    | out    | real    | $Q_S$                                      |
  | `bddep(i0:i1,j0:j1,0:NBD)` | out    | real    | Discontinuity boundary depths \[km\]       |

 By creating a Fortran subroutine that
returns the medium parameters `rho, lam, mu, qp` and `qs` at locations
given in the input of the subroutines `xc, yc`, and `zc`, it is easy to
add a new velocity model.

The topography and bathymetry are automatically investigated in the
`m_medium` module after calling the `vmodel_ast` routine. To make this
investigation work properly, the medium parameter `mu` must be zero in
the air and ocean columns and `lam` must be zero in the air column.

The variables `bddep(:,:,0)` are assumed to be the topography, and are
used for the snapshot output. The other values of `bddep(:,:,1:NBD)` are
used to fit the source and/or station location to the discontinuity
depths. Providing dummy values of these functions is not necessary.

## Defining Your own Source Time Function


The source time function is called by the `source__momentrate` `Fortran`
function in `m_source.F90` based on the choice of `stftype`. The
definitions of the source time functions are given in
`share/m_fdtool.F90`. It is easy to add a new source time function here
and to add the call to the new function in the `m_source` module.

All of the pre-defined source time functions take two time parameters,
`tbeg` and `trise`. In the source code, they are stored in the array
variable `srcprm(:)`. If the new source time function requires more than
three parameters, the user can expand the array `srcprm(:)` to store
them.


## Appending New Control Parameters

In many `Fortran` modules, the first set-up is performed by subroutines
called `(modulename)__setup` during the first computation. Some of the
setup modules read parameters from the input parameter file. These
parameters are read by the subroutine `readini`, which is defined in
`shared/m_readini.F90`.