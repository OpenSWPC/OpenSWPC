# Compilation and Execution


## make

The directories `src/swpc_3d`, `src/swpc_psv`, `src/swpc_sh`, and `src/tools` contain makefile. Execute the `make` command in each directory to
generate the executable binaries. An executable file (with a `*.x`
extension) will be stored in the `bin` directory.

## Specifying Compiler Options


In the makefiles, the following variables must be specified according to
the environment:


| variable | description |
| ----- | ------ |
| `FC`  | compiler name |
| `FFLAGS` | compiler option |
| `NCFLAG` | `NetCDF` flag |
| `NCLIB` | location of the  `NetCDF` library directory |
| `NCINC`  | location of the  `NetCDF` header file directory|
| `NETCDF` | linker option for `NetCDF`|


If `NCFLAG = -D_NETCDF` is specified, the `make` command will try to compile `OpenSWPC` with `NetCDF` library.

A set of the above variables under different computer environments is
defined in `src/shared/makefile.arch` and
`src/shared/makefile-tools.arch`. The former is for the compilation of
FDM codes, and the latter is for the compilation of misc tools. The user
can specify the `arch` option in the `make` command as in the following example: 

```bash
make arch=mac-gfortran
```

In the above case, the appropriate compiler options for the architecture `mac-gfortran` is automatically set. In some environments, the `debug` option is predefined. This option is used in the following way: 

```bash
make arch=eic debug=true
```

The list of pre-defined architecture (`arch`) options is described in the following table.

| `arch` name | target  |`NetCDF` location |
| ------------------ | ---------------------------------------- |  ----------------------- |
|  mac-intel     | Mac OSX + Intel Compiler + OpenMPI  |    `${HOME}/local` |
|  mac-gfortran  | Mac OSX + gfortran + Open MPI   |     `/usr/local` |
|  eic           | EIC (ERI, UTokyo) with the Intel Compiler  |    `${HOME}/local` |
|  fx            | Fujitsu FX10, FX100 and the K-computer | `${HOME}/xlocal` |
|  es3           | The Earth Simulator 3 (obsolete)            |   Provided by the system  |
|  es4           | The Earth Simulator 4                       |   Provided by the system  |
|  ubuntu-gfortran   | Ubuntu 16.04LTS + gfortran + Open MPI  |  Installation by `apt` |
| ofp (or oak) | Oakforest-PACS of the University of Tokyo (obsolete) | automatically specified by the `module` command |
| obcx | Oakbridge-CX of the University of Tokyo (obsolete) | automatically specified by the `module` command |
| bdec-o |Wisteria/BDEC-01 (Odyssey)  of the University of Tokyo  | automatically specified by the `module` command |
| mac-m1  | macOS + gfortran (Apple Silicon (M1/M2) + Homebrew) | `/opt/homebrew/` |


## More about the `NetCDF` library

The `NetCDF` library consists of the following items:

- `libnetcdf.*`:   `NetCDF` library file
- `libnetcdff.*`:  `NetCDF` Fortran library file (only for the `NetCDF` version 4 or later)
- `netcdf.mod`:    Fortran module information file

The extension of the library files may be `*.a` (static library) or
`*.so` (dynamic library), depending on the installation. All these files
are necessary for successful compilation with `NetCDF`. In particular,
the `netcdf.mod` file must be created by the same Fortran compiler as
`OpenSWPC`. If `NetCDF` is installed using packaging tools such as
`yum, apt`, or `homebrew`, the use of `gfortran` is implicitly assumed.

## Parameters embedded in the source code

Although most of the controlling parameter of the `OpenSWPC` are given in the input parameter file at the time of execution, a few parameters shown below are embedded in the source code for improving computational performance. All of these parameters are defined in `m_global.F90`. If the users edit these parameters, re-compilation (`make clean; make`) is necessary. 


!!! Info "Parameters"
    **`UC`**
    : Unit conversion coefficient from the computation unit (a mixed unit system using km for distance, g/cm$^3$ for mass density, and km/s for seismic wave velocity) to the SI unit. If the user changes the unit of input, this parameter should be changed. The default values are `1e-15` for the 3D code and `1e-12` for the 2D codes. 
    
    **`MP`**
    : Preicsion constant for the finite difference calculation. If `MP=DP` (this is the default), the key calculation of the finite difference will be performed in double precision, while the other calculations are in single precision. If one changes it to `MP=SP`, all calculation will be done in single precision. In this case, the memory used in the 3D simulation is 2/3 as much as in the `MP=DP` case, and compuational speed also is faster than the default. On the other hand, single-precision calculation may lead numerical instability around the earthquake source after long-term integration. The values of `DP` and `SP` may differ at different compiler, but usually they take `8` and `4`. 

    **`NM`**
    : Number of viscoelastic bodies adopted in the generalized Zener body. If this number is greater than `1`, anelastic attenuation is automatically adjusted so that Q-values are approximately frequency independent in a frequency range specified by the paramreters `fq_*`. If this parameter is zero, the simulation is performed under the assumption of perfectly elastic body, and thus no intrinsic attenuation is taken into account. For 3D simulations, this parameter significantly affects the computation time and required memory size. 
