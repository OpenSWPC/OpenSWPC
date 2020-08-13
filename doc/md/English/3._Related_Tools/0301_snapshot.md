# Snapshot data handling

## `read_snp.x` 

Snapshot files in both `NetCDF` and the originally defined binary format
can be extracted or visualized by the program `read_snp.x`.


``` bash
read_snp.x -i snapshotfile [-h] [-ppm|-bmp] [-pall] 
              [-mul var | -mul1 var -mul2 var ...] [-bin] [-asc] [-skip n]
```

  `-h`
  : Print the header information defined in the snapshot, as in the
    following example.

    ``` txt
     $ ../bin/read_snp.x -i swpc_3d.xz.ps.snp -h

     [binary type]   : STREAMIO
     [code type]     : SWPC_3D
     [header version]:          3
     [title]         : swpc_3d
     [date generated]: 1408015126
                       2014-08-14T11-18-46
     [coordinate]    : xz
     [data type]     : ps
     [ns1]           :        256
     [ns2]           :        256
     [beg1]          :       -63.87500
     [beg2]          :        -9.87500
     [ds1]           :         0.25000
     [ds2]           :         0.25000
     [dt]            :         0.05000
     [na1]           :         20
     [na2]           :         20
     [nmed]          :          3
     [nsnp]          :          2
     [clon]          :       143.50000
     [clat]          :        42.00000
    ```

  `-ppm`, `-bmp`
  : Visualize and export the image files in ppm or bmp format. The `ppm`
    or `bmp` directory will be automatically created in the current
    directory and image files with sequential numbers will be stored
    there. If the snapshot file is displacement or velocity, the
    absolute values of the vertical and horizontal amplitudes will be
    colored red and green, respectively. For the `PS` file, the absolute
    values of the divergence and rotation vector will be colored
    similarly. If the absolute value option is specified, the
    black-red-yellow-white color palette (similar to the "hot" color
    palette in GMT) will be adopted. For cross sections along the
    surface (`ob`, `fs`), the topography color map will be overlaid. For
    other cross sections, the velocity structure in the section will be
    overlaid.

!!! bug "Limitation of the BMP format"
    Output to `bmp` format occasionally fails as this format has a restriction to the image size. `ppm` format is recommended. 
      
  `-pall`
  : Visualize including the absorbing boundary region. This option works
    only if it is used with . By default, the absorbing boundary region
    will be clipped.

  `-mul`
  : Multiply `var` by the amplitude for visualization. Adjust the
    visualized color by changing this value. Optionally, by specifying
    `-mul1` or `-mul2`, for example, one may change the weight of the
    amplitude by component.

  `-abs`
  : Visualize the absolute value of the vector. This only works with the
    velocity or displacement snapshots.

  `-bin`, `-asc`
  : Export the snapshot data to binary (`-bin`) or ascii (`-asc`) files.
    The data file will be created in the automatically created `bin` or
    `asc` directories. The binary formatted data can be directly used in
    `GMT` with the `xyz2grd` module by appending the `-bis` option.

  `-skip n`
  : Skip the first $n$ snapshots for visualization or data exports.

  `-notim` (after v5.1)
  : Do *not* plot the elapsed time in the snapshot figures.


## `diff_snp.x`

This program takes the difference between two snapshots and exports it
to another snapshot file.

``` bash
$ diff_snp.x snap1 snap2 diffile
```

The output file format (`NetCDF` or binary) depends on the input file
format.


## `fs2grd.x`  

**New in v5.1**

Although snapshot data from `OpenSWPC` along the ground surface or ocean bottom follows the NetCDF format, they cannot be used in the GMT's grdimage command because they are not evenly-spaced grid data along longitude and latitude. 

The utility program `fs2grd.x` resamples the `OpenSWPC`'s output dataset in longitude/latitude coordinate to convert it to GMT-friendly grd-format dataset. 

``` bash
$ fs2grd.x -i input.nc -v variable_name 
           -R region -dlon delta_lon -dlat delta_lat 
```


  `-i`
  : Specify the output file of `OpenSWPC` in NetCDF format to be resampled. The snapshot must be on the free surface (`fs`), ocean bottom (`ob`) or `xy` coordinates. 

  `-R` `lon0/lon1/lat0/lat1`
  : The region of resampling. Minimum(`lon0`) and maximum(`lon1`) longigude, minimum(`lat0`) and maximum(`lat1`) latitude. The formatting is similar to that in the GMT, but a blank space is necessary between `-R` and arguments. 
  
  `-dlon` delta_lon, `-dlat` delta_lat
  : Grid spacings in longitudal and latitudal directions. 

  `-v` variable_name
  : Variable name to be resampled. Any horizontal-space 2D variables can be specified in the snapshot. One may confirm the list of the variables in the NetCDF file by `ncdump -h` command. If the specified variable is time-dependent, such as `Vx, Vy, Vz` or `div, rot_x, rot_y, rot_z`, the `fs2grd.x` will export time-sequential (many) files. 