"""Accompanying python module for OpenSWPC

Copyright 2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.
"""

import sys

def prm_new():

    """ A standard parameterset of the OpenSWPC

    Parameters
    ----------
    None

    Return
    ------
    prm: dict
        Contains parameters of OpenSWPC equivalent to the example/input.inf
    """
    
    prm = {
        'title'           : 'swpc', 
        'odir'            : './out', 
        'ntdec_r'         : 50, 
        'strict_mode'     : False,
        'nproc_x'         : 2,
        'nproc_y'         : 2,
        'nx'              : 384,
        'ny'              : 384,
        'nz'              : 384,
        'nt'              : 1000,
        'dx'              : 0.5,
        'dy'              : 0.5,
        'dz'              : 0.5,
        'dt'              : 0.02,
        'vcut'            : 1.5,
        'xbeg'            : -96.0,
        'ybeg'            : -96.0,
        'zbeg'            : -10.0,
        'tbeg'            : 0.0,
        'clon'            : 139.7604,
        'clat'            : 35.7182,
        'phi'             : 0.0,
        'fq_min'          : 0.02,
        'fq_max'          : 2.00,
        'fq_ref'          : 1.0,  
        'snp_format'      : 'netcdf',
        'xy_ps%sw'        : False,
        'xz_ps%sw'        : True,
        'yz_ps%sw'        : False,
        'fs_ps%sw'        : False,
        'ob_ps%sw'        : True,
        'xy_v%sw'         : False,
        'xz_v%sw'         : True,
        'yz_v%sw'         : False,
        'fs_v%sw'         : False,
        'ob_v%sw'         : True,
        'xy_u%sw'         : False,
        'xz_u%sw'         : True,
        'yz_u%sw'         : False,
        'fs_u%sw'         : False,
        'ob_u%sw'         : True,
        'z0_xy'           :  7.0,
        'x0_yz'           :  0.0,
        'y0_xz'           :  0.0,
        'ntdec_s'         : 5,
        'idec'            : 2,
        'jdec'            : 2,
        'kdec'            : 2,
        'sw_wav_v'        : True,
        'sw_wav_u'        : False,
        'sw_wav_stress'   : False,
        'sw_wav_strain'   : False,
        'ntdec_w'         : 5,
        'ntdec_w_prg'     : 0, 
        'st_format'       : 'xy',
        'fn_stloc'        : './example/stloc.xy',
        'wav_format'      : 'sac',
        'stf_format'      : 'xym0ij',
        'stftype'         : 'kupper',
        'fn_stf'          : './example/source.dat',
        'sdep_fit'        : 'asis',
        'bf_mode'         : False,
        'pw_mode'         : False,
        'pw_ztop'         : 100.,
        'pw_zlen'         : 30.,
        'pw_ps'           : 'p',
        'pw_strike'       : 0.0,
        'pw_dip'          : 0.0,
        'pw_rake'         : 0.0,
        'abc_type'        : 'pml',
        'na'              : 20,
        'stabilize_pml'   : False,
        'fullspace_mode'  : False, 
        'vmodel_type'     : 'lhm',
        'is_ocean'        : True,
        'topo_flatten'    : False,
        'munk_profile'    : True,
        'earth_flattening': False,
        'vp0'             : 5.0,
        'vs0'             : 3.0,
        'rho0'            : 2.7,
        'qp0'             : 200,
        'qs0'             : 200,
        'topo0'           : 0,
        'dir_grd'         : 'dataset/vmodel/ejivsm/',
        'fn_grdlst'       : './example/grd.lst',
        'node_grd'        :  0,
        'fn_lhm'          : 'example/lhm.dat',
        'dir_rmed'        : './in/',
        'fn_grdlst_rmed'  : './example/grd.lst',
        'rhomin'          :  1.0,
        'fn_rmed0'        : 'dummy.ns',
        'green_mode'      : False,
        'green_stnm'      : 'st01',
        'green_cmp'       : 'z',
        'green_trise'     : 1.0,
        'green_bforce'    : False,
        'green_maxdist'   : 550.,
        'green_fmt'       : 'llz',
        'fn_glst'         : 'example/green.lst',
        'stopwatch_mode'  : False,
        'benchmark_mode'  : False,
        'ipad'            : 0,
        'jpad'            : 0,
        'kpad'            : 0}
    
    return prm


def prm_read(fn): 

    """Read OpenSWPC parameter file
    
    Parameters
    ----------
    fn: str
        parameter filename (with path)
    
    Return
    ------
    p: dict
        A dictionary which contains control parameters
    
    """
    p = prm_new()
    with open(fn, 'r') as f: 
        for line in f.readlines():
            s = line.rstrip().split()

            if s == []:
                continue
            if s[0] == '': 
                continue
            if s[0][0:1] == '!':
                continue
            k = s[0]
            v = s[2]
            
            if k in p:
                t = type(p[k])
                if isinstance(p[k], bool):
                    if v == '.true.' or v == '.True.' or v == '.TRUE.':
                        p[k] = True
                    else:
                        p[k] = False                        
                elif isinstance(p[k], str):
                    p[k] = v.replace("'", "").replace('"', '')
                elif isinstance(p[k], int):
                    p[k] = int(v)
                elif isinstance(p[k], float):
                    p[k] = float(v)
                else:
                    p[k] = v
            else:
                pass
    return p


def prm_print(prm, io=sys.stdout):

    """ export control parameters in OpenSWPC's format

    Parameters
    ----------
    prm: dict
        A dictionary which contains control parameters    
    io: file handler (Optional)
        Output file handler

    Return
    ------
    None

    Example
    -------
    # read a parameter file
    prm = prm_read('./example/input.inf')
    # modify it
    prm['title'] = 'new simulation'
    # output it to a new file
    with open('./in/input_new.inf', 'w') as fp:
        prm_print(prm, io=fp)
    """
    
    for k, v in prm.items():
        if type(v) is str:
            print(f"{k.ljust(18)} =   '{v}'", file=io)
        elif type(v) is bool:
            if v:
                print(f"{k.ljust(18)} =   .true.", file=io)
            else:
                print(f"{k.ljust(18)} =   .false.", file=io)
        else:
            print(f"{k.ljust(18)} =   {v}", file=io)


def prm_save(prm, filename): 

    """ save control parameters of OpenSWPC into file

    Parameters
    ----------
    prm: dict
        A dictionary which contains control parameters    
    filename: str
        Output filename    
    """

    with open(filename, 'w') as f:
        prm_print(prm, io=f)


def chkregion_2d(prm): 

    """ Visualize physical model size and MPI node partitioning

    Requires matplotlib. 

    Parameters
    ----------
    prm: dict
        A dictionary which contains control parameters

    Return
    ------
    fig: figure object of matplotlib.pyplot    
    """

    import matplotlib.pyplot as plt

    xbeg = prm['xbeg']
    zbeg = prm['zbeg']
    nx = prm['nx']
    nz = prm['nz']
    dx = prm['dx']
    dz = prm['dz']
    nproc_x = prm['nproc_x']

    col_node = (0., 0.3, 1.0)

    xend = xbeg + (nx-1) * dx
    zend = zbeg + (nz-1) * dz

    
    fig, ax = plt.subplots()
    
    ax.set(xlim   = [xbeg, xend], ylim   = [zbeg, zend], 
           xticks = [xbeg, xend], yticks = [zbeg, zend], 
           xlabel = 'x [km]',     ylabel = 'z [km]')
    ax.invert_yaxis()
    ax.grid()
    ax.set_aspect(1.0)

    for i in range(1, nproc_x):
        x = xbeg + i * (xend - xbeg) / nproc_x
        ax.plot([x, x], [zbeg, zend], color=col_node, alpha=0.2, lw=2)

    return fig


def chkregion_3d(prm):

    """ Visualize physical model size and MPI node partitioning

    Requires matplotlib. 

    Parameters
    ----------
    prm: dict
        A dictionary which contains control parameters

    Return
    ------
    fig: figure object of matplotlib.pyplot    
    """    

    import numpy as np
    import matplotlib.pyplot as plt
    
    xbeg = prm['xbeg']
    ybeg = prm['ybeg']
    zbeg = prm['zbeg']
    nx = prm['nx']
    ny = prm['ny']
    nz = prm['nz']
    dx = prm['dx']
    dy = prm['dy']
    dz = prm['dz']
    nproc_x = prm['nproc_x']
    nproc_y = prm['nproc_y']
    phi = prm['phi']

    xend = xbeg + (nx - 1) * dx
    yend = ybeg + (ny - 1) * dy
    zend = zbeg + (nz - 1) * dz

    Lx = xend - xbeg
    Ly = yend - ybeg
    Lz = zend - zbeg
    
    fig = plt.figure(figsize=(8, 8))

    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[Ly, Lz],
        height_ratios=[Lx, Lz],
        wspace=0.1,
        hspace=0.125,
    )

    if np.isclose(phi, 0): 
        xlabel = rf"$x$ (Northing) [km]"
        ylabel = rf"$y$ (Easting) [km]"
    elif np.isclose(phi, 90): 
        xlabel = rf"$x$ (Easting) [km]"
        ylabel = rf"$y$ (Southing) [km]"
    elif np.isclose(phi, 180): 
        xlabel = rf"$x$ (Southing) [km]"
        ylabel = rf"$y$ (Westing) [km]"
    elif np.isclose(phi, 270): 
        xlabel = rf"$x$ (Westing) [km]"
        ylabel = rf"$y$ (Northing) [km]"        
    else:
        xlabel = rf"$x$ (N{round(phi)}$\degree$E) [km]"
        ylabel = rf"$y$ (N{round(phi+90)}$\degree$E) [km]"
    zlabel = "$z$ (Depth) [km]"
    
    # XY
    ax_xy = fig.add_subplot(gs[0, 0])
    ax_xy.set(xlim   = [ybeg, yend], ylim   = [xbeg, xend], 
              xticks = [ybeg, yend], yticks = [xbeg, xend], 
              xlabel = ylabel,       ylabel = xlabel)
    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.set_anchor("E")
    ax_xy.tick_params(top    = True,  labeltop    = True, 
                      bottom = False, labelbottom = False,
                      left   = True,  labelleft   = True, 
                      right  = False, labelright  = False)
    ax_xy.xaxis.set_label_position('top')

    for j in range(1, nproc_y):
        y = ybeg + j * (yend - ybeg) / nproc_y
        ax_xy.plot([y, y], [xbeg, xend], color=(0., 0.3, 1.0), alpha=0.2, lw=2)

    for i in range(1, nproc_x):
        x = xbeg + i * (xend - xbeg) / nproc_x
        ax_xy.plot([ybeg, yend], [x, x],  color=(1.0, 0.3, 0.), alpha=0.2, lw=2)

    
    # YZ
    ax_yz = fig.add_subplot(gs[1, 0])
    ax_yz.set(xlim   = [ybeg, yend], ylim   = [zbeg, zend], 
              xticks = [ybeg, yend], yticks = [zbeg, zend], 
              xlabel = ylabel,       ylabel = zlabel)
    ax_yz.invert_yaxis()
    ax_yz.set_aspect("equal", adjustable="box")
    ax_yz.set_anchor("E")    # 右寄せ

    for j in range(1, nproc_y):
        y = ybeg + j * (yend - ybeg) / nproc_y
        ax_yz.plot([y, y], [zbeg, zend], color=(0., 0.3, 1.0), alpha=0.2, lw=2)


    # XZ
    ax_xz = fig.add_subplot(gs[0, 1])
    ax_xz.set(xlim   = [zbeg, zend], ylim   = [xbeg, xend], 
              xticks = [zbeg, zend], yticks = [xbeg, xend], 
              xlabel = zlabel,       ylabel = xlabel)
    ax_xz.set_aspect("equal", adjustable="box")
    ax_xz.yaxis.set_label_position('right')
    ax_xz.set_anchor("W")    # 左寄せ
    ax_xz.tick_params(top=False, bottom=True, left=False, right=True, 
                      labeltop=False, labelbottom=True, labelleft=False, labelright=True)

    for i in range(1, nproc_x):
        x = xbeg + i * (xend - xbeg) / nproc_x
        ax_xz.plot([zbeg, zend], [x, x],  color=(1.0, 0.3, 0.), alpha=0.2, lw=2)
    
    return fig


### Architecture Dependent Function

def bdec_rankmap(nproc_x, nproc_y, io=sys.stdout):

    for j in range(nproc_y):
        for i in range(nproc_x):
            print(f"({i//2},{j//2})", file=io)
