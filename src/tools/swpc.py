"""Accompanying python module for OpenSWPC"""

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
        'strict_mode'     : True,
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
        'xz_ps%sw'        : False,
        'yz_ps%sw'        : False,
        'fs_ps%sw'        : False,
        'ob_ps%sw'        : False,
        'xy_v%sw'         : False,
        'xz_v%sw'         : False,
        'yz_v%sw'         : False,
        'fs_v%sw'         : False,
        'ob_v%sw'         : False,
        'xy_u%sw'         : False,
        'xz_u%sw'         : False,
        'yz_u%sw'         : False,
        'fs_u%sw'         : False,
        'ob_u%sw'         : False,
        'z0_xy'           :  7.0,
        'x0_yz'           :  0.0,
        'y0_xz'           :  0.0,
        'ntdec_s'         : 5,
        'idec'            : 2,
        'jdec'            : 2,
        'kdec'            : 2,
        'sw_wav_v'        : False,
        'sw_wav_u'        : False,
        'sw_wav_stress'   : False,
        'sw_wav_strain'   : False,
        'ntdec_w'         : 5,
        'st_format'       : 'xy',
        'fn_stloc'        : './example/stloc.xy',
        'wav_format'      : 'sac',
        'wav_calc_dist'   : False,
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
        'dir_grd'         : '${DATASET}/vmodel/ejivsm/',
        'fn_grdlst'       : './example/grd.lst',
        'node_grd'        :  0,
        'fn_lhm'          : 'example/lhm.dat',
        'dir_rmed'        : './in/',
        'fn_grdlst_rmed'  : './example/grd.lst',
        'rhomin'          :  1.0,
        'fn_rmed0'        : 'dummy.ns',
        'is_ckp'          : False,
        'ckpdir'          : './out/ckp',
        'ckp_interval'    : 1000000,
        'ckp_time'        : 1000000.,
        'ckp_seq'         : True,
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
                    if v == '.true.' or '.True.' or '.TRUE.':
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

    Parameter
    ---------
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


### Architecture Dependent Function

def bdec_rankmap(nproc_x, nproc_y, io=sys.stdout):

    for j in range(nproc_y):
        for i in range(nproc_x):
            print(f"({i//2},{j//2})", file=io)
