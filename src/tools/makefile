## ------------------------------------------------------------------------- ##
##
## makefile for tools in SWPC
##
## @copyright
##   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
##
## --

all: read_snp.x diff_snp.x mapregion.x fdmcond.x qmodel_tau.x grdsnp.x \
	gen_rmed2d.x gen_rmed3d.x ll2xy.x xy2ll.x wav2sac.x fs2grd.x wvconv.x

## file search path

VPATH = ../shared/:../../bin/

##
## Compiler & Compile Option
##
include ../shared/makefile-tools.arch

.F90.o:
	$(FC) $(NCINC) $(FFLAGS) $(NCFLAG) -c $< -o $@ $(NETCDF)

## file dependency

m_std.o     : m_std.F90
m_system.o  : m_system.F90 m_std.o
m_getopt.o  : m_getopt.F90 m_std.o m_system.o
m_daytim.o  : m_daytim.F90 m_std.o
m_gk.o      : m_gk.F90   m_std.o
m_stamp.o   : m_stamp.F90 m_std.o m_pnm.o
m_pnm.o     : m_pnm.F90    m_std.o m_debug.o
m_fdsnap.o  : m_fdsnap.F90 m_std.o m_daytim.o m_geomap.o
m_rmedia.o  : m_rmedia.F90 m_std.o m_fk.o m_gammaf.o
m_filt2d.o  : m_filt2d.F90 m_std.o m_fk.o
m_fdtool.o  : m_fdtool.F90 m_std.o m_gk.o
m_readini.o : m_readini.F90 m_std.o m_system.o m_debug.o
m_geomap.o  : m_geomap.F90 m_gk.o
m_bmp.o     : m_bmp.F90 m_std.o m_debug.o
m_debug.o   : m_debug.F90 m_std.o
m_bicubic.o : m_bicubic.F90 m_std.o
m_color.o   : m_color.F90 m_std.o
m_fk.o      : m_fk.F90 m_std.o m_rfft.o
m_rfft.o    : m_rfft.F90 m_std.o
m_gammaf.o  : m_gammaf.F90 m_std.o
m_sac.o    : m_sac.F90 m_std.o

.SUFFIXES:
.SUFFIXES: .F90 .o

../../bin/gen_rmed2d.x: gen_rmed2d.F90 m_std.o m_system.o m_getopt.o m_rmedia.o m_daytim.o m_fk.o m_rfft.o m_gammaf.o  m_debug.o
	$(FC) $(NCLIB) $(NCINC) $(NCFLAG) $(FFLAGS) -o $@ $^ -D_DEBUG $(NETCDF)

../../bin/gen_rmed3d.x: gen_rmed3d.F90 m_std.o m_system.o m_getopt.o m_rmedia.o m_fk.o m_rfft.o m_gammaf.o m_debug.o
	$(FC) $(NCLIB) $(NCINC) $(NCFLAG) $(FFLAGS) -o $@ $^ -D_DEBUG $(NETCDF)

../../bin/qmodel_tau.x: qmodel_tau.F90 m_std.o m_system.o m_getopt.o m_fdtool.o m_readini.o m_debug.o
	$(FC) $(FFLAGS) -o $@ $^

../../bin/grdsnp.x: grdsnp.F90 m_std.o m_system.o m_getopt.o m_fdtool.o m_readini.o m_bicubic.o m_geomap.o m_gk.o m_debug.o
	$(FC) $(NCLIB) $(NCINC) $(NCFLAG) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/read_snp.x: read_snp.F90 m_std.o m_system.o m_getopt.o m_pnm.o m_stamp.o m_fdsnap.o m_daytim.o m_color.o m_fk.o m_rfft.o m_filt2d.o m_geomap.o m_gk.o m_bmp.o m_debug.o
	$(FC) $(NCLIB) $(NCINC) $(NCFLAG) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/diff_snp.x: diff_snp.F90 m_std.o m_system.o m_fdsnap.o m_daytim.o m_gk.o m_debug.o
	$(FC) $(NCLIB) $(NCINC) $(NCFLAG) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/mapregion.x: mapregion.F90 m_std.o m_system.o m_fdtool.o m_gk.o m_readini.o m_getopt.o m_geomap.o m_debug.o
	$(FC) $(FFLAGS) -o $@ $^

../../bin/fdmcond.x: fdmcond.F90 m_std.o m_getopt.o m_system.o m_debug.o
	$(FC) $(FFLAGS) -o $@ $^

../../bin/ll2xy.x: ll2xy.F90 m_std.o m_geomap.o m_system.o m_gk.o m_debug.o
	$(FC) $(FFLAGS) -o $@ $^ 

../../bin/xy2ll.x: xy2ll.F90 m_std.o m_geomap.o m_system.o m_gk.o m_debug.o
	$(FC) $(FFLAGS) -o $@ $^ 

../../bin/wav2sac.x: wav2sac.F90 m_std.o m_system.o m_daytim.o m_sac.o 
	$(FC) $(FFLAGS) -o $@ $^ 

../../bin/fs2grd.x: fs2grd.F90 m_std.o m_system.o m_gk.o m_geomap.o m_getopt.o m_bicubic.o m_fdsnap.o m_daytim.o
	$(FC) $(NCLIB) $(NCINC) $(NCFLAG) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/wvconv.x: wvconv.F90 m_std.o m_rfft.o m_fk.o m_sac.o 
	$(FC) $(FFLAGS) -o $@ $^ 

clean:
	/bin/rm -f *.o *.mod *~ ../../bin/read_snp.x ../../bin/diff_snp.x ../../bin/mapregion.x ../../bin/fdmcond.x ../../bin/qmodel_tau.x ../../bin/grdsnp.x ../../bin/gen_rmed2d.x ../../bin/gen_rmed3d.x ../../bin/ll2xy.x ../../bin/xy2ll.x
