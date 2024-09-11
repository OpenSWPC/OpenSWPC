## ------------------------------------------------------------------------- ##
##
## makefile for tools in SWPC
##
## @copyright
##   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
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

.f90.o:
	$(FC) $(NCINC) $(FFLAGS) -c $< -o $@ $(NETCDF)

## file dependency

m_std.o     : m_std.f90
m_system.o  : m_system.f90 m_std.o
m_getopt.o  : m_getopt.f90 m_std.o m_system.o
m_daytim.o  : m_daytim.f90 m_std.o
m_gk.o      : m_gk.f90   m_std.o
m_stamp.o   : m_stamp.f90 m_std.o m_pnm.o
m_pnm.o     : m_pnm.f90    m_std.o m_debug.o
m_fdsnap.o  : m_fdsnap.f90 m_std.o m_daytim.o m_geomap.o
m_rmedia.o  : m_rmedia.f90 m_std.o m_fk.o m_gammaf.o
m_filt2d.o  : m_filt2d.f90 m_std.o m_fk.o
m_fdtool.o  : m_fdtool.f90 m_std.o m_gk.o
m_readini.o : m_readini.f90 m_std.o m_system.o m_debug.o
m_geomap.o  : m_geomap.f90 m_gk.o
m_bmp.o     : m_bmp.f90 m_std.o m_debug.o
m_debug.o   : m_debug.f90 m_std.o
m_bicubic.o : m_bicubic.f90 m_std.o
m_color.o   : m_color.f90 m_std.o
m_fk.o      : m_fk.f90 m_std.o m_rfft.o
m_rfft.o    : m_rfft.f90 m_std.o
m_gammaf.o  : m_gammaf.f90 m_std.o
m_sac.o    : m_sac.f90 m_std.o
m_version.o: m_version.f90

.SUFFIXES:
.SUFFIXES: .f90 .f90 .o

../../bin/gen_rmed2d.x: gen_rmed2d.f90 m_std.o m_system.o m_getopt.o m_rmedia.o m_daytim.o m_fk.o m_rfft.o m_gammaf.o  m_debug.o m_version.o
	$(FC) $(NCLIB) $(NCINC) $(FFLAGS) -o $@ $^ -D_DEBUG $(NETCDF)

../../bin/gen_rmed3d.x: gen_rmed3d.f90 m_std.o m_system.o m_getopt.o m_rmedia.o m_fk.o m_rfft.o m_gammaf.o m_debug.o m_version.o
	$(FC) $(NCLIB) $(NCINC) $(FFLAGS) -o $@ $^ -D_DEBUG $(NETCDF)

../../bin/qmodel_tau.x: qmodel_tau.f90 m_std.o m_system.o m_getopt.o m_fdtool.o m_readini.o m_debug.o m_version.o
	$(FC) $(FFLAGS) -o $@ $^

../../bin/grdsnp.x: grdsnp.f90 m_std.o m_system.o m_getopt.o m_fdtool.o m_readini.o m_bicubic.o m_geomap.o m_gk.o m_debug.o m_version.o
	$(FC) $(NCLIB) $(NCINC) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/read_snp.x: read_snp.f90 m_std.o m_system.o m_getopt.o m_pnm.o m_stamp.o m_fdsnap.o m_daytim.o m_color.o m_fk.o m_rfft.o m_filt2d.o m_geomap.o m_gk.o m_bmp.o m_debug.o m_version.o
	$(FC) $(NCLIB) $(NCINC) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/diff_snp.x: diff_snp.f90 m_std.o m_system.o m_fdsnap.o m_daytim.o m_gk.o m_debug.o m_version.o
	$(FC) $(NCLIB) $(NCINC) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/mapregion.x: mapregion.f90 m_std.o m_system.o m_fdtool.o m_gk.o m_readini.o m_getopt.o m_geomap.o m_debug.o m_version.o
	$(FC) $(FFLAGS) -o $@ $^

../../bin/fdmcond.x: fdmcond.f90 m_std.o m_getopt.o m_system.o m_debug.o m_version.o
	$(FC) $(FFLAGS) -o $@ $^

../../bin/ll2xy.x: ll2xy.f90 m_std.o m_geomap.o m_system.o m_gk.o m_debug.o m_version.o m_getopt.o
	$(FC) $(FFLAGS) -o $@ $^ 

../../bin/xy2ll.x: xy2ll.f90 m_std.o m_geomap.o m_system.o m_gk.o m_debug.o m_version.o  m_getopt.o
	$(FC) $(FFLAGS) -o $@ $^ 

../../bin/wav2sac.x: wav2sac.f90 m_std.o m_system.o m_daytim.o m_sac.o m_version.o m_getopt.o
	$(FC) $(FFLAGS) -o $@ $^ 

../../bin/fs2grd.x: fs2grd.f90 m_std.o m_system.o m_gk.o m_geomap.o m_getopt.o m_bicubic.o m_fdsnap.o m_daytim.o m_version.o
	$(FC) $(NCLIB) $(NCINC) $(FFLAGS) -o $@ $^ $(NETCDF)

../../bin/wvconv.x: wvconv.f90 m_std.o m_rfft.o m_fk.o m_sac.o m_version.o m_getopt.o m_system.o
	$(FC) $(FFLAGS) -o $@ $^ 

clean:
	/bin/rm -f *.o *.mod *~ ../../bin/read_snp.x ../../bin/diff_snp.x ../../bin/mapregion.x ../../bin/fdmcond.x ../../bin/qmodel_tau.x ../../bin/grdsnp.x ../../bin/gen_rmed2d.x ../../bin/gen_rmed3d.x ../../bin/ll2xy.x ../../bin/xy2ll.x
