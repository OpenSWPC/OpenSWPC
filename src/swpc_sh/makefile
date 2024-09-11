## ------------------------------------------------------------------------------------------------------------------------------ ##
##
## makefile for SWPC_SH
##
## @copyright
##   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
##
## ------------------------------------------------------------------------------------------------------------------------------ ##

##
## Target
##
all: ../../bin/swpc_sh.x

##
## Settings
##
VPATH = ../shared/:../../bin/       # file search path
.SUFFIXES:                          # neglect .mod files
.SUFFIXES: .f90 .f90 .o                  # only uses .f90 and .o file in makefile


##
## Compiler & Compile Option
##
include ../shared/makefile.arch


##
## Source Code Specification
##



## kernel module
s_global  = m_global.f90
s_kernel  = m_kernel.f90
s_source  = m_source.f90
s_fdtool  = m_fdtool.f90
s_snap  = m_snap.f90
s_wav  = m_wav.f90
s_report  = m_report.f90
s_medium  = m_medium.f90
s_absorb_c  = m_absorb_c.f90
s_absorb_p  = m_absorb_p.f90
s_absorb  = m_absorb.f90
s_vmodel0 = m_vmodel_user.f90
s_vmodel1 = m_vmodel_uni.f90
s_vmodel2 = m_vmodel_grd.f90
s_vmodel3 = m_vmodel_lhm.f90
s_vmodel1r = m_vmodel_uni_rmed.f90
s_vmodel2r = m_vmodel_grd_rmed.f90
s_vmodel3r = m_vmodel_lhm_rmed.f90
s_version = m_version.f90

## common library modules
s_system  = m_system.f90
s_std     = m_std.f90
s_getopt  = m_getopt.f90
s_pwatch  = m_pwatch.f90
s_gk      = m_gk.f90
s_bicubic = m_bicubic.f90
s_daytim  = m_daytim.f90
s_sac    = m_sac.f90
s_readini = m_readini.f90
s_geomap  = m_geomap.f90
s_pdebug  = m_pdebug.f90
s_rdrmed  = m_rdrmed.f90

## main program
s_main      = main.f90


##
## Objectd Filenames
##

o_std     = $(s_std:.f90=.o)
o_pdebug  = $(s_pdebug:.f90=.o)
o_getopt  = $(s_getopt:.f90=.o)
o_system  = $(s_system:.f90=.o)
o_pwatch  = $(s_pwatch:.f90=.o)
o_gk      = $(s_gk:.f90=.o)
o_bicubic = $(s_bicubic:.f90=.o)
o_daytim  = $(s_daytim:.f90=.o)
o_sac    = $(s_sac:.f90=.o)
o_readini = $(s_readini:.f90=.o)
o_geomap  = $(s_geomap:.f90=.o)
o_rdrmed  = $(s_rdrmed:.f90=.o)
o_version = $(s_version:.f90=.o)

o_global  = $(s_global:.f90=.o)
o_kernel  = $(s_kernel:.f90=.o)
o_source  = $(s_source:.f90=.o)
o_fdtool  = $(s_fdtool:.f90=.o)
o_snap  = $(s_snap:.f90=.o)
o_wav = $(s_wav:.f90=.o)
o_report  = $(s_report:.f90=.o)
o_medium  = $(s_medium:.f90=.o)
o_absorb  = $(s_absorb:.f90=.o)
o_absorb_c  = $(s_absorb_c:.f90=.o)
o_absorb_p  = $(s_absorb_p:.f90=.o)
o_vmodel0 = $(s_vmodel0:.f90=.o)
o_vmodel1 = $(s_vmodel1:.f90=.o)
o_vmodel2 = $(s_vmodel2:.f90=.o)
o_vmodel3 = $(s_vmodel3:.f90=.o)
o_vmodel1r = $(s_vmodel1r:.f90=.o)
o_vmodel2r = $(s_vmodel2r:.f90=.o)
o_vmodel3r = $(s_vmodel3r:.f90=.o)

##
## Module File Code Dependencies
##

$(o_std):     $(s_std)
$(o_getopt):  $(s_getopt)  $(o_std) $(o_system)
$(o_system):  $(s_system)  $(o_std)
$(o_pdebug):  $(s_pdebug)  $(o_std)
$(o_pwatch):  $(s_pwatch)  $(o_std)
$(o_gk):      $(s_gk)    $(o_std)
$(o_bicubic): $(s_bicubic) $(o_std)
$(o_daytim):  $(s_daytim)  $(o_std)
$(o_sac):    $(s_sac)     $(o_std) 
$(o_readini): $(s_readini) $(o_std) $(o_system) $(o_pdebug)
$(o_geomap):  $(s_geomap)  $(o_std) $(o_gk)
$(o_rdrmed):  $(s_rdrmed)  $(o_std) $(o_pdebug)
$(o_version): $(s_version) 

$(o_snap):  $(s_snap)  $(o_std) $(o_global) $(o_pwatch) $(o_fdtool) $(o_readini) $(o_geomap) $(o_pdebug)
$(o_wav):     $(s_wav)     $(o_std) $(o_global) $(o_pwatch) $(o_readini) $(o_geomap) $(o_pdebug) $(o_sac) $(o_system)
$(o_global):  $(s_global)  $(o_std) $(o_fdtool) $(o_pwatch) $(o_system) $(o_pwatch) $(o_daytim) $(o_readini) $(o_pdebug)
$(o_fdtool):  $(s_fdtool)  $(o_std) $(o_pdebug)
$(o_source):  $(s_source)  $(o_std) $(o_global) $(o_fdtool) $(o_pwatch) $(o_readini) $(o_geomap) $(o_pdebug)
$(o_kernel):  $(s_kernel)  $(o_std) $(o_global) $(o_pwatch) $(o_medium) $(o_readini) $(o_pdebug)
$(o_report):  $(s_report)  $(o_std) $(o_global) $(o_kernel) $(o_pwatch) $(o_readini) $(o_pdebug) $(o_version)
$(o_medium):  $(s_medium)  $(o_std) $(o_global) $(o_pwatch) $(o_readini) $(o_vmodel1) $(o_vmodel2) $(o_vmodel3) $(o_vmodel4) $(o_pdebug) $(o_vmodel1r) $(o_vmodel3r) $(o_vmodel2r) $(o_vmodel0)
$(o_absorb):  $(s_absorb)  $(o_std) $(o_global) $(o_pwatch) $(o_readini) $(o_absorb_p) $(o_absorb_c)  $(o_pdebug)
$(o_absorb_c):  $(s_absorb_c)  $(o_std) $(o_fdtool) $(o_global) $(o_pdebug)
$(o_absorb_p):  $(s_absorb_p)  $(o_std) $(o_fdtool) $(o_global) $(o_readini) $(o_pdebug)
$(o_vmodel0): $(s_vmodel0) $(o_std) $(o_system) $(o_global) $(o_fdtool) $(o_bicubic) $(o_readini) $(o_geomap)
$(o_vmodel1): $(s_vmodel1) $(o_std) $(o_global) $(o_readini)  $(o_pdebug)
$(o_vmodel2): $(s_vmodel2) $(o_std) $(o_system) $(o_global) $(o_fdtool) $(o_bicubic) $(o_readini) $(o_geomap) $(o_pdebug)
$(o_vmodel3): $(s_vmodel3) $(o_std) $(o_global) $(o_readini)  $(o_pdebug)

$(o_vmodel1r): $(s_vmodel1r) $(o_std) $(o_global) $(o_readini)  $(o_pdebug) $(o_rdrmed)
$(o_vmodel2r): $(s_vmodel2r) $(o_std) $(o_global) $(o_system) $(o_fdtool) $(o_bicubic) $(o_readini) $(o_geomap) $(o_pdebug) $(o_rdrmed) $(o_fdtool)
$(o_vmodel3r): $(s_vmodel3r) $(o_std) $(o_global) $(o_readini)  $(o_pdebug) $(o_rdrmed) $(o_fdtool)

##
## Module Compile Rule
##
.f90.o:
	$(FC) $(NCINC) $(NCFLAG) $(FFLAGS)  -o $@	  -c $<


##
## Targets
##

../../bin/swpc_sh.x: $(s_main) $(o_std) $(o_system) $(o_fdtool) $(o_global) $(o_kernel) $(o_getopt) $(o_source) $(o_medium) $(o_report) $(o_pwatch) $(o_snap) $(o_absorb)  $(o_gk) $(o_bicubic) $(o_daytim) $(o_sac) $(o_readini) $(o_geomap) $(o_absorb_p) $(o_absorb_c) $(o_vmodel1) $(o_vmodel2) $(o_vmodel3)  $(o_pdebug) $(o_rdrmed) $(o_vmodel1r) $(o_rdrmed) $(o_vmodel3r) $(o_vmodel2r) $(o_vmodel0) $(o_version) $(o_wav)
	$(FC) $(NCINC)  $(NCLIB)  $(FFLAGS) -o $@ $^ $(NETCDF)

clean:
	/bin/rm -f $(o_std) $(o_system) $(o_fdtool) $(o_global) $(o_kernel) $(o_getopt) $(o_source) $(o_medium) $(o_report) $(o_pwatch) $(o_snap) $(o_absorb) $(o_gk) $(o_bicubic) $(o_daytim) $(o_sac) $(o_readini) $(o_geomap)  $(o_absorb_p) $(o_absorb_c) $(o_vmodel1) $(o_vmodel2) $(o_vmodel3) $(o_pdebug) $(o_rdrmed) $(o_vmodel1r) $(o_rdrmed) $(o_vmodel3r) $(o_vmodel2r)  $(o_vmodel0) $(o_version) $(o_wav) ../../bin/swpc_sh.x
	/bin/rm -f *.mod
	/bin/rm -f *~
	/bin/rm -f *.L

clean2:
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *~
	/bin/rm -f *.L
