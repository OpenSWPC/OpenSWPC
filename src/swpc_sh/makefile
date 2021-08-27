## ------------------------------------------------------------------------------------------------------------------------------ ##
##
## makefile for SWPC_SH
##
## @copyright
##   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
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
.SUFFIXES: .F90 .o                  # only uses .F90 and .o file in makefile


##
## Compiler & Compile Option
##
include ../shared/makefile.arch


##
## Source Code Specification
##



## kernel module
s_global  = m_global.F90
s_kernel  = m_kernel.F90
s_source  = m_source.F90
s_fdtool  = m_fdtool.F90
s_output  = m_output.F90
s_report  = m_report.F90
s_medium  = m_medium.F90
s_absorb_c  = m_absorb_c.F90
s_absorb_p  = m_absorb_p.F90
s_absorb  = m_absorb.F90
s_ckprst  = m_ckprst.F90
s_vmodel0 = m_vmodel_user.F90
s_vmodel1 = m_vmodel_uni.F90
s_vmodel2 = m_vmodel_grd.F90
s_vmodel3 = m_vmodel_lhm.F90
s_vmodel1r = m_vmodel_uni_rmed.F90
s_vmodel2r = m_vmodel_grd_rmed.F90
s_vmodel3r = m_vmodel_lhm_rmed.F90

## common library modules
s_system  = m_system.F90
s_std     = m_std.F90
s_getopt  = m_getopt.F90
s_pwatch  = m_pwatch.F90
s_gk      = m_gk.F90
s_bicubic = m_bicubic.F90
s_daytim  = m_daytim.F90
s_sac    = m_sac.F90
s_readini = m_readini.F90
s_geomap  = m_geomap.F90
s_pdebug  = m_pdebug.F90
s_rdrmed  = m_rdrmed.F90

## main program
s_main      = main.F90


##
## Objectd Filenames
##

o_std     = $(s_std:.F90=.o)
o_pdebug  = $(s_pdebug:.F90=.o)
o_getopt  = $(s_getopt:.F90=.o)
o_system  = $(s_system:.F90=.o)
o_pwatch  = $(s_pwatch:.F90=.o)
o_gk      = $(s_gk:.F90=.o)
o_bicubic = $(s_bicubic:.F90=.o)
o_daytim  = $(s_daytim:.F90=.o)
o_sac    = $(s_sac:.F90=.o)
o_readini = $(s_readini:.F90=.o)
o_geomap  = $(s_geomap:.F90=.o)
o_rdrmed  = $(s_rdrmed:.F90=.o)

o_global  = $(s_global:.F90=.o)
o_kernel  = $(s_kernel:.F90=.o)
o_source  = $(s_source:.F90=.o)
o_fdtool  = $(s_fdtool:.F90=.o)
o_output  = $(s_output:.F90=.o)
o_report  = $(s_report:.F90=.o)
o_medium  = $(s_medium:.F90=.o)
o_absorb  = $(s_absorb:.F90=.o)
o_absorb_c  = $(s_absorb_c:.F90=.o)
o_absorb_p  = $(s_absorb_p:.F90=.o)
o_ckprst  = $(s_ckprst:.F90=.o)
o_vmodel0 = $(s_vmodel0:.F90=.o)
o_vmodel1 = $(s_vmodel1:.F90=.o)
o_vmodel2 = $(s_vmodel2:.F90=.o)
o_vmodel3 = $(s_vmodel3:.F90=.o)
o_vmodel1r = $(s_vmodel1r:.F90=.o)
o_vmodel2r = $(s_vmodel2r:.F90=.o)
o_vmodel3r = $(s_vmodel3r:.F90=.o)

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

$(o_output):  $(s_output)  $(o_std) $(o_global) $(o_pwatch) $(o_fdtool) $(o_sac) $(o_readini) $(o_geomap) $(o_pdebug)
$(o_global):  $(s_global)  $(o_std) $(o_fdtool) $(o_pwatch) $(o_system) $(o_pwatch) $(o_daytim) $(o_readini) $(o_pdebug)
$(o_fdtool):  $(s_fdtool)  $(o_std) $(o_pdebug)
$(o_source):  $(s_source)  $(o_std) $(o_global) $(o_fdtool) $(o_pwatch) $(o_readini) $(o_geomap) $(o_pdebug)
$(o_kernel):  $(s_kernel)  $(o_std) $(o_global) $(o_pwatch) $(o_medium) $(o_readini) $(o_pdebug)
$(o_report):  $(s_report)  $(o_std) $(o_global) $(o_kernel) $(o_pwatch) $(o_readini) $(o_pdebug)
$(o_medium):  $(s_medium)  $(o_std) $(o_global) $(o_pwatch) $(o_readini) $(o_vmodel1) $(o_vmodel2) $(o_vmodel3) $(o_vmodel4) $(o_pdebug) $(o_vmodel1r) $(o_vmodel3r) $(o_vmodel2r) $(o_vmodel0)
$(o_absorb):  $(s_absorb)  $(o_std) $(o_global) $(o_pwatch) $(o_readini) $(o_absorb_p) $(o_absorb_c)  $(o_pdebug)
$(o_absorb_c):  $(s_absorb_c)  $(o_std) $(o_fdtool) $(o_global) $(o_pdebug)
$(o_absorb_p):  $(s_absorb_p)  $(o_std) $(o_fdtool) $(o_global) $(o_readini) $(o_pdebug)
$(o_ckprst):  $(s_ckprst)  $(o_std) $(o_system) $(o_global) $(o_kernel) $(o_medium) $(o_absorb) $(o_report) $(o_source) $(o_pwatch) $(o_readini) $(o_output) $(o_pdebug)
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
.F90.o:
	$(FC) $(NCINC) $(NCFLAG) $(FFLAGS)  -o $@	  -c $<


##
## Targets
##

../../bin/swpc_sh.x: $(s_main) $(o_std) $(o_system) $(o_fdtool) $(o_global) $(o_kernel) $(o_getopt) $(o_source) $(o_medium) $(o_report) $(o_pwatch) $(o_output) $(o_absorb)  $(o_gk) $(o_bicubic) $(o_daytim) $(o_sac) $(o_ckprst) $(o_readini) $(o_geomap) $(o_absorb_p) $(o_absorb_c) $(o_vmodel1) $(o_vmodel2) $(o_vmodel3)  $(o_pdebug) $(o_rdrmed) $(o_vmodel1r) $(o_rdrmed) $(o_vmodel3r) $(o_vmodel2r) $(o_vmodel0)
	$(FC) $(NCINC)  $(NCLIB)  $(FFLAGS) -o $@ $^ $(NETCDF)

clean:
	/bin/rm -f $(o_std) $(o_system) $(o_fdtool) $(o_global) $(o_kernel) $(o_getopt) $(o_source) $(o_medium) $(o_report) $(o_pwatch) $(o_output) $(o_absorb) $(o_gk) $(o_bicubic) $(o_daytim) $(o_sac) $(o_ckprst) $(o_readini) $(o_geomap)  $(o_absorb_p) $(o_absorb_c) $(o_vmodel1) $(o_vmodel2) $(o_vmodel3) $(o_pdebug) $(o_rdrmed) $(o_vmodel1r) $(o_rdrmed) $(o_vmodel3r) $(o_vmodel2r)  $(o_vmodel0) ../../bin/swpc_sh.x
	/bin/rm -f *.mod
	/bin/rm -f *~
	/bin/rm -f *.L

clean2:
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *~
	/bin/rm -f *.L
