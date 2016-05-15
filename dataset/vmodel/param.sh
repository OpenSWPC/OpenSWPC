#!/bin/bash

#
# parameter settings for generating velocity model files
#

FC=gfortran

region=129/147/30/47

# grid spacing
#dlon=0.0125
#dlat=0.00833
dlon=0.005
dlat=0.005

# Topography grd file
topo=./ETOPO1_Bed_g_gmt4.grd
