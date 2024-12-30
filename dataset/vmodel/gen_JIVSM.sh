#!/bin/bash

# generate NetCDF (GMT grd) files for Japan Integrated Velocity Structure Model
# This script assumes that the original distribution archive
#   * lp2012nankai-e_str.zip
#   * lp2012nankai-w_str.zip
# are downloaded and stored in the current directly.
#

function usage {
  cat << EOF
Usage: gen_JIVSM.sh [OPTION]
  -h         Display help
  -d dh      set grid spacing to dh [0.005]
  -r region  specify region in lon0/lon1/lat0/lat1 (GMT's -R) format [129/147/30/47]
  -p npara   number of parallel processings [8]
  -f fc      specify fortran compiler name [gfortran]
  -c         classic dataset (obsolete, for archiving purpose)
EOF
exit
}

# default parameters
dlon=0.005
dlat=0.005
region=129/147/30/47
NPARA=8 # number of parallel processing
FC=gfortran # fortran compiler for extrapolation
CLASSIC_MODE=false

while getopts ":d:r:p:ch" optkey; do
case "$optkey" in 
  d)
    dlon=${OPTARG}
    dlat=${OPTARG}
    ;;
  r)
    region=${OPTARG}
    ;;
  p)
    NPARA=${OPTARG}
    ;;
  c)
    CLASSIC_MODE=true
    ;;
  f) 
    FC=${OPTARG}
    ;;
  h) 
    usage
    ;;
esac
done

echo "REGION       $region"
echo "dlon/dlat    $dlon/$dlat"
echo "NPARA        $NPARA"
echo "CLASSIC MODE $CLASSIC_MODE"
echo "FC           $FC"


topo=ETOPO1_Bed_g_gmt4.grd


#
# File download
#

# Please un-comment the following section to newly download the dataset

# JIVSM original data

[ ! -d jivsm_data ] && mkdir jivsm_data

if "${CLASSIC_MODE}"; then 
  # JIVSM for shallower structure
  archive_jivsm_e=lp2012nankai-e_str.zip
  archive_jivsm_w=lp2012nankai-w_str.zip

  if [ ! -e ${archive_jivsm_e} ]; then
    curl  -o ${archive_jivsm_e} https://www.jishin.go.jp/main/chousa/12_choshuki/dat/nankai/lp2012nankai-e_str.zip
  fi
  unzip -o ${archive_jivsm_e} &
  if [ ! -e ${archive_jivsm_w} ]; then
    curl  -o ${archive_jivsm_w} https://www.jishin.go.jp/main/chousa/12_choshuki/dat/nankai/lp2012nankai-w_str.zip
  fi
  unzip -o ${archive_jivsm_w} &
  wait

  edata=Ejapan_path20111110.dat
  wdata=Wjapan_path20111110.dat
  mv -f $edata ./jivsm_data/
  mv -f $wdata ./jivsm_data/  
else

  archive_jivsm=JIVSM_V1R.zip
  if [ ! -e ${archive_jivsm} ]; then
    curl -o ${archive_jivsm} http://taro.eri.u-tokyo.ac.jp/saigai/models/JIVSM_V1R.zip
  fi
  unzip -o ${archive_jivsm}
  edata=JIVSMeast_V1R.dat
  wdata=JIVSMwest_V1R.dat
  mv -f $edata ./jivsm_data/
  mv -f $wdata ./jivsm_data/
  rm -f JIVSMproperties.xlsx
  rm -f comparison_of_versions.docx
  rm -f Tectono2009.pdf WCEE2008.pdf WCEE2012.pdf sample_JIVSM.f README.txt
  rm -f JIVSMeast_V1R.prm JIVSMwest_V1R.prm
fi

# Topography & bathymetry data for extrapotation

if [ ! -e ${topo} ]; then
  curl -o ETOPO1_Bed_g_gmt4.grd.gz  https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz
  gunzip ETOPO1_Bed_g_gmt4.grd.gz
fi


#
# velocity structure
#

echo "1  1.7 0.35 1.80 119  70
2  1.8 0.5  1.95 170 100
3  2.0 0.6  2.00 204 120
4  2.1 0.7  2.05 238 140
5  2.2 0.8  2.07 272 160
6  2.3 0.9  2.10 306 180
7  2.4 1.0  2.15 340 200
8  2.7 1.3  2.20 442 260
9  3.0 1.5  2.25 510 300
10 3.2 1.7  2.30 578 340
11 3.5 2.0  2.35 680 400
12 4.2 2.4  2.45 680 400
13 5.0 2.9  2.60 680 400
14 5.5 3.2  2.65 680 400
15 5.8 3.4  2.70 680 400
16 6.4 3.8  2.80 680 400
17 7.5 4.5  3.20 850 500
18 5.0 2.9  2.40 340 200
19 6.8 4.0  2.90 510 300
20 8.0 4.7  3.20 850 500
21 5.4 2.8  2.60 340 200
22 6.5 3.5  2.80 510 300
23 8.1 4.6  3.40 850 500" > jivsm_layer.dat

#
# split original data file
#


(( p = 0 ))
for (( i=1; i<=23; i++ ));do

    ii=`echo $i | awk '{printf("%.2d",$1)}'`
    ic=`echo $i | awk '{printf("%.2d",$1+2)}'` # column in the datafile

    echo "Extract JIVSM layer $ii"
    awk '{printf("%.5f %.5f %12.5f\n",$1,$2,-$'$ic')}' \
        jivsm_data/$wdata > jivsm_data/sw.$ii.dat &
    awk '{printf("%.5f %.5f %12.5f\n",$1,$2,-$'$ic')}' \
        jivsm_data/$edata> jivsm_data/ne.$ii.dat &
    (( p += 2 )) 
    if (( p == NPARA )); then
      (( p = 0 ))
      wait
    fi
done

wait


#
# generate grd data
#

lyr_name=( "DUMMY"    "GSURFACE" "BASEMENT" "BASEMENT" "BASEMENT" \
           "BASEMENT" "BASEMENT" "BASEMENT" "BASEMENT" "BASEMENT" \
           "BASEMENT" "BASEMENT" "BASEMENT" "BASEMENT" "UPCRUST1" \
           "UPCRUST2" "LOWCRUST" "C-MANTLE" "PHS-LYR2" "PHS-LYR3" \
           "PHS-MNTL" "PAC-LYR2" "PAC-LYR3" "PAC-MNTL" )

if "${CLASSIC_MODE}"; then 
  dir_jivsm=jivsm-o
else
  dir_jivsm=jivsm
fi
[ ! -d ${dir_jivsm} ] && mkdir ${dir_jivsm}

rm -f jivsm.lst

(( p = 0 ))
for (( i=1; i<=23; i++ )); do

    ii=`echo $i | awk '{printf("%.2d",$1)}'`
    echo "Convert layer $ii"

    vp=`awk '$1=='$i'{printf("%.2f",$2)}' jivsm_layer.dat`
    vs=`awk '$1=='$i'{printf("%.2f",$3)}' jivsm_layer.dat`
    ro=`awk '$1=='$i'{printf("%.2f",$4)}' jivsm_layer.dat`
    qp=`awk '$1=='$i'{printf("%.4d",$5)}' jivsm_layer.dat`
    qs=`awk '$1=='$i'{printf("%.4d",$6)}' jivsm_layer.dat`

    grdfile=JIVSM_${ii}_${lyr_name[$i]}_${ro}_${vp}_${vs}_${qp}_${qs}_.grd

    cat jivsm_data/sw.$ii.dat jivsm_data/ne.$ii.dat | \
        gmt nearneighbor -R${region} -S2k -I${dlon}/${dlat} -G${dir_jivsm}/${grdfile} &
    (( p += 1 ))

    if (( p == NPARA )); then 
      (( p = 0 ))
      wait
    fi

done
wait
rm -rif jivsm_data

#
# list file for FDM input
#

rm -f jivsm.lst
for (( i=1; i<=23; i++ )); do

    ii=`echo $i | awk '{printf("%.2d",$1)}'`
    vp=`awk '$1=='$i'{printf("%.2f",$2)}' jivsm_layer.dat`
    vs=`awk '$1=='$i'{printf("%.2f",$3)}' jivsm_layer.dat`
    ro=`awk '$1=='$i'{printf("%.2f",$4)}' jivsm_layer.dat`
    qp=`awk '$1=='$i'{printf("%.4d",$5)}' jivsm_layer.dat`
    qs=`awk '$1=='$i'{printf("%.4d",$6)}' jivsm_layer.dat`

    grdfile=JIVSM_${ii}_${lyr_name[$i]}_${ro}_${vp}_${vs}_${qp}_${qs}_.grd

    if [ $ii = "18" ]; then
      ilyr=1
    elif [ $ii = "21" ]; then
      ilyr=2
    else
      ilyr=0
    fi
    echo "'"$grdfile"'" "  $ro  $vp  $vs  $qp  $qs  $ilyr" >> jivsm.lst
done
rm -f jivsm_layer.dat


##
## ejivsm part
##

#
# Generate NetCDF (GMT grd) files for extended-version of jivsm
# This script requires to have jivsm grd data are already generated.
# It requires wider area topography data in grd format, which is defined
# in the first part of this file
#


if "${CLASSIC_MODE}"; then 
  dir_ejivsm=ejivsm-o
else
  dir_ejivsm=ejivsm
fi
[ ! -e ${dir_ejivsm} ] && mkdir ${dir_ejivsm}

gmt grdcut -R$region $topo -Gtopo.japan.grd
topo_org=`/bin/ls ${dir_jivsm}/*_01_*.grd`
topo_new=${dir_ejivsm}/e`basename $topo_org`
echo "TOPO     = $topo_org"
echo "TOPO_NEW = $topo_new"

# change topography data to be positive downward
gmt grdmath topo.japan.grd -1 MUL = topo.japan2.grd
gmt grdsample topo.japan2.grd -R$region -I$dlon/$dlat -G$topo_new

#
# differential depth between internal discontinuity and topography
# for shallow structure
#
cd code
${FC} m_std.f90 m_system.f90 extrap.f90 -o extrap.x
cd ..

for (( i=2; i<=13; i++ ));
do
    ii=`echo $i | awk '{printf("%.2d",$1)}'`

    grd=`/bin/ls ${dir_jivsm}/*_${ii}_*.grd`
    out=${dir_ejivsm}/e`basename $grd`

    jj=`echo $i | awk '{printf("%.2d",$1-1)}'`
    grd_up=`/bin/ls ${dir_jivsm}/*_${jj}_*.grd`
    out_up=${dir_ejivsm}/e`basename $grd_up`

    echo $ii $grd

    # first subtract the topography
    gmt grdmath $grd $topo_org SUB = tmp.${i}.grd
    gmt grd2xyz tmp.${i}.grd -bo > tmp.${i}.dat

    # extrapolation
    ./code/extrap.x tmp.${i}.dat tmp2.${i}.dat $dlon $dlat 129 147 30 47 1
    gmt surface -bi tmp2.${i}.dat -R$region -I$dlon/$dlat -Gtmp2.${i}.grd

    # add the (new) topography
    gmt grdmath tmp2.${i}.grd ${topo_new} ADD ${out_up} MAX = $out


done

for (( i=14; i<=23; i++ ));
do
    ii=`echo $i | awk '{printf("%.2d",$1)}'`
    grd=`/bin/ls ${dir_jivsm}/*_${ii}_*.grd`
    echo "GRD = $grd"
    out=${dir_ejivsm}/e`basename $grd`
    jj=`echo $i | awk '{printf("%.2d",$1-1)}'`
    grd_up=`/bin/ls ${dir_jivsm}/*_${jj}_*.grd`
    out_up=${dir_ejivsm}/e`basename $grd_up`

    # simple extrapolation in deeper structure
    if (( i>=18 && i<=20 )); then
      imode=2
    else
      imode=1
    fi
    gmt grd2xyz $grd -bo > tmp.$i.dat
    ./code/extrap.x tmp.${i}.dat tmp2.${i}.dat $dlon $dlat 129 147 30 47 $imode
    gmt surface -bi tmp2.$i.dat -R$region -I$dlon/$dlat -Gtmp.grd
    gmt grdmath tmp.grd ${out_up} MAX = $out

done
rm -f tmp*dat topo*grd tmp*grd



# velocity list file
rm -f ejivsm.lst
for f in ${dir_ejivsm}/*.grd
do
    g=` basename $f `
    ii=` echo $f | awk -F_ '{print $2}' `
    vp=` echo $f | awk -F_ '{print $5}' `
    vs=` echo $f | awk -F_ '{print $6}' `
    ro=` echo $f | awk -F_ '{print $4}' `
    qp=` echo $f | awk -F_ '{print $7}' `
    qs=` echo $f | awk -F_ '{print $8}' `

    if [ $ii = "18" ]; then
      ilyr=1
    elif [ $ii = "21" ]; then
      ilyr=2
    else
      ilyr=0
    fi
    echo "'"$g"'" "  $ro  $vp  $vs  $qp  $qs  $ilyr" >> ejivsm.lst
done

rm -f ${topo}
rm -f gmt.history
