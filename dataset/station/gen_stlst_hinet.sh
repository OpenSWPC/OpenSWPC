#!/bin/bash

# generate Hi-net station location list file for OpenSWPC

# Before executing this script, please obtain following station data file
# "NIED_SeismicStation_YYYYMMDD.csv"
#    from http://www.hinet.bosai.go.jp/st_info/detail/dataset.php
#

lst=`/bin/ls -t NIED_SeismicStation_*csv | head -1`
awk -F, 'NR>1{printf("%10.4f %10.4f %10.4f    %s    obb\n", $9, $8, -$13/1000.,$3)}' $lst | sort | uniq > hinet.stlst
