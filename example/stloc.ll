#                                                    -*- mode:sh -*-
# stloc.ll
#
# station location data by geographical format. 
# lines starting from '#' and blank lines are omitted. 
#
# zsw: controls station depth
#      'dep': use the depth
#      'fsb': locate one-grid below from the free surface/sea surface
#      'obb': locate one-grid below from the groud surface/seafloor
#      'oba': locate one-grid above from the groud surface/seafloor
#      'bd{i}' (i=0,...,9) i-th boundary interface
#
#     lon       lat    dep     stnm   zsw
# --------------------------------------
 139.7602   35.7183      0     st01    'obb'
 139.8602   35.7283      0     st02    'obb'
 139.8602   34.7283      0     st03    'obb'
