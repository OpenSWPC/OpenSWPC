## ------------------------------------------------------------------------------------------------------------------------------ ##
##
## compile all SWPC codes
##
## @copyright
##   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
##
## ------------------------------------------------------------------------------------------------------------------------------ ##

# legacy option
ifdef build
arch=$(build)
endif

all:
	cd swpc_3d;  make arch=$(arch)
	cd swpc_sh;  make arch=$(arch)
	cd swpc_psv; make arch=$(arch)
	cd tools;    make arch=$(arch)

clean:
	cd swpc_3d;  make clean
	cd swpc_sh;  make clean
	cd swpc_psv; make clean
	cd tools;    make clean
