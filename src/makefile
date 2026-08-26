## ------------------------------------------------------------------------------------------------------------------------------ ##
##
## compile all SWPC codes
##
## @copyright
##   Copyright 2013-2026 Takuto Maeda. All rights reserved. This project is released under the MIT license.
##
## ------------------------------------------------------------------------------------------------------------------------------ ##

# legacy option
ifdef build
arch=$(build)
endif

.PHONY: all swpc_3d swpc_sh swpc_psv tools clean
export arch

all: swpc_3d swpc_sh swpc_psv tools

swpc_3d: 
	$(MAKE) -C swpc_3d arch=$(arch)

swpc_psv:
	$(MAKE) -C swpc_psv arch=$(arch)

swpc_sh:
	$(MAKE) -C swpc_sh arch=$(arch)

tools:
	$(MAKE) -C tools arch=$(arch)

clean: clean_3d clean_sh clean_psv clean_tools
	rm -f shared/*.o shared/*.mod

clean_3d:
	$(MAKE) -C swpc_3d clean

clean_psv:
	$(MAKE) -C swpc_psv clean

clean_sh:
	$(MAKE) -C swpc_sh clean


clean_tools:
	$(MAKE) -C tools clean
