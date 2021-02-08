# ecRad Makefile - read the README file before editing

#############################
### --- CONFIGURATION --- ###
#############################

# Use the nf-config utility, if available, to set the NETCDF_INCLUDE
# and NETCDF_LIB flags
HAVE_NFCONFIG := $(shell nf-config --version 2> /dev/null)
ifdef HAVE_NFCONFIG
$(info *** Using nf-config to obtain NetCDF flags)
NETCDF_INCLUDE = $(shell nf-config --fflags)
NETCDF_LIB     = $(shell nf-config --flibs)
ifeq ($(shell nf-config --has-nc4),yes)
NETCDF4        = 1
endif
else
$(info *** nf-config not found)
endif

# make can be invoked using "make PROFILE=<prof>" in which case your
# local configuration parameters will be obtained from
# Makefile_include.<prof>
ifndef PROFILE
$(info *** No "PROFILE" variable provided, assuming "gfortran")
PROFILE = gfortran
endif

# Include a platform-specific makefile that defines FC, FCFLAGS and
# LIBS
include	Makefile_include.$(PROFILE)

# Check for presence of the NETCDF_INCLUDE and NETCDF_LIB flags
ifndef NETCDF_INCLUDE
$(info *** You may need to set NETCDF_INCLUDE manually)
endif
ifndef NETCDF_LIB
$(info *** You may need to set NETCDF_LIB manually)
endif

# Add single-precision flag if SINGLE_PRECISION=1 was given on the
# "make" command line
ifdef SINGLE_PRECISION
CPPFLAGS += -DSINGLE_PRECISION
endif

# Use shortwave reflectance-transmittance routine from ecRAD,
# which uses the same equations but should be faster
#ifdef USE_RTE_REFTRANS_SW
CPPFLAGS += -DUSE_RTE_REFTRANS_SW
#endif

# ------------- NEW FOR OPT -------------- 
ifdef FAST_EXPONENTIAL
CPPFLAGS += -DFAST_EXPONENTIAL
endif

ifdef BLOCK_DERIVED_TYPES
CPPFLAGS += -DBLOCK_DERIVED_TYPES
endif

# Use the GPTL timing library?
ifeq ($(GPTL_TIMING),1)
CPPFLAGS += -DUSE_TIMING
TIMING_INCLUDE += -I$(TIME_DIR)/include
LDFLAGS_TIME += -L$(TIME_DIR)/lib -Wl,-rpath=$(TIME_DIR)/lib
LIBS_TIMING += -lgptl  -rdynamic  
# Use the GPTL timing library with PAPI enabled to estimate computational intensity
else ifeq ($(GPTL_TIMING),2)
# GPTL + PAPI 
CPPFLAGS += -DUSE_TIMING  -DUSE_PAPI
TIMING_INCLUDE += -I$(TIME_DIR)/include 
LDFLAGS_TIME = -L$(TIME_DIR)/lib -Wl,-rpath=$(TIME_DIR)/lib
#LIBS_TIMING += -L$(TIME_DIR)/lib -Wl,-rpath=$(TIME_DIR)/lib -lgptl -rdynamic -lpapi
LIBS_TIMING +=  -lpthread -lgptl -rdynamic -lpapi
endif

# Use LIBDXSMM for small matrix-matrix multiplications in SPARTACUS
ifdef USE_LIBXSMM
LIBXSMM_INCLUDE = -I$(LIBXSMM_DIR)/include
CPPFLAGS += -DUSE_LIBXSMM
LDFLAGS_LIBXSMM = $(LIBXSMM_DIR)/lib/libxsmmf.a $(LIBXSMM_DIR)/lib/libxsmm.a -lpthread -lrt -ldl -lm 
endif
# ------------- NEW FOR OPT -------------- 
# If PRINT_ENTRAPMENT_DATA=1 was given on the "make" command line
# then the SPARTACUS shortwave solver will write data to fort.101 and
# fort.102
ifdef PRINT_ENTRAPMENT_DATA
CPPFLAGS += -DPRINT_ENTRAPMENT_DATA 
endif
# For backwards compatibility we allow the following as well
ifdef PRINT_ENCROACHMENT_DATA
CPPFLAGS += -DPRINT_ENTRAPMENT_DATA 
endif
# Allow the capability to write NetCDF4/HDF5 files, provided the code
# is compiled against the NetCDF4 library
ifdef NETCDF4
$(info *** Building with NetCDF4/HDF5 support)
CPPFLAGS += -DNC_NETCDF4
endif

# Consolidate flags
export FC
# ------------- NEW FOR OPT -------------- 
export FCFLAGS = $(WARNFLAGS) $(BASICFLAGS) $(CPPFLAGS) -I../include \
	$(OPTFLAGS) $(DEBUGFLAGS) $(LIBXSMM_INCLUDE) $(NETCDF_INCLUDE) $(TIMING_INCLUDE)  $(OMPFLAG)
export LIBS    = $(LDFLAGS) -L../lib  $(LDFLAGS_TIME) -lradsurf -lradiation -lutilities \
	-lifsrrtm  -ldrhook -lifsaux  -lrrtmgp -lrte -lneural -lnetcdff -lnetcdf $(FCLIBS) $(NETCDF_LIB) $(LIBS_TIMING) $(LDFLAGS_LIBXSMM) $(OMPFLAG) 
# ------------- NEW FOR OPT -------------- 
#export FCFLAGS = $(WARNFLAGS) $(BASICFLAGS) $(CPPFLAGS) -I../include \
#	$(OPTFLAGS) $(DEBUGFLAGS) $(NETCDF_INCLUDE) $(OMPFLAG)
#export LIBS    = $(LDFLAGS) -L../lib -lradsurf -lradiation -lutilities \
#	-lifsrrtm -ldrhook -lifsaux $(FCLIBS) $(NETCDF_LIB) $(OMPFLAG)
ifdef DR_HOOK
LIBS += -ldl -lrt
export CFLAGS = -g -O2
endif



#############################
### --- BUILD TARGETS --- ###
#############################

all: build

help:
	@echo "Usage:"
	@echo "  make PROFILE=<prof>"
	@echo "where <prof> is one of gfortran, pgi, intel or cray (see Makefile_include.<prof>)"
	@echo "Other arguments to make are:"
	@echo "  DEBUG=1              Compile with debug settings on and optimizations off"
	@echo "  SINGLE_PRECISION=1   Compile with single precision"
	@echo "  DR_HOOK=1            Compile with the Dr Hook profiling system"
	@echo "  test                 Run test cases in test directory"
	@echo "  clean                Remove all compiled files"

ifdef DR_HOOK
build: directories libifsaux librrtmgp libdrhook libutilities libifsrrtm libradiation libradsurf driver symlinks
else
build: directories libifsaux librrtmgp libdummydrhook libutilities libifsrrtm libradiation libradsurf driver symlinks
endif

# git cannot store empty directories so they may need to be created 
directories: mod lib
mod:
	mkdir -p mod
lib:
	mkdir -p lib

deps: clean-deps
	cd ifsaux && $(MAKE) deps
	cd ifsrrtm && $(MAKE) deps

clean-deps:
	rm -f include/*.intfb.h

libifsaux:
	cd ifsaux && $(MAKE)

librrtmgp:
	cd rte-rrtmgp-nn && $(MAKE)

libdrhook:
	cd drhook && $(MAKE)

libdummydrhook:
	cd drhook && $(MAKE) dummy

libutilities:
	cd utilities && $(MAKE)

libifsrrtm:
	cd ifsrrtm && $(MAKE)

libradiation:
	cd radiation && $(MAKE)

libradsurf:
	cd radsurf && $(MAKE)

driver:
	cd driver && $(MAKE)

symlinks: clean-symlinks
	cd practical && ln -s ../bin/ecrad
	cd practical && ln -s ../data

test: test_ifs test_i3rc

test_ifs:
	cd test/ifs && $(MAKE) test

test_i3rc:
	cd test/i3rc && $(MAKE) test

test_surface:
	cd test/surface && $(MAKE) test

clean: clean-tests clean-toplevel clean-utilities clean-mods clean-symlinks

clean-tests:
	cd test/ifs && $(MAKE) clean
	cd test/i3rc && $(MAKE) clean
	cd test/surface && $(MAKE) clean

clean-toplevel:
	cd rte-rrtmgp-nn && $(MAKE) clean
	cd radiation && $(MAKE) clean
	cd radsurf && $(MAKE) clean
	cd driver && $(MAKE) clean

clean-utilities:
	cd ifsaux && $(MAKE) clean
	cd utilities && $(MAKE) clean
	cd ifsrrtm && $(MAKE) clean
	cd drhook && $(MAKE) clean

clean-mods:
	rm -f mod/*.mod mod/*__genmod.f90

clean-symlinks:
	rm -f practical/ecrad practical/data

clean-autosaves:
	rm -f *~ .gitignore~ */*~ */*/*~

.PHONY: all build help deps clean-deps libifsaux libdrhook libutilities libifsrrtm \
	libradiation libradsurf driver symlinks clean clean-toplevel test test_ifs \
	test_i3rc test_surface clean-tests clean-utilities clean-mods clean-symlinks
