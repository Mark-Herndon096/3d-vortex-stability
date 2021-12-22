# Makefile for vortex stability solver
# Written by: Mark Herndon
# Lehigh University: Department of Mechanical
# Engineering and Mechanics

# Common source directories
SRC_DIR = src
GSL_DIR = GSL_INTERFACE
# Need to compile GNU scientific library serparately and set path here
GSL_ROOT = /opt/GSL
#GSL_ROOT = /custom_builds/GSL
export GSL_ROOT
# Compiler (ifort, gfortran)
FC = ifort
CC = icc

export FC
export CC

MODDIR = .mod

ifneq ($(MODDIR),)
  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
  FCFLAGS+= -module $(MODDIR)
endif

# GSL include and library flags
GSL_INC   = -I${GSL_ROOT}/include
GSL_LIBS  = -L${GSL_ROOT}/lib 
GSL_LD    = -lgsl -lgslcblas -lm
GSL_FLAGS = ${GSL_INC} ${GSL_LIBS} ${GSL_LD}

# MKL LIBS (set MKL by sourcing intel OneAPI env prior to running Makefile)
MKL_INC   = -I${MKLROOT}/include/intel64/ilp64 -I"${MKLROOT}/include"
MKL_LIBS  = -L${MKLROOT}/lib
MKL_LD    = -lmkl_blas95_ilp64 -lmkl_lapack95_ilp64 -lmkl_intel_ilp64 \
	    -lmkl_sequential -lmkl_core -lpthread -lm -ldl
MKL_FLAGS = ${MKL_INC} ${MKL_LIBS} ${MKL_LD}

# Include directories
# Libraries 

LIBRARIES = ${GSL_FLAGS} ${MKL_FLAGS}

COMMONFLAGS = -r8 -traceback -qopenmp
PRODFLAGS = -O3


COMPFLAGS = ${COMMONFLAGS} ${PRODFLAGS}

export COMPFLAGS
export LIBRARIES

# Executable name
EXEC_NAME = vortex_solver.exe

# Object list
OBJECTS = $(GSL_DIR)/special_function_wrapper.o   \
	  $(GSL_DIR)/special_function_interface.o \
	  $(SRC_DIR)/mod_global.o                 \
	  $(SRC_DIR)/mod_numerical_routines.o     \
	  $(SRC_DIR)/mod_io.o                \
	  $(SRC_DIR)/main.o
solver:
	$(MAKE) -C $(GSL_DIR) gsl_objs
	$(MAKE) -C $(SRC_DIR) src_objs
	$(FC) $(COMPFLAGS) $(FCFLAGS) -o $(EXEC_NAME) $(OBJECTS) $(LIBRARIES)

clean:
	$(MAKE) -C $(GSL_DIR) clean
	$(MAKE) -C $(SRC_DIR) clean
	rm -rf $(EXEC_NAME) .mod
