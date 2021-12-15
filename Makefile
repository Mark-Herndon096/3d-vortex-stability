# Common source directories
SRC_DIR = src
GSL_DIR = GSL_INTERFACE
# Need to compile GNU scientific library serparately and set path here
GSL_ROOT = /custom_builds/GSL
# Compiler (ifort, gfortran)
FC = ifort
CC = icc

export FC
export CC

# GSL include and library flags
GSL_INC   = -I${GSL_ROOT}/include
GSL_LIBS  = -L${GSL_ROOT}/lib 
GSL_LD    = -lgsl -lgslcblas -lm
GSL_FLAGS = ${GSL_INC} ${GSL_LIBS} ${GSL_LD}

# MKL LIBS (set MKL by sourcing intel OneAPI env prior to running Makefile)
MKL_INC   = -I${MKLROOT}/include/intel64/ilp64 -I"${MKLROOT}/include"
MKL_LIBS  = -L${MKL_ROOT}/lib/intel64
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
	$(FC) $(COMPFLAGS) -o $(EXEC_NAME) $(OBJECTS) $(LIBRARIES)

clean:
	$(MAKE) -C $(GSL_DIR) clean
	$(MAKE) -C $(SRC_DIR) clean
	rm $(EXEC_NAME)
