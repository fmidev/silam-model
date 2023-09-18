ARCH = linux_gnu

F90C = gfortran
# Flags for optimization, etc.
OPTIMIZATION = -O3 -fopenmp -ffast-math
DEBUG = -g  -fbacktrace  #-fcheck=all -fbacktrace   -finit-real=snan  #-DDEBUG
FIXED = -ffixed-form
PREPROCESS = -x f95-cpp-input  -D WITH_BZIP2 -DUSE_PROJ6
FFLAGS = $(OPTIMIZATION) $(DEBUG) $(TESTOPTIONS) $(INCLUDE) -ffree-line-length-none -mno-avx
INCLUDE = -I$(OBJDIR) -I/usr/include $(NETCDF4_INCLUDE) $(ECCODES_INCLUDE)

# Libraries given as linker options
LFLAGS = -fopenmp  -lbz2  $(NETCDF4_LIB) $(ECCODES_LIB) $(proj_LIB) $(OPENBLAS_LIB)
