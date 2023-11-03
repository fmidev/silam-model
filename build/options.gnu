ARCH = linux_gnu

F90C = gfortran
#TESTOPTIONS = -D ONES_TEST # Init massmap from cellmass before first adv step 
TESTOPTIONS = 
# Flags for optimization, etc.
OPTIMIZATION = -O3 -fopenmp -ffast-math
OPENMP = 
DEBUG = -g -fbacktrace # -finit-real=snan   -fcheck=all      
FIXED = -ffixed-form
PREPROCESS = -x f95-cpp-input  -D WITH_BZIP2 -DUSE_PROJ6
FFLAGS = $(OPTIMIZATION) $(DEBUG) $(TESTOPTIONS) $(INCLUDE) -ffree-line-length-none
#FFLAGS = $(OPTIMIZATION) $(DEBUG) $(TESTOPTIONS)  $(INCLUDE) -ffree-line-length-none -fbacktrace -msse4.2 #-fdefault-real-8
INCLUDE = -I$(OBJDIR) -I/usr/include

# netcdf install directory
NCDIR =
# grib api install directory
GRDIR =

# Libraries passed to ifort
SILAM_LIBS = 
# Libraries given as linker options
LFLAGS = -fopenmp -llapack -lblas -llapack -leccodes  -leccodes_f90 -lnetcdff -lnetcdf  -lbz2 -lproj
