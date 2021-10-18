#################################################################
#
# Makefile for for silam version 4.2
#
# author: J. Vira / FMI
#
#################################################################

# The executable

-include Makefile.branch
PROGRAM ?= $(BRANCH)
# The setup variable: if defined, compiler options will be read from
# options.$(SETUP), and the executable will be named as $(PROGRAM).$(SETUP)
SETUP = 

SILAM_ROOT = ../
# directory with makefiles
BUILDDIR = $(SILAM_ROOT)build/

ENDING = $(if $(SETUP),.$(SETUP))$(if $(SUFFIX), .$(SUFFIX))
OPTIONS = $(BUILDDIR)/options/$(ENDING)
include $(BUILDDIR)/$(join options, $(ENDING))
PROGRAM := $(PROGRAM)$(ENDING)

# source code.
SRCDIR = $(SILAM_ROOT)source/
# object and .mod files
OBJDIR = $(SILAM_ROOT)obj/$(SETUP)/
# directory for the executable
BINDIR = $(SILAM_ROOT)bin/

# directory for a distribution package
DISTRDIR=$(SILAM_ROOT)/distr

# Build rules: if this is set to simple, every module is assumed to depend
# on its own source only. Caveat emptor.
RULES = full

# The location of main program for the model
MODEL_MAIN = silam_main
# The  location of main program for the standalone 3dvar
ASSIMILATOR_MAIN = assimilator
# The location of main program for the box model
BOX_MAIN = box_model

# The main program for grads-grib converter
GRADS2GRIB_MAIN = grads_2_grib2

# The main program for is4fires
IS4FIRES_MAIN = frp_main_v2_0

# build targets
all : silam

export

# Only the actual fortran compilation shall be parallel - otherwise
# we'd be finding the dependencies simultaneously with compilation, etc.
.NOTPARALLEL : 

silam : expand last.revision config 
	mkdir -p $(OBJDIR) $(BINDIR)	
	$(MAKE) -f$(BUILDDIR)/makefile.model

silam_only : 
	mkdir -p $(OBJDIR) $(BINDIR)	
	$(MAKE) -f$(BUILDDIR)/makefile.model

box : expand config_box last.revision
	mkdir -p $(OBJDIR) $(BINDIR)	
	$(MAKE) -f$(BUILDDIR)/makefile.box

grads2grib : expand config_grads2grib last.revision
	mkdir -p $(OBJDIR) $(BINDIR)	
	$(MAKE) -f$(BUILDDIR)/makefile.grads2grib 

is4fires : expand config_is4fires last.revision
	mkdir -p $(OBJDIR) $(BINDIR)	
	$(MAKE) -f$(BUILDDIR)/makefile.is4fires


clean:
	rm -f $(OBJDIR)*.o $(OBJDIR)*.mod $(BINDIR)$(PROGRAM)

config: 	
	ls -1 $(SRCDIR)*.f90 > $(BUILDDIR)/sourcelist
	perl $(BUILDDIR)/config.pl $(BUILDDIR)/sourcelist $(MODEL_MAIN) $(ARCH) \
$(if $(MPI_CHOICE), $(MPI_CHOICE), no_mpi)
	rm $(BUILDDIR)/sourcelist
	@echo Config OK.

config_box:
	ls -1 $(SRCDIR)*.f90 > $(BUILDDIR)/sourcelist
	perl $(BUILDDIR)/config.pl $(BUILDDIR)/sourcelist $(BOX_MAIN) $(ARCH) no_mpi
	rm $(BUILDDIR)/sourcelist
	@echo Config OK.

config_grads2grib:
	ls -1 $(SRCDIR)*.f90 > $(BUILDDIR)/sourcelist
	perl $(BUILDDIR)/config.pl $(BUILDDIR)/sourcelist $(GRADS2GRIB_MAIN) $(ARCH) \
$(if $(MPI_CHOICE), $(MPI_CHOICE), no_mpi)
	rm $(BUILDDIR)/sourcelist
	@echo Config OK.

config_is4fires:
	ls -1 $(SRCDIR)*.f90 > $(BUILDDIR)/sourcelist
	perl $(BUILDDIR)/config.pl $(BUILDDIR)/sourcelist $(IS4FIRES_MAIN) $(ARCH) \
$(if $(MPI_CHOICE), $(MPI_CHOICE), no_mpi)
	rm $(BUILDDIR)/sourcelist
	@echo Config OK.



tags: $(SRCDIR)*.f90
	cd $(SRCDIR) && $(MAKE) -f makefile.tags

lib : expand config
	$(MAKE) -f $(BUILDDIR)/makefile.lib

dot : expand config
#	dot -Gratio=1.4 -Granksep=0.4 -Gnodesep=0.3 -Nfontsize=18 -Tsvg silam.dot > $(PROGRAM).svg
	dot -Granksep=0.4 -Gnodesep=0.3 -Nfontsize=18 -Tsvg $(MODEL_MAIN).dot > $(PROGRAM).svg
	@rsvg-convert -f pdf -o $(PROGRAM).pdf $(PROGRAM).svg || echo 'Sorry, failed to make pdf out of svg!'

# If no f90 file has changed since last running check_revision, no need to do anything.
# If something has been changed, check_revision will run and generate svndiff.mod.f90. 
# revision.mod.f90 is only generated if revision is changed.
last.revision : $(filter-out %svndiff.mod.f90, $(wildcard $(SRCDIR)/*.f90))
# if the filter-out is removed, svndiff will be always regenerated because it's newer than last.revision

#	which svn && svn info ../source > last.revision.tmp && mv last.revision.tmp last.revision || rm -f last.revision.tmp
	python $(BUILDDIR)/check_revision $(SRCDIR)


compact :
	rm -rf $(SILAM_ROOT)/obj/*
#	cd $(SILAM_ROOT) && tar cvf source.tar source && gzip source.tar
	cd $(SILAM_ROOT) && zip -r source.zip source && rm -rf source

expand :
	cd $(SILAM_ROOT) && { if [ -e source.zip ]; then unzip source.zip && rm source.zip; fi }
#	cd $(SILAM_ROOT) && tar xvf source.tar.gz

up : expand
	cd $(SRCDIR) && svn up

KPP_DIR=$(shell readlink -f ../kpp/)
KPP=$(KPP_DIR)/bin/kpp
# Create chemistries using KPP. KPP_DIR must be defined!
kpp : $(KPP)
$(KPP) : 
	cd $(KPP_DIR) && $(MAKE)

cbm4 : 
	KPP_HOME=$(KPP_DIR) python kpp2silam.py $(KPP2SILAM_FLAGS) --kpp-dir=$(KPP_DIR)/cbm4 \
	--silam-dir=$(SRCDIR) cbm4 # if start using photolysis for cbm4, use the line below
#	--silam-dir=$(SRCDIR) -u photolysis --photoarr cbm4

cbm4_SOA : 
	KPP_HOME=$(KPP_DIR) python kpp2silam.py --transf-id=5011 $(KPP2SILAM_FLAGS) --kpp-dir=$(KPP_DIR)/cbm4_SOA \
	--silam-dir=$(SRCDIR) cbm4_SOA

cbm42_strato :
	KPP_HOME=$(KPP_DIR) python kpp2silam.py --transf-id=5008 $(KPP2SILAM_FLAGS) --kpp-dir=$(KPP_DIR)/cbm42_strato \
	--silam-dir=$(SRCDIR) -u photolysis cbm42_strato

cbm42_strato_SOA : 
	KPP_HOME=$(KPP_DIR) python kpp2silam.py --transf-id=5012 $(KPP2SILAM_FLAGS) --kpp-dir=$(KPP_DIR)/cbm42_strato_SOA \
	--silam-dir=$(SRCDIR) -u photolysis cbm42_strato_SOA


# Create a source distribution package. Libraries are not
# included. The check_revision script is replace with a dummy, since
# the .svn directory is not included either.

DISTRO = $(DISTRDIR)/$(PROGRAM)
distr : config
	python $(BUILDDIR)/check_revision --distr $(SRCDIR)
	rm -rf  $(DISTRO)
	mkdir -p $(DISTRO) && cd $(DISTRO) && mkdir -p source bin build obj ini kpp
	cp -vp makefile makefile.model config.pl last.revision library.objects options.*[!~] $(DISTRO)/build
	touch $(DISTRO)/build/check_revision
	cp -vp ../source/Makefile ../source/*.f90 $(DISTRO)/source
	cp -vrp ../ini/*  $(DISTRO)/ini
	cp -vrp ../kpp/*  $(DISTRO)/kpp
	cd  $(DISTRO)/kpp && 	make distclean

	cp README $(DISTRO)
	cd $(DISTRDIR) && tar cvjf $(PROGRAM).tar.bz2 $(PROGRAM)
#	rm ../source/svndiff.mod.f90 # Do not keep distros svndiff.mod.f90
	touch ../source/silam_main.f90 #Touch something to force regeneration of svndiff.mod.f90
	echo "module svndiff" >  ../source/svndiff.mod.f90 # Make config happy
	echo "use globals" >>  ../source/svndiff.mod.f90 # Must be present, so config works correctly


# Make documentation. So far just the data assimilation "readme" document.
docs : 
	cd ../doc && $(MAKE)

