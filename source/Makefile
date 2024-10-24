
#Default setup. Overriden by DefaultSetup file
SE=gnu
DefaultSetup=DefaultSetup

# Include if exists
-include ../build/Makefile.branch

PROGRAM=$(BRANCH)

include $(DefaultSetup)



#Fastest make by defult
all:   $(DefaultSetup)
	echo Using RULES=simple SETUP=$(SE)
#	cd ../build && $(MAKE) RULES=simple 
	cd ../build && $(MAKE) RULES=simple SETUP=$(SE)
	rm -f ../bin/$(PROGRAM)
	ln -s $(PROGRAM).$(SE) ../bin/$(PROGRAM)

profiletest:
	cd ../build && $(MAKE) -f makefile.generic SETUP=craygnu_debug PROGRAM=profiletest

gnumpi:
	echo "SE=gnumpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

craygnu_mpi:
	echo "SE=craygnu_mpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

craygnu_debug:
	echo "SE=craygnu_debug" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

craygnu:
	echo "SE=craygnu" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

craygnu_mpi_debug:
	echo "SE=craygnu_mpi_debug" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

gnu_debug:
	echo "SE=gnu_debug" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`
gnumpi_debug:
	echo "SE=gnumpi_debug" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

gnu:
	echo "SE=gnu" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

newcray:
	echo "SE=newcray" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

newcray_mpi:
	echo "SE=newcray_mpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

lumimpi:
	echo "SE=lumimpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

lumi:
	echo "SE=lumi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

mahti_mpi:
	echo "SE=mahti_mpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

mahti:
	echo "SE=mahti" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

puhti_mpi:
	echo "SE=puhti_mpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

puhti:
	echo "SE=puhti" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

tuimai:
	echo "SE=tuimai" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

tuima_impi:
	echo "SE=tuima_impi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

tuima_mpi:
	echo "SE=tuima_mpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

tuima:
	echo "SE=tuima" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

ec:
	echo "SE=ec" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

ecmpi:
	echo "SE=ecmpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

$(DefaultSetup):
	echo "SE=$(SE)" > $(DefaultSetup)	
	echo Switching make to `cat $(DefaultSetup)`


conf:
	cd ../build && make config

full :
	echo Using -j 8 SETUP=$(SE)
#	cd ../build && $(MAKE) -j 8
	cd ../build && $(MAKE) -j 8 SETUP=$(SE)
	rm -f ../bin/$(PROGRAM)
	ln -s $(PROGRAM).$(SE) ../bin/$(PROGRAM)

box : 
	cd ../build && $(MAKE) box -j 8 SETUP=$(SE)

tags : *.f90
	etags --ignore-case-regex='/ *interface *.*/' --ignore-case-regex='/ *integer,.*public.*/' --ignore-case-regex='/ *real,.*public.*/' \
--ignore-case-regex='/ *type [a-z_]*/' *.f90

ctags : *.f90
	ctags -R --fortran-kinds=+i

up : 
	cd ../build/ && $(MAKE) up

clean:
	cd ../build && $(MAKE) clean  SETUP=$(SE)

what:
	cat DefaultSetup

distr:
	cd ../build/ && $(MAKE) distr

grads2grib:
	cd ../build/ && $(MAKE) grads2grib -j 8  SETUP=$(SE)

ascii2nc:
	cd ../build/ && $(MAKE) ascii2nc -j 8  SETUP=$(SE)

is4fires :
	cd ../build/ && $(MAKE) is4fires SETUP=$(SE)

dot:
	cd ../build/ && $(MAKE) dot  SETUP=$(SE)

cbm% : 
	cd ../build/ && $(MAKE) $@  SETUP=$(SE)
