
#Default setup. Overriden by DefaultSetup file
SE=gnu
DefaultSetup=DefaultSetup

PROGRAM ?= silam_v5_7

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

puhti_mpi:
	echo "SE=puhti_mpi" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

puhti:
	echo "SE=puhti" > $(DefaultSetup)
	echo Switching make to `cat $(DefaultSetup)`

DefaultSetup:
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

is4fires :
	cd ../build/ && $(MAKE) is4fires SETUP=$(SE)

dot:
	cd ../build/ && $(MAKE) dot  SETUP=$(SE)
