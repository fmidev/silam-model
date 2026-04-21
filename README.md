# Silam public distribution 

This is a full-featured source code of the Silam model
with striped revision history

## BUILDING SILAM v6_1 in Ubuntu  20.04, 22.04, 24.04

Get the source code

`$ git clone https://github.com/fmidev/silam-model.git`

Install needed packages

`$ sudo apt install make python3 gfortran libeccodes-dev libnetcdf-dev libnetcdff-dev liblapack-dev libblas-dev libbz2-dev libproj-dev`

Compile the binary

`$ cd silam-model/source/`

`$ make gnu`

`$ make`

The latter command might require a second try. Finally, it should create a binary in ../bin,
that can be launched as:

`$ ../bin/silam_v6_1pub.gnu`

Silam should run and complain about missing silam.ini.
Then the binary is ready to use! You can test it with a test case from
https://github.com/fmidev/silam-toypoint



## Troubleshooting

In Ubuntu one might get an error message like

`grib_api.mod not found`

There are several versions of gfortran available in 20.04. Some libraries
have headers in a location specific for "module version", where gfortran can't find them:

`/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/`

Others have them in gfortran-version-specific location

`/usr/lib/gcc/x86_64-linux-gnu/10/finclude/`

or 

`/usr/lib/gcc/x86_64-linux-gnu/19/finclude/`

A workaround would be to either explicitly call FORTRAN with `-I
/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15` (by adding this option to
FFLAGS  in build/options.gnu).

or add symlinks to the needed .mod files to your gfortran include
directory (as root).

`# ln -s /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/grib_api.mod /usr/lib/gcc/x86_64-linux-gnu/10/finclude/`

or

`# ln -s /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/grib_api.mod /usr/lib/gcc/x86_64-linux-gnu/9/finclude/`

The exact command depends on your gcc version, which can be checked with 
`$ gfortran -v`

The issue has been reported at
https://bugs.launchpad.net/ubuntu/+source/gcc-defaults/+bug/1883855.
Please consider confirming the bug if it affects you.



  
