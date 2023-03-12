# Silam public distribution 

This is a full-featured source code of the Silam model
with striped revision history

## BUILDING SILAM v5.8 in Ubuntu 18.04, 20.04, 22.04

Get the source code

`$ git clone https://github.com/fmidev/silam-model.git`

Install needed packages

`$ sudo apt install make python3 gfortran libeccodes-dev libnetcdf-dev libnetcdff-dev liblapack-dev libblas-dev libbz2-dev libproj-dev`

Compile the binary

`$ cd silam-model/source/`

`$ make gnu`

`$ make`

The latter command should create a binary in ../bin

Now run the binary

`$ ../bin/silam_v5_7pub.gnu`

Silam should run and complain about missing silam.ini.
Then the binary is ready to use! You can test it with a test case from
https://github.com/fmidev/silam-toypoint

In Ubuntu 20.04 or 22.04 one might get an error message like

`grib_api.mod not found`

There are several versions of gfortran available in 20.04 and 22.04. Some libraries
have headers in a location specific for "module version", where gfortran can't find them:

`/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/`

others have them in gfortran-version-specific locarion

`/usr/lib/gcc/x86_64-linux-gnu/10/finclude/`

or 

`/usr/lib/gcc/x86_64-linux-gnu/19/finclude/`

A workaround would be to either explicitly call fortran with `-I
/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15` (by adding this option to
FFLAGS in the file build/options.gnu).

or add symlinks to the needed .mod files to your gfortran include
directory (as root).

`# ln -s /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/grib_api.mod /usr/lib/gcc/x86_64-linux-gnu/10/finclude/`

or

`# ln -s /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/grib_api.mod /usr/lib/gcc/x86_64-linux-gnu/9/finclude/`

or wherever your fortran can find them...
The issue has been reported at
https://bugs.launchpad.net/ubuntu/+source/gcc-defaults/+bug/1883855 .
Please consider confirming the bug if it affects you.





## Instructions for Centos

(contributed by Lars Örtegren <leo@apertum.se>)

### For build server

Enable Repo EPEL

`$ dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm -y`

Enable repo PowerTools

`$ yum install dnf-plugins-core`
`$ yum config-manager --set-enabled PowerTools`

Install libs

`$ dnf install zlib zlib-devel netcdf netcdf-devel netcdf-fortran-openmpi-devel jasper-libs jasper-devel lapack lapack-devel eccodes eccodes-devel`

Add to options.gnu:

`INCLUDE = ... -I/usr/lib64/gfortran/modules -I/usr/lib64/gfortran/modules/openmpi`

Add path to NETCDF shared libraries (as root):

`$ echo "/usr/lib64/openmpi/lib" > /etc/ld.so.conf/netcdf.conf`
`$ lsconfig -v`


In source directory

`$ make gnu`
`$ make`

### For run server 

Enable Repo EPEL

`$ dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm -y`

Enable repo PowerTools

`$ dnf install dnf-plugins-core`
`$ dnf config-manager --set-enabled PowerTools`

Install packages

`$ dnf install lapack netcdf netcdf-fortran netcdf-fortran-openmpi eccodes`





## Old RAEDME file (kept for historical reasons)

==============Seems to be somewhat outdated======================
* BUILDING SILAM v5.1
 
** Prerequiseites

*** Compilers
    
The following compilers and platforms have been tested on Linux:
- Intel Fortran v. 12, x86-64. Works with SILAM v5.0. Runtime issues 
  with SILAM 4.5.5.
- Intel Fortran v. 11, ia32 and x86-64 platforms. Recently compilation 
  issues with SILAM 5.2.
- Intel Fortran v. 10, ia32, ia64 and x86-64. Both were used with SILAM
  v4.5. Seems to have sporadic compilation issues in the
  ini_boundary_conditions module.
- gfortran < 4.4.5, various platforms. Compiles, but the executables
  do not work. Do not use. Early versions will not even compile the source
  because they lack the ISO C-binding module.
- gfortran 4.4.5, ia32. Seems to work.
- gfortran 4.5.3 and 4.6.1, x86-64. Works.
- gfortran 4.8.2, x86-64. Works.
- Pathscale, some version, x86-64. Works but has memory issues at least
  with the standard settings.
- Portland Group, some version, x86-64. Does not compile, has various issues.

Gfortran 4.5.0 on 64-bit Mac OS seems to work.

The build system requires GNU Make, Perl, and if used with version
control, Subversion and Python.

*** Libraries

The following freely available I/O libraries are required:
- NetCDF 3.6.1 or later, download from UCAR or possibly install with a
  package manager.
- GRIB API 1.8.0 or later, download from ECMWF, or possibly install with a
  package manager.
- Jasper. Required by GRIB API, usually available from package
  managers.
- A LAPACK implementation. Typically already installed.

In Ubuntu (at least in 14.04) these libraries can be found from package manager:
libgrib-api-dev libnetcdf-dev 

These libraries are installed following the standard procedure
(./configure, make, make install) but make sure to use the same
compiler and version (at least for fortran).

** Compiling SILAM

The makefile referes the setup files, which are provided for the Intel
(options.intel), Intel on Cray (options.cray) and for the GNU
compilers (options.gnu). The options file needs to be modified with
the correct file paths, including
- the -I flags for the header files (variable INCLUDE)
- the locations of the libraries or the -L and -l flags (variables
  NCDIR, LIBDIR, GRDIR, SILAM_LIBS, LFLAGS).

If called without arguments, make will look for file named simply
'options'. This can be a link to any of the specific option
files, or otherwise type e.g.

make SETUP=gnu

to use the file 'options.gnu'. When the setup is not the default, its
name will be appended to executable name. The different setups have
separate directories for the intermediate files.

*** Other compilation options

The options files in this directory are for Intel release (with
full optimization), Intel debug and for Gfortran. 

In addition to the SETUP variable, the make command can be given the
RULES variable which can be either 'full' (default) or 'simple'. The
former implies normal compilation resolving all inter-module
dependencies. The latter defines each module to only depend on its
source code. In this case, only the modules that changed are compiled,
but not the dependent modules. 

Additional build targets are 
- silam: build the model without first scanning the dependencies
- revision: update the revision.f90 module (in non-svn distributions
  this is replaced with a dummy)
- lib: create a library of the objects listed in lib_codes
- dot: visualize the module dependencies using the graphviz software
  (dot). The ouput will be in svg format.
- distr: create a .tar.gz file containing the sources and the build
  files.



  
