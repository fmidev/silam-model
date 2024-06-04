MODULE max_sizes_and_limits

  ! Description:
  ! This module contains maximum sizes of certain arrays, which are
  ! not handled dynamically. Also other kind of maximum limits
  ! regarding memory or cpu usage are found here.
  !
  ! Author: Mika Salonoja, FMI email Mika.Salonoja@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90

  IMPLICIT NONE

  !----------------------------------------------------------------------
  !
  ! Limits species array sizes
  !
  integer, public, parameter :: maxSizeModes = 128   ! Maximum number of modes allowed
  integer, public, parameter :: max_species = 500 ! max number of species allowed

  integer, public, parameter :: max_descriptors_in_source = 100 ! max descriptors in area source
  !----------------------------------------------------------------------

  INTEGER, PARAMETER, PUBLIC :: max_met_srcs = 10 ! max number of
  ! diffferent data sources used in single model run
 
  INTEGER, PARAMETER, PUBLIC :: max_met_files = 10 ! max number of
  ! diffferent files containing data for a single time 
 
  INTEGER, PARAMETER, PUBLIC :: max_levels = 200 ! the maximum
  ! number of level-value-pairs or fields in a field_3d

  INTEGER, PARAMETER, PUBLIC :: max_quantities = 200 ! max number of
  ! different quantities in memory and also max number of quantities
  ! given in shopping list when requesting data from external source
  
  INTEGER, PARAMETER, PUBLIC :: max_variables = max_quantities + 5*max_species ! max number of
  !  variables (quantity_species)  in grads file
  
  integer, parameter, public :: max_q_4_derived = 50 ! Max nbr of
  ! quantities needed to create one derived quantity

  INTEGER, PARAMETER, PUBLIC :: max_2d_fields = 3000 ! max number of
  ! 2d fields per time step. Note that 3d field contains nLevs 2d fields,
  ! so max_2d_fields < max_levels * max_quantities - not all
  ! quantities are alowed to be 3-dimensional.

  integer, PARAMETER, PUBLIC :: max_divisions = 130
  
  INTEGER, PARAMETER, PUBLIC :: max_times = 5000 

  INTEGER, PARAMETER, PUBLIC :: worksize = 10000000 ! size of some work
  ! -vectors used when the actual sizes are yet unknown. Should be bigger 
  ! than number of gridpoints in largest 2d nwp data grid!
    
  INTEGER, PARAMETER, PUBLIC :: worksize_2dx = 2500 ! size of some work
  ! -vectors used when the actual sizes are yet unknown
    
  INTEGER, PARAMETER, PUBLIC :: worksize_2dy = 1500 ! size of some work
  ! -vectors in weatherserver used, when the actual sizes are yet
  ! unknown
    
  INTEGER, PARAMETER, PUBLIC :: worksize_string = 2000
  ! size of a work string

  INTEGER, PARAMETER, PUBLIC :: max_work_arrays = 1600 ! max number
  ! of 1d and 2d work arrays in use simultaneously. Their
  ! sizes also given here.
  ! in project_a_src_second_grd max_work_arrays needed is Maximum of
  ! species per source times nThreads
  integer, parameter, public :: max_threads = 1024


  ! ----------------------------------------------------
  !
  ! The following parameters are related to the actual separate dose
  ! assessment:

  ! Maximum number of nuclides (e.g. in the source term and cocktails):

  INTEGER, PARAMETER, PUBLIC :: max_nuclides = 84


  ! Note: It was noticed on March 19, 1997, when doing the INEX-2 preliminary
  ! exercise, that after changing max_matrices the program had to be compiled
  ! twice in a row, both times with 'touch *.f90 & make'. Otherwise, it would
  ! not function correctly. The same might be true for other parameters as
  ! well.

  ! 'max_nuclides_used' is the dimension of some small arrays in module
  ! 'nuclides' (nuc_array is allocated dynamically). It is normally equal to
  ! the number of nuclides in the SILJA nuclide database.

  INTEGER, PARAMETER, PUBLIC :: max_nuclides_used = 496

  ! 'max_nuclide_energies' is the dimension of energy and emission probabilities
  ! as read from the nuclide namelist. Should be changed if/when someone adds
  ! more data to nuclide namelist
  integer, parameter, public :: max_nuclide_energies = 50
  

  ! --------------------------------------------------------------------------
  ! The following parameters are related to the coding system of radioactive
  ! decay chains:

  INTEGER, PARAMETER, PUBLIC :: max_chains = 23

  INTEGER, PARAMETER, PUBLIC :: max_nuclides_in_chain = 10
  INTEGER, PARAMETER, PUBLIC :: max_daughters = 3

  ! ----------------------------------------------------

  ! Same thing for nuclide names (e.g. 'TE-131M')

  INTEGER, PARAMETER, PUBLIC :: nuc_name_len = 7

  INTEGER, PARAMETER, PUBLIC :: max_trajectory_steps = 10000 
      ! Trajectory will be cut after that many steps



END MODULE max_sizes_and_limits







