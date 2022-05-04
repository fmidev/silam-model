MODULE source_terms_bomb

  ! This module contains the general description nuclear bomb source-term of SILAM
  ! Source-term describes the spatial- and time- distribution of the 
  ! release to atmosphere, and the amount of chemical/radioactive 
  ! materials released.
  !
  ! NOTE. Co-ordinates of the bomb source (position) are always stored
  !       in the geographical grid as in ini file. 
  !
  ! Currently the module contains description of BOMB SOURCE
  !
  ! All units: NOT NECESSARILY SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Code owner: Mikhail Sofiev, FMI
  ! Edited 2018: Sebastian Heinonen
  ! 
  USE source_terms_time_params
  use chemical_setup

  IMPLICIT NONE
  private

  ! The public functions and subroutines available in this module:
  
  ! Functions universal for all types of the sources
  !
  public reserve_bomb_source
  public add_input_needs
  public link_source_to_species
  public add_source_species_b_src
  public defined
  PUBLIC report
  PUBLIC fu_name
  public fu_sector
  PUBLIC fu_start_position
  PUBLIC fu_start_time
  PUBLIC fu_end_time
  PUBLIC fu_duration
  public total_bomb_src_species_unit
  public fu_source_id_nbr
  public fu_source_nbr
  public fill_b_src_from_namelist
  public source_2_map_bomb_source
  public source_2_second_grid
  public create_source_containing_grid
  public inject_emission_euler_b_src
  public inject_emission_lagr_b_src

  ! Nuclear bomb source specific functions
  PUBLIC fu_bomb_source_yield

  ! The private functions and subroutines for bomb source
  private uranium_ind_stuk
  private plutonium_stuk
  private add_input_needs_b_src
  private fu_bomb_source_defined
  PRIVATE fu_bomb_source_position
  PRIVATE fu_bomb_source_start_time
  PRIVATE fu_bomb_source_end_time
  PRIVATE fu_bomb_source_duration
  private fu_name_b_src
  private fu_sector_b_src
  private init_mushroom_b_src
  private initialize_vertical_fraction
  private fu_source_id_nbr_of_b_src
  private fu_source_nbr_of_b_src
  private link_b_src_to_species
  private project_b_src_second_grd  ! projects to the given grid and stores new position
  private create_src_cont_grd_b_src
  private report_bomb_source
  private fu_venting_fraction
  private fu_fireball_buried_fraction
  private fu_fission_yield
  private fu_fireball_radius
  private get_horizontal_coverage_and_mom_4_corners
  private set_cocktail_bomb_source
  private set_two_modes_bomb
  private init_bomb_derived_params
  private get_bomb_bins_from_2_lognorm_modes

  ! Generic names and operator-interfaces of some functions:

  interface add_input_needs
    module procedure add_input_needs_b_src
  end interface

  INTERFACE report
    MODULE PROCEDURE report_bomb_source
  END INTERFACE

  interface defined
    module procedure fu_bomb_source_defined
  end interface

  INTERFACE fu_name
    module procedure fu_name_b_src
  END INTERFACE

  INTERFACE fu_sector
    module procedure fu_sector_b_src
  END INTERFACE

  INTERFACE fu_start_position
    MODULE PROCEDURE fu_bomb_source_position
  END INTERFACE

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_bomb_source_start_time
  END INTERFACE

  INTERFACE fu_end_time
    MODULE PROCEDURE fu_bomb_source_end_time
  END INTERFACE

  INTERFACE fu_duration
    MODULE PROCEDURE fu_bomb_source_duration
  END INTERFACE

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_b_src
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_b_src
  end interface

  interface link_source_to_species
    module procedure link_b_src_to_species
  end interface


  interface source_2_second_grid
    module procedure project_b_src_second_grd
  end interface

  interface create_source_containing_grid  
    module procedure create_src_cont_grd_b_src
  end interface

  !
  ! The bomb source
  !
  TYPE silam_bomb_source
    PRIVATE
    CHARACTER(len = clen) :: src_nm, sector_nm, bomb_type, dist_type ! Source name, sector name, bomb type
    integer :: src_nbr, id_nbr  ! source and id number in the WHOLE source list
    integer :: nDescriptors, location_switch, nSpecies
    TYPE(silam_grid_position) :: position ! contains also start time
    type(silam_vertical) :: vert          ! this is stupid: stored and needed only once
    type(silam_species), dimension(:), pointer :: species
    ! activities are in same order as in cocktail.
    ! height_frac in same order as vertical
    real, dimension(:), pointer :: activities, height_frac
    type(chemical_adaptor) :: adaptor2Trn
    type(Taerosol) :: aerosolSrc       ! storage place for the requested bins
    REAL :: total_yield, blast_height_m, fission_fraction, venting_fraction, &
          & fireball_buried_fraction, fireball_radius, fission_yield
    real :: fxDispGrid, fyDispGrid  ! coordinates of the source in the dispersion_grid
    logical :: ifOverWater, if_inside_domain
    type(silja_logical) :: defined
  END TYPE silam_bomb_source

  type b_src_ptr
    TYPE(silam_bomb_source) :: b_src
  end type b_src_ptr
  public b_src_ptr

  ! Types of the blast allowed here
  integer, private, parameter :: blast_underground = 6600
  integer, private, parameter :: blast_on_surface = 6601
  integer, private, parameter :: blast_in_air = 6602
  
  ! Supplementary strucure for handling the release from the blast
  type Tblast_nuclide
    character(len=nuc_name_len) :: nuc_name = ''
    real :: bq_per_kt = real_missing
  end type Tblast_nuclide
  private Tblast_nuclide
  
  
CONTAINS



  !**************************************************************************

  subroutine reserve_bomb_source(b_src, &        ! Src to initialise
                               & iSrcNbr, &      ! Src number in the grand list
                               & iSrcIdNbr)    ! SrcID number
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_bomb_source), intent(inout) :: b_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    b_src%src_nm = ''
    b_src%sector_nm = ''
    nullify(b_src%species)
!    nullify(b_src%nSpeciesInDescr)
!    nullify(b_src%fDescr2SpeciesUnit)
!    nullify(b_src%pEmisSpeciesMapping)
    nullify(b_src%activities)
    nullify(b_src%height_frac)
!    allocate(b_src%cocktail_descr(1))
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    b_src%src_nbr = iSrcNbr
    b_src%id_nbr = iSrcIdNbr
    b_src%nDescriptors = 1
    !
    ! Finally, mark the source as incomplete
    !
    b_src%defined = silja_false

  end subroutine reserve_bomb_source


  !*************************************************************************

  subroutine add_input_needs_b_src(b_src, q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_bomb_source), intent(in) :: b_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static
    !Internal
    integer :: iTmp
    !
    ! Land mask needed
    !
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_static)
    
    if(fu_fails(fu_merge_integer_to_array(height_flag, q_met_dynamic) > 0,'Failed height_flag', &
                                                              & 'add_input_needs_bomb_source'))return
    !
  end subroutine add_input_needs_b_src


  !*****************************************************************

  subroutine link_b_src_to_species(species_list, b_src)
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and cocktail. The link goes via descr%iEmisCocktSpeciesMapping 
    ! and  descr%factor_to_basic_unit.
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    ! That has to happen in two steps. Firstly, we establish these links using the 
    ! single descriptor per source. Then, these connections are distributed to each time slot
    implicit none
    ! Imported parameters
    type(silam_bomb_source), intent(inout) :: b_src
    type(silam_species), dimension(:), pointer :: species_list
    !biovoc source does this similarly
    !call link_src_species_to_given_list(1, &   !b_src%nDescriptors, &
    !                                  & b_src%nSpeciesInDescr, &
    !                                  & b_src%cocktail_descr, &
    !                                  & b_src%fDescr2SpeciesUnit, &
    !                                  & b_src%pEmisSpeciesMapping, &
    !                                  & species_list)
!    call create_adaptor(fu_species(b_src%cocktail_descr(1)), species_list, b_src%adaptor2Trn)
    call create_adaptor(b_src%species, species_list, b_src%adaptor2Trn)
  end subroutine link_b_src_to_species


  !*******************************************************************

  subroutine add_source_species_b_src(b_src, speciesLst, nSpecies)
    !
    ! Fills-in the given list with the own species. Checks for the duplicates
    ! 
    implicit none

    ! Improted parameters
    type(silam_bomb_source), intent(in) :: b_src
    type(silam_species), dimension(:), pointer :: speciesLst
    integer, intent(inout) :: nSpecies

    call addSpecies(speciesLst, nSpecies, b_src%species, b_src%nSpecies)

  end subroutine add_source_species_b_src


  !*****************************************************************

  subroutine total_bomb_src_species_unit(b_src, species, nspecies, amounts, start, duration, layer_)
    !
    ! Returns the amount of the released material starting from 
    ! start during the duration time interval. For bomb source it actually 
    ! means that the interval either covers or does not cover the explosion
    ! moment. So, the amount is either whole mass or zero.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_bomb_source), intent(in) :: b_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nspecies
    real, dimension(:) :: amounts
    type(silja_time), intent(in), optional :: start
    type(silja_interval), intent(in), optional :: duration
    type(silja_level), intent(in), optional :: layer_  ! if defined, include only fraction emitted in it

    ! Local variables
    real :: amount, fOverlap
    type(silja_level) :: layer
    type(chemical_adaptor) :: adaptor
    integer :: iSpecies, nSpeciesSrc
    type(silam_species), dimension(:), pointer :: speciesSrc
    !
    ! The blast is assumed to last one minute. Check the overlap
    !
    nspecies = 0
    
    if(present(start))then
      if(fu_time(b_src%position) + one_minute < start)return
      if(present(duration))then
        if(fu_time(b_src%position) > start + duration)return
        fOverlap = (fu_earliest_time((/fu_time(b_src%position)+one_minute, start+duration/)) - &
                  & fu_latest_time((/fu_time(b_src%position), start/))) / one_minute
        if(fOverlap < 1e-5)return
      else
        fOverlap = 1.0  ! no duration and bomb blasts after the start
      endif   ! present duration
    else 
      fOverlap = 1.0  ! no time boundaries, always whole emission
    endif  ! present start
    !
    ! Get the relative amount in layer (if specified).
    !
    if (present(layer_)) then
      layer = layer_
      amount = b_src%height_frac(nint(fu_level_index(layer,b_src%vert)))
      if(error)return
    else
      amount = 1.0
    end if

    call msg('Emission for bomb source: ' + fu_name(b_src) + ':', amount)

    if(amount == 0.) return
    !
    ! Conversion into the species units is comparatively straightforward.
    !
    nullify(species)
    nSpecies = 0
    call add_source_species_b_src(b_src, species, nSpecies)
    amounts(1:nSpecies) = 0.0
    !
    ! Now explore the descriptors
    !
    if (fOverlap > 0.0) then
       do iSpecies = 1, nSpecies
          !call set_error('Scaling from yield to total emission is not properly introduced', &
          !     & 'inject_emission_bomb_source')
          !call unset_error('inject_emission_bomb_source')
          amounts(iSpecies) = amount * &               !Fraction of activity in level (or 1.0 if level not given)
                            & b_src%activities(iSpecies) * &!Total activity of species 
                            & fOverlap                      !Fraction of overlap of time and bomb (0 or 1)
       end do
    end if
  end subroutine total_bomb_src_species_unit

  
  !*****************************************************************
  
  subroutine initialize_vertical_fraction(b_src, vertical)
    !
    ! Stores fractions of total activity to bomb source for each 
    ! level in the vertical.
    ! The hat and stem are divided into three equally tall pieces
    ! and relative activities of them are distributed similarly as
    ! in KDFOC3. The total of the fractions of these 6 pieces is 1.
    !
    implicit none
    ! Imported parameters
    type(silam_bomb_source), intent(inout) :: b_src
    type(silam_vertical), intent(in) :: vertical
    !Local variables
    real :: stem_bottom, stem_top, hat_top, stem_radius, hat_radius
    real :: base_surge_radius, tot_frac, fraction, forced_act, dummy_sum
    real :: hatFractions(3) = (/0.3, 0.3, 0.18/) !KDFOC3
    real :: stemFractions(3) = (/0.02, 0.05, 0.15/)
    !Stem and Hat thirds are the layers that consist of the activity 
    !fractions in the same order as fractions above
    type(silja_level) :: levBaseSurge, levStemThirds(3), levHatThirds(3), level
    integer :: i, iLev, nlevels
    !Vertical stored
    character(len=*), parameter :: sub_name = 'initialize_vertical_fraction'
    b_src%vert = vertical
    !
    ! Stupidity check
    !
    if (.not. defined(vertical)) call set_error('Undefined vertical given',sub_name)
    if (error) return
    !
    !Check that lower bound of lowest level is 0.0. 
    !Error if not (will have to tweak this if SILAM sometimes has different lowest bound).
    !
    if (fu_leveltype(vertical) /= layer_btw_2_height) then
          call msg("Vertical for bomb initialize_vertical_fraction")
          call report(vertical,.True.) 
          call msg("bottom value", fu_bottom_of_layer_value(fu_level(vertical,1)))
          call set_error("Bomb source works on z vertical only..", sub_name)
    endif
    if (fu_bottom_of_layer_value(fu_level(vertical,1)) /= 0.0) call set_error('Lowest bound not 0.0 meters???',sub_name)
    
    if (error) return
    !
    !Get cloud measures in meters above surface
    !
    call init_mushroom_b_src(b_src, & 
                           & stem_bottom, stem_top, hat_top, &
                           & stem_radius, hat_radius, &
                           & base_surge_radius)
    !
    !Set layers for base surge AND hat and stem divided into three equally tall parts
    !
    levBaseSurge = fu_set_layer_between_two(fu_set_constant_height_level(0.0), &
                                          & fu_set_constant_height_level(stem_bottom))
    do i = 1,3
       levStemThirds(i) = fu_set_layer_between_two(fu_set_constant_height_level((stem_top-stem_bottom)*((real(i)-1.0)/3.0) + stem_bottom), &
                                                 & fu_set_constant_height_level((stem_top-stem_bottom)*(real(i)/3.0) + stem_bottom))
       levHatThirds(i) = fu_set_layer_between_two(fu_set_constant_height_level((hat_top-stem_top)*((real(i)-1.0)/3.0) + stem_top), &
                                                & fu_set_constant_height_level((hat_top-stem_top)*(real(i)/3.0) + stem_top))
       !stem_top == hat_bottom
    end do
    !
    !Loop through levels and find out if all bomb layers are represented by the SILAM layers
    !Simultaneously get height distribution of the levels that are represented.
    !
    nlevels = fu_NbrOfLevels(vertical)
    call enlarge_array(b_src%height_frac,nlevels)
    b_src%height_frac = 0.0
    do iLev = 1,nlevels
      level = fu_level(vertical, iLev)
      fraction = 0.0
      select case(b_src%location_switch)
        case (blast_in_air)                       !Air burst: only hat present
          tot_frac = sum(hatFractions)
          do i = 1,3
             fraction = fraction + ((hatFractions(i)/tot_frac) * fu_vert_overlap_fraction(levHatThirds(i), level))
          end do
          b_src%height_frac(iLev) = b_src%height_frac(iLev) + fraction
        case (blast_underground)                  !Underground burst: hat, stem and base surge
          tot_frac = (1.0 - b_src%venting_fraction) + 1.0
          do i = 1,3
             fraction = fraction + ((stemFractions(i)/tot_frac) * fu_vert_overlap_fraction(levStemThirds(i), level)) + &
                                 & ((hatFractions(i) /tot_frac) * fu_vert_overlap_fraction(levHatThirds(i),  level))
          end do
          fraction = fraction + (((1.0 - b_src%venting_fraction)/tot_frac) * fu_vert_overlap_fraction(levBaseSurge, level))
          b_src%height_frac(iLev) = b_src%height_frac(iLev) + fraction
        case (blast_on_surface)                    !Surface burst: hat and stem
          do i = 1,3
             fraction = fraction + (stemFractions(i) * fu_vert_overlap_fraction(levStemThirds(i), level)) + &
                                 & (hatFractions(i)  * fu_vert_overlap_fraction(levHatThirds(i),  level))
          end do
          b_src%height_frac(iLev) = b_src%height_frac(iLev) + fraction
        case default
          call set_error('Unknown blast location switch:' + fu_str(b_src%location_switch), sub_name)
      end select
      if(fu_fails(b_src%height_frac(iLev) >= 0.0, 'About to inject negative mass to a level',sub_name))return
    end do  ! iLev
    
    !
    ! Stuff beyond the projection
    !

    forced_act = sum(b_src%height_frac) !! Total of fractions numerics can lead to negatives
    if (abs(1.0 - forced_act) < 2e-7) then 
      forced_act = 0
    else
      forced_act = 1.0 - forced_act
    endif
    !
    ! Error if activity forced too much, warning if quite a lot
    !
    if (forced_act > 0.01) then
       call msg_warning('Some of the activity is above SILAM top and is forced to the highest SILAM level',sub_name)
       call msg('The fraction of stuff forced from bwyond the domain',forced_act)
       if (forced_act > 0.5) then
         call msg("Mushroom heights: stem_bottom, stem_top, hat_top", (/stem_bottom, stem_top, hat_top/))
         call msg("Dispersion vertical:")
         call report(vertical, .True.)
         call set_error('Too much stuff emitted beyond the dispersion vertical', sub_name)

       endif
    end if
    !
    ! Force the remaining fraction of the bomb the the highest SILAM level
    !
    b_src%height_frac(nlevels) = b_src%height_frac(nlevels) + forced_act

    !
    if (.not. abs(sum(b_src%height_frac) - 1.0) < 0.001)  then
        call msg("bomb fraction per levels", b_src%height_frac(i))
        call set_error('Fraction total not 1',sub_name)

        return
    endif

    if (any(b_src%height_frac(:) < 0.)) then
      call set_error("Negative Fraction", sub_name)
    endif
       
  end subroutine initialize_vertical_fraction

  !*****************************************************************

  integer function fu_source_id_nbr_of_b_src(b_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_bomb_source), intent(in) :: b_src

    ! Stupidity check
    if(.not. (b_src%defined == silja_true))then
      call set_error('Undefined source given','fu_source_nbr_of_bomb_source')
      return
    endif
    fu_source_id_nbr_of_b_src = b_src%id_nbr

  end function fu_source_id_nbr_of_b_src


  !*****************************************************************

  integer function fu_source_nbr_of_b_src(b_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_bomb_source), intent(in) :: b_src

    ! Stupidity check
    if(.not. (b_src%defined == silja_false))then
      fu_source_nbr_of_b_src = b_src%src_nbr
    else
      fu_source_nbr_of_b_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_bomb_source')
      return
    endif

  end function fu_source_nbr_of_b_src


  !=================================================================

  function  fu_name_b_src(b_src) result(nm)
    implicit none
    character(len=clen) :: nm
    type(silam_bomb_source), intent(in) :: b_src
    nm = b_src%src_nm
  end function fu_name_b_src

  !=================================================================

  function  fu_sector_b_src(b_src) result(sect)
    implicit none
    character(len=clen) :: sect
    type(silam_bomb_source), intent(in) :: b_src
    sect = b_src%sector_nm
  end function fu_sector_b_src


  ! ***************************************************************

  subroutine fill_b_src_from_namelist(nlSrc, b_src, chBombSrcFileVersion)
    !
    ! Reads and sets one bomb source term from an external file.
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_bomb_source), intent(inout) :: b_src
    character(len=*), intent(in) :: chBombSrcFileVersion
    type(Tsilam_namelist), pointer :: nlSrc

    ! Local declarations:
    REAL :: fLat, fLon
    integer :: iLocal, iTmp
    character(len=10) :: chUnit
    character(len=*), parameter :: sub_name = 'fill_b_src_from_namelist'

    ! A bit of preparations
    !
    b_src%defined = silja_false

    if (chBombSrcFileVersion /= '4') then
      call set_error("only BOMB_SOURCE_4 is supported", sub_name)
      return
    endif


    
    !
    !  Name of the source should come
    !
    b_src%src_nm = fu_content(nlSrc,'source_name')
    !
    ! Lat & lon
    !
    fLat = fu_content_real(nlSrc,'source_latitude')
    if (fLat == real_missing) call set_error('Failed to read source',sub_name)
    fLon = fu_content_real(nlSrc,'source_longitude')
    if (fLon == real_missing) call set_error('Failed to read source',sub_name)
    IF (error) RETURN

    !
    ! There is only one time parameter - actual time of explosion. Duration is one minute
    !
    b_src%position = fu_set_pos_in_geo_global_grid(fLon, fLat, 92500., &
                                         & fu_io_string_to_time(fu_content(nlSrc,'explosion_time')))
    if(error)return
    !
    ! The bomb yield
    !
    b_src%total_yield = fu_set_named_value(fu_content(nlSrc,'bomb_yield'))
    if(fu_fails(.not.((b_src%total_yield .eps. real_missing) .or. b_src%total_yield < 0.), &
                                    & 'Failed to read yield', sub_name))return
    b_src%total_yield = b_src%total_yield * 1.0e-6   ! set_named_value puts all to SI, kg in this case
    !
    !  Height of the blast, negative means underground and surface blast if not given
    !
    b_src%blast_height_m = fu_set_named_value(fu_content(nlSrc,'blast_height'))
    if(fu_fails(.not.(b_src%blast_height_m == real_missing),'Failed to read blast_height', &
                                      & sub_name))return
    !
    !  Bomb type, options: URANIUM, PLUTONIUM
    !
    b_src%bomb_type = fu_content(nlSrc,'bomb_type')
    if (fu_fails(b_src%bomb_type == 'URANIUM' &
             & .or. b_src%bomb_type == 'PLUTONIUM' &
             & .or. b_src%bomb_type == 'PASSIVE', &
          & 'Missing or wrong bomb_type (can be URANIUM or PLUTONIUM or PASSIVE)',sub_name))return
    !
    !  Lognormal distribution based on, options: KDFOC3, BAKER
    !  KDFOC3: A Nuclear Fallout Assessment Capability Harvey et al. 1992, UCRL-TM-222788
    !  BAKER: Implications of Atmospheric Test Fallout Data for Nuclear Winter, Baker 1987
    b_src%dist_type = fu_content(nlSrc,'particle_size_distribution')
    if(fu_fails(b_src%dist_type == 'KDFOC3' .or. b_src%dist_type == 'BAKER', &
                        & 'Missing or wrong particle_size_distribution (can be KDFOC3 or BAKER)',sub_name))return
    !
    ! Fission fraction
    ! -The worst case scenario is that fission fraction is 1
    !  i.e. all of the bombs energy comes from fission and not fusion.
    !  Largest pure fission bomb detonated ever was 500kt.
    !  A rule of thumb says that 500kt or larger will have at least 0.5 fusion fraction.
    !  We will use a linear dependence around this value so that 300kt is still pure fission
    !  and 700kt and larger must be at least half fusion.
    ! -Fission fraction and yield determine the amount of activity and hopefully
    !  this doesn't underestimate it.
    ! -Would be unrealistic to assume that user would have any idea what the fission
    !  fraction is in a situation; thus it is not a parameter that a user could give.
    if (b_src%total_yield > 700.0) then
       b_src%fission_fraction = 0.5
    else if (b_src%total_yield < 300.0) then
       b_src%fission_fraction = 1.0
    else
       b_src%fission_fraction = 0.5 + 0.5*((b_src%total_yield-300.0)/(700.0-300.0))
    end if
    !
    ! Aerosol modes are in the namelist, get them
    !
    call set_aerosol(nlSrc, b_src%aerosolSrc)
    if(error)return
    !
    ! Call function to fill in b_src derivec attributes
    !
    call init_bomb_derived_params(b_src, .false.)
    if(error)return
    !
    !NOTE: derived values already initialized, but they are misleading if
    ! surface isnt normal land. They are initialized more accurately later
    ! and these are only used for reporting.
    !
    b_src%defined = silja_true
    !
    ! Set the cocktail. It depends on the blast height 
    ! Needs defined source
    !

    call set_cocktail_bomb_source(b_src)
    if(error)return
    !
    ! Things have been done but do NOT destroy the namelist - it was made above
    !
    call msg('Bomb source name:' + b_src%src_nm)

  END subroutine fill_b_src_from_namelist


  !******************************************************************

  subroutine init_mushroom_b_src(b_src, & 
                               & stem_bottom, stem_top, hat_top, &
                               & stem_radius, hat_radius, &
                               & base_surge_radius)
    !
    ! Converts the source yield in kilotones to parameters of the mushroom
    !
    implicit none

    ! Imported parameters
    type(silam_bomb_source) :: b_src
    real, intent(out) :: stem_bottom, stem_top, hat_top
    real, intent(out) :: base_surge_radius, stem_radius, hat_radius
    
    ! Calculate stem and cloud dimensions:
    select case(b_src%location_switch)
      case (blast_in_air)                       !Air burst: only hat present
        !Only hat created
        stem_bottom = 0.0
        !stem top defines the hat bottom
        IF (b_src%total_yield <= 20) THEN
          stem_top = 1670.0*b_src%total_yield**0.477
        ELSE
          stem_top = 4450.1*b_src%total_yield**0.159
        END IF
        ! Currently bomb detonation height doesn't affect the height of the hat for air bursts.
        ! This model is incapable of simulating a bomb cloud that is detonated at high altitude.
        ! So if this parametrization leads to lower hat bottom as detonation height error is invoked.
        if (stem_top < b_src%blast_height_m) call set_error('Too high detonation height set','init_mushroom_b_src')
        !
        IF (b_src%total_yield <= 20) THEN
          hat_top = 3365.0*b_src%total_yield**0.38
        ELSE
          hat_top = 6370.3*b_src%total_yield**0.177 !realistic-ish even for tsar bomba
        END IF
        hat_radius = 970.0*b_src%total_yield**0.42
        stem_radius = 0.0
        base_surge_radius = 0.0
        
      case(blast_on_surface)                        !Surface burst
        !Hat and stem are created
        stem_bottom = 0.0
        IF (b_src%total_yield <= 20) THEN
          stem_top = 1670.0*b_src%total_yield**0.477
        ELSE
          stem_top = 4450.1*b_src%total_yield**0.159
        END IF
       
        IF (b_src%total_yield <= 20) THEN
          hat_top = 3365.0*b_src%total_yield**0.38
        ELSE
          hat_top = 6370.3*b_src%total_yield**0.177
        END IF
        hat_radius = 970.0*b_src%total_yield**0.42
        ! Radius of stem is calculated by linear interpolation from known stem radii
        ! for 20 kt and 1 Mt. For small explosions, this may lead to greater radius
        ! of stem than that of cloud, which is avoided by forcing stem radius to that
        ! of cloud:
        ! EDIT: changed so that stem can be half the radius of hat at maximum
        stem_radius = max(b_src%fireball_radius*3.0, 1705.0 + (2647.5 - 1705.0) * (b_src%total_yield - 20.0)/(1000.0 - 20.0))
        base_surge_radius = 0.0
        
      case(blast_underground)                        !bomb detonated underwater or underground
        !Argueably the amount of escaped material is larger when bomb is detonated underground than underwater
        ! but underwater explosions vaporize water that condensates as droplets and fall relatively faster.
        !The unknowns in both cases disable our ability to determine even if either of these detonations 
        ! lead to higher effective yield as other.
        !Nevertheless the detonation depth is already taken into account in the venting fraction calculation.
        !Therefore here underwater and underground detonations are treated the same.
        IF (b_src%total_yield <= 20) THEN
          stem_top = 1670.0*b_src%total_yield**0.477
        ELSE
          stem_top = 4450.1*b_src%total_yield**0.159
        END IF
       
        IF (b_src%total_yield <= 20) THEN
          hat_top = 3365.0*b_src%total_yield**0.38
        ELSE
          hat_top = 6370.3*b_src%total_yield**0.177
        END IF
        hat_radius = 970.0*b_src%total_yield**0.42
        stem_radius = max(b_src%fireball_radius*3.0, 1705.0 + (2647.5 - 1705.0) * (b_src%total_yield - 20.0)/(1000.0 - 20.0))
        !Base surge height is thought to reach stem top height in the limit of 0 venting fraction
        ! -limit of veting fraction = 1 base surge height is 0
       
        stem_bottom = (1.0 - b_src%venting_fraction)*stem_top
        !Base surge radius is thought to reach 2*hat radius in limit of venting fraction = 0
        ! -limit of veting fraction = 1 base surge radius is stem radius
        base_surge_radius = (1.0 - b_src%venting_fraction)*2.0*stem_radius + stem_radius
      case default
       !None of the logicals, just an extra error setter
       call set_error('Impossible bomb','init_mushroom_b_src')
    end select

  end subroutine init_mushroom_b_src


  ! ***************************************************************

  SUBROUTINE report_bomb_source(bs)

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_bomb_source), INTENT(in) :: bs

    call msg('')
    call msg('=========================== BOMB SOURCE ===============================')
    call msg('Nuclear bomb source:' + bs%src_nm)
    call msg('start point: ' +  fu_geographical_latlon_string(bs%position,geo_global_grid))
    call msg('time of release: ' + fu_str(fu_time(bs%position)))
    call msg('bomb yield:',bs%total_yield)
    call msg('=======================================================================')

  END SUBROUTINE report_bomb_source


  !****************************************************************

  subroutine source_2_map_bomb_source(bs, dataPtr, id)
    !
    ! Projects the point source to the map, having the given id as a 
    ! template: the id deterines the grid, level and substance name
    !
    implicit none

    ! Imported parameters
    type(silam_bomb_source), intent(in) :: bs
    real, dimension(:), intent(inout) :: dataPtr
    type(silja_field_id), intent(in) :: id

    ! Local variables
    real :: xOut, yOut !, fRate, fTmp
    type(silam_vertical) :: vertSrc
!    type(Tcocktail_descr) :: cocktail_descr
    integer :: iOut, nx, ny, nSpecies !, iDescr
    real, dimension(:), pointer :: fractions, amounts
    type(silam_species), dimension(:), pointer :: species
!    integer, dimension(:), pointer :: species_index
!    logical :: ifFound
    !
    ! Stupidity check
    !
    if(.not.(bs%defined == silja_true))then
      call set_error('Undefined bomb source given','source_2_map_bomb_source')
      return
    endif
    if(.not.defined(id))then
      call set_error('Undefined id given','source_2_map_bomb_source')
      return
    endif
    if(.not.defined(fu_grid(id)))then
      call set_error('Undefined grid given','source_2_map_bomb_source')
      return
    endif
    if(.not.fu_number_of_gridpoints(fu_grid(id)) > size(dataPtr))then
      call set_error('Too small data array given','source_2_map_bomb_source')
      return
    endif
    !
    ! Find the index and mass fraction of the requested name. 
    !
    fractions => fu_work_array()
    
    !
    ! Find the grid location of the source. A trick: the point source is always
    ! in geographical co-ordinates, while the output grid can be whatever
    !
    call project_point_to_grid(geo_global_grid, fu_x(bs%position), fu_y(bs%position), &
                             & fu_grid(id), xOut, yOut)
    call grid_dimensions(fu_grid(id), nx, ny)
    if(xOut < 0.5 .or. xOut > nx+0.5)then
      call msg('Bomb source is outside the x-limits.x=', xOut)
      call msg_warning('Bomb source is outside the x-limits.Skipping','source_2_map_bomb_source')
      return
    endif
    if(yOut < 0.5 .or. yOut > ny+0.5)then
      call msg('Bomb source is outside the y-limits;y=', yOut)
      call msg_warning('Bomb is outside the y-limits.Skipping','source_2_map_bomb_source')
      return
    endif
    iOut = int(xOut+0.5) + (int(yOut+0.5)-1)*nx
    !
    !
    ! Linear variation goes between time slots. Let's read the start time
    !
    amounts => fu_work_array()
    !
    !
    nullify(species)
    nSpecies = 0
    call total_bomb_src_species_unit(bs, &
                                   & species, nSpecies, amounts, &
                                   & fu_accumulation_start_time(id), &
                                   & fu_accumulation_length(id), &
                                   & fu_level(id))
    if(error)return
    if(amounts(1) < 1e-10)then  ! quite arbitrary but small anyway
      call free_work_array(amounts)
      return
    endif
!    nullify(species)
!    nSpecies = 0
!    !
!    ! Once the full inventory is returned, have to select the specific species for each descriptor.
!    ! In fact, this is almost the worst possible solution but we cannot afford keeping all the
!    ! maps for all the sources
!    !
!    species_index => fu_work_int_array()
!    if(error)return
!    ifFound = .false.
!    do iDescr = 1, bs%nDescriptors
!      call get_inventory(bs%cocktail_descr(1), species, nSpecies)
!      species_index(iDescr) = fu_index(fu_species(id), species, nSpecies)
!      if(species_index(iDescr) >= 1 .and. species_index(iDescr) <= nSpecies) ifFound = .true.
!    end do
!    species_index = fu_index(fu_species(id), species, nSpecies)
!    if(.not. ifFound)then
!      call free_work_array(amounts)
!      call free_work_array(species_index)
!      return
!    endif
!    !
!    ! Preparatory work is over
!    !
!    do iDescr = 1, bs%nDescriptors
!      if(species_index(iDescr) < 1)cycle    ! Some descriptors may not have the needed species
!      dataPtr(iOut) = dataPtr(iOut) + amounts(species_index(iDescr)) * &
!                                    & bs%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
!    end do

    nx = fu_index(fu_species(id), species, nSpecies)
    if (nx /= int_missing) dataPtr(iOut) = dataPtr(iOut) + amounts(nx)
    
    call free_work_array(amounts)
!    call free_work_array(species_index)

!    call set_error('Scaling of the emission to yield is not introduced','source_2_map_bomb_source')
!    call unset_error('source_2_map_bomb_source')

  end subroutine source_2_map_bomb_source


 !*****************************************************************

  subroutine project_b_src_second_grd(bs, grid)
    !
    ! Just finds the grid-position in the given grid and stores it to the
    ! position_disp_grd
    !
    implicit none

    ! Imported parameters
    type(silam_bomb_source), intent(inout) :: bs
    type(silja_grid), intent(in) :: grid

    ! Local variables
    real :: xOut, yOut
    integer :: nx, ny

    !
    ! Just create new position from the original one with modified horizontal coordinates
    !
    call project_point_to_grid(geo_global_grid, fu_x(bs%position), fu_y(bs%position), &
                             & grid, bs%fxDispGrid, bs%fyDispGrid)
    call grid_dimensions(grid,nx,ny)
    !
    ! Check that it is inside the grid
    !
    if(bs%fxDispGrid < 0.5 .or. bs%fxDispGrid > nx + 0.5 .or. &
     & bs%fyDispGrid < 0.5 .or. bs%fyDispGrid > ny + 0.5) then
      bs%if_inside_domain = .false.
    else
      bs%if_inside_domain = .true.
    endif

  end subroutine project_b_src_second_grd


  !****************************************************************

  subroutine create_src_cont_grd_b_src(b_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! This routine creates the grid similar to the given grid_template
    ! covering all sources of the em_source
    !
    implicit none

    ! Imported parameters
    type(silam_bomb_source), intent(in) :: b_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! Local variables
    integer :: iSrc, nx, ny, iCell
    real :: x, y

    if(.not. (b_src%defined == silja_true))then
      call set_error('Undefiend source','create_src_cont_grd_b_src')
      return
    endif

    if(defined(grid_template))then
      if(ifMinimal)then
        !
        ! The requested grid is to be minimal-size just covering the current source
        ! We ignore the current size of the grid and its resolution but still keep the 
        ! projection and grid type
        !
        call project_point_to_grid(fu_x(b_src%position), fu_y(b_src%position), grid_template, x, y)
        if(error)return
        
        call make_minimal_grid(grid_template, x, y)
        
      else
        !
        ! For each source we scan its central points and extending the template_grid 
        ! if these points are outside. For point and bomb sources it is only one point,
        ! while area source has plenty of them.
        !
        call grid_dimensions(grid_template, nx,ny)

        call project_point_to_grid(geo_global_grid, &
                                 & fu_x(b_src%position), &
                                 & fu_y(b_src%position), &
                                 & grid_template, x, y)
        if(x<1 .or. x>nx .or. y<1 .or. y>ny) then
          ifExtended = .true.
          if(ifVerbose)then
            call msg('Extending the grid for the bomb source:',iSrc)
            call report(b_src)
          endif
          call extend_grid_to_coordinates(grid_template, x, y)
        else
          ifExtended = .false.
        endif
      endif ! ifMinimal
    else
      !
      ! Grid_template is undefined. Create it using this source as the starting point
      !
      grid_template = fu_set_grid('', lonlat, pole_geographical, &
                                & fu_x(b_src%position)-0.01, &
                                & fu_y(b_src%position)-0.01, &
                                & 3, 3, 0.01, 0.01)   ! nx, ny, dx, dy
      if(error)return
    endif  ! if grid_template is defined

  end subroutine create_src_cont_grd_b_src


  !****************************************************************

  subroutine inject_emission_euler_b_src(b_src, &
                                       & mapEmis, mapCoordX, mapCoordY, mapCoordZ, &
                                       & met_buf, disp_buf, &
                                       & now, timestep, &
                                       & ifSpeciesMoment, &
                                       & fMassInjected, &
                                       & pHorizInterpMet2DispStruct, &
                                       & ifHorizInterp)
    !
    ! Adds the emission flux to the concentration map and stores the position
    ! of the injected mass to map*Coord. All maps may contain some data before,
    ! so no overwriting - just adding.
    ! A simplification: the explosion is supposed to last less than a model time step,
    ! therefore no time integration, just checking the right moment
    !
    implicit none

    ! Imported parameters
    TYPE(silam_bomb_source), intent(inout) :: b_src
    type(Tmass_map), intent(inout):: mapEmis, mapCoordX, mapCoordY, mapCoordZ
    TYPE(Tfield_buffer), intent(in) :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real(8), dimension(:), intent(inout) :: fMassInjected
    logical, intent(in) :: ifSpeciesMoment, ifHorizInterp
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct

    ! Local declarations
    real :: stem_bottom, stem_top, hat_top, land, momin, momout, cnc
    real :: stem_radius, hat_radius, base_surge_radius
    real :: base_surge_stem_r, stem_hat_r
    real :: fDist, xShift, yShift, fTmp, fCellArea, fDX, fDY
    real :: fStemBottomLev, fStemTopLev, fHatTopLev, activityInLevel, trueFissionYield
    real, dimension(:,:,:), pointer :: cell_fraction_in_cloud, momX, momY
    integer :: ix, iy, iTmp, iDistCells, ixSrc, iySrc, iLev, iSpecies, nSpecies, iEmisSpecies, iC
    real, dimension(:), pointer :: amounts
    type(chemical_adaptor) :: adaptor
    type(silam_species), dimension(:), pointer :: species
    real, dimension(5) :: allRadii
    integer :: i, s
    
    if(.not. b_src%if_inside_domain)return  ! either not this MPI subdomain or out of whole domain
    
    amounts => fu_work_array()
    if(error)return
    call total_bomb_src_species_unit(b_src, species, nspecies, amounts, now, timestep)
    !
    if(sum(amounts(1:nSpecies)) < 1e-10)then
      call free_work_array(amounts)
      return
    endif
    call create_adaptor(species, mapEmis%species, adaptor)
    if(error)return
    !
    ixSrc = nint(b_src%fxDispGrid)
    iySrc = nint(b_src%fyDispGrid)
    !
    ! Bomb over land or water? Use physiography-maintained land fraction field
    !
    land = fu_get_value(fraction_of_land_fld, nx_meteo, ixSrc, iySrc, &
                      & pHorizInterpMet2DispStruct, .not. (meteo_grid == dispersion_grid))
    
    if(fu_fails(land<=1.00001 .and. land >= 0., 'strange land fraction:' + fu_str(land),'inject_emission_euler_b_src'))return

    land = max(land, 1.0)

    if (land < 0.1) then
      ! 
      ! Water
      !
      call msg('Bomb over water, land fraction:',land)
      !
      ! Init params again using density weighted detonation height if underwater
      !
      call init_bomb_derived_params(b_src, .true.)
      !
      amounts => fu_work_array()
      if(error)return
      call total_bomb_src_species_unit(b_src, species, nspecies, amounts, now, timestep)
      !
      ! Remove surface and buried activation effect
      !
      trueFissionYield = b_src%total_yield * (b_src%fission_fraction + &
                                            & (1.0 - b_src%fission_fraction) * 0.02)
      amounts(:) = amounts(:)*trueFissionYield/b_src%fission_yield
      !
      ! Remove induced activity if there is any
      !
      do i = 1,nspecies
         if (if_activation_isotope(fu_name(species(i)%material))) then
            amounts(i) = 0.0
         end if
      end do
    else 
      !
      ! Land
      !
      call msg('Bomb over land, land fraction:',land)
      !
      call init_bomb_derived_params(b_src, .false.)
      amounts => fu_work_array()
      if(error)return
      call total_bomb_src_species_unit(b_src, species, nspecies, amounts, now, timestep)
    end if  ! land or water ?
    !
    ! Get the parameters of the cloud and density of emission for stem and the hat
    !
    call init_mushroom_b_src(b_src, stem_bottom, stem_top, hat_top, stem_radius, hat_radius, &
                           & base_surge_radius)
    call initialize_vertical_fraction(b_src, mapEmis%vertTemplate)
    if(error)return
    !
    ! The cloud can be tens of km in diameter, so have to carefully count the number of 
    ! dispersion cells affected by the emission. Since the resolution will anyway be crude,
    ! do the following: for each grid cell check its corners and attribute 0-25-50-75-100%
    ! of cell area to be covered by the cloud depending on 0-1-2-3-4 corners to be
    ! within the radius distance from the centre of the explosion.
    !
    xShift = b_src%fxDispGrid - ixSrc  ! shift from centre of the cell
    yShift = b_src%fyDispGrid - iySrc
    fCellArea = fu_cell_size(dispersion_grid,ixSrc,iySrc)
    fDX = fu_dx_cell_m(dispersion_grid,ixSrc,iySrc)
    fDY = fu_dy_cell_m(dispersion_grid,ixSrc,iySrc)
    iDistCells = nint(max(base_surge_radius, hat_radius) / sqrt(fCellArea)) + 1
    
    !Last dimension is for different radii
    allocate(cell_fraction_in_cloud(-iDistCells:iDistCells,-iDistCells:iDistCells,5), &
                             & momX(-iDistCells:iDistCells,-iDistCells:iDistCells,5), &
                             & momY(-iDistCells:iDistCells,-iDistCells:iDistCells,5), &
                             & stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed to allocate temporary for grid cell coverage', &
                         & 'inject_emission_b_src'))return
    
    !For levels between two radii average radius is implemented
    base_surge_stem_r = (base_surge_radius + stem_radius)/2.0
    stem_hat_r = (stem_radius + hat_radius)/2.0
    allRadii = (/base_surge_radius, &
               & base_surge_stem_r, &
               & stem_radius, &
               & stem_hat_r, &
               & hat_radius/)
    do iC = 1,5
       !For each radii get the coverage and moments, each radius is referred now with iC
       call get_horizontal_coverage_and_mom_4_corners(cell_fraction_in_cloud(:,:,iC), momX(:,:,iC), momY(:,:,iC), &
                                                    & allRadii(iC),iDistcells, fDX, fDY, xShift, yShift)
    end do
    ! This mass has to be distributed along the vertical. All what we need is the top of the stem
    ! and the top of the cloud. Both of them are given in meters and have to be projected to 
    ! the dispersion vertical
    !
    fStemBottomLev = fu_level_index( &
                          & fu_level_to_vertical_crude( &
                                 & fu_set_constant_height_level(stem_bottom), mapEmis%vertTemplate), &
                          & mapEmis%vertTemplate)
    fStemTopLev = fu_level_index( &
                          & fu_level_to_vertical_crude( &
                                 & fu_set_constant_height_level(stem_top), mapEmis%vertTemplate), &
                          & mapEmis%vertTemplate)
    fHatTopLev = fu_level_index( &
                          & fu_level_to_vertical_crude( &
                                 & fu_set_constant_height_level(hat_top), mapEmis%vertTemplate), &
                          & mapEmis%vertTemplate)
    !
    do iLev = 1, fu_NbrOfLevels(mapEmis%vertTemplate)
       activityInLevel = b_src%height_frac(iLev) !fraction
       !
       ! What part of the cloud is this level
       !
       if (activityInLevel == 0.0) cycle !Not in cloud
       if (activityInLevel < 0.0) call set_error('About to inject negative mass','inject_emission_euler_b_src')
       if (error) return
       if ((base_surge_radius > 0.0) .and. (iLev+0.5 <= fStemBottomLev)) then !Base surge cloud
          iC = 1
       else if ((iLev-0.5 >= fStemBottomLev) .and. (iLev+0.5 <= fStemTopLev)) then !Stem
          iC = 3
       else if (iLev-0.5 >= fStemTopLev) then !Hat
          iC = 5
       else !Level is belong in two different sections of the cloud
          if ((iLev+0.5 > fStemBottomLev) .and. (iLev-0.5 < fStemBottomLev)) then !Combination of base surge and stem
             iC = 2
          else if ((iLev+0.5 > fStemTopLev) .and. (iLev-0.5 < fStemTopLev)) then !Combination of stem and hat
             iC = 4
          end if
       end if
       if (iLev == fu_NbrOfLevels(mapEmis%vertTemplate)) iC = 5 !Force highest level as hat
       if (sum(cell_fraction_in_cloud(:,:,iC)) == 0.0) call set_error('If there is activity then it should be in some grid cells', &
                                                                    & 'inject_emission_euler_b_src')
       if (error) return
       !
       ! Inject the specific substances directly into mapEmis using pEmisCocktailMapping, which is
       ! the coding rule for species in the map. Total released stuff goes into the coordinate
       ! map as a scaling factor for momentum.
       !
       Y : do iy = -iDistCells, iDistCells
          do ix = -iDistCells, iDistCells
             !
             ! Cycle with zero activity and when trying to reach out of bounds
             !
             if ((cell_fraction_in_cloud(ix,iy,iC) == 0.0) .or. (activityInLevel == 0.0)) cycle
             if (lbound(mapEmis%arM,dim=4) > ix+ixSrc) cycle
             if (lbound(mapEmis%arM,dim=5) > iy+iySrc) cycle Y
             if (ubound(mapEmis%arM,dim=4) < ix+ixSrc) cycle
             if (ubound(mapEmis%arM,dim=5) < iy+iySrc) cycle Y
             !
             !fCellTotal = 0.
             do iSpecies = 1, nSpecies
                if (amounts(iSpecies) <= 0.0) cycle
                fTmp = amounts(iSpecies) * &           !total(species)
                     & activityInLevel * &             !relative fraction in level
                     & cell_fraction_in_cloud(ix,iy,iC)!relative fraction in cell
                iEmisSpecies = adaptor%iSp(iSpecies)
                mapEmis%arM(iEmisSpecies, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) = &
                     & mapEmis%arM(iEmisSpecies, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) + &
                     & fTmp
                !fCellTotal = fCellTotal + fTmp
                fMassInjected(iEmisSpecies) = fMassInjected(iEmisSpecies) + fTmp
                !Below message commented, not quite smart to print emission on
                !each grid cell in each level for each species.
!call msg('Plus this:',mapEmis%arM(iEmisSpecies, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc))
                if (ifSpeciesMoment) then
                   mapCoordX%arm(iEmisSpecies, b_src%id_nbr, ilev, ix+ixsrc, iy+iysrc) = &
                        & mapCoordX%arm(iEmisSpecies, b_src%id_nbr, ilev, ix+ixsrc, iy+iysrc) + &
                        & fTmp * momX(ix,iy,iC)
                   mapCoordY%arm(iEmisSpecies, b_src%id_nbr, ilev, ix+ixsrc, iy+iysrc) = &
                        & mapCoordY%arm(iEmisSpecies, b_src%id_nbr, ilev, ix+ixsrc, iy+iysrc) + &
                        & fTmp * momY(ix,iy,iC)
                   mapCoordZ%arM(iEmisSpecies, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) = &
                        & mapCoordZ%arM(iEmisSpecies, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) + &
                        & fTmp * 0.0
                end if
             end do   ! iSpecies
             
             !
             ! Note that we store coordinates inside each grid cell, relative for horizontal grid
             ! and absolute for vertical one
             !
             fTmp = sum(amounts(1:nSpecies)) * activityInLevel * cell_fraction_in_cloud(ix,iy,iC)
             if (.not. ifSpeciesMoment) then
                mapCoordX%arM(1, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) = &
                     & fTmp * momX(ix,iy,iC)
                mapCoordY%arM(1, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) = &
                     & fTmp * momY(ix,iy,iC)
                !
                ! Here we simply assume that all emitted is set into the centre of the layer
                !
                mapCoordZ%arM(1, b_src%id_nbr, iLev, ix+ixSrc, iy+iySrc) = 0.0
             end if
             mapEmis%ifColumnValid(b_src%id_nbr, ix+ixSrc, iy+iySrc) = .true.
             mapEmis%ifGridValid(iLev, b_src%id_nbr) = .true.
          end do   ! ix
       end do Y    ! iy
    end do  ! iLev
!call report(mapEmis)
    call msg('=============================================================')
    call msg('| Nuclear cloud injected. Dimensions:  Radius       Bottom       Top')
    if (base_surge_radius > 0.0) call msg(' | Base surge cloud                  :', (/base_surge_radius, 0.0, stem_bottom/))
    call msg(' | Stem                              :', (/stem_radius, stem_bottom, stem_top/))
    call msg(' | Hat                               :', (/hat_radius, stem_top, hat_top/))
    call msg('=============================================================')

    deallocate(cell_fraction_in_cloud, stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed deallocating temporary grid cell coverage', &
              & 'inject_emission_bomb_source')) call unset_error('inject_emission_b_src')
    call free_work_array(amounts)
    
  end subroutine inject_emission_euler_b_src

  
  !**************************************************************************

  subroutine inject_emission_lagr_b_src(b_src, &
                                      & lpSet, arParticleMass, & ! Lagrangian
                                      & ChemRunSetup, &  ! translate emission species to transport
                                      & met_buf, disp_buf, &
                                      & now, timestep, &
                                      & fMassInjected)
    !
    ! Adds the emission flux to Lagrangian structure by starting new particles.
    ! Note that particles fly in the meteorological grid to utilise max of available 
    ! dynamic information and also to have a simpler connection to pressure and omega-wind.
    !
    ! Therefore, all source variables with "dispersion-grid"-meaning here mean "meteo grid"
    ! Refer to the project_b_src_2_grids, where the selection is done.
    !
    implicit none

    TYPE(silam_bomb_source), INTENT(in) :: b_src
    type(Tlagrange_particles_set), INTENT(inout), target :: lpSet
    real, dimension(:) :: arParticleMass
    type(TchemicalRunSetup), pointer :: ChemRunSetup
    TYPE(Tfield_buffer), pointer :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real(8), dimension(:), intent(inout) :: fMassInjected

    ! Local variables
    real :: stem_bottom, stem_top, hat_top, stem_radius, hat_radius, fOverlap, base_surge_radius
    integer :: iLev, ix, iSrc, nSpEmis, iSpDescr, nP, iParticle, iP, ind_height, ind_pressure
    integer :: iSpEmis, iSpTo, iSpTransp, nSpeciesSrc, iSpSrc
    real :: fTimeScale, fWeightPastSrc, factor, timestep_sec, fDx, fDy, stem_mass_fract
    real :: fStemBottomPressure, fStemTopPressure, fHatTopPressure
    real, dimension(:), pointer ::  xSize, ySize
    real, dimension(max_species) :: amounts
    real, dimension(max_levels) :: meteo_heights, meteo_pressure
    integer, dimension(:), pointer :: nPartInLev
    logical :: ifFound
    type(silam_species), dimension(:), pointer :: speciesSrc
    type(TspeciesReference), dimension(:), pointer :: references
    integer :: iSlotStart, iSlotEnd
    type(TVertInterpStruct), pointer ::  interpCoefVert_void
    type(THorizInterpStruct), pointer ::  interpCoefHoriz_void
    type(field_4d_data_ptr), pointer :: fldPressure

    real, dimension(:,:), pointer :: arDyn, arMass
    integer, dimension(:), pointer :: arStatus
    !
    ! Is there anything happening?
    !
    if(.not. fu_between_times(fu_time(b_src%position), now, now+timestep, .true.))return
    !
    ! The blast is assumed to last one minute. Check the overlap
    !
    fOverlap = (fu_earliest_time((/fu_time(b_src%position)+one_minute, now+timestep/)) - &
              & fu_latest_time((/fu_time(b_src%position), now/))) / one_minute
    if(fOverlap < 1e-5) return
    !
    ! Get the parameters of the cloud
    !
    call init_mushroom_b_src(b_src, stem_bottom, stem_top, hat_top, stem_radius, hat_radius, base_surge_radius)
    if(error)return
    xSize => fu_grid_data(meteo_cell_x_size_fld)
    ySize => fu_grid_data(meteo_cell_y_size_fld)
    ix = nint(b_src%fxDispGrid) + nint(b_src%fyDispGrid-1) * nx_meteo
    fDx = 1./xSize(ix)
    fDy = 1./ySize(ix)

    ind_height = int_missing
    nullify(fldPressure)
    do ix = 1, size(met_buf%buffer_quantities)
      if(met_buf%buffer_quantities(ix) == int_missing)exit
      if(met_buf%buffer_quantities(ix) == height_flag)then
        ind_height = ix
      elseif(met_buf%buffer_quantities(ix) == pressure_flag)then
        ind_pressure = ix
        fldPressure => met_buf%p4d(ix)
      endif
    enddo
    if(fu_fails(ind_height /= int_missing, 'height field is not found in meteo buffer', &
                                         & 'inject_emis_lagr_bomb_src'))return
    if(fu_fails(associated(fldPressure), 'Pressure is not found in meteo buffer', &
                                       & 'inject_emis_lagr_bomb_src'))return
!call msg('1')

    nullify(speciesSrc)
    references => chemRunSetup%refEmis2Transp_mass
    nSpEmis = size(references)
    !
    ! General parameters of the release
    !

!call msg('2')
    !
    ! Rate has to be integrated over time. We have such function:
    ! total_from_p_src_species_unit(a_src, start, duration, layer, ifRateOnly, &
    !                                 species, nSpecies, amounts, z_moment)
    !
    call total_bomb_src_species_unit(b_src, speciesSrc, nSpeciesSrc, amounts, now, timestep)
    if(error)return

!call msg('3')

    !
    ! Lagrangian injection goes particle by particle. At least one particle is
    ! always released, so the particle masses are NOT identical - but an effort is 
    ! made to have them close. Within one release step they are identical.
    !
    ! Mass of a single particle and number of released particles: take the 
    ! low-mass threshold for each emission species, compare it to the rate of
    ! emission for this species, and take the largest ratio as the number of particles.
    ! Other species are distributed among these particles
    !
    nP = 1
    do iSpSrc = 1, nSpeciesSrc
      nP = max(nP, int(amounts(iSpSrc) / arParticleMass(b_src%adaptor2Trn%iSp(iSpSrc))))
!call msg('fMassTmp(iSpTransp), arParticleMass(iSpTransp)', fMassTmp(iSpTransp), arParticleMass(iSpTransp))
!call msg('ratio:',fMassTmp(iSpTransp) / arParticleMass(iSpTransp), nP)
    end do
    !
    ! Get the meteo level heights and project there the stem and hat limits, then turn to pressure
    ! Note however that final vertical location is still relative. Pressure is just to get the 
    ! distribution right
    !
    call column_from_buffer(met_buf, ind_height, nint(b_src%fxDispGrid), nint(b_src%fyDispGrid), &
                          & nx_meteo, meteo_heights, &
                          & interpCoefHoriz_void, interpCoefVert_void, &  ! whatever
                          & .false., .false., met_buf%weight_past)

    call column_from_buffer(met_buf, ind_pressure, nint(b_src%fxDispGrid), nint(b_src%fyDispGrid), &
                          & nx_meteo, meteo_pressure, &
                          & interpCoefHoriz_void, interpCoefVert_void, &  ! whatever
                          & .false., .false., met_buf%weight_past)
    if(error)return
    fStemBottomPressure = fu_4d_interpolation(fldPressure, &
                                            & b_src%fxDispGrid, b_src%fyDispGrid, &
                                            & max(0.50001, fu_value_index_in_array(stem_bottom, &
                                                                  & meteo_heights, nz_meteo) - 0.5), &
                                            & nx_meteo, ny_meteo, nz_meteo, &
                                            & met_buf%weight_past, &
                                            & linear, linear, notallowed)
    fStemTopPressure = fu_4d_interpolation(fldPressure, &
                                         & b_src%fxDispGrid, b_src%fyDispGrid, &
                                         & max(0.50001, fu_value_index_in_array(stem_top, &
                                                                  & meteo_heights, nz_meteo) - 0.5), &
                                         & nx_meteo, ny_meteo, nz_meteo, &
                                         & met_buf%weight_past, &
                                         & linear, linear, notallowed)
    fHatTopPressure = fu_4d_interpolation(fldPressure, &
                                        & b_src%fxDispGrid, b_src%fyDispGrid, &
                                        & max(0.50001, fu_value_index_in_array(hat_top, &
                                                                  & meteo_heights, nz_meteo) - 0.5), &
                                        & nx_meteo, ny_meteo, nz_meteo, &
                                        & met_buf%weight_past, &
                                        & linear, linear, notallowed)
    !
    ! The rest is simple: particles are distributed between the stem and hat and spread homogeneously 
    ! in each, then pressure is converted to relative index
    !
    iParticle = lpset%iFirstEmptyParticle
    arDyn    => lpset%lpDyn
    arMass   => lpset%lpMassTrn
    arStatus => lpset%lpStatus
    do iP = 1, nP
      !
      ! Find free space
      !
!call msg('11, iLev,iP',iLev,iP)
      do while(arStatus(iParticle) /= int_missing)
        iParticle = iParticle + 1

        if(iParticle >= lpset%nop)then
call msg('Enlarging the number of particles, 1:',  lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          call enlarge_lagrange_particles_set(lpset, lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          exit
        endif
      end do
      !
      ! Put the masses
      ! Emit the mass going species by species. Note multi-step re-indexing: from species index in 
      ! descriptor to species index in emission list, then to species index in transport list.
      ! Simialrly, MassTimeCommon is mass of a descriptor cocktail, which needs to be explored
      !
      do iSpSrc = 1, nSpeciesSrc
        arMass(b_src%adaptor2trn%iSp(iSpSrc), iParticle) = amounts(iSpSrc) / real(nP)
      end do
!call msg('12')

      !
      ! Position, turbulent movement, starting size.
      ! Note that vertical has to be converted to relative
      !
      if(iP <= nP * stem_mass_fract)then
        arDyn(lp_x, iParticle) = fu_random_number_center(b_src%fxDispGrid, min(0.4999,stem_radius*fDx))
        arDyn(lp_y, iParticle) = fu_random_number_center(b_src%fyDispGrid, min(0.4999,stem_radius*fDy))
        arDyn(lp_z, iParticle) = min(nz_meteo+0.49999, max(0.50001, &
                                   & fu_value_index_in_array( &
                                            & fu_random_number_boundaries(fStemBottomPressure, &
                                                                        & fStemTopPressure), &
                                            & meteo_pressure, nz_meteo) - 0.5))
        arDyn(lp_dx, iParticle) = stem_radius    ! metres
        arDyn(lp_dy, iParticle) = stem_radius    ! metres
        arDyn(lp_dz, iParticle) = stem_top - stem_bottom  ! metres
      else
        arDyn(lp_x, iParticle) = fu_random_number_center(b_src%fxDispGrid, min(0.4999,hat_radius*fDx))
        arDyn(lp_y, iParticle) = fu_random_number_center(b_src%fyDispGrid, min(0.4999,hat_radius*fDx))
        arDyn(lp_z, iParticle) = fu_random_number_boundaries(fStemTopPressure, fHatTopPressure)
        arDyn(lp_z, iParticle) = min(nz_meteo+0.49999, max(0.50001, &
                                   & fu_value_index_in_array( &
                                            & fu_random_number_boundaries(fStemTopPressure, &
                                                                        & fHatTopPressure), &
                                            & meteo_pressure, nz_meteo) - 0.5))
        arDyn(lp_dx, iParticle) = hat_radius    ! metres
        arDyn(lp_dy, iParticle) = hat_radius    ! metres
        arDyn(lp_dz, iParticle) = hat_top - stem_top  ! metres
      endif
      arStatus(iParticle) = b_src%id_nbr
      arDyn(lp_uT:lp_wT, iParticle) = 0.0    ! turbulent-wind motion

      !Keep track on fresh particles
      lpSet%nNewPart=lpSet%nNewPart+1
      lpSet%NewPartList(lpSet%nNewPart) = iParticle

      iParticle = iParticle + 1

!call msg('14')

      if(iParticle >= lpset%nop)then
call msg('Enlarging the number of particles, 2:',  lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          call enlarge_lagrange_particles_set(lpset, lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
      endif
    end do  ! iP

!call msg('15')

    !
    ! Fix the starting particle. Not exactly (place can be occupied) but better than starting from 1
    !
    lpset%iFirstEmptyParticle = iParticle

    
  end subroutine inject_emission_lagr_b_src
                                          
                                          
  !==========================================================================
  !==========================================================================
  !==========================================================================
  !
  !  Some encapsulation stuff
  !
  !==========================================================================
  !==========================================================================
  !==========================================================================

  logical function fu_bomb_source_defined(b_src)
    implicit none
    type(silam_bomb_source), intent(in) :: b_src
    fu_bomb_source_defined = b_src%defined == silja_true
  end function fu_bomb_source_defined

  !==========================================================================
  REAL FUNCTION fu_bomb_source_yield(b_src)
    IMPLICIT NONE
    TYPE(silam_bomb_source), INTENT(in) :: b_src
    fu_bomb_source_yield = b_src%total_yield
  END FUNCTION fu_bomb_source_yield

  !==========================================================================
  FUNCTION fu_bomb_source_position(b_src)
    IMPLICIT NONE
    TYPE(silam_grid_position) :: fu_bomb_source_position
    TYPE(silam_bomb_source), INTENT(in) :: b_src
    fu_bomb_source_position = b_src%position
  END FUNCTION fu_bomb_source_position

  !==========================================================================
  FUNCTION fu_bomb_source_start_time(b_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_bomb_source_start_time
    TYPE(silam_bomb_source), INTENT(in) :: b_src
    fu_bomb_source_start_time = fu_time(b_src%position)
  END FUNCTION fu_bomb_source_start_time

  !==========================================================================
  FUNCTION fu_bomb_source_end_time(b_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_bomb_source_end_time
    TYPE(silam_bomb_source), INTENT(in) :: b_src
    fu_bomb_source_end_time = fu_time(b_src%position) + one_minute ! + b_src%duration
  END FUNCTION fu_bomb_source_end_time

  !==========================================================================
  FUNCTION fu_bomb_source_duration(b_src)
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_bomb_source_duration
    TYPE(silam_bomb_source), INTENT(in) :: b_src
    fu_bomb_source_duration = one_minute   ! One minute is assumed
  END FUNCTION fu_bomb_source_duration
  !==========================================================================
  function fu_venting_fraction(yield, height) result(vf)
    implicit none
    real :: yield, height, sdob, vf
    if (height >= 0.0) then
       vf = 1.0
    else !underground detonation
       sdob = -height*3.2808399/(yield**0.294) !height transfered to feet and positive
       if (sdob > 200.0) then !when scaled depth of burial is deeper than 200 feet bomb has no effect
          vf = 0.0
          call set_error('Given bomb size and depth buries the bomb totally, no emission or emission would be highly uncertain', &
                       & 'fu_venting_fraction')
       else
          !function fitted to match KDFOC venting fraction figure:
          !small burial doesn't give much effect on the escaping amount,
          !but then at some point it starts quite rapidly to be significant
          vf = -0.6831437 + (1.000321 + 0.6831437)/(1 + (sdob/183.9856)**4.245358)
          if (vf < 0.0) call set_error('Given bomb size and depth buries the bomb totally, no emission or emission would be highly uncertain', &
                                     & 'fu_venting_fraction')
          if (vf > 0.97) vf = 0.97
       end if
    end if
  end function fu_venting_fraction
  
  !==========================================================================
  function fu_fission_yield(yield, height, fission_fraction) result(fission_yield)
    !Fusion neutron activation products induce 'fission' yield (KDFOC3)
    ! in addition to the actual fission yield.
    implicit none
    real :: yield, height, fission_fraction, &
          & fission_yield, induced_buried, induced_surface, induced_bomb
    !2% more for bomb material activation
    induced_bomb = 0.02
    !8% more for soil activation in surface bursts
    induced_surface = 0.08
    !8% more for soil activation for buried bursts
    if (height < 0.0) then !buried
       induced_buried = 0.08
    else
       induced_buried = 0.0
    end if
    fission_yield = yield*(fission_fraction + &
                         & (1.0 - fission_fraction) * &
                         & (induced_bomb + induced_surface + induced_buried) &
                         & )
  end function fu_fission_yield
  
  !==========================================================================
  function fu_fireball_buried_fraction(yield, height) result(ff)
    !Is one when height is 0 and is zero when height is
    !fireball_radius = 55.0*(yield**0.3) !DELFIC
    !fireball_buried_fraction = 0.45345**(0.05048*height*(yield**(1.0/3.0))) !DELFIC
    implicit none
    real :: yield, height, fireball_radius, ff
    !Fireball radius is calculated assuming activation products from fusion neutron
    ! activation from surface and bomb material is maxed: overestimates the fireball
    ! size for bombs above surface in order not to underestimate the actual activation products
    ! calculated later.
    ! NOTE: fireball radius is different here than in other uses (KDFOC3)
    fireball_radius = 55.0*(yield**0.4)
    if (height < 0.0) then !underground
       ff = 1.0
    else if (height >= fireball_radius) then !too high, fireball doesn't touch the ground
       ff = 0.0
    else !not surface detonation, but fireball reaches surface
       !ff is the volume of the half sphere that is underground divided
       !   by the volume of the half sphere
       ff = (fireball_radius - height)**2.0 * (2.0*fireball_radius + height) / & ! underground volume
          & (2.0 * fireball_radius**3.0)         ! half-sphere volume, Pi/3 factors cancelled
    end if
  end function fu_fireball_buried_fraction
  
  !==========================================================================
  function fu_fireball_radius(yield) result(fr)
    implicit none
    real :: yield, fr
    fr = 30.0*(yield**0.333) !KDFOC3
  end function fu_fireball_radius
  
  !==========================================================================
  subroutine get_horizontal_coverage_and_mom_4_corners(cell_fraction_in_cloud, momX, momY, &
                                              & radius, iDistcells, fDX, fDY, xShift, yShift)
    !
    ! Projects circle to grid: Check if centers of the refined  grid are within the circle
    ! 

    implicit none
    !out
    real, dimension(-iDistCells:iDistCells,-iDistCells:iDistCells), &
         & intent(out) :: cell_fraction_in_cloud, momX, momY
    !in
    integer, intent(in) :: iDistcells
    real, intent(in) :: radius, fDX, fDY, xShift, yShift
    !local
    integer :: ix, iy
    real :: fDist2, R2,total
    
    cell_fraction_in_cloud = 0.0
    momX = 0.0
    momY = 0.0
    
    R2 = radius*radius
    do iy = -iDistCells, iDistCells
       do ix = -iDistCells, iDistCells
          
          ! Corner subcell 1 (-x,-y)
          fDist2 = (((ix-0.25 - xShift)*fDX)**2 + ((iy-0.25 - yShift)*fDY)**2.0)
          if (fDist2 < R2) then
             cell_fraction_in_cloud(ix,iy) = cell_fraction_in_cloud(ix,iy) + 0.25
             momX(ix,iy) = momX(ix,iy) - 0.125
             momY(ix,iy) = momY(ix,iy) - 0.125
          end if
          ! Corner subcell 2 (+x,-y)
          fDist2 = (((ix+0.25 - xShift)*fDX)**2 + ((iy-0.25 - yShift)*fDY)**2.0)
          if (fDist2 < R2) then
             cell_fraction_in_cloud(ix,iy) = cell_fraction_in_cloud(ix,iy) + 0.25
             momX(ix,iy) = momX(ix,iy) + 0.125
             momY(ix,iy) = momY(ix,iy) - 0.125
          end if
          ! Corner subcell 3 (-x,+y)
          fDist2 = (((ix-0.25 - xShift)*fDX)**2 + ((iy+0.25 - yShift)*fDY)**2.0)
          if (fDist2 < R2) then
             cell_fraction_in_cloud(ix,iy) = cell_fraction_in_cloud(ix,iy) + 0.25
             momX(ix,iy) = momX(ix,iy) - 0.125
             momY(ix,iy) = momY(ix,iy) + 0.125
          end if
          ! Corner subcell 4 (+x,+y)
          fDist2 = (((ix+0.25 - xShift)*fDX)**2 + ((iy+0.25 - yShift)*fDX)**2)
          if (fDist2 < R2) then
             cell_fraction_in_cloud(ix,iy) = cell_fraction_in_cloud(ix,iy) + 0.25
             momX(ix,iy) = momX(ix,iy) + 0.125
             momY(ix,iy) = momY(ix,iy) + 0.125
          end if
       end do   ! ix +- iDistCells
    end do  ! iy +- iDistCells

    total = sum(cell_fraction_in_cloud)

    if (total == 0.) then !! None of the cells hit: circle within 
      cell_fraction_in_cloud(0,0) = 1.
      momX(0,0) = xShift
      momY(0,0) = yShift
    else  !! Just  normalize
      cell_fraction_in_cloud = cell_fraction_in_cloud / total
    endif


  end subroutine get_horizontal_coverage_and_mom_4_corners

  !==========================================================================
  
  subroutine uranium_ind_stuk(blast_nuclides)
    !
    !Returns list of nuclides and corresponding activities in becquerel/kiloton
    !compiled by STUK and based mostly on two papers by Kraus and Foster and
    !England and Rider both from 2013. For a possible improvised nuclear device
    !with 94% U_235 and 6% U_238 as device fuel. Detailed description by STUK 
    !of origin and considerations hopefully saved somewhere.
    !Nuclide-activity/kt pairs listed parallel for easy value checking in the 
    !case that activities need to be later revised.
    !
    !Ps. list gives untrimmed nuclide names -> make sure to trim them before use
    !
    implicit none
    
    ! Imported parameters
    type(Tblast_nuclide), dimension(:), allocatable, intent(out) :: blast_nuclides
    
    ! Local parameters and variables
    integer, parameter :: number_of_nuclides = 56
    
    allocate(blast_nuclides(number_of_nuclides))

    blast_nuclides(1:number_of_nuclides) = &
              & (/Tblast_nuclide('KR_85M' ,  8.24E+16), &
                & Tblast_nuclide('KR_85  ',  8.66E+09), &
                & Tblast_nuclide('KR_87  ',  5.45E+17), &
                & Tblast_nuclide('KR_88  ',  3.27E+17), &
                & Tblast_nuclide('XE_133M',  2.15E+13), &
                & Tblast_nuclide('XE_133 ',  3.10E+12), &
                & Tblast_nuclide('XE_135M',  1.93E+17), &
                & Tblast_nuclide('XE_135 ',  3.48E+15), &
                & Tblast_nuclide('XE_137 ',  2.63E+19), &
                & Tblast_nuclide('XE_138 ',  6.93E+18), &
                & Tblast_nuclide('BA_140 ',  5.37E+15), &
                & Tblast_nuclide('BA_141 ',  5.31E+18), &
                & Tblast_nuclide('BA_142 ',  8.48E+18), &
                & Tblast_nuclide('CE_141 ',  1.08E+10), &
                & Tblast_nuclide('CE_143 ',  1.11E+13), &
                & Tblast_nuclide('CE_144 ',  2.10E+14), &
                & Tblast_nuclide('CO_58  ',  4.54E+13), &
                & Tblast_nuclide('CO_58M ',  8.43E+15), &
                & Tblast_nuclide('CS_137 ',  2.38E+11), &
                & Tblast_nuclide('CS_138 ',  2.60E+17), &
                & Tblast_nuclide('I_131  ',  1.47E+12), &
                & Tblast_nuclide('I_132  ',  1.32E+15), &
                & Tblast_nuclide('I_133  ',  5.24E+15), &
                & Tblast_nuclide('I_134  ',  2.27E+17), &
                & Tblast_nuclide('I_135  ',  2.62E+17), &
                & Tblast_nuclide('LA_141 ',  3.16E+14), &
                & Tblast_nuclide('LA_142 ',  5.15E+15), &
                & Tblast_nuclide('LA_143 ',  6.55E+18), &
                & Tblast_nuclide('MN_54  ',  8.45E+12), &
                & Tblast_nuclide('MN_56  ',  7.65E+17), &
                & Tblast_nuclide('MO_99  ',  2.52E+16), &
                & Tblast_nuclide('MO_101 ',  5.83E+18), &
                & Tblast_nuclide('RU_103 ',  9.96E+14), &
                & Tblast_nuclide('RU_106 ',  2.09E+13), &
                & Tblast_nuclide('SB_128 ',  1.27E+14), &
                & Tblast_nuclide('SB_129 ',  5.37E+16), &
                & Tblast_nuclide('SB_130 ',  1.27E+17), &
                & Tblast_nuclide('SB_131 ',  2.10E+18), &
                & Tblast_nuclide('SN_128 ',  1.44E+17), &
                & Tblast_nuclide('SR_89  ',  9.88E+14), &
                & Tblast_nuclide('SR_90  ',  5.84E+12), &
                & Tblast_nuclide('SR_91  ',  1.63E+17), &
                & Tblast_nuclide('SR_92  ',  5.88E+17), &
                & Tblast_nuclide('TC_104 ',  2.04E+18), &
                & Tblast_nuclide('TE_131 ',  5.32E+16), &
                & Tblast_nuclide('TE_131M',  2.23E+15), &
                & Tblast_nuclide('TE_132 ',  1.67E+16), &
                & Tblast_nuclide('TE_133 ',  3.39E+18), &
                & Tblast_nuclide('TE_133M',  8.46E+17), &
                & Tblast_nuclide('TE_134 ',  2.57E+18), &
                & Tblast_nuclide('Y_92   ',  2.16E+14), &
                & Tblast_nuclide('Y_93   ',  1.69E+17), &
                & Tblast_nuclide('Y_94   ',  5.37E+18), &
                & Tblast_nuclide('Y_95   ',  9.76E+18), &
                & Tblast_nuclide('ZR_95  ',  2.54E+12), &
                & Tblast_nuclide('ZR_97  ',  9.81E+16) /)
  end subroutine uranium_ind_stuk

!============================================================================

  subroutine plutonium_stuk(blast_nuclides)
    !Does same as uranium_ind_stuk but for a plutonium device.
    implicit none
    
    type(Tblast_nuclide), dimension(:), allocatable, intent(out) :: blast_nuclides

    integer, parameter :: number_of_nuclides = 56
    
    allocate(blast_nuclides(number_of_nuclides))

    blast_nuclides(1:number_of_nuclides) = &
              &(/ Tblast_nuclide('KR_85M ',  3.71E+16), &
                & Tblast_nuclide('KR_85  ',  2.81E+10), &
                & Tblast_nuclide('KR_87  ',  2.29E+17), &
                & Tblast_nuclide('KR_88  ',  1.27E+17), &
                & Tblast_nuclide('XE_133M',  2.47E+14), &
                & Tblast_nuclide('XE_133 ',  3.51E+13), &
                & Tblast_nuclide('XE_135M',  9.26E+17), &
                & Tblast_nuclide('XE_135 ',  1.88E+16), &
                & Tblast_nuclide('XE_137 ',  2.45E+19), &
                & Tblast_nuclide('XE_138 ',  5.60E+18), &
                & Tblast_nuclide('BA_140 ',  4.85E+15), &
                & Tblast_nuclide('BA_141 ',  4.65E+18), &
                & Tblast_nuclide('BA_142 ',  7.01E+18), &
                & Tblast_nuclide('CE_141 ',  8.81E+10), &
                & Tblast_nuclide('CE_143 ',  2.44E+14), &
                & Tblast_nuclide('CE_144 ',  1.50E+14), &
                & Tblast_nuclide('CO_58  ',  4.54E+13), &
                & Tblast_nuclide('CO_58M ',  8.43E+15), &
                & Tblast_nuclide('CS_137 ',  1.06E+12), &
                & Tblast_nuclide('CS_138 ',  3.55E+17), &
                & Tblast_nuclide('I_131  ',  2.91E+13), &
                & Tblast_nuclide('I_132  ',  2.13E+16), &
                & Tblast_nuclide('I_133  ',  1.74E+16), &
                & Tblast_nuclide('I_134  ',  3.82E+17), &
                & Tblast_nuclide('I_135  ',  2.58E+17), &
                & Tblast_nuclide('LA_141 ',  4.93E+15), &
                & Tblast_nuclide('LA_142 ',  4.93E+16), &
                & Tblast_nuclide('LA_143 ',  5.12E+18), &
                & Tblast_nuclide('MN_54  ',  8.45E+12), &
                & Tblast_nuclide('MN_56  ',  7.65E+12), &
                & Tblast_nuclide('MO_99  ',  2.53E+16), &
                & Tblast_nuclide('MO_101 ',  7.65E+18), &
                & Tblast_nuclide('RU_103 ',  2.02E+15), &
                & Tblast_nuclide('RU_106 ',  1.36E+14), &
                & Tblast_nuclide('SB_128 ',  1.14E+15), &
                & Tblast_nuclide('SB_129 ',  9.20E+16), &
                & Tblast_nuclide('SB_130 ',  2.43E+17), &
                & Tblast_nuclide('SB_131 ',  2.10E+18), &
                & Tblast_nuclide('SN_128 ',  2.24E+17), &
                & Tblast_nuclide('SR_89  ',  3.98E+14), &
                & Tblast_nuclide('SR_90  ',  2.25E+12), &
                & Tblast_nuclide('SR_91  ',  7.38E+16), &
                & Tblast_nuclide('SR_92  ',  3.11E+17), &
                & Tblast_nuclide('TC_104 ',  6.05E+18), &
                & Tblast_nuclide('TE_131 ',  1.68E+17), &
                & Tblast_nuclide('TE_131M',  6.15E+15), &
                & Tblast_nuclide('TE_132 ',  1.84E+16), &
                & Tblast_nuclide('TE_133 ',  1.66E+18), &
                & Tblast_nuclide('TE_133M',  1.11E+18), &
                & Tblast_nuclide('TE_134 ',  1.91E+18), &
                & Tblast_nuclide('Y_92   ',  1.59E+15), &
                & Tblast_nuclide('Y_93   ',  1.05E+17), &
                & Tblast_nuclide('Y_94   ',  3.78E+18), &
                & Tblast_nuclide('Y_95   ',  7.45E+18), &
                & Tblast_nuclide('ZR_95  ',  1.69E+13), &
                & Tblast_nuclide('ZR_97  ',  8.76E+16) /)
    
  end subroutine plutonium_stuk

!============================================================================

  logical function if_activation_isotope(name)
    !
    !True if isotope is formed via activation of ground.
    !
    !Input
    character(len=*) :: name
    !Internal
    select case(trim(name))
      case('CO_58', 'CO_58M', 'MN_54', 'MN_56')
        if_activation_isotope = .true.
      case default
        if_activation_isotope = .false.
    end select
  
  end function if_activation_isotope
  
!============================================================================
  
  subroutine set_cocktail_bomb_source(b_src)
    !
    ! Routine to set cocktail description with the list of nuclides and
    ! particle size bins.
    !
    implicit none
    
    ! Input parameter
    type(silam_bomb_source), intent(inout) :: b_src
    
    ! Local variables
    integer :: n, i, m
    type(silam_species) :: species_tmp
    real, dimension(max_species) :: act_fractTmp
    real :: totalFraction
    type(silam_material), pointer :: tmp_material
    real, dimension(2) :: act_frac_bomb_modes
    integer :: iBombMode
    type(Tblast_nuclide), dimension(:), allocatable :: blast_nuclides
    type(Taerosol_mode), dimension(2) :: aerModesBomb
    real :: fTmp
    !
    if(fu_fails(fu_true(b_src%defined),'Source not defined','set_cocktail_bomb_source'))return
    !
    ! Get nuclides and activities
    !
    if (b_src%bomb_type == 'URANIUM') then
      call uranium_ind_stuk(blast_nuclides)
    else if (b_src%bomb_type == 'PLUTONIUM') then
      call plutonium_stuk(blast_nuclides)
    else if (b_src%bomb_type == 'PASSIVE') then
       allocate(blast_nuclides(2))
       blast_nuclides(1:2) = &
              &(/ Tblast_nuclide('passive',  1e10), &
                & Tblast_nuclide('PM  ',     1e10) /)
    else
      call set_error('Invalid bomb type given','set_cocktail_bomb_source')
    end if
    if (error) return
    !
    ! Get the lognormal distributions for the given elevation of the blast
    !
    call set_two_modes_bomb(b_src, aerModesBomb, act_frac_bomb_modes)
    if(error)return
    !
    ! scan the nuclides released by the blast, explore the aerosol bins for the particulates
    ! We shall keep the budget and dump the un-allocated activity fraction to the log file
    !
    nullify(b_src%species)
    b_src%nSpecies = 0
    do n = 1, size(blast_nuclides)
      !
      ! Nuclides can come from the boms or from the contaminated soil lifted in the air
      ! Activation species depend on the amount of activation of the ground
      !
      if (if_activation_isotope(blast_nuclides(n)%nuc_name))then
        if (b_src%fireball_buried_fraction == 0.0) cycle  ! no soil lifted to the air by the fireball
        blast_nuclides(n)%bq_per_kt = blast_nuclides(n)%bq_per_kt * b_src%fireball_buried_fraction
       end if
      tmp_material => fu_get_material_ptr(trim(blast_nuclides(n)%nuc_name))
      ! gas or aerosol?
      if (fu_true(fu_if_gas(fu_get_material_ptr(trim(blast_nuclides(n)%nuc_name))))) then 
        !set gas 
          call set_species(species_tmp, tmp_material, in_gas_phase)
          call addSpecies(b_src%species, b_src%nSpecies, (/species_tmp/), 1)
        act_fractTmp(b_src%nSpecies) = blast_nuclides(n)%bq_per_kt * b_src%fission_yield ! activity is determined by fission yield, not total yield
      else 
        ! if it can be aerosol set aerosol(s) 
        ! evidently wrong for multi-phase species but for now this is the solution
        ! only a fraction of total amount that falls within specified modes goes to SILAM
        totalFraction = 0.0   ! will accumulate amounts for the modes  of the material
        do i = 1,(b_src%aerosolSrc%n_modes)      ! Loop over bins of the required aerosol
          ! add the species
          call set_species(species_tmp, tmp_material,b_src%aerosolSrc%modes(i))
          call addSpecies(b_src%species, b_src%nSpecies, (/species_tmp/), 1)

          act_fractTmp(b_src%nSpecies) = 0
          do iBombMode = 1, 2
            fTmp = act_frac_bomb_modes(iBombMode) 
            if (fTmp > 0) act_fractTmp(b_src%nSpecies) = act_fractTmp(b_src%nSpecies) + &
                                         & fTmp * fu_integrate_volume(fu_min_d(b_src%aerosolSrc%modes(i)), &
                                                             & fu_max_d(b_src%aerosolSrc%modes(i)), &
                                                             & aerModesBomb(iBombMode)) 
          end do
           ! collect what included
          totalFraction = totalFraction + act_fractTmp(b_src%nSpecies)
          ! ... and its activity
          act_fractTmp(b_src%nSpecies) = act_fractTmp(b_src%nSpecies) * &
                                       & blast_nuclides(n)%bq_per_kt * b_src%fission_yield
        end do  ! b_src%aerosolSrc bins

        ! If something is left behind, dump it to the log file
        if(totalFraction < 1.0)then
          call msg('>> Accounted fraction and amount for:' + trim(blast_nuclides(n)%nuc_name),&
                 & totalFraction, &
                 & totalFraction * blast_nuclides(n)%bq_per_kt * b_src%fission_yield)
        endif
      end if  ! gas or aerosol
    end do  ! sizeof blast_nuclides
    !
    ! Allocate the main activity split and copy them from the temporary
    !
    allocate(b_src%activities(b_src%nSpecies))
    b_src%activities(1:b_src%nSpecies) = act_fractTmp(1:b_src%nSpecies)
    deallocate(blast_nuclides)

  end subroutine set_cocktail_bomb_source


  !============================================================================
  
  subroutine set_two_modes_bomb(b_src, aerModesBomb, act_frac)
    !
    ! Determines correct aerosol activity distribution based on bomb source.
    ! Is based on KDFOC3 model.
    ! Values are in micrometers.
    !
    implicit none
    !Input
    type(silam_bomb_source), intent(in) :: b_src
    !Output
    real, dimension(2), intent(out) :: act_frac !1 is for smaller and 2 is for larger diameter
    type(Taerosol_mode), dimension(2), intent(out) :: aerModesBomb
    
    ! Local parameters
    ! from KDFOC3: A Nuclear Fallout Assessment Capability 1993 p.26-28
    real, dimension(2), parameter :: KDFOC_surf_r = (/14.44, 151.41/)
    real, dimension(2), parameter :: KDFOC_buri_r = (/90.02, 298.87/)
    real, dimension(2), parameter :: KDFOC_surf_s = (/4.01, 2.69/)
    real, dimension(2), parameter :: KDFOC_buri_s = (/2.01, 1.82/)
    real, dimension(2), parameter :: KDFOC_uL = (/0.23, 0.65/) !surface, buried
    ! MATCH for airburst
    real, parameter :: MATCH_rA = 4.4, MATCH_sA = 1.5
    
    ! Implications of Atmospheric Test Fallout Data for Nuclear Winter, 
    ! Baker, 1987 and LANL report: User Guide fro the Air Force Nuclear Weapons Center Dust 
    ! Cloud Calculator Version 1.0, St Ledger, John W. 2015
    real, dimension(2), parameter :: BAKER_r = (/0.42, 9.34/)  ! radius, um
    real, dimension(2), parameter :: BAKER_s = (/2.0, 4.0/)  
    real, parameter :: BAKER_uL = 0.75
    
    ! Local variables
    real :: fTmp
    real, dimension(2) :: modes_r, stds
    integer :: iMode
    
    modes_r = real_missing
    stds = real_missing
    act_frac = real_missing
    !
    ! Have two lognormal modes, scan both. Explcit coefs for activity is not an accident
    !
    if (b_src%dist_type == 'KDFOC3') then
      select case(b_src%location_switch)
        case(blast_underground)
          if (b_src%blast_height_m <= -15.0) then !deep underground
            modes_r(:) = KDFOC_buri_r(:)
            stds(:) = KDFOC_buri_s(:)
            act_frac(2) = KDFOC_uL(2)
            act_frac(1) = 1.0 - act_frac(2)
          else   ! (b_src%blast_height_m <= 0.0) then !shallow underground
            ! Linear change from surface burst mean and s to deep one in shallow bursts
            fTmp = b_src%blast_height_m / (-15.)
            modes_r(:) = KDFOC_surf_r(:) - ((KDFOC_buri_r(:) + KDFOC_surf_r(:)) * fTmp)
            stds(:) = KDFOC_buri_s(:) - ((KDFOC_surf_s(:) - KDFOC_buri_s(:)) * (1.0 - fTmp))
            act_frac(2) = KDFOC_uL(1) + ((KDFOC_uL(2) - KDFOC_uL(1)) * fTmp)
            act_frac(1) = 1.0 - act_frac(2)
          endif
          
        case(blast_in_air)  !Air burst, fireball doesn't touch the ground
          !KDFOC3 description doesn't handle this case, but MATCH model that is based partly
          !on it makes bins from the given lognormal distribution. When handling air bursts it
          !concentrates all of the activity to the two smallest ones r = 2.2 and 4.4(micrometers).
          !(4.4 only minuscule). Thus we will base our guess on MATCH's airburst configuration,
          !but make the system continuous from surface bursts.
          !-> In air bursts the larger lognormal distribution will disappear and the smaller
          !is centered in d = 4.4(r=2.2) with standard deviation 1.0 (resulting in a reasonable distribution).
          modes_r(1) = MATCH_rA
          modes_r(2) = 0.0
          stds(1) = MATCH_sA
          stds(2) = 0.0
          act_frac(1) = 1.0
          act_frac(2) = 0.0

        case(blast_on_surface)
          ! Fireball touches the ground, linear interpolation between surface and air burst, based on fireball fraction
          modes_r(:) = MATCH_rA + (KDFOC_surf_r(:) - MATCH_rA) * b_src%fireball_buried_fraction
          stds(1) = MATCH_sA + (KDFOC_surf_s(1) - 1.0)*b_src%fireball_buried_fraction
          stds(2) = KDFOC_surf_s(2)
          act_frac(1) = KDFOC_uL(1) + ((1.0 - KDFOC_uL(1))*(1.0 - b_src%fireball_buried_fraction))
          act_frac(2) = 1.0 - act_frac(1)
        case default
          call set_error('Unknown blast height / type','set_two_modes_bomb')
          return
      end select
        
    else if (b_src%dist_type == 'BAKER') then
      
      select case(b_src%location_switch)
        case(blast_in_air)
          modes_r(1) = BAKER_r(1)
          modes_r(2) = 0.0
          stds(1) = BAKER_s(1)
          stds(2) = 0.0
          act_frac(1) = 1.0
          act_frac(2) = 0.0
        case default  ! surface and underground
          modes_r(1:2) = BAKER_r(1:2)
          stds(1:2) = BAKER_s(1:2)
          act_frac(1) = b_src%fireball_buried_fraction*(1.0 - BAKER_uL)
          act_frac(2) = 1.0 - act_frac(1)
      end select
    end if  ! types of parameterizations

    do iMode = 1, 2
      if (modes_r(iMode) > 0) then
        fTmp = 2*modes_r(iMode) * 1e-6 ! diameter im meters
        call set_aerosol_mode(aerModesBomb(iMode), &
                   &'BombMode'//trim(b_src%dist_type)//trim(fu_str(iMode)), &     ! mode, chNm, 
                   & fTmp, &      ! fp1
                   & stds(iMode), &         ! fp2  -- dimensionless for lognormal
                   & fTmp, 2.5e3, &      ! mass_mean_d, dens, 
                   & lognormal_flag, 1)     ! distr_type, solubil
        if(error)return
      else
        aerModesBomb(iMode) = aerosol_mode_missing
      endif
      call msg("Bomb mode", iMode)
      call report(aerModesBomb(iMode))
    end do
    
    if (sum(act_frac) > 1.00001) call set_error('Mass fraction sum illegal','set_two_modes_bomb')
    if (any(abs(act_frac) > 1.00001)) call set_error('Fraction cant be over unity','set_two_modes_bomb')
    if (error) return
    call msg('=============================================================')
    if (b_src%location_switch == blast_in_air) call msg('Air blast, no, Distribution type: '//trim(b_src%dist_type))
    if (b_src%location_switch == blast_on_surface) call msg('Surface blast, Distribution type: '//trim(b_src%dist_type))
    if (b_src%location_switch == blast_underground) call msg('Underground blast, Distribution type: '//trim(b_src%dist_type))
    call msg('|LOGNORMAL DISTRIBUTIONS CREATED. Small        Larger')
    call msg(' |Radii (in micrometers)         :',(/modes_r(1), modes_r(2)/))
    call msg(' |Standard deviations            :',(/stds(1), stds(2)/))
    call msg(' |Activity fractions             :',(/act_frac(1), act_frac(2)/))
    call msg('=============================================================')
    
  end subroutine set_two_modes_bomb
 

  !============================================================================
  
  subroutine init_bomb_derived_params(b_src, ifWater)
    implicit none
    type(silam_bomb_source) :: b_src
    logical, intent(in) :: ifWater
    real :: blast_height
    real :: j
    integer :: i
    !
    ! Normal density of land/soil is 1600kg/m3 and water is 1000kg/m3
    ! the ratio of densities affect the depth of detonation and venting
    ! so that lighter burial material needs deeper burial for same degree
    ! of venting and vice versa.
    !
    ! Currently only water used land material information would enable using similar 
    ! logic for other medias of burial.
    !
!    ifWater = .false.
    if (ifWater .and. (b_src%blast_height_m < 0.0)) then
       blast_height = b_src%blast_height_m*(1.0/1.6)
    else
       blast_height = b_src%blast_height_m
    end if
    !
    ! Venting fraction determined based on blast height and yield
    !
    b_src%venting_fraction = fu_venting_fraction(b_src%total_yield, blast_height)
    !
    ! Fireball fraction determined based on blast height and yield
    !
    b_src%fireball_buried_fraction = fu_fireball_buried_fraction(b_src%total_yield, &
                                                               & b_src%blast_height_m)
    !
    ! Parameter fission_yield determines the activities induced by the fission
    ! part of the bomb. The total_yield depicts the power of the bomb; reach of cloud etc.
    !
    b_src%fission_yield = fu_fission_yield(b_src%total_yield, &
                                         & b_src%blast_height_m, &
                                         & b_src%fission_fraction)
    ! Fireball radius
    !
    b_src%fireball_radius = fu_fireball_radius(b_src%total_yield)
    !
    ! Three logicals to easily check which type of detonation
    !
    if (b_src%blast_height_m < 0.0) then
      b_src%location_switch = blast_underground
    else if (b_src%fireball_buried_fraction == 0.0) then
      b_src%location_switch = blast_in_air
    else
      b_src%location_switch = blast_on_surface
    end if
    
  end subroutine init_bomb_derived_params

  
  !============================================================================

  subroutine get_bomb_bins_from_2_lognorm_modes(nbins, min_d, max_d, mu1, mu2, s1, s2, f1, f2, &
                                              & diameters, total_fractions)
    !
    ! Creates number of bins requested from two lognormal modes.
    ! - Larger than max_d will be thrown out/instant deposition.
    ! - If nbins is 10, 9 bins will be created and the tenth 'bin' will be
    !   values higher than max_d and not included because SILAM
    !   can't handle too large aerosols (FIXME Really???).
    ! - Smaller than min_d will be included in smallest bin.
    ! - Routine returns the representing bins and corresponding activity fractions.
    ! - The ratio of the sizes of parallel bins is constant.
    ! - Takes diameters, deals internally with radii and returns diameters.
    !
    implicit none
    !Input
    real, intent(in) :: min_d, max_d !Diameters
    integer, intent(in) :: nbins
    real, intent(in) :: mu1, mu2, s1, s2, f1, f2 !1 -> smaller mode, 2 -> larger
    !Output
    real, dimension(nbins), intent(out) :: total_fractions
    real, dimension(nbins-1), intent(out) :: diameters
    
    !Internal
    real :: factor, min_r, max_r, fractions1, fractions2, radius
    real(r8k) :: erf_input1, erf_input2
    integer :: i, j
    real, dimension(nbins) :: bounds
    ! Handle with radii inside the subroutine
    min_r = min_d*0.5
    max_r = max_d*0.5
    ! Ratio between parallel bins
    factor = (max_r/min_r)**(1.0/(nbins-1))
    ! Bounds based on ratio
    do i = 1,nbins
       bounds(i) = min_r*(factor**(i-1))
    end do
    ! From two lognormal distributions get fractions for each bin
    ! Representing diameter is the average within a bin
    do i = 1,(nbins-1)
       ! Total from both with given fractions f1 and f2 of the modes (not bins)
       total_fractions(i)  = 0.
       if (f1 > 0) then
         ! Error function input has to be real 8 kind
         erf_input1 = (log(bounds(i+1)) - log(mu1))/(log(s1)*(2.0**(1.0/2.0)))
         ! Cumulative lognormal distribution (fraction of total) under the given radius
         fractions1 = (1.0 + ERF(erf_input1))/(2.0)
         total_fractions(i) = total_fractions(i) + fractions1*f1
       endif
       if (f2 > 0) then
         erf_input2 = (log(bounds(i+1)) - log(mu2))/(log(s2)*(2.0**(1.0/2.0)))
         fractions2 = (1.0 + ERF(erf_input2))/(2.0)
         total_fractions(i) = total_fractions(i) + fractions2*f2
       endif

       ! Not the cumulative
       if (i > 1) then
          total_fractions(i) = total_fractions(i) - sum(total_fractions(:(i-1)))
       end if
       ! Radius as average
       radius = (bounds(i) + bounds(i+1))/2.0
       ! Diameter outputted
       diameters(i) = radius * 2.0
    end do
    ! 'Thrash' bin gets rest of the fraction
    total_fractions(nbins) = 1.0 - sum(total_fractions(:(nbins - 1)))
  end subroutine get_bomb_bins_from_2_lognorm_modes
  

END MODULE source_terms_bomb
