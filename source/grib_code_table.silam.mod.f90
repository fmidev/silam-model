MODULE grib_code_table

  ! This module contains all necessary information about GRIB code tables
  ! used/allowed to use in SILAM.
  ! It underpins the nwpm_administrations module and uses the 
  ! names_of_quantities. The purpose of this module is to establish 
  ! the 2-way link between the SILAM quantities and their coding in GRIB
  ! in various code tables 2.
  ! 
  !
  ! Author: Mikhail Sofiev, FMI mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  !
!  USE globals
!  USE names_of_quantities
  !use levels
  use stacks  !grads_io

  implicit none
  
  private

  public read_GRIB_code_table_v5
  public get_silam_params_for_grib_1
  public get_silam_params_for_grib_2
  public get_GRIB_codes

!  public defined
  public report

  private detect_SILAM_level
  private CTE_defined
  private report_code_table

!  interface defined
!    module procedure CTE_defined
!  end interface

  interface report
    module procedure report_code_table
  end interface

  !
  ! Baseline information dataset for GRIB coding procedure for GRIB 1 and GRIB 2 formats.
  ! Note that this format should be sufficient for everything: centre is there and it is 
  ! believed that one centre will not use more than one code table.
  !
  type Tgrib_code_table_entry
    private
    integer, dimension(10) :: centre, sub_centre ! GRIB 1 & 2
    integer, dimension(10) :: grib1_table_version            ! GRIB 1
    integer, dimension(10) :: grib1_model_id  ! GRIB 1
    integer :: discipline, paramCategory, paramNumber    ! GRIB 2
    character(len=clen) :: name                  ! GRIB 1 & 2
    character(len=substNmLen) :: short_name      ! GRIB 1 & 2
    character(len=unitNmLen)  :: unit            ! GRIB 1 & 2
    integer :: typeOfFirstFixedSurface                   ! GRIB 2
    real    :: ValueOfFirstFixedSurface                  ! GRIB 2
    integer :: parameterId                    ! GRIB 1
    integer :: level_type                     ! GRIB 1
    real    :: level_value                    ! GRIB 1
    integer :: silam_quantity                             ! SILAM
    type(silja_level) :: silam_level                      ! SILAM
    character(len=substNmLen) :: silam_species_io_str     ! SILAM
    type(silja_logical) :: defined
  end type Tgrib_code_table_entry

  type Tgrib_code_table
    private
    integer :: nEntries
    type(Tgrib_code_table_entry), dimension(:), pointer :: pCte
    type(silja_logical) :: defined
  end type Tgrib_code_table

  type(Tgrib_code_table), private, save :: grib_code_table_main




CONTAINS



  !************************************************************************************

  subroutine read_GRIB_code_table_v5(chFNm)
    !
    ! Reads all code tables available in the given file and fills-in the table structure
    ! Made for namelist-type file, which corresponds to code table v4 and further
    ! Also, for the first time, reads the species names
    !
    implicit none

    ! Imported parameter
    character(len=*), intent(in) :: chFNm

    ! Local variables
    type(Tsilam_namelist_group), pointer :: nlGrp
    type(Tsilam_namelist), pointer :: nlTmp
    integer :: iUnit, io_status, nEntries, iEntry, nItems, iItem
    character(len=clen) :: line
    character(len=20) :: chTmp
    logical :: eof
    type(Tgrib_code_table_entry), pointer :: pE
    type(silam_sp) :: sp

    !
    ! Open the file and get the namelist group
    !
    iUnit = fu_next_free_unit()
    open(iUnit, file=chFNm, status='old', iostat=io_status)
    if(io_status /= 0)then
      call set_error('Cannot open the file:' + chFNm,'read_GRIB_code_tables_v5')
      return
    endif

    nlGrp => fu_read_namelist_group(iUnit, .false., 'END_GRIB_CODE_TABLE_V5')
    if(error)return

    close(iUnit)

    nEntries = fu_nbr_of_namelists(nlGrp)
    if(error .or. nEntries < 1 .or. nEntries > 10000)then
      call set_error('Problems with the number of entries in the code table:' + fu_str(nEntries), &
                   & 'read_GRIB_code_tables_v5')
      return
    endif
    
    sp%sp => fu_work_string()
    if(error)return
    !
    ! Allocate the main table
    !
    allocate(grib_code_table_main%pCte(nEntries), stat=io_status)
    if(io_status /= 0)then
      call set_error('Failed to allocate the grib code table','read_GRIB_code_tables_v5')
      return
    endif
    grib_code_table_main%nEntries = nEntries
    !
    ! Go entry by entry putting them into the table
    !
    do iEntry = 1, nEntries
      
      nlTmp => fu_namelist(nlGrp, iEntry)
      pE => grib_code_table_main%pCte(iEntry)
      pE%centre(:) = int_missing
      pE%sub_centre(:) = int_missing
      pE%grib1_model_id(:) = int_missing
      !
      ! There can be several centres sharing the same definition
      !
      call get_int_array_from_items(nlTmp, 'centre', .false., pE%centre, 10)
!      call get_items(nlTmp, 'centre', pItems, nItems)
!      if(nItems > 10)then
!        call set_error('Too many centre items:'+fu_str(nItems),'read_GRIB_code_tables_v5')
!        return
!      endif
!      do iItem = 1, nItems
!        pE%centre(iItem) = fu_content_int(pItems(iItem))
!      end do
!      call destroy_items(pItems)
      !
      ! There can be several sub-centres sharing the same definition
      !
      call  get_int_array_from_items(nlTmp, 'sub_centre', .true., pE%sub_centre, 10)
!      call get_items(nlTmp, 'sub_centre', pItems, nItems)
!      if(nItems > 10)then
!        call set_error('Too many sub_centre items:'+fu_str(nItems),'read_GRIB_code_tables_v5')
!        return
!      elseif(nItems == 0)then
!          pE%sub_centre(1) = accept_all
!      elseif(nItems == 1)then
!        if(index(fu_content(pItems(1)),'*')>0)then
!          pE%sub_centre(1) = accept_all
!        else
!          pE%sub_centre(1) = fu_content_int(pItems(1))
!        endif
!      else
!        do iItem = 1, nItems
!          pE%sub_centre(iItem) = fu_content_int(pItems(iItem))
!        end do
!      endif
!      call destroy_items(pItems)
      !
      ! Same for models: more than one can use the definition
      !

      call get_int_array_from_items(nlTmp, 'grib1_model_identification', &
           & .true., pE%grib1_model_id, 10)
!      call get_items(nlTmp, 'grib1_model_identification', pItems, nItems)
!      if(nItems > 10)then
!        call set_error('Too many grib1_model_identification items:'+fu_str(nItems), &
!                     & 'read_GRIB_code_tables_v5')
!        return
!      endif
!      if(nItems == 1)then
!        if(index(fu_content(pItems(1)),'*')>0)then
!          pE%grib1_model_id(1) = accept_all
!        else
!          pE%grib1_model_id(1) = fu_content_int(pItems(1))
!        endif
!      else
!        do iItem = 1, nItems
!          pE%grib1_model_id(iItem) = fu_content_int(pItems(iItem))
!        end do
!      endif
!      call destroy_items(pItems)

!      pE%grib1_table_version = fu_content_int(nlTmp,'grib1_table_version')
      call get_int_array_from_items(nlTmp, 'grib1_table_version', &
           & .true., pE%grib1_table_version, 10)

      pE%discipline = fu_content_int(nlTmp,'discipline')
      pE%paramCategory = fu_content_int(nlTmp,'parameterCategory')
      pE%paramNumber = fu_content_int(nlTmp,'parameterNumber')
      pE%name = fu_content(nlTmp,'name')
      pE%short_name = fu_content(nlTmp,'short_name')
      pE%unit = fu_content(nlTmp,'unit')
      pE%typeOfFirstFixedSurface = fu_content_int(nlTmp,'typeOfFirstFixedSurface')
      if(pE%typeOfFirstFixedSurface == int_missing)pE%typeOfFirstFixedSurface = accept_all
      pE%ValueOfFirstFixedSurface = fu_content_real(nlTmp,'ValueOfFirstFixedSurface')
      if(pE%ValueOfFirstFixedSurface .eps. real_missing)pE%ValueOfFirstFixedSurface = accept_all
      pE%parameterId = fu_content_int(nlTmp,'parameterId')
      !
      ! SILAM stuff
      !
      ! There can be several levels meaning the same quantity
      !
      sp%sp = fu_content(nlTmp,'level_type')
      nItems = fu_nbrOfWords(sp%sp)
      do iItem = 1, nItems
        read(unit=sp%sp, fmt=*) chTmp
        if(chTmp == '*')then
          pE%level_type = accept_all
        else
          read(unit=chTmp, fmt=*) pE%level_type
        endif
      end do
      sp%sp = fu_content(nlTmp,'level_value')
      do iItem = 1, nItems
        read(unit=sp%sp, fmt=*) chTmp
        if(chTmp == '*')then
          pE%level_value = accept_all
        else
          read(unit=chTmp, fmt=*) pE%level_value
        endif
      end do
      !
      ! Finally, SILAM variable: quantity and species io string
      !
      pE%silam_quantity = fu_get_silam_quantity(fu_content(nlTmp,'silam_quantity_name'))
      if (error) then
         call report(nlTmp)
      endif
      if(fu_str_u_case(fu_content(nlTmp,'silam_level')) == 'UNDEFINED')then
        pE%silam_level = level_missing
      else
        call set_named_level_with_fract(fu_get_item(nlTmp,'silam_level'), pE%silam_level)
      endif
      pE%silam_species_io_str = fu_content(nlTmp,'silam_species_io_string')
      if(error)then
        pE%defined = silja_false
        return
      else
        pE%defined = silja_true
      endif
    end do  ! over the table entries

    if(error)then
      grib_code_table_main%defined = silja_false
    else
      grib_code_table_main%defined = silja_true
    endif

    call free_work_array(sp%sp)

    contains
    subroutine get_int_array_from_items(nlTmp, item_name, ifAcceptAllAllowed, arr, max_items)
      type(Tsilam_namelist), pointer :: nlTmp
      character(len=*) :: item_name
      logical, intent(in) :: ifAcceptAllAllowed
      integer, dimension(:), intent(out) :: arr
      integer, intent(in) :: max_items

      ! local
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems

      nullify(pItems)

      arr(1:max_items) = int_missing
      call get_items(nlTmp, item_name, pItems, nItems)
      if(nItems > max_items)then
        call set_error('Too many "'+item_name+'" items:'+fu_str(nItems),'read_GRIB_code_tables_v5')
      elseif(nItems == 0 .and. ifAcceptAllAllowed)then
          arr(1) = accept_all
      elseif(nItems == 1)then
        if(index(fu_content(pItems(1)),'*')>0 .and. ifAcceptAllAllowed)then
          arr(1) = accept_all
        else
          arr(1) = fu_content_int(pItems(1))
        endif
      else
        do iItem = 1, nItems
          arr(iItem) = fu_content_int(pItems(iItem))
        end do
      endif
      call destroy_items(pItems)
    end subroutine get_int_array_from_items


  end subroutine read_GRIB_code_table_v5


  !*************************************************************************

  subroutine detect_SILAM_level(chSILAMLevType,fLevVal, level)
    !
    ! Having the level type and level value, creates the SILAM level
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chSILAMLevType
    real, intent(in) :: fLevVal
    type(silja_level), intent(out) :: level

    if(chSILAMLevType == '*')then
      level = level_missing

    elseif(trim(chSILAMLevType) == 'surface')then
      level = surface_level

    elseif(trim(chSILAMLevType) == 'top_of_the_atmosphere')then
      level = top_atmosphere_level

    elseif(trim(chSILAMLevType) == 'constant_pressure')then
      if(fLevVal .eps. real_missing)then
        call set_error('No pressure value for pressure level','read_code_tables')
        return
      else
        level = fu_set_pressure_level(fLevVal)
      endif

    elseif(trim(chSILAMLevType) == 'mean_sea')then
      level = mean_sea_level

    elseif(trim(chSILAMLevType) == 'constant_altitude')then
      if(fLevVal .eps. real_missing)then
        call set_error('No altitude value for altitude level','read_code_tables')
        return
      else
        level = fu_set_constant_altitude_level(fLevVal)
      endif

    elseif(trim(chSILAMLevType) == 'constant_height')then
      if(fLevVal .eps. real_missing)then
        call set_error('No height value for height level','read_code_tables')
        return
      else
        level = fu_set_constant_height_level(fLevVal)
      endif

    elseif(trim(chSILAMLevType) == 'sigma_level')then
      if(fLevVal .eps. real_missing)then
        call set_error('No sigma value for sigma level','read_code_tables')
        return
      else
        level = fu_set_sigma_level(fLevVal)
      endif

    elseif(trim(chSILAMLevType) == 'entire_atmosphere_single_layer')then
      level = entire_atmosphere_mean_level

    else
      call set_error('Non-supported level type:' + chSILAMLevType, 'read_code_tables')
      return
    end if

  end subroutine detect_SILAM_level


  !***************************************************************

  subroutine get_silam_params_for_grib_1(tabVersion,tabCentre,tabSubCentre,tabModel, &
                                       & paramId, &
                                       & levType, levVal, &
                                       & silam_quantity, silam_species_string, silam_level)
    !
    ! Returns the SILAM quantity for the given set of GRIB parameters or
    ! int_missing if nothing is found. Since it is possible that some
    ! GRIB fields are not used in SILAM, no error is set if SILAM variable
    ! is not found.
    ! A complexity is that there can be "general-level" quantities and "specific-level"
    ! ones. For instance, temperature and temperature at 2m. The first one is 3D, the second 
    ! one is at specific level. In some cases that can be expressed via the same GRIB codes
    ! and only levels will differ.
    ! Algorithm:
    ! Each time check all entries ensuring that the entry with more coinsiding parameters
    ! prevails over the one with wider acceptance conditions.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: tabCentre,tabSubCentre, tabVersion,tabModel, paramId, levType
    real, intent(in) :: levVal
    integer, intent(out) :: silam_quantity
    type(silja_level), intent(out) :: silam_level
    character(len=*), intent(out) :: silam_species_string

    ! Local variables
    integer :: iEntry, iPrevEntry
    logical :: ifFullMatch
    !
    ! Search the whole database. As soon as some entry is found to be valid, store it
    ! and continue until the whole database is searched - may be, the more accurate match
    ! will be found.
    !
    iPrevEntry = -1

    do iEntry = 1, grib_code_table_main%nEntries

!      call check_entry_decode_grib_1(grib_code_table_main%pCte(iEntry), &
      call check_entry_alternate(iEntry, iPrevEntry, ifFullMatch)
      if(ifFullMatch)return

    end do  ! Cycle over the whole database
    !
    ! Nothing found
    !
    if(iPrevEntry < 1)then
!      call msg_warning('Failed GRIB-1 decoding for parameter nbr:' + fu_str(paramId),'get_GRIB_codes')
      silam_quantity = int_missing
      silam_level = level_missing
      silam_species_string = ''
    endif

    CONTAINS
    
    !================================================================================
    subroutine check_entry_alternate(iEntry, iPrevEntry, ifFullMatch)
      implicit none
      integer, intent(in) :: iEntry
      integer, intent(inout) :: iPrevEntry
      logical, intent(out) :: ifFullMatch

      type(Tgrib_code_table_entry), dimension(:), pointer :: pCte
      logical :: accept

      pCte => grib_code_table_main%pCte

      ! Doest this match
      ! No -> return
      ! Yes
      ! - > Is there a previous match
      !     No -> Accept so far
      !     Yes -> Is this better?
      !            Yes -> Accept and set as prev
      !            No ->  Keep the prev
      ifFullMatch = .false.
      if(any(tabVersion == pCte(iEntry)%grib1_table_version) .and. &
       & (pCTE(iEntry)%centre(1) == accept_all .or. any(pCTE(iEntry)%centre(:) == tabCentre)) .and. &
       & (pCTE(iEntry)%sub_centre(1) == accept_all .or. any(pCTE(iEntry)%sub_centre(:) == tabSubCentre)) .and. &
       & (pCTE(iEntry)%grib1_model_id(1) == accept_all .or. any(pCTE(iEntry)%grib1_model_id(:) == tabModel)) .and. &
       & (paramId == pCTE(iEntry)%parameterId))then     ! minimum match 
        
        if (pCTE(iEntry)%level_type == levType)then
          if(pCTE(iEntry)%level_value .eps. levVal)then
            ! A full match - all done.
            silam_quantity = pCTE(iEntry)%silam_quantity
            silam_level = pCTE(iEntry)%silam_level
            silam_species_string = pCTE(iEntry)%silam_species_io_str
            ifFullMatch = .true.                     ! No further search needed.
            iPrevEntry = iEntry

          else if (pCTE(iEntry)%level_value .eps. real(accept_all)) then
            ! A partial match, but better than nothing, or accepting all level types.
            accept = iPrevEntry < 1 ! the fun of not having short-circuit logical operators...
            if (.not. accept) accept = pCTE(iPrevEntry)%level_type == accept_all 
            if (accept) then
              silam_quantity = pCTE(iEntry)%silam_quantity
              silam_level = pCTE(iEntry)%silam_level
              silam_species_string = pCTE(iEntry)%silam_species_io_str
              iPrevEntry = iEntry
            end if
          end if
        end if
        
        if (pCTE(iEntry)%level_type == accept_all .and. iPrevEntry < 1) then
          ! A worst possible match, no need to check for being better than previous.
          silam_quantity = pCTE(iEntry)%silam_quantity
          silam_level = pCTE(iEntry)%silam_level
          silam_species_string = pCTE(iEntry)%silam_species_io_str
          iPrevEntry = iEntry
        end if
        
      end if
    end subroutine check_entry_alternate

!!$    subroutine check_entry_decode_grib_1(pCTE, pCTEPrev, ifFullMatch)
!!$      !
!!$      ! Checks whether the given entry is the right one for the given quantity and
!!$      ! returns the decoded SILAM parameter.
!!$      ! The problem: 
!!$      ! Several entries may suit for the decoding. In this case, we shall accept the
!!$      ! entry that is more restrictive, i.e. the entry, which accepts all levels is weaker
!!$      ! than the entry, which requires specific level. Then we will use the one with
!!$      ! specific level.
!!$      !
!!$      implicit none
!!$      type(Tgrib_code_table_entry), intent(in) :: pCte, pCTEPrev
!!$      logical, intent(out) :: ifFullMatch
!!$
!!$      !
!!$      ! First of all, check that the entry meets the GRIB definitions.
!!$      ! Do not miss the case of complete match of all input parameters.
!!$      !
!!$      ifFullMatch = .false.
!!$
!!$      if((tabVersion == pCTE%grib1_table_version) .and. &
!!$       & (pCTE%centre(1) == accept_all .or. any(pCTE%centre(:) == tabCentre)) .and. &
!!$       & (pCTE%grib1_model_id(1) == accept_all .or. any(pCTE%grib1_model_id(:) == tabModel)) .and. &
!!$       & (paramId == pCTE%parameterId))then     ! minimum match 
!!$        !
!!$        ! Minimum match reached. The entry can be used unless levels mismatch. 
!!$        ! The previous entry is either the same or at least matches at the same level.
!!$        ! Let's check whether the new entry is any better than the previous one
!!$        !
!!$        if(pCTE%level_type == levType)then  ! levType ==
!!$          if(pCTE%level_value .eps. levVal)then  ! level ==
!!$              silam_quantity = pCTE%silam_quantity
!!$              silam_level = pCTE%silam_level
!!$              silam_species_string = pCTE%silam_species_io_str
!!$              ifFullMatch = .true.                     ! No further search needed.
!!$          else
!!$            !
!!$            ! levVal /=. Previous entry can then be used if it matches or accepts all
!!$            !
!!$            if((pCTEPrev%level_value .eps. levVal) .or. &
!!$             & (pCTEPrev%level_value .eps. real(accept_all)))then
!!$              !
!!$              ! Previous entry (or the same one given second time) has the same or stronger fit - use it
!!$              !
!!$              if(.not. ifPrevEntryValid)then
!!$                !
!!$                ! Actualy, previous entry is the same, just given twice
!!$                !
!!$                silam_quantity = pCTE%silam_quantity
!!$                silam_level = pCTE%silam_level
!!$                silam_species_string = pCTE%silam_species_io_str
!!$                ifPrevEntryValid = .true.
!!$              endif
!!$            endif
!!$          endif  ! level value fit
!!$        else   
!!$          !
!!$          ! levType /=. Means also that level value is of no use. Previous entry can be used if matches.
!!$          !
!!$          if(pCTEPrev%level_type == accept_all)then
!!$            !
!!$            ! Previous entry has the same or stronger fit - use it
!!$            !
!!$            if(.not. ifPrevEntryValid)then
!!$              !
!!$              ! Actualy, previous entry is the same, just given twice
!!$              !
!!$              silam_quantity = pCTE%silam_quantity
!!$              silam_level = pCTE%silam_level
!!$              silam_species_string = pCTE%silam_species_io_str
!!$              ifPrevEntryValid = .true.
!!$            endif
!!$          elseif(pCTEPrev%level_type == levType)then  
!!$            !
!!$            ! level type of prev entry defined and matches. Value?
!!$            !
!!$            if((pCTEPrev%level_value .eps. levVal) .or. &
!!$             & (pCTEPrev%level_value .eps. real(accept_all)))then
!!$              !
!!$              ! Previous entry has the same or stronger fit - use it
!!$              !
!!$              if(.not. ifPrevEntryValid)then
!!$                !
!!$                ! Actualy, previous entry is the same, just given twice
!!$                !
!!$                silam_quantity = pCTE%silam_quantity
!!$                silam_level = pCTE%silam_level
!!$                silam_species_string = pCTE%silam_species_io_str
!!$                ifPrevEntryValid = .true.
!!$              endif
!!$            endif
!!$          endif  ! previous-entry level type
!!$        endif  ! level type fit
!!$      endif ! basic table parameters fit
!!$
!!$    end subroutine check_entry_decode_grib_1

  end subroutine get_silam_params_for_grib_1


  !***************************************************************

  subroutine get_silam_params_for_grib_2(discipline, paramCategory, paramNbr, &
                                       & levType, levVal, &
                                       & silam_quantity, silam_species_string, silam_level)
    !
    ! Returns the SILAM quantity for the given set of GRIB parameters or
    ! int_missing if nothing is found. Since it is possible that some
    ! GRIB fields are not used in SILAM, no error is set if SILAM variable
    ! is not found.
    ! A complexity is that there can be "general-level" quantities and "specific-level"
    ! ones. For instance, temperature and temperature at 2m. The first one is 3D, the second 
    ! one is at specific level. In some cases that can be expressed via the same GRIB codes
    ! and only levels will differ.
    ! Algorithm:
    ! Each time check all entries ensuring that the entry with more coinsiding parameters
    ! prevails over the one with wider acceptance conditions.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: discipline, paramCategory, paramNbr, levType
    real, intent(in) :: levVal
    integer, intent(out) :: silam_quantity
    type(silja_level), intent(out) :: silam_level
    character(len=*), intent(out) :: silam_species_string

    ! Local variables
    integer :: iEntry, iPrevEntry
    logical :: ifPrevEntryValid, ifFullMatch

    !
    ! Search the whole database. As soon as some entry is found to be valid, store it
    ! and continue until the whole database is searched - may be, the more accurate match
    ! will be found.
    !
    iPrevEntry = -1

    do iEntry = 1, grib_code_table_main%nEntries
!      call check_entry_decode_grib_2(grib_code_table_main%pCte(iEntry), &
      call check_entry_alternate(iEntry, iPrevEntry, ifFullMatch)
      if(ifFullMatch)return
    end do  ! Cycle over the whole database
    !
    ! Nothing found
    !
    if(iPrevEntry < 1)then
      if (debug_level > 1) call msg_warning('Failed GRIB-2 decoding for discipline:' + fu_str(discipline) + &
                                          & ', category:' + fu_str(paramCategory) + &
                                          & ', parameter nbr:' + fu_str(paramNbr),'get_GRIB_codes')
      silam_quantity = int_missing
      silam_level = level_missing
      silam_species_string = ''
    endif

  CONTAINS

    subroutine check_entry_alternate(iEntry, iPrevEntry, ifFullMatch)
      implicit none
      integer, intent(in) :: iEntry
      integer, intent(inout) :: iPrevEntry
      logical, intent(out) :: ifFullMatch
      
      type(Tgrib_code_table_entry), dimension(:), pointer :: pCte
      logical :: accept

      pCte => grib_code_table_main%pCte

      ! Doest this match
      ! No -> return
      ! Yes
      ! - > Is there a previous match
      !     No -> Accept so far
      !     Yes -> Is this better?
      !            Yes -> Accept and set as prev
      !            No ->  Keep the prev
      ifFullMatch = .false.
      if((discipline == pCTE(iEntry)%discipline) .and. &
       & (paramCategory == pCTE(iEntry)%paramCategory) .and. &
       & (paramNbr == pCTE(iEntry)%paramNumber))then     ! minimum match 
        
        if (pCte(iEntry)%typeOfFirstFixedSurface == levType)then
          if(pCte(iEntry)%ValueOfFirstFixedSurface .eps. levVal)then
            ! A full match - all done.
            silam_quantity = pCte(iEntry)%silam_quantity
            silam_level = pCte(iEntry)%silam_level
            silam_species_string = pCte(iEntry)%silam_species_io_str
            ifFullMatch = .true.                     ! No further search needed.
            iPrevEntry = iEntry

          else if (pCte(iEntry)%ValueOfFirstFixedSurface .eps. real(accept_all)) then
            ! A partial match, but better than nothing, or accepting all level types.
            accept = iPrevEntry < 1
            if (.not. accept) accept = pCte(iPrevEntry)%typeOfFirstFixedSurface == accept_all 
            if (accept) then
              silam_quantity = pCte(iEntry)%silam_quantity
              silam_level = pCte(iEntry)%silam_level
              silam_species_string = pCte(iEntry)%silam_species_io_str
              iPrevEntry = iEntry
            end if
          end if

        end if
        
        if (pCte(iEntry)%typeOfFirstFixedSurface == accept_all .and. iPrevEntry < 1) then
          ! A worst possible match, no need to check for being better than previous.
          silam_quantity = pCte(iEntry)%silam_quantity
          silam_level = pCte(iEntry)%silam_level
          silam_species_string = pCte(iEntry)%silam_species_io_str
          iPrevEntry = iEntry
        end if
        
      end if
    end subroutine check_entry_alternate

    subroutine check_entry_decode_grib_2(pCTE, pCTEPrev, ifFullMatch)
      !
      ! Checks whether the given entry is the right one for the given quantity and
      ! returns the GRIB-2 encoding
      !
      implicit none
      type(Tgrib_code_table_entry), intent(in) :: pCte, pCTEPrev
      logical, intent(out) :: ifFullMatch
      !
      ! Start from the minimum match check
      !
      if(discipline == pCTE%discipline .and. &
       & paramCategory == pCTE%paramCategory .and. &
       & paramNbr == pCTE%paramNumber)then     ! minimum match 
        !
        ! Minimum match reached. The entry can be used unless levels mismatch. 
        ! The previous entry is either the same or at least matches at the same level.
        ! Let's check whether the new entry is any better than the previous one
        !
        if(pCTE%typeOfFirstFixedSurface == levType)then  ! levType ==
          if(pCTE%ValueOfFirstFixedSurface .eps. levVal)then  ! level ==
              silam_quantity = pCTE%silam_quantity
              silam_level = pCTE%silam_level
              silam_species_string = pCTE%silam_species_io_str
              ifFullMatch = .true.                     ! No further search needed.
          else
            !
            ! levVal /=. Previous entry can then be used if it matches or accepts all
            !
            if((pCTEPrev%ValueOfFirstFixedSurface .eps. levVal) .or. &
             & (pCTEPrev%ValueOfFirstFixedSurface .eps. real(accept_all)))then
              !
              ! Previous entry has the same or stronger fit - use it
              !
              if(.not. ifPrevEntryValid)then
                !
                ! Actualy, previous entry is the same, just given twice
                !
                silam_quantity = pCTE%silam_quantity
                silam_level = pCTE%silam_level
                silam_species_string = pCTE%silam_species_io_str
                ifPrevEntryValid = .true.
              endif
            endif
          endif  ! level value fit
        else   
          !
          ! levType /=. Means also that level value is of no use. Previous entry can be used if matches.
          !
          if(pCTEPrev%typeOfFirstFixedSurface == accept_all)then
            !
            ! Previous entry has the same or stronger fit - use it
            !
            if(.not. ifPrevEntryValid)then
              !
              ! Actualy, previous entry is the same, just given twice
              !
              silam_quantity = pCTE%silam_quantity
              silam_level = pCTE%silam_level
              silam_species_string = pCTE%silam_species_io_str
              ifPrevEntryValid = .true.
            endif
          elseif(pCTEPrev%typeOfFirstFixedSurface == levType)then  
            !
            ! level type of prev entry defined and matches. Value?
            !
            if((pCTEPrev%ValueOfFirstFixedSurface .eps. levVal) .or. &
             & (pCTEPrev%ValueOfFirstFixedSurface .eps. real(accept_all)))then
              !
              ! Previous entry has the same or stronger fit - use it
              !
              if(.not. ifPrevEntryValid)then
                !
                ! Actualy, previous entry is the same, just given twice
                !
                silam_quantity = pCTE%silam_quantity
                silam_level = pCTE%silam_level
                silam_species_string = pCTE%silam_species_io_str
                ifPrevEntryValid = .true.
              endif
            endif
          endif  ! previous-entry level type
        endif  ! level type fit

      endif   ! Basic parameters match

    end subroutine check_entry_decode_grib_2

  end subroutine get_silam_params_for_grib_2


  ! ***************************************************************

  subroutine get_GRIB_codes(SILAM_quantity, SILAM_species_io_str, &  ! SILAM variable  IN
                          & iGribType, &                             ! 1 or 2      IN
                          & tabCentre, tabSubCentre, tabModel, &     ! universal   IN
                          & tabVersion, &                            ! for GRIB 1  INOUT
                          & iParamId, &                                  ! GRIB 1    OUT
                          & iDiscipline, iParamCategory, iParamNbr, &    ! GRIB 2    OUT
                          & iLevelType_required, fLevelValueRequired) ! univ, OUT
    !
    ! Returns the GRIB code of the SILAM quantity from the code_table. 
    ! If the SILAM quantity does not exist the integer_missing is returned.
    ! Distinction is made between the GRIB 1 and GRIB 2. Some quantities require
    ! level for being understandable in specific table version. Should this happen,
    ! this level will be returned, otherwise it will be level_missing
    !
    IMPLICIT NONE

    ! Imported parameters:
    INTEGER, INTENT(IN) :: SILAM_quantity, iGribType, tabCentre, tabSubCentre, tabModel
    integer, intent(inout) :: tabVersion
    character(len=*), intent(in) :: SILAM_species_io_str
    integer, intent(out) :: iParamId, iDiscipline, iParamCategory, iParamNbr, iLevelType_required
    real, intent(out) :: fLevelValueRequired

    ! Local variables:
    INTEGER :: iEntry
    integer, save :: iPrevEntry=0  ! Previous-time entry

    !
    ! First, check the previous-time entry. May be, it is just what we need. Can be a good 
    ! speedup for multi-level fields. If not, search the whole list looking for the exact fit.
    ! Note that the rules for == check are different for GRIB-1 and GRIB-2
    !
    ! Is the previous-time entry OK?    
    !
    if(iPrevEntry > 0)then
      if(iGribType == 1)then
        if(fu_if_entry_OK_encode_grib_1(grib_code_table_main%pCte(iPrevEntry))) return
      elseif(iGribType == 2)then
        if(fu_if_entry_OK_encode_grib_2(grib_code_table_main%pCte(iPrevEntry)))return
      else
        call set_error('E.1 Unknown GRIB type:' + fu_str(iGribType),'get_GRIB_codes')
        return
      endif
    endif  ! if previous entry index > 0

    !
    ! If not previous-time entry, search the whole database
    !
    do iEntry = 1, grib_code_table_main%nEntries
      if(iGribType == 1)then
        if(fu_if_entry_OK_encode_grib_1(grib_code_table_main%pCte(iEntry)))then
          iPrevEntry = iEntry
          return
        endif
      elseif(iGribType == 2)then
        if(fu_if_entry_OK_encode_grib_2(grib_code_table_main%pCte(iEntry)))then
          iPrevEntry = iEntry
          return
        endif
      else
        call set_error('E.2 Unknown GRIB type:' + fu_str(iGribType),'get_GRIB_codes')
        return
      endif    ! Grib 1 or 2
    end do  ! Cycle over the whole database
    !
    ! Nothing found
    !
    call set_error('Failed GRIB coding for:'+fu_quantity_string(SILAM_quantity) + ',' + &
                 & SILAM_species_io_str,'get_GRIB_codes')
    iParamId = int_missing
    iDiscipline = int_missing
    iParamCategory = int_missing
    iParamNbr = int_missing
    iLevelType_required = int_missing
    fLevelValueRequired = real_missing


    CONTAINS
    
    !================================================================================

    logical function fu_if_entry_OK_encode_grib_1(pCTE)
      !
      ! Checks whether the given entry is the right one for the given quantity and
      ! returns the GRIB-1 encoding
      !
      implicit none
      type(Tgrib_code_table_entry), intent(in) :: pCte

      if(SILAM_quantity == pCTE%silam_quantity                                     .and. &
       & (tabVersion == accept_all .or. pCTE%grib1_table_version(1) == accept_all .or. & 
                                     & any(pCTE%grib1_table_version(:) == tabVersion))  .and. &
       & (tabCentre == accept_all .or. pCTE%centre(1) == accept_all .or. & 
                                     & any(pCTE%centre(:) == tabCentre))           .and. &
       & (tabSubCentre == accept_all .or. pCTE%sub_centre(1) == accept_all .or. & 
                                     & any(pCTE%sub_centre(:) == tabSubCentre))    .and. &
       & (tabModel == accept_all .or. pCTE%grib1_model_id(1) == accept_all .or. & 
                                    & any(pCTE%grib1_model_id(:) == tabModel))     .and. &
       & fu_str_l_case(SILAM_species_io_str) == fu_str_l_case(pCTE%silam_species_io_str))then
          iParamId = pCTE%parameterId
          tabVersion = pCTE%grib1_table_version(1)
          iLevelType_required = pCTE%level_type
          fLevelValueRequired = pCTE%level_value
          fu_if_entry_OK_encode_grib_1 = .true.
          return
      endif

      fu_if_entry_OK_encode_grib_1 = .false.

    end function fu_if_entry_OK_encode_grib_1

    !===============================================================================

    logical function fu_if_entry_OK_encode_grib_2(pCTE)
      !
      ! Checks whether the given entry is the right one for the given quantity and
      ! returns the GRIB-2 encoding
      !
      implicit none
      type(Tgrib_code_table_entry), intent(in) :: pCte

      if(SILAM_quantity == pCTE%silam_quantity                                     .and. &
       & (tabCentre == accept_all .or. pCTE%centre(1) == accept_all .or. & 
                                     & any(pCTE%centre(:) == tabCentre))           .and. &
       & (tabSubCentre == accept_all .or. pCTE%sub_centre(1) == accept_all .or. & 
                                     & any(pCTE%sub_centre(:) == tabSubCentre))    .and. &
       & (tabModel == accept_all .or. pCTE%grib1_model_id(1) == accept_all .or. & 
                                    & any(pCTE%grib1_model_id(:) == tabModel))     .and. &
       & trim(SILAM_species_io_str) == trim(pCTE%silam_species_io_str))then
          iDiscipline = pCTE%discipline
          iParamCategory = pCTE%paramCategory
          iParamNbr = pCTE%paramNumber
          tabVersion = pCTE%grib1_table_version(1)
          iLevelType_required = pCTE%typeOfFirstFixedSurface
          fLevelValueRequired = pCTE%ValueOfFirstFixedSurface
          fu_if_entry_OK_encode_grib_2 = .true.
          return
      endif

      fu_if_entry_OK_encode_grib_2 = .false.

    end function fu_if_entry_OK_encode_grib_2

  END subroutine get_GRIB_codes


  !*****************************************************************************
  !
  !  PRIVATE functions of this module
  !
  !*****************************************************************************

  logical function CTE_defined(CTE)
    implicit none
    type(Tgrib_code_table_entry), intent(in) :: CTE
    CTE_defined = (CTE%defined == silja_true)
  end function CTE_defined

  !****************************************************************

  subroutine report_code_table()
    !
    ! Prints all the content read for the code table
    !
    implicit none

    ! Local variables
    integer :: iEntry, iTmp
    type(silam_sp) :: sp

    if(grib_code_table_main%nEntries < 1)then
      call msg('GRIB code table has no entries')
      return
    endif

    sp%sp => fu_work_string()
    if(error)return

    call msg('**************** GRIB CODE TABLE content *******************')

    do iEntry = 1, grib_code_table_main%nEntries
      if(.not. CTE_defined(grib_code_table_main%pCTE(iEntry)))exit

      do iTmp = 1, size(grib_code_table_main%pCTE(iEntry)%centre)
        if(grib_code_table_main%pCTE(iEntry)%centre(iTmp) == int_missing)exit
        call msg('centre = ',grib_code_table_main%pCTE(iEntry)%centre(iTmp))
      end do

      do iTmp = 1, size(grib_code_table_main%pCTE(iEntry)%sub_centre)
        if(grib_code_table_main%pCTE(iEntry)%sub_centre(iTmp) == int_missing)exit
        call msg('sub_centre = ',grib_code_table_main%pCTE(iEntry)%sub_centre(iTmp))
      end do

      do iTmp = 1, size(grib_code_table_main%pCTE(iEntry)%grib1_table_version)
        if(grib_code_table_main%pCTE(iEntry)%grib1_table_version(iTmp) == int_missing)exit
        call msg('grib1_table_version = ',grib_code_table_main%pCTE(iEntry)%grib1_table_version(iTmp))
      end do

      do iTmp = 1, size(grib_code_table_main%pCTE(iEntry)%grib1_model_id)
        if(grib_code_table_main%pCTE(iEntry)%grib1_model_id(iTmp) == int_missing)exit
        call msg('grib1_model_identification = ', &
                                     & grib_code_table_main%pCTE(iEntry)%grib1_model_id(iTmp))
      end do

      call msg('discipline = ', grib_code_table_main%pCTE(iEntry)%discipline)
      call msg('parameterCategory = ', grib_code_table_main%pCTE(iEntry)%paramCategory)
      call msg('parameterNumber = ', grib_code_table_main%pCTE(iEntry)%paramNumber)
      call msg('name = ' // trim(grib_code_table_main%pCTE(iEntry)%name))
      call msg('short_name = ' // trim(grib_code_table_main%pCTE(iEntry)%short_name))
      call msg('unit = ' // trim(grib_code_table_main%pCTE(iEntry)%unit))
      call msg('typeOfFirstFixedSurface = ', &
                                & grib_code_table_main%pCTE(iEntry)%typeOfFirstFixedSurface)
      call msg('ValueOfFirstFixedSurface = ', &
                                & grib_code_table_main%pCTE(iEntry)%ValueOfFirstFixedSurface)
      call msg('parameterId = ', grib_code_table_main%pCTE(iEntry)%parameterId)
      call msg('level_type = ', grib_code_table_main%pCTE(iEntry)%level_type)
      call msg('level_value = ', grib_code_table_main%pCTE(iEntry)%level_value)
      call msg('silam_quantity_name = ' // &
                         & trim(fu_quantity_string(grib_code_table_main%pCTE(iEntry)%silam_quantity)))
      call level_to_short_string(grib_code_table_main%pCTE(iEntry)%silam_level, sp%sp)
      call msg('silam_level = ' // trim(sp%sp))
      call msg('silam_species_io_string = ' // &
                         & trim(grib_code_table_main%pCTE(iEntry)%silam_species_io_str))

    end do  ! table entries

    call msg('************* End of GRIB CODE TABLE content ***************')

    call free_work_array(sp%sp)

  end subroutine report_code_table


END MODULE grib_code_table
