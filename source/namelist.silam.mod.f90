MODULE silam_namelist
  !
  ! Current module supports the Tsilam_namelist class, which has a few nice features
  ! compare to usual FORTRAN namelist. In particular, there is no need to define
  ! the list in advance - it is a completely dynamic object.
  ! Second, the list is, in fact, an array of namelists, each with own name. So, 
  ! there is one more level of hierarchy.
  ! Third, the content of one item can be whatever, including several words
  ! separated by space.
  ! A payment for universality is - Tsilam_namelist has only character content,
  ! so its actual meaning has to be handled by the object calling the namelist.
  !
  ! Namelist item is, by definition, anything with non-empty title and some, 
  ! possibly, empty, content.
  !
  ! Namelist may have several items with the same title - just use appropriate
  ! functions for handling. A single-content function picks the first item with
  ! given name adn ignores the others. 
  !
  ! Language: ANSI FORTRAN-90
  !
  ! Code owner Mikhail Sofiev
  !
  use toolbox

  implicit none

  public fu_create_namelist_group
  public fu_read_namelist_group  ! creates the namelist group
  public destroy_namelist_group  ! eliminates the namelist group and frees the memory
  public reset_namelist_group  ! cleans the namelist group , so it can be further used
  public write_namelist_group  ! output in the form of namelist group
  public fu_add_namelist   ! adds namelist and returns a pointer to it
  public add_namelist_item !allocates one more namelist item and fills it, if needed
  public replace_namelist_item ! replaces the requested item with new name and content
  public fu_read_namelist  ! creates the namelist
  public destroy_namelist  ! eliminates the namelist and frees the memory
  public reset_namelist    ! cleans the namelist, so it can be further used
  public fu_make_missing_namelist  ! Instead of namelist_missing
  public write_namelist      ! output in the form of namelist
  public write_namelist_item
  public fu_create_namelist  ! create an empty namelist
  public next_named_line_from_file  ! Searches for the next named line 
  public fu_namelist       ! get namelist from the group
  public get_items         ! get all items with given name
  public destroy_items     ! and destroy them
  public fu_get_item
  public fu_name           ! name of the namelist item
  public fu_content        ! Content of the namelist item
  public fu_content_real   ! Content of the namelist item, must be real number
  public fu_content_int    ! Content of the namelist item, must be integer number
  public content_real_with_unit ! A value with a given unit
  public content_int_with_unit ! A value with a given unit
  public clean  ! items, namelist, group
  public report
  public defined
  public empty
  public fu_nbr_of_items ! in the namelist or nlGroup
  public fu_nbr_of_namelists  ! in the namelist group
  
  !
  ! Namelist operations within nl-goup
  !
  private fu_add_new_namelist
  private fu_add_existing_namelist
  private increase_ptr_array_nl_group

  !
  ! Extracting all items with given names
  !
  private items_from_namelist ! all items with given name
  private items_from_nl_group ! all items with given namelist and item names
  private destroy_namelist_item_list  ! and destroys them
  private add_existing_namelist_item
  private add_namelist_item_basic
  private replace_namelist_item_by_name
  private replace_namelist_item_content
  !
  ! Getting the content 
  !
  private fu_item_name
  private fu_itemPtr_name
  private fu_item_content
  private fu_itemPtr_content
  private fu_item_content_real
  private fu_itemPtr_content_real
  private fu_item_content_int
  private fu_itemPtr_content_int
  private print_nlItem
  private print_nlItemPtr
  private fu_content_from_namelist
  private fu_content_from_namelist_real
  private fu_content_from_namelist_int
  private fu_content_from_nl_group
  private fu_content_from_nl_group_real
  private fu_content_from_nl_group_int
  private report_namelist
  private write_namelist_item_char
  private write_namelist_item_int
  private write_namelist_item_real
  private fu_namelist_defined
  private fu_namelist_empty
  private fu_namelist_group_empty
  private fu_namelist_by_name
  private fu_namelist_by_index
  private fu_nbr_of_items_of_namelist
  private fu_nbr_of_items_of_nl_group
  private fu_get_item_by_index
  private fu_get_item_by_name
  private clean_item    ! Cleaning without destroying memory
  private clean_itemPtr
  private clean_namelist
  private clean_namelist_group
  private item_cnt_real_with_unit
  private item_cnt_int_with_unit
  private itemptr_cnt_real_with_unit
  private itemptr_cnt_int_with_unit
  private cnt_nl_real_with_unit
  private cnt_nl_int_with_unit
  private cnt_nl_group_real_with_unit
  private cnt_nl_grp_int_with_unit

  interface  get_items
    module procedure items_from_namelist ! all items with given name
    module procedure items_from_nl_group ! all items with given namelist and item names
  end interface

  interface destroy_items
    module procedure destroy_namelist_item_list
  end interface

  interface fu_get_item
    module procedure fu_get_item_by_index
    module procedure fu_get_item_by_name
  end interface

  interface add_namelist_item
    module procedure add_existing_namelist_item
    module procedure add_namelist_item_basic
  end interface

  interface replace_namelist_item
    module procedure replace_namelist_item_by_name
    module procedure replace_namelist_item_content
  end interface

  interface fu_name
    module procedure fu_namelist_name
    module procedure fu_item_name
    module procedure fu_itemPtr_name
  end interface

  interface fu_content
    module procedure fu_item_content    ! from item 
    module procedure fu_itemPtr_content    ! from item 
    module procedure fu_content_from_namelist  ! from namelist, having item name
    module procedure fu_content_from_nl_group  ! from nl-group, having namelist and item names
  end interface

  interface fu_content_real
    module procedure fu_item_content_real
    module procedure fu_itemPtr_content_real
    module procedure fu_content_from_namelist_real
    module procedure fu_content_from_nl_group_real
  end interface

  interface fu_content_int
    module procedure fu_item_content_int
    module procedure fu_itemPtr_content_int
    module procedure fu_content_from_namelist_int
    module procedure fu_content_from_nl_group_int
  end interface

  interface content_int_with_unit
    module procedure item_cnt_int_with_unit
    module procedure itemPtr_cnt_int_with_unit
    module procedure cnt_nl_int_with_unit
    module procedure cnt_nl_grp_int_with_unit
  end interface

  interface content_real_with_unit
    module procedure item_cnt_real_with_unit
    module procedure itemPtr_cnt_real_with_unit
    module procedure cnt_nl_real_with_unit
    module procedure cnt_nl_group_real_with_unit
  end interface

  interface fu_add_namelist
    module procedure fu_add_new_namelist
    module procedure fu_add_existing_namelist
  end interface

  interface clean 
    module procedure clean_item    ! Cleaning without destroying memory
    module procedure clean_itemPtr
    module procedure clean_namelist
    module procedure clean_namelist_group
  end interface

  interface report
    module procedure report_namelist
    module procedure print_nlItem
    module procedure print_nlItemPtr
  end interface

  interface write_namelist_item
    module procedure write_namelist_item_char
    module procedure write_namelist_item_int
    module procedure write_namelist_item_real
  end interface

  interface fu_namelist
    module procedure fu_namelist_by_name
    module procedure fu_namelist_by_index
  end interface

  interface empty
    module procedure fu_namelist_empty
    module procedure fu_namelist_group_empty
  end interface

  interface defined
    module procedure fu_namelist_defined
  end interface

  interface fu_nbr_of_items
    module procedure fu_nbr_of_items_of_namelist
    module procedure fu_nbr_of_items_of_nl_group
  end interface

  !
  ! Namelist item type
  !
  type Tsilam_namelist_item
    private
    character(len=clen) :: chName
!    character(len=180) :: chContent
    character(len=fnlen) :: chContent
  end type Tsilam_namelist_item

  type Tsilam_nl_item_ptr
    private
    type(Tsilam_namelist_item) :: nli
  end type Tsilam_nl_item_ptr

  !
  ! Namelist itself
  !
  type Tsilam_namelist
    private
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    integer :: iLastFilledItem, iTotalSize, nAllocated
    type(silja_logical) :: defined
    character(len=clen) :: chName
  end type Tsilam_namelist

  type Tsilam_namelist_ptr
!    private
    type(Tsilam_namelist), pointer :: nl
  end type Tsilam_namelist_ptr

  !
  ! Next level of aggregation - namelist group
  !
  type Tsilam_namelist_group
    private
    type(Tsilam_namelist_ptr), dimension(:), pointer :: lists
    integer :: iLastFilledNL, iTotalSize, nAllocated
    character(len=clen) :: chName
    type(silja_logical) :: defined
  end type Tsilam_namelist_group

!  !
!  ! A global namelist describing certain internal multi-run options
!  !
!  type(Tsilam_namelist), public, pointer, save :: multi_run_namelist


CONTAINS

  !*****************************************************************

  function fu_create_namelist_group(chNm) result(nlGrpPtr)
    !
    ! Creates an empty new namelist and returns a pointer on it
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in), optional :: chNm

    ! Return value
    type(Tsilam_namelist_group), pointer :: nlGrpPtr

    ! Internal variables
    integer :: iStat

    allocate(nlGrpPtr, stat = iStat)
    if(iStat /= 0)then
      call set_error('Failed to allocate namelist','fu_create_namelist_group')
      nullify(nlGrpPtr)
      return
    endif

    nullify(nlGrpPtr%lists)
    nlGrpPtr%iLastFilledNL = 0
    nlGrpPtr%iTotalSize = 0
    nlGrpPtr%nAllocated = 0
    if(present(chNm))then
      nlGrpPtr%chName = trim(chNm)
    else
      nlGrpPtr%chName = ''
    endif
    nlGrpPtr%defined = silja_true

  end function fu_create_namelist_group


  !******************************************************************************

  function fu_read_namelist_group (iUnit, ifEnv, chStopLine) result(nlGroup)
    !
    ! Creates the namelists one-by-one and reads them from the file pointed by iUnit
    ! The reading is stopped at the chStopLine, which is also 
    ! consumed from the file. If this line has a namelist format - it is added
    ! to the namelist. If not - skipped
    ! In general, all named lines are included, all non-named (without '=' sign)
    ! are skipped
    !
    implicit none

    ! Return value
    type(Tsilam_namelist_group), pointer :: nlGroup

    ! Imported parameters
    integer, intent(in) :: iUnit
    character(len=*), intent(in), optional :: chStopLine
    logical, intent(in) :: ifEnv

    ! Local variables
    integer :: i !, iList
    logical :: eof
    character(len=fnlen) :: line
!    type(Tsilam_namelist_group), pointer :: nlPtr
    type(Tsilam_namelist), pointer :: nlPtr

    ! Stupidity check
    !
    if(present(chStopLine))then
      if(len_trim(chStopLine) == 0)then
        call set_error('Empty stop line','fu_read_namelist_group')
        return
      endif
    endif

    !
    ! Create an empty namelist group
    !
    allocate(nlGroup, stat=i)
    if(i /= 0)then
      call msg_test('Strange allocation status: ',i)
      call set_error('Failed to allocate memory for namelist','fu_read_namelist_group')
      return
    endif
    nlGroup%iTotalSize=0
    nlGroup%iLastFilledNL=0
    nlGroup%nAllocated=0
    nlGroup%chName = ''
!    iList = 0
    nullify(nlPtr)
!    nlPtr => nlGroup
    
    !
    ! Simply read the file through until eof or stop line
    !
    eof=.false.
    do while(.not.(eof.or.error))
      call next_line_from_input_file(iUnit, line, eof)
      if(error)return
      if(eof)then
        if(present(chStopLine)) then 
          call msg_warning('EOF reached prior to stopline', 'fu_read_namelist_group')
          call msg("Stop line was '"//trim(chStopLine)//"'")
        endif
        if(nlGroup%iLastFilledNL > 0) then
          nlGroup%defined = silja_true
        else
          nlGroup%defined = silja_false
        endif
        return
      endif
      !
      ! Namelist group is over ?
      !
      if(present(chStopLine))then
        if(trim(line) == trim(chStopLine)) then
          nlGroup%defined = silja_true
          return
        endif
      endif

      !
      ! Split the line
      !
      i=index(line,'=')
      if(i == 0)then
        ! Stop line is usually not a named line either.
        if(index(line, trim(chStopLine)) == 0 .and. index(chStopLine,trim(line)) == 0) & 
                          & call msg_warning('Not a named line:' + line, 'fu_read_namelist_group')
      else
        !
        ! Named line - store it somewhere
        !
!        if(iList == 0)then 
        if(.not. associated(nlPtr))then 

          if(trim(fu_str_u_case(line(1:i-1))) == 'LIST_GROUP')then
            !
            ! This is just a name of list group , store it
            !
            nlGroup%chName = fu_str_u_case(adjustl(line(i+1:len_trim(line))))

          elseif(trim(fu_str_u_case(line(1:i-1))) == 'END_LIST_GROUP')then
            !
            ! Namelist group has ended
            !
            nlGroup%defined = silja_true
            return

          else
            !
            ! Pointer to the namelist is undefined - start a new namelist
            !
            nlPtr => fu_add_namelist(nlGroup) ! nlPtr now points to the new list
            if(error)return

            if(trim(fu_str_u_case(line(1:i-1))) == 'LIST')then
              !
              ! New list has a name, store it
              !
              nlPtr%chName = fu_str_u_case(adjustl(line(i+1:len_trim(line))))

            else
              !
              ! Some ordinary line => namelist does not have any name, just fill it
              !
              nlPtr%chName = ''
              nlPtr%items(nlPtr%iLastFilledItem+1)%nli%chName = adjustl(line(1:i-1))
              if (ifEnv) then
                 nlPtr%items(nlPtr%iLastFilledItem+1)%nli%chContent = adjustl(fu_expand_environment(line(i+1:len_trim(line))))
              else
                 nlPtr%items(nlPtr%iLastFilledItem+1)%nli%chContent = adjustl(line(i+1:len_trim(line)))
              endif
              nlPtr%iLastFilledItem = nlPtr%iLastFilledItem + 1
              nlPtr%defined = silja_true
              call add_namelist_item(nlPtr) ! More space needed ?
              if(error)return
!              nlGroup%iLastFilledNL = nlGroup%iLastFilledNL + 1

            endif
          endif

        else
          !
          ! iList pointer is OK => just continue the namelist filling
          !
          if(fu_str_u_case(adjustl(line(1:i-1))) == 'END_LIST')then
            !
            ! List is over => close it and set the pointer to undefined
            !
            nullify(nlPtr)   !iList = 0

          else
            !
            ! Ordinary line => add it to the list
            !
            nlPtr%items(nlPtr%iLastFilledItem+1)%nli%chName = adjustl(line(1:i-1))
            if (ifEnv) then
               nlPtr%items(nlPtr%iLastFilledItem+1)%nli%chContent = adjustl(fu_expand_environment(line(i+1:len_trim(line))))
            else
               nlPtr%items(nlPtr%iLastFilledItem+1)%nli%chContent = adjustl(line(i+1:len_trim(line)))
            endif
            nlPtr%iLastFilledItem = nlPtr%iLastFilledItem + 1
            call add_namelist_item(nlPtr) ! More space needed ?
            if(error)return
!            nlGroup%iLastFilledNL = nlGroup%iLastFilledNL + 1

          endif ! if line ends the namelist

        endif ! iList==0

      endif  ! if a line is named line

    end do  ! through the input file

    nlGroup%defined = silja_true

  end function fu_read_namelist_group


  !*******************************************************************************

  subroutine destroy_namelist_group(nlGrp)
    !
    ! Frees the memory allocated for the namelist nl
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), pointer :: nlGrp

    ! Local variables
    integer :: i

    if(associated(nlGrp))then
      if(nlGrp%defined == silja_true)then
        do i=1, nlGrp%iTotalSize
          call destroy_namelist(nlGrp%lists(i)%nl)
        end do
      end if
      deallocate(nlGrp)
      nullify(nlGrp)
    endif

  end subroutine destroy_namelist_group


  !********************************************************************************

  subroutine reset_namelist_group(nlGrp)
    !
    ! cleans the namelist, so it can be further used
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), pointer :: nlGrp

    ! Local variables
    integer :: i

    if(associated(nlGrp))then
      if(nlGrp%defined == silja_true)then
        do i=1,nlGrp%iLastFilledNL
          call reset_namelist(nlGrp%lists(i)%nl)
        end do
        nlGrp%defined = silja_false
      end if
    endif

  end subroutine reset_namelist_group


  !*********************************************************************************

  subroutine write_namelist_group(iUnit, nlGrp, ifPublic) 
    !
    ! output in the form of namelist
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), pointer :: nlGrp
    integer, intent(in) :: iUnit
    logical, intent(in), optional :: ifPublic

    ! Local variables
    integer :: i

    if(associated(nlGrp))then
      if(nlGrp%defined == silja_true)then
        do i=1, nlGrp%iLastFilledNL
          if(present(ifPublic))then
              call write_namelist(iUnit, nlGrp%lists(i)%nl, .true., ifPublic) ! header public
          else
            call write_namelist(iUnit, nlGrp%lists(i)%nl, .true., .false.)
          endif
        end do
      else
        call msg('Undefined namelist group')
      end if
    endif
  end subroutine write_namelist_group


  !*****************************************************************************

  function fu_namelist_by_name(nlGrp, chNm)
    !
    ! Searches the namelist in the group and returns the Tsilam_namelist
    ! pointer
    !
    implicit none

    ! Return value
    type(Tsilam_namelist), pointer :: fu_namelist_by_name

    ! Imported parameters
    type(Tsilam_namelist_group), pointer :: nlGrp
    character(len=*), intent(in) :: chNm

    ! Local variables
    integer :: i

!call msg('Start selecting. Name:' + chNm)
!call msg('Total namelists:',nlGrp%iLastFilledNL)
!call write_namelist_group(run_log_funit, nlGrp, .true.)

    do i=1,nlGrp%iLastFilledNL

!call msg('i',i)
!call msg('nlGrp%lists(i)%nl%chName' + nlGrp%lists(i)%nl%chName)

      if(fu_str_u_case(nlGrp%lists(i)%nl%chName) == fu_str_u_case(chNm))then
        fu_namelist_by_name => nlGrp%lists(i)%nl
        return
      end if
    end do
    nullify(fu_namelist_by_name)

  end function fu_namelist_by_name


  !*****************************************************************************

  function fu_namelist_by_index(nlGrp, indexNl)
    !
    ! Searches the namelist in the group and returns the Tsilam_namelist
    ! pointer
    !
    implicit none

    ! Return value
    type(Tsilam_namelist), pointer :: fu_namelist_by_index

    ! Imported parameters
    type(Tsilam_namelist_group), intent(in) :: nlGrp
    integer, intent(in) :: indexNl

    ! Local variables
    integer :: i

    if(indexNl < 1 .or. indexNl > nlGrp%iLastFilledNL)then
      call msg('Namelist index is outside the limits:', indexNl)
      call msg('Number of existing namelists:',nlGrp%iLastFilledNL)
      call set_error('Namelist index is outside the limits','fu_namelist_by_index')
      nullify(fu_namelist_by_index)
    else
      fu_namelist_by_index => nlGrp%lists(indexNl)%nl
    end if

  end function fu_namelist_by_index


  !******************************************************************

  logical function fu_namelist_group_empty(nlGrp) 
    !
    ! Checks whether the namelist is reasonable
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist_group), intent(in) :: nlGrp

    ! Local variables 
    integer :: i

    if(nlGrp%iLastFilledNl <= 0)then
      fu_namelist_group_empty = .true.
    else
      do i=1, nlGrp%iLastFilledNl
        if(.not. fu_namelist_empty(nlGrp%lists(i)%nl)) then
          fu_namelist_group_empty = .false.
          return
        endif
      end do
      fu_namelist_group_empty = .true.
    endif

  end function fu_namelist_group_empty


  !********************************************************************

  integer function fu_nbr_of_items_of_nl_group(nlGrp)
    !
    ! Returns the number filled items
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist_group), intent(in) :: nlGrp

    ! Local variables 
    integer :: i

    fu_nbr_of_items_of_nl_group = 0
    if(nlGrp%iLastFilledNl <= 0)then
      return
    else
      do i=1, nlGrp%iLastFilledNl
        if(.not. fu_namelist_empty(nlGrp%lists(i)%nl)) then
          fu_nbr_of_items_of_nl_group = fu_nbr_of_items_of_nl_group + &
                                      & nlGrp%lists(i)%nl%iLastFilledItem
        endif
      end do
    endif
  end function fu_nbr_of_items_of_nl_group


  !********************************************************************

  integer function fu_nbr_of_namelists(nlGrp)
    !
    ! Returns the number of filled namelists
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist_group), intent(in) :: nlGrp

    fu_nbr_of_namelists = nlGrp%iLastFilledNl

  end function fu_nbr_of_namelists


  !==============================================================================
  !==============================================================================
  !
  !  Namelist 
  !
  !==============================================================================
  !==============================================================================

  !******************************************************************************

  function fu_read_namelist(iUnit, ifEnv, chStopLine, chIgnoreItem, nItems, chNotNamedLine) result(nl)
    !
    ! Creates the namelist nl by reading it from the file pointed by iUnit
    ! The namelist reading is stopped at the chStopLine, which is also 
    ! consumed from the file. If this line has a namelist format - it is added
    ! to the namelist. If not - skipped
    ! In general, all named ilnes are included, all non-named (without '=' sign)
    ! are skipped
    !
    implicit none

    ! Return value
    type(Tsilam_namelist), pointer :: nl

    ! Imported parameters
    integer, intent(in) :: iUnit
    logical, intent(in) :: ifEnv
    character(len=*), intent(in), optional :: chStopLine, chIgnoreItem, chNotNamedLine
    integer, intent(in), optional :: nItems

    ! Local variables
    integer :: i
    logical :: eof
    character(len=fnlen) :: line

    ! Stupidity check
    !
    if(present(chStopLine))then
      if(len_trim(chStopLine) == 0)then
        call set_error('Empty stop line','fu_read_namelist')
        return
      endif
    endif

    !
    ! Create the empty namelist
    !
    allocate(nl, stat=i)
    if(i /= 0)then
      call msg_test('Strange allocation status: ',i)
      call set_error('Failed to allcoate memory for namelist','fu_read_namelist')
      return
    endif
    nl%iTotalSize=0
    nl%iLastFilledItem=0
    nl%chName = ''
    
    !
    ! Add a few items to start with
    !
    if(present(nItems))then
      if(nItems > 50 .and. nItems < 1000000)then
        call add_namelist_item(nl,nItemsTotal=nItems)
      else
        call add_namelist_item(nl)
      endif
    else
      call add_namelist_item(nl)
    endif
    if(error)return
    nl%defined = silja_true

    !
    ! Simply read the file through until eof or stop line
    !
    eof=.false.
    do while(.not.(eof.or.error))
      call next_line_from_input_file(iUnit, line, eof)
      if(error)return
      if(eof)then
        if(present(chStopLine)) then 
                call msg_warning('EOF reached prior to stopline', &
                                               & 'fu_read_namelist')
                call msg("Stop line was '"//trim(chStopLine)//"'")
        endif
        return
      endif
      !
      ! Namelist is over ?
      !
      if(present(chStopLine))then
!        call msg(fu_connect_strings(trim(line),'<>',trim(chStopLine)))
        if(trim(line) == trim(chStopLine)) then
!          call msg('GOT the end line')
          return
        endif
      endif

      !
      ! Split the line
      !
      nl%items(nl%iLastFilledItem+1)%nli%chName = ''
      nl%items(nl%iLastFilledItem+1)%nli%chContent = ''
      i=index(line,'=')
      if(i == 0)then
        if(present(chNotNamedLine))then
          if(index(line,chNotNamedLine) == 0) call msg_warning('Not named line:'+line,'fu_read_namelist')
        else
          call msg_warning('Not a named line:'+line,'fu_read_namelist')
        endif
        
      elseif(trim(fu_str_u_case(line(1:i-1))) == 'LIST')then
        !
        ! New list has a name, store it
        !
        nl%chName = fu_str_u_case(adjustl(line(i+1:len_trim(line))))

      elseif(trim(fu_str_u_case(line(1:i-1))) == 'END_LIST')then  ! this list is over
        if(nl%chName == fu_str_u_case(adjustl(line(i+1:len_trim(line)))))return

      else
        nl%items(nl%iLastFilledItem+1)%nli%chName = adjustl(line(1:i-1))
        if (ifEnv) then
           nl%items(nl%iLastFilledItem+1)%nli%chContent = adjustl(fu_expand_environment(line(i+1:len_trim(line))))
        else
           nl%items(nl%iLastFilledItem+1)%nli%chContent = adjustl(line(i+1:len_trim(line)))
        endif
        !
        ! If there is chIgnoreItem - skip the corresponding items
        !
        if(present(chIgnoreItem))then
          if(chIgnoreItem == nl%items(nl%iLastFilledItem+1)%nli%chName) then
            nl%items(nl%iLastFilledItem+1)%nli%chName = ''
            nl%items(nl%iLastFilledItem+1)%nli%chContent = ''
            cycle
          endif
        endif
        nl%iLastFilledItem = nl%iLastFilledItem + 1
        call add_namelist_item(nl) ! More space needed ?
        if(error)return
      endif
!      if(mod(nl%iLastFilledItem ,100) == 0) call msg('Namelist:',nl%iLastFilledItem )
    end do  ! through the input file

  end function fu_read_namelist



  !*******************************************************************************

  subroutine destroy_namelist(nl)
    !
    ! Frees the memory allocated for the namelist nl
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl

    ! Local variables
    integer :: i
    type(Tsilam_namelist), pointer :: nl__

    nl__ => nl

    if(associated(nl))then
      if(nl%defined == silja_true)then
!        do i=1, nl%nAllocated
!          deallocate(nl%items(i)%nli)
!          nullify(nl%items(i)%nli)
!        end do
        if (nl%nAllocated > 0)  deallocate(nl%items)
      end if
      deallocate(nl)
      nullify(nl)
    endif

  end subroutine destroy_namelist


  !********************************************************************************

  subroutine reset_namelist(nl)
    !
    ! cleans the namelist, so it can be further used
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl

    ! Local variables
    integer :: i

    if(associated(nl))then
      if(nl%defined == silja_true)then
        do i=1,nl%iLastFilledItem
          nl%items(i)%nli%chName = ''
          nl%items(i)%nli%chContent = ''
        end do
      end if
      nl%iLastFilledItem = 0
    endif

  end subroutine reset_namelist  


  !********************************************************************

  function fu_make_missing_namelist() result(nl)
    !
    ! Creates the missing namelist nl 
    !
    implicit none

    ! Return value
    type(Tsilam_namelist), pointer :: nl

    ! Local variables
    integer :: i

    !
    ! Create the empty namelist
    !
    allocate(nl, stat=i)
    if(i /= 0)then
      call msg_test('Strange allocation status: ',i)
      call set_error('Failed to allcoate memory for namelist','fu_make_missing_namelist')
      return
    endif
    nl%iTotalSize=0
    nl%iLastFilledItem=0
    nl%nAllocated=0 
    nl%chName = ''
    !
    ! Add an empty item - probably, not needed
    !
    call add_namelist_item(nl)
    if(error)return

  end function fu_make_missing_namelist


  !*********************************************************************************

  subroutine write_namelist(iUnit, nl, ifName, ifPublic) 
    !
    ! output in the form of namelist
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in):: nl
    integer, intent(in) :: iUnit
    logical, intent(in) :: ifName, ifPublic

    ! Local variables
    integer :: i, j
    integer, dimension(2) :: fUnits

    fUnits(1:2) = (/iUnit, 6/)

    if(nl%defined == silja_true)then
      do j = 1,2
        if(ifName)then
          write(fUnits(j),*)
          write(fUnits(j),fmt='(2A)')'LIST = ',trim(nl%chName)
        endif
        do i=1, nl%iLastFilledItem
          write(fUnits(j),fmt='(3A)') trim(nl%items(i)%nli%chName),' = ', &
                                & trim(nl%items(i)%nli%chContent)
        end do
        if(ifName)then
          write(fUnits(j),fmt='(2A)')'END_LIST = ',trim(nl%chName)
        endif
        if (smpi_global_rank /= 0) exit !Be quiet 
        if (.not. ifPublic) exit 
      enddo
    end if

  end subroutine write_namelist


  !*****************************************************************
  
  subroutine write_namelist_item_char(chName, chContent, iUnit)
    !
    ! Writes the given namelist item to the given file. If unit is not given,
    ! namelist goes into the log file
    !
    implicit none
    
    ! Imported parameters
    character(len=*), intent(in) :: chName, chContent
    integer, intent(in), optional :: iUnit
    
    ! Local variables
    integer :: iStat
    type(silam_sp) :: sp
    
    if(present(iUnit))then
      write(unit=iUnit,fmt='(3A)',iostat=iStat) trim(chName),' = ',trim(chContent)
    else
      sp%sp => fu_work_string()
      if(error)return
      write(unit=sp%sp,fmt='(3A)',iostat=iStat) trim(chName),' = ',trim(chContent)
      call msg(sp%sp)
      call free_work_array(sp%sp)
    endif
    if(iStat /= 0) call set_error('Failed writing the namelist item','write_namelist_item_char')
    
  end subroutine write_namelist_item_char


  !*****************************************************************
  
  subroutine write_namelist_item_int(chName, iContent, iUnit)
    !
    ! Writes the given namelist item to the given file. If unit is not given,
    ! namelist goes into the log file
    !
    implicit none
    
    ! Imported parameters
    character(len=*), intent(in) :: chName
    integer, intent(in) :: iContent
    integer, intent(in), optional :: iUnit
    
    ! Local variables
    integer :: iStat
    type(silam_sp) :: sp
    
    if(present(iUnit))then
      write(unit=iUnit,fmt='(2A,I15)',iostat=iStat) trim(chName),' = ',iContent
    else
      sp%sp => fu_work_string()
      if(error)return
      write(unit=sp%sp,fmt='(2A,I15)',iostat=iStat) trim(chName),' = ',iContent
      call msg(sp%sp)
      call free_work_array(sp%sp)
    endif
    if(iStat /= 0) call set_error('Failed writing the namelist item','write_namelist_item_int')
    
  end subroutine write_namelist_item_int


  !*****************************************************************
  
  subroutine write_namelist_item_real(chName, fContent, iUnit)
    !
    ! Writes the given namelist item to the given file. If unit is not given,
    ! namelist goes into the log file
    !
    implicit none
    
    ! Imported parameters
    character(len=*), intent(in) :: chName
    real, intent(in) :: fContent
    integer, intent(in), optional :: iUnit
    
    ! Local variables
    integer :: iStat
    type(silam_sp) :: sp
    
    if(present(iUnit))then
      if(abs(fContent) < 1.e5 .and. abs(fContent) > 1.0)then    ! from 1 to 99999
        write(unit=iUnit,fmt='(2A,F10.4)',iostat=iStat) trim(chName),' = ',fContent
      elseif(abs(fContent) <= 1.0 .and. abs(fContent) > 1.e-5)then  ! from 0.00001 to 1.0
        write(unit=iUnit,fmt='(2A,F9.6)',iostat=iStat) trim(chName),' = ',fContent
      else
        write(unit=iUnit,fmt='(2A,E12.6)',iostat=iStat) trim(chName),' = ',fContent
      endif
    else
      sp%sp => fu_work_string()
      if(error)return
      if(abs(fContent) < 1.e5 .and. abs(fContent) > 1.0)then
        write(unit=sp%sp,fmt='(2A,F10.4)',iostat=iStat) trim(chName),' = ',fContent
      elseif(abs(fContent) <= 1.0 .and. abs(fContent) > 1.e-5)then  ! from 0.00001 to 1.0
        write(unit=sp%sp,fmt='(2A,F9.6)',iostat=iStat) trim(chName),' = ',fContent
      else
        write(unit=sp%sp,fmt='(2A,E12.6)',iostat=iStat) trim(chName),' = ',fContent
      endif
      call msg(sp%sp)
      call free_work_array(sp%sp)
    endif
    if(iStat /= 0) call set_error('Failed writing the namelist item','write_namelist_item_real')
    
  end subroutine write_namelist_item_real


  !*****************************************************************

  function fu_create_namelist(chNm) result(nlPtr)
    !
    ! Creates an empty new namelist and returns a pointer on it
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in), optional :: chNm

    ! Return value
    type(Tsilam_namelist), pointer :: nlPtr

    ! Internal variables
    integer :: iStat

    allocate(nlPtr, stat = iStat)
    if(iStat /= 0)then
      call set_error('Failed to allocate namelist','fu_create_namelist')
      nullify(nlPtr)
      return
    endif

    nullify(nlPtr%items)
    nlPtr%iLastFilledItem = 0
    nlPtr%iTotalSize = 0
    nlPtr%nAllocated = 0
    if(present(chNm))then
      nlPtr%chName = trim(chNm)
    else
      nlPtr%chName = ''
    endif
    nlPtr%defined = silja_true

  end function fu_create_namelist


  ! ***************************************************************

  subroutine next_named_line_from_file(FUnit, ifMust, ifEnv, chName, chContent, eof)
    ! 
    ! Returns next non-comment line fron an open input file (source
    ! term file, control parameter file or trajectory input file).
    ! Also if there's comment in the end of line, it is removed.
    ! All lines beginning with # or ! are considered to be comments.
    ! Line is split to name and the content, which are to be separated by the 
    ! left-most '=' sign. So, the line name must NEVER contain this sign,
    ! while the content may.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: FUnit
    logical, intent(in) :: ifMust ! if the line must be named line
    logical, intent(in) :: ifEnv ! if expand the environment

    ! Imported parameters with intent OUT:
    CHARACTER (LEN=*), INTENT(out) ::  chName, chContent
    LOGICAL, INTENT(out) :: eof

    ! Local declarations:
    integer :: i
    character(len=fnlen) :: line

    !
    ! Read a string
    !
    call next_line_from_input_file(FUnit, line, eof)
    if(error.or.eof)return
    !
    ! Split the line
    !
    chName = ''; chContent = ''
    i=index(line,'=')
    if(i == 0)then
      if(ifMust)then
        call set_error(fu_connect_strings('Not a named line:',line),'next_named_line_from_file')
      else
        chName = trim(line)
      endif
    else
      chName(1:i-1) = line(1:i-1)
      chName = adjustl(chName)
      if (ifEnv) then
         chContent = adjustl(fu_expand_environment(line(i+1:len_trim(line))))
      else
         chContent = adjustl(line(i+1:len_trim(line)))
      endif
    endif

  end subroutine next_named_line_from_file


  !===========================================================================
  !===========================================================================
  !
  !    Private routines
  !
  !===========================================================================
  !===========================================================================


  function fu_add_existing_namelist(nlGrp, nl) result(nlPtr)
    !
    ! Adds an existing namelist into the group, returns its index
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), pointer :: nlGrp
    type(Tsilam_namelist), pointer :: nl

    ! return value
    type(Tsilam_namelist), pointer :: nlPtr

    !
    ! Do we have free pointers ?
    !
    if(nlGrp%iLastFilledNL == nlGrp%iTotalSize)then
      !
      ! Have to increase the pointer array. 
      !
      call increase_ptr_array_nl_group(nlGrp)
      if(error)return
    endif

    nlGrp%iLastFilledNL = nlGrp%iLastFilledNL + 1
    nlGrp%lists(nlGrp%iLastFilledNL)%nl => nl
    nlGrp%nAllocated = nlGrp%iLastFilledNL

    nlPtr => nlGrp%lists(nlGrp%iLastFilledNL)%nl

  end function fu_add_existing_namelist


  !***********************************************************************

  function fu_add_new_namelist(nlGrp, chNamelistName) result(nlPtr)
    !
    ! Allocates space for one more namelist in the group.
    ! Idea: the namelist itself contains only pointers to the items,
    ! which can be allocated one-by-one with further setting the
    ! address.
    ! What I will have to do is to check the number of pointers
    ! in the namelist is adequate.
    ! The principle is: do NOT keep empty items, they are too big
    ! While keeping empty pointers is OK - they are small.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), pointer :: nlGrp
!    integer, intent(out) :: iPtr
    character(len=*),intent(in), optional :: chNamelistName

    ! Return variable
    type(Tsilam_namelist), pointer :: nlPtr
    
    ! Local variable
    integer :: i

    !
    ! Do we have free pointers ?
    !
    if(nlGrp%iLastFilledNL == nlGrp%iTotalSize)then
      !
      ! Have to increase the pointer array. 
      !
      call increase_ptr_array_nl_group(nlGrp)
      if(error)return
    endif
    !
    ! Just add one more list and put a few items in it
    !
    allocate(nlGrp%lists(nlGrp%iLastFilledNL+1)%nl, stat=i)
    if(i /= 0)then
      call msg_test('Strange allocation status: ',i)
      call set_error('Failed to allocate memory for namelist item','add_new_namelist')
      return
    endif

    nlGrp%iLastFilledNL = nlGrp%iLastFilledNL + 1  ! Increase the number of filled-in namelists

    nlGrp%nAllocated = nlGrp%iLastFilledNL
    nlGrp%lists(nlGrp%iLastFilledNL)%nl%iTotalSize = 0
    nlGrp%lists(nlGrp%iLastFilledNL)%nl%iLastFilledItem = 0
    nlGrp%lists(nlGrp%iLastFilledNL)%nl%nAllocated = 0
    call add_namelist_item(nlGrp%lists(nlGrp%iLastFilledNL)%nl)
    if(error)return
    nlGrp%lists(nlGrp%iLastFilledNL)%nl%defined = silja_true
    if(present(chNamelistName))then
      nlGrp%lists(nlGrp%iLastFilledNL)%nl%chName = chNamelistName
    else
      nlGrp%lists(nlGrp%iLastFilledNL)%nl%chName = ''
    endif

    nlPtr => nlGrp%lists(nlGrp%iLastFilledNL)%nl

  end function fu_add_new_namelist


  !**********************************************************************

  subroutine increase_ptr_array_nl_group(nlGrp)

    implicit none

    ! Imported parameter
    type(Tsilam_namelist_group), pointer :: nlGrp

    ! Local variables
    type(Tsilam_namelist_ptr), dimension(:), pointer :: nlPtr
    integer :: i

    if(nlGrp%iTotalSize > 0)then
      !
      ! Create temporary place, store the 
      ! existing pointers to it and then re-allocate space 
      !
      allocate(nlPtr(nlGrp%iTotalSize), stat=i)
      if(i /= 0)then
        call msg_test('Strange allocation status: ',i)
        call set_error('Failed to allocate memory for temporary pointers', &
                     & 'increase_ptr_array_nl_group')
        return
      endif
      !
      ! Store pointers to temporary place
      !
      do i=1, nlGrp%iTotalSize
        nlPtr(i)%nl => nlGrp%lists(i)%nl
        nullify(nlGrp%lists(i)%nl)  ! Break the connection 
      end do

      deallocate(nlGrp%lists) ! Only pointers !

    end if  ! Something to store temporarily
    !
    ! Re-allocate the main array of pointers
    !
    allocate(nlGrp%lists(nlGrp%iTotalSize+2),stat=i)
    if(i /= 0)then
      call msg_test('Strange allocation status: ',i)
      call set_error('Failed to re-allocate memory for namelists', &
                   & 'increase_ptr_array_nl_group')
      return
    endif

    if(nlGrp%iTotalSize > 0)then
      !
      ! Return the pointers to the right place
      !
      do i=1,nlGrp%iTotalSize
        nlGrp%lists(i)%nl => nlPtr(i)%nl
      end do
      do i=nlGrp%iTotalSize+1, size(nlGrp%lists)
        nullify(nlGrp%lists(i)%nl)
      end do
      deallocate(nlPtr)
    endif

    nlGrp%iTotalSize = nlGrp%iTotalSize + 2

  end subroutine increase_ptr_array_nl_group


  !***********************************************************************
  
  subroutine add_existing_namelist_item(nl, item, nItemsTotal)
    !
    ! Adds an existing namelist item into the namelist. In fact, just calls the below sub
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl
    type(Tsilam_nl_item_ptr), pointer :: item
    integer, intent(in), optional :: nItemsTotal
    
    if(present(nItemsTotal))then
      call add_namelist_item_basic(nl, item%nli%chName, item%nli%chContent, nItemsTotal)
    else
      call add_namelist_item_basic(nl, item%nli%chName, item%nli%chContent)
    endif
  end subroutine add_existing_namelist_item


  !***********************************************************************

  subroutine add_namelist_item_basic(nl, chNm, chContent, nItemsTotal)  
    !
    ! Allocates one more namelist item.
    ! Idea: the namelist itself contains only pointers to the items,
    ! which can be allocated one-by-one with further setting the
    ! address.
    ! What I will have to do is to check the number of pointers
    ! in the namelist is adequate.
    ! The principle is: do NOT keep empty items, they are too big
    ! While keeping empty pointers is OK - they are small.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl
    character(len=*), intent(in), optional :: chNm, chContent
    integer, intent(in), optional :: nItemsTotal

    ! Local variables
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: arNlPtr
    integer :: i, nNewSize

    if(.not.associated(nl))then
      call set_error('Empty namelist pointer','add_namelist_item_basic')
      return
    endif
    !
    ! Size recommendation may be given
    !
    if(present(nItemsTotal))then
      if(nItemsTotal > 50 .and. nItemsTotal < 1000000)then
        nNewSize = nItemsTotal
      else
        nNewSize = 100
      endif
    else
        nNewSize = 100
    endif
    !
    ! Do we have free pointers ?
    !
    if(nl%iLastFilledItem == nl%iTotalSize)then
      !
      ! Have to increase the pointer array. 
      !
      if(nl%iTotalSize > 0)then
        !
        ! If there is already something in the list, create temporary place, store the 
        ! existing item pointers to it and then re-allocate space 
        !
        allocate(arNlPtr(nl%iTotalSize), stat=i)
        if(i /= 0)then
          call msg_test('Strange allocation status: ',i)
          call set_error('Failed to allocate memory for temporary item pointers', &
                       & 'add_namelist_item_basic')
          return
        endif
        !
        ! Store pointers to temporary place
        !
        do i=1, nl%iTotalSize
          arNlPtr(i)%nli%chName = nl%items(i)%nli%chName
          arNlPtr(i)%nli%chContent = nl%items(i)%nli%chContent

!          nullify(nl%items(i)%nli)  ! Break the connection 
        end do
        !
        ! Re-allocate the main array of pointers
        !
!        call msg('10')
        deallocate(nl%items) ! does not destroy items themselves
!        call msg('20')
      endif ! List is not empty
      !
      ! Allocate proper space for the list items
      !
      allocate(nl%items(max(int(nl%iTotalSize*1.5),nNewSize)),stat=i)
!      allocate(nl%items(nl%iTotalSize+500),stat=i)
      if(i /= 0)then
        call msg_test('Strange allocation status: ',i)
        call set_error('Failed to re-allocate memory for namelist items', &
                     & 'add_namelist_item_basic')
        return
      endif

      if(nl%iTotalSize > 0)then
        !
        ! Return the pointers to the right place
        !
        do i=1,nl%iTotalSize
          nl%items(i)%nli%chName = arNlPtr(i)%nli%chName
          nl%items(i)%nli%chContent = arNlPtr(i)%nli%chContent
        end do
!        do i=nl%iTotalSize+1, size(nl%items)
!          nullify(nl%items(i)%nli)
!        end do

        deallocate(arNLPtr)
      endif

      nl%iTotalSize = size(nl%items)

    endif  ! If spare pointers are available

    !
    ! Just add one more item
    !
!    allocate(nl%items(nl%iLastFilledItem+1)%nli, stat=i)
!    if(i /= 0)then
!      call msg_test('Strange allocation status: ',i)
!      call set_error('Failed to allocate memory for namelist item','add_namelist_item')
!      return
!    endif
    nl%nAllocated = nl%iLastFilledItem+1
    !
    ! Should some meaningfull information is given - write it down to the created item
    !
    if(present(chNm))then
      nl%items(nl%iLastFilledItem+1)%nli%chName = chNm
      if(present(chContent))then
        nl%items(nl%iLastFilledItem+1)%nli%chContent = chContent
      else
        nl%items(nl%iLastFilledItem+1)%nli%chContent = ''
      endif
      nl%iLastFilledItem = nl%iLastFilledItem + 1
      nl%defined = silja_true
    else
      nl%items(nl%iLastFilledItem+1)%nli%chName = ''
    endif

  end subroutine add_namelist_item_basic


  !********************************************************************
  
  subroutine replace_namelist_item_by_name(nl, chOldNm, chNewNm, chNewContent)
  
    implicit none
    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl
    character(len=*), intent(in) :: chOldNm, chNewNm, chNewContent

    ! Local variables
    integer :: iTmp, iItemCount
    !
    ! The task is not trivial: the namelist may have several items with the same name.
    ! Set error in this case
    !
    iItemCount=0
    do iTmp=1, nl%iLastFilledItem
      if(nl%items(iTmp)%nli%chName == chOldNm)then
        if(iItemCount == 0)then
          nl%items(iTmp)%nli%chName = chNewNm
          nl%items(iTmp)%nli%chContent = chNewContent
          iItemCount = iItemCount+1
        else
          call set_error('Namelist item >>' + chOldNm + '<< is met for the second time', &
                       & 'replace_namelist_item_by_name')
          return
        endif
      endif
    end do

  end subroutine replace_namelist_item_by_name


  !********************************************************************
  
  subroutine replace_namelist_item_content(nlItem, chOldNm, chNewNm, chNewContent)
  
    implicit none
    ! Imported parameters
    type(Tsilam_nl_item_ptr), pointer :: nlItem
    character(len=*), intent(in) :: chOldNm, chNewNm, chNewContent

    nlItem%nli%chName = chNewNm
    nlItem%nli%chContent = chNewContent

  end subroutine replace_namelist_item_content


  !*********************************************************************

  subroutine items_from_namelist(nl, chNm, pItems, nItems) ! all items with given name
    !
    ! Searches the whole namelist and returns an array of item pointers
    ! pointing to the found items.
    ! Important: here we do whatever we need with the pointer array, but if
    ! it happened to point to somewhere, those data are not touched. Pointers
    ! are just redirected, leaving the old-pointed memory as it was.
    ! So, NEVER use the pItems for allocation of real namelist items
    !
    implicit none
    !
    ! Imported parameters
    !
    type(Tsilam_namelist), intent(in) :: nl
    character(len=*), intent(in) :: chNm ! name of the items to search
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems ! output array
    integer, intent(out) :: nItems

    ! Local variables
    !
    integer :: iItemCount, iTmp, iSize

    if(associated(pItems))then
      iSize = size(pItems)
    else
      iSize = 0
    endif
    !
    ! Starting preparations. 
    !
    iItemCount=0
    do iTmp=1, nl%iLastFilledItem
      if(nl%items(iTmp)%nli%chName == chNm)iItemCount=iItemCount+1
    end do
    nItems = iItemCount

    if(iItemCount == 0)then ! No such items, destroy the pointer array
      if(iSize > 0) deallocate(pItems)
      nullify(pItems)
      return
    elseif(iItemCount /= iSize)then ! Have to re-allocate memory for pointers
      if(associated(pItems) .and. iSize > 0 .and. iSize < 1000000) then
        deallocate(pItems, stat=iTmp)
        if(iTmp /= 0)then
          call msg('status, iItemCount, array iSize:' + fu_str(iTmp), iItemCount, iSize)
          call set_error('Failed to deallocate item pointer array','items_from_namelist')
          call unset_error('items_from_namelist')
          return
        endif
      endif
      allocate(pItems(iItemCount), stat=iTmp)
      if(iTmp /= 0)then
        call set_error('Failed to allocate item pointer array','items_from_namelist')
        return
      endif
    else        ! No need to reallocate memory, size matches exactly
    endif
    !
    ! Store the pointers. Note: the old-pointed addresses are just dropped !!
    !
    iItemCount=1
    do iTmp=1, nl%iLastFilledItem
      if(nl%items(iTmp)%nli%chName == chNm)then
        pItems(iItemCount)%nli%chName = nl%items(iTmp)%nli%chName
        pItems(iItemCount)%nli%chContent = nl%items(iTmp)%nli%chContent
        iItemCount = iItemCount+1
      endif
    end do

  end subroutine items_from_namelist


  !**********************************************************************

  subroutine destroy_namelist_item_list(pItems)
    !
    ! Deallocates the items array. 
    ! ATTENTION. If items are reverted to an array of pointers, they also MUST 
    ! be deallocated or some other treatment must be taken
    !
    implicit none

    ! Imported parameter
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems ! output array

    if(associated(pItems))then
      if(size(pItems) > 0 .and. size(pItems) < worksize*10)then
        deallocate(pItems)
      endif
    endif
    nullify(pItems)

  end subroutine destroy_namelist_item_list


  !*********************************************************************

  function fu_get_item_by_index(nl, indexNl)result(item)
    !
    ! Returns an item, which is pointed by the index in the namelist.
    ! CAREFUL! Nobody promised any order inthe namelist. Use this function only if
    ! the namelist is created inside the programme, so that the order of items is under 
    ! full control
    !
    implicit none

    ! return value
    type(Tsilam_nl_item_ptr),pointer :: item
    
    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl
    integer, intent(in) :: indexNl

    ! local variables
    type(Tsilam_namelist), pointer :: nlp

    nlp => nl

    if(indexNl < 1 .or. indexNl > nl%iLastFilledItem)then
      call msg('Strange index:',indexNl)
      call set_error('Strange index','fu_get_item_by_index')
      nullify(item)
      return
    endif

    item => nl%items(indexNl)

  end function fu_get_item_by_index


  !*********************************************************************

  function fu_get_item_by_name(nl, chItemName)result(item)
    !
    ! Returns an item, which is pointed by the index in the namelist.
    ! CAREFUL! Nobody promised any order inthe namelist. Use this function only if
    ! the namelist is created inside the programme, so that the order of items is under 
    ! full control
    !
    implicit none

    ! return value
    type(Tsilam_nl_item_ptr),pointer :: item
    
    ! Imported parameters
    type(Tsilam_namelist), pointer :: nl
    character(len=*), intent(in) :: chItemName

    ! Local variables
    integer :: i

    do i=1,nl%iLastFilledItem
      if(fu_str_u_case(nl%items(i)%nli%chName) == fu_str_u_case(chItemName))then
        item => nl%items(i)
        return
      endif
    end do

  end function fu_get_item_by_name


  !*********************************************************************

  subroutine items_from_nl_group(nlGrp, chNlName, chItemNm, pItems, nItems)
    !
    ! Searches the whole namelist group and returns an array of item pointers
    ! pointing to the found items. The search is done until the first namelist 
    ! with the given name is met. 
    !
    implicit none
    !
    ! Imported parameters
    !
    type(Tsilam_namelist_group), pointer :: nlGrp
    character(len=*), intent(in) :: chNlName, chItemNm  ! name of namelist and items
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems ! output array
    integer, intent(out) :: nItems

    ! local vars
    integer :: iList

    do iList=1,nlGrp%iLastFilledNL
      if(fu_str_u_case(chNlName) == fu_str_u_case(nlGrp%lists(iList)%nl%chName)) then
        call items_from_namelist(nlGrp%lists(iList)%nl, chItemNm, pItems, nItems)
        return
      endif
    end do

  end subroutine items_from_nl_group


  !******************************************************************
  !******************************************************************
  !
  !  Extraction of various contents from items, namelists and namelist groups
  !
  !******************************************************************
  !******************************************************************

  function fu_namelist_name(nl) result(name)
    implicit none
    character(len=clen) :: name
    type(Tsilam_namelist) :: nl
    name = nl%chName
  end function fu_namelist_name


  !*********************************************************************

  function fu_item_name(nlItem) result(chName)
    !
    ! Actually, encapsulation
    !
    implicit none

    character(len=clen) :: chName
    type(Tsilam_namelist_item), intent(in) :: nlItem

    chName = nlItem%chName

  end function fu_item_name



  !*********************************************************************

  function fu_itemPtr_name(nlItem) result(chName)
    !
    ! Actually, encapsulation
    !
    implicit none

    character(len=clen) :: chName
    type(Tsilam_nl_item_ptr), intent(in) :: nlItem

    chName = nlItem%nli%chName

  end function fu_itemPtr_name



  !*********************************************************************

  function fu_item_content(nlItem) result(chContent)
    !
    ! Actually, encapsulation
    !
    implicit none

    character(len=fnlen) :: chContent
    type(Tsilam_namelist_item), intent(in) :: nlItem

    chContent = nlItem%chContent

  end function fu_item_content


  !*********************************************************************

  real function fu_item_content_real (nlItem) result(fContent)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_namelist_item), intent(in) :: nlItem
    integer :: status

    read(unit=nlItem%chContent,fmt=*,iostat=status) fContent
    if(status /= 0)then
      call set_error(fu_connect_strings('Cannot get real content from:', &
                                      & nlItem%chContent), &
                   & 'fu_item_content_real')
      fContent = real_missing
    endif

  end function fu_item_content_real


  !*********************************************************************

  integer function fu_item_content_int(nlItem) result(iContent)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_namelist_item), intent(in) :: nlItem
    integer :: status

    read(unit=nlItem%chContent,fmt=*,iostat=status) iContent
    if(status /= 0)then
      call set_error(fu_connect_strings('Cannot get int content value from:', &
                                      & nlItem%chContent), &
                   & 'fu_item_content_int')
      iContent = int_missing
    endif

  end function fu_item_content_int


  !*********************************************************************

  function fu_itemPtr_content(nlItemPtr) result(chContent)
    !
    ! Actually, encapsulation
    !
    implicit none

    character(len=fnlen) :: chContent
    type(Tsilam_nl_item_ptr), intent(in) :: nlItemPtr

    chContent = nlItemPtr%nli%chContent

  end function fu_itemPtr_content


  !*********************************************************************

  real function fu_itemPtr_content_real (nlItemPtr) result(fContent)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_nl_item_ptr), intent(in) :: nlItemPtr
    integer :: status

    read(unit=nlItemPtr%nli%chContent,fmt=*,iostat=status) fContent
    if(status /= 0)then
      call set_error(fu_connect_strings('Cannot get real content from:', &
                                      & nlItemPtr%nli%chContent), &
                   & 'fu_item_content_real')
      fContent = real_missing
    endif

  end function fu_itemPtr_content_real


  !*********************************************************************

  integer function fu_itemPtr_content_int(nlItemPtr) result(iContent)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_nl_item_ptr), intent(in) :: nlItemPtr
    integer :: status

    read(unit=nlItemPtr%nli%chContent,fmt=*,iostat=status) iContent
    if(status /= 0)then
      call set_error(fu_connect_strings('Cannot get int content value from:', &
                                      & nlItemPtr%nli%chContent), &
                   & 'fu_item_content_int')
      iContent = int_missing
    endif

  end function fu_itemPtr_content_int


  !**********************************************************************

  function fu_content_from_namelist(nl, chItemName) result(chContent)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Return value
    character(len=fnlen) :: chContent

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nl
    character(len=*), intent(in) :: chItemName

    ! local vars
    integer :: i

    do i=1,nl%iLastFilledItem
      if(fu_str_u_case(nl%items(i)%nli%chName) == fu_str_u_case(chItemName))then
        chContent = nl%items(i)%nli%chContent
        return
      endif
    end do
    chContent = ''

  end function fu_content_from_namelist


  !**********************************************************************

  real function fu_content_from_namelist_real(nl, chItemName) result(fContent)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nl
    character(len=*), intent(in) :: chItemName

    ! local vars
    integer :: i, status

    do i=1,nl%iLastFilledItem
      if(fu_str_u_case(nl%items(i)%nli%chName) == fu_str_u_case(chItemName))then
        read(unit=nl%items(i)%nli%chContent,fmt=*,iostat=status) fContent
        if(status /= 0)then
          call set_error('Cannot get real content value from:' + nl%items(i)%nli%chContent, &
                       & 'fu_content_from_namelist_real')
          fContent = real_missing
        endif
        return
      endif
    end do
    fContent = real_missing ! If nothing found...

  end function fu_content_from_namelist_real


  !**********************************************************************

  integer function fu_content_from_namelist_int(nl, chItemName) result(iContent)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nl
    character(len=*), intent(in) :: chItemName

    ! local vars
    integer :: i, status

    do i=1,nl%iLastFilledItem
      if(fu_str_u_case(nl%items(i)%nli%chName) == fu_str_u_case(chItemName))then
        read(unit=nl%items(i)%nli%chContent,fmt=*,iostat=status) iContent
        if(status /= 0)then
          call set_error('Cannot get int content value from:' + nl%items(i)%nli%chContent, &
                       & 'fu_content_from_namelist_int')
          iContent = int_missing
        endif
        return
      endif
    end do
    iContent = int_missing  ! If nothing found...

  end function fu_content_from_namelist_int


  !**********************************************************************

  function fu_content_from_nl_group(nlGrp, chItemName, chListNm) result(chContent)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Return value
    character(len=fnlen) :: chContent

    ! Imported parameters
    type(Tsilam_namelist_group), intent(in) :: nlGrp
    character(len=*), intent(in) :: chItemName
    character(len=*), intent(in), optional :: chListNm

    ! local vars
    integer :: iList, i

    do iList=1,nlGrp%iLastFilledNL
      if(present(chListNm))then   ! Skip the list if name is different
        if(fu_str_u_case(chListNm) /= fu_str_u_case(nlGrp%lists(iList)%nl%chName)) continue 
      endif
      do i=1,nlGrp%lists(iList)%nl%iLastFilledItem
        if(fu_str_u_case(nlGrp%lists(iList)%nl%items(i)%nli%chName) == &
         & fu_str_u_case(chItemName))then
          chContent = nlGrp%lists(iList)%nl%items(i)%nli%chContent
          return
        endif
      end do
    end do
    chContent = ''

  end function fu_content_from_nl_group


  !**********************************************************************

  real function fu_content_from_nl_group_real(nlGrp, chItemName, chListNm) result(fContent)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), intent(in) :: nlGrp
    character(len=*), intent(in) :: chItemName
    character(len=*), intent(in), optional :: chListNm

    ! local vars
    integer :: iList, i, status

    do iList=1, nlGrp%iLastFilledNL
      if(present(chListNm))then   ! Skip the list if name is different
        if(fu_str_u_case(chListNm) /= fu_str_u_case(nlGrp%lists(iList)%nl%chName)) continue 
      endif
      do i=1,nlGrp%lists(iList)%nl%iLastFilledItem
        if(fu_str_u_case(nlGrp%lists(iList)%nl%items(i)%nli%chName) == &
         & fu_str_u_case(chItemName))then
          read(unit=nlGrp%lists(iList)%nl%items(i)%nli%chContent,fmt=*,iostat=status) fContent
          if(status /= 0)then
            call set_error(fu_connect_strings('Cannot get real content value from:', &
                                            & nlGrp%lists(iList)%nl%items(i)%nli%chContent), &
                         & 'fu_content_from_nl_group_real')
            fContent = real_missing
          endif
          return
        endif
      end do
    end do
    fContent = real_missing ! If nothing found...

  end function fu_content_from_nl_group_real


  !**********************************************************************

  integer function fu_content_from_nl_group_int(nlGrp, chItemNm, chListNm) result(iContent)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), intent(in) :: nlGrp
    character(len=*), intent(in) :: chItemNm
    character(len=*), intent(in), optional :: chListNm

    ! local vars
    integer :: iList, i, status

    do iList=1,nlGrp%iLastFilledNL
      if(present(chListNm))then   ! Skip the list if name is different
        if(fu_str_u_case(chListNm) /= fu_str_u_case(nlGrp%lists(iList)%nl%chName)) continue 
      endif
      do i=1,nlGrp%lists(iList)%nl%iLastFilledItem
        if(fu_str_u_case(nlGrp%lists(iList)%nl%items(i)%nli%chName) == fu_str_u_case(chItemNm))then
          read(unit=nlGrp%lists(iList)%nl%items(i)%nli%chContent,fmt=*,iostat=status) iContent
          if(status /= 0)then
            call set_error('Cannot get int content value from:' + &
                                            & nlGrp%lists(iList)%nl%items(i)%nli%chContent, &
                         & 'fu_content_from_namelist_group_int')
            iContent = int_missing
          endif
          return
        endif
      end do
    end do
    iContent = int_missing  ! If nothing found...

  end function fu_content_from_nl_group_int


  !*********************************************************************

  subroutine item_cnt_real_with_unit(nlItem, fContent, chUnit)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_namelist_item), intent(in) :: nlItem
    real, intent(out) :: fContent
    character(len=*), intent(out) :: chUnit
    integer :: status

    read(unit=nlItem%chContent,fmt=*,iostat=status) fContent, chUnit
    if(status /= 0)then
      call set_error('Cannot get real content from:' + nlItem%chContent, 'fu_item_content_real')
      fContent = real_missing
      chUnit = ''
    endif

  end subroutine item_cnt_real_with_unit


  !*********************************************************************

  subroutine item_cnt_int_with_unit(nlItem, iContent, chUnit)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_namelist_item), intent(in) :: nlItem
    integer, intent(out) :: iContent
    character(len=*), intent(out) :: chUnit
    integer :: status

    read(unit=nlItem%chContent,fmt=*,iostat=status) iContent, chUnit
    if(status /= 0)then
      call set_error('Cannot get int content value from:' + nlItem%chContent, 'fu_item_content_int')
      iContent = int_missing
      chUnit = ''
    endif

  end subroutine item_cnt_int_with_unit

  !*********************************************************************

  subroutine itemPtr_cnt_real_with_unit(nlItemPtr, fContent, chUnit)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_nl_item_ptr), intent(in) :: nlItemPtr
    real, intent(out) :: fContent
    character(len=*), intent(out) :: chUnit
    integer :: status

    read(unit=nlItemPtr%nli%chContent,fmt=*,iostat=status) fContent,chUnit
    if(status /= 0)then
      call set_error('Cannot get real content from:' + nlItemPtr%nli%chContent, 'fu_item_content_real')
      fContent = real_missing
      chUnit = ''
    endif

  end subroutine itemPtr_cnt_real_with_unit

  !*********************************************************************

  subroutine itemPtr_cnt_int_with_unit(nlItemPtr, iContent, chUnit)
    !
    ! Actually, encapsulation
    !
    implicit none

    type(Tsilam_nl_item_ptr), intent(in) :: nlItemPtr
    integer, intent(out) :: iContent
    character(len=*), intent(out) :: chUnit
    integer :: status

    read(unit=nlItemPtr%nli%chContent,fmt=*,iostat=status) iContent, chUnit
    if(status /= 0)then
      call set_error('Cannot get int content value from:' + nlItemPtr%nli%chContent, &
                   & 'fu_item_content_int')
      iContent = int_missing
      chUnit = ''
    endif

  end subroutine itemPtr_cnt_int_with_unit

  !**********************************************************************

  subroutine cnt_nl_real_with_unit(nl, chItemName, fContent, chUnit)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nl
    character(len=*), intent(in) :: chItemName
    real, intent(out) :: fContent
    character(len=*), intent(out) :: chUnit

    ! local vars
    integer :: i, status

    do i=1,nl%iLastFilledItem
      if(fu_str_u_case(nl%items(i)%nli%chName) == fu_str_u_case(chItemName))then
        !
        ! Careful: unit may be split by / sign, which is treated by read as a delimiter
        !
        if(index(nl%items(i)%nli%chContent,' ') == 0)then
          call set_error('No space delimiter found in:' + nl%items(i)%nli%chContent, &
                       & 'fu_content_from_namelist_real')
        else
          read(unit=nl%items(i)%nli%chContent,fmt=*,iostat=status) fContent
          if(status /= 0)then
            call set_error('Cannot get real content value from:' + nl%items(i)%nli%chContent, &
                         & 'fu_content_from_namelist_real')
            fContent = real_missing
            chUnit = ''
          endif
          chUnit = nl%items(i)%nli%chContent(index(nl%items(i)%nli%chContent,' ')+1:)
          chUnit = adjustl(chUnit)
        endif
        return
      endif
    end do
    fContent = real_missing ! If nothing found...
    chUnit = ''

  end subroutine cnt_nl_real_with_unit


  !**********************************************************************

  subroutine cnt_nl_int_with_unit(nl, chItemName, iContent, chUnit)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nl
    character(len=*), intent(in) :: chItemName
    integer, intent(out) :: iContent
    character(len=*), intent(out) :: chUnit

    ! local vars
    integer :: i, status

    do i=1,nl%iLastFilledItem
      if(fu_str_u_case(nl%items(i)%nli%chName) == fu_str_u_case(chItemName))then
        read(unit=nl%items(i)%nli%chContent,fmt=*,iostat=status) iContent, chUnit
        if(status /= 0)then
          call set_error('Cannot get int content value from:' + nl%items(i)%nli%chContent, &
                       & 'fu_content_from_namelist_int')
          iContent = int_missing
          chUnit = ''
        endif
        return
      endif
    end do
    iContent = int_missing  ! If nothing found...
    chUnit = ''

  end subroutine cnt_nl_int_with_unit

  !**********************************************************************

  subroutine cnt_nl_group_real_with_unit(nlGrp, chItemName, fContent, chUnit, chListNm)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), intent(in) :: nlGrp
    character(len=*), intent(in) :: chItemName
    character(len=*), intent(in), optional :: chListNm
    real, intent(out) :: fContent
    character(len=*), intent(out) :: chUnit

    ! local vars
    integer :: iList, i, status

    do iList=1, nlGrp%iLastFilledNL
      if(present(chListNm))then   ! Skip the list if name is different
        if(fu_str_u_case(chListNm) /= fu_str_u_case(nlGrp%lists(iList)%nl%chName)) continue 
      endif
      do i=1,nlGrp%lists(iList)%nl%iLastFilledItem
        if(fu_str_u_case(nlGrp%lists(iList)%nl%items(i)%nli%chName) == &
         & fu_str_u_case(chItemName))then
          read(unit=nlGrp%lists(iList)%nl%items(i)%nli%chContent, &
             & fmt=*, &
             & iostat=status) fContent, chUnit
          if(status /= 0)then
            call set_error('Cannot get real content value from:' + &
                                            & nlGrp%lists(iList)%nl%items(i)%nli%chContent, &
                         & 'fu_content_from_nl_group_real')
            fContent = real_missing
            chUnit = ''
          endif
          return
        endif
      end do
    end do
    fContent = real_missing ! If nothing found...
    chUnit = ''

  end subroutine cnt_nl_group_real_with_unit


  !**********************************************************************
  subroutine cnt_nl_grp_int_with_unit(nlGrp, chItemNm, chListNm, iContent, chUnit)
    !
    ! Searches the item name and returns the content of the 
    ! appropriate item
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist_group), intent(in) :: nlGrp
    character(len=*), intent(in) :: chItemNm
    character(len=*), intent(in), optional :: chListNm
    integer, intent(out) :: iContent
    character(len=*), intent(out) :: chUnit

    ! local vars
    integer :: iList, i, status

    do iList=1,nlGrp%iLastFilledNL
      if(present(chListNm))then   ! Skip the list if name is different
        if(fu_str_u_case(chListNm) /= fu_str_u_case(nlGrp%lists(iList)%nl%chName)) continue 
      endif
      do i=1,nlGrp%lists(iList)%nl%iLastFilledItem
        if(fu_str_u_case(nlGrp%lists(iList)%nl%items(i)%nli%chName) == fu_str_u_case(chItemNm))then
          read(unit=nlGrp%lists(iList)%nl%items(i)%nli%chContent, &
             & fmt=*, &
             & iostat=status) iContent, chUnit
          if(status /= 0)then
            call set_error('Cannot get int content value from:' + &
                                            & nlGrp%lists(iList)%nl%items(i)%nli%chContent, &
                         & 'fu_content_from_namelist_group_int')
            iContent = int_missing
            chUnit = ''
          endif
          return
        endif
      end do
    end do
    iContent = int_missing  ! If nothing found...
    chUnit = ''

  end subroutine cnt_nl_grp_int_with_unit


  !***********************************************************************

  subroutine print_nlItem(nlItem)
    !
    ! Overloading the report sub
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist_item), intent(inout) :: nlItem

    call msg(fu_connect_strings(nlItem%chName,' = ', nlItem%chContent))

  end subroutine


  !***********************************************************************

  subroutine print_nlItemPtr(nlItemPtr)
    !
    ! Overloading the report sub
    !
    implicit none

    ! Imported parameter
    type(Tsilam_nl_item_ptr), intent(inout) :: nlItemPtr

    call msg(fu_connect_strings(nlItemPtr%nli%chName,' = ', nlItemPtr%nli%chContent))

  end subroutine


  !******************************************************************
  !******************************************************************
  !
  !  Cleaning routines
  !
  !******************************************************************
  !******************************************************************


  subroutine clean_item(nlItem)    
    !
    ! Cleaning without destroying memory
    !
    implicit none

    ! imported parameter
    type(Tsilam_namelist_item), intent(inout) :: nlItem

    nlItem%chContent = ''
    nlItem%chName = ''

  end subroutine clean_item


  !******************************************************************

  subroutine clean_itemPtr(nlItemPtr)
    !
    ! Cleaning without destroying memory
    !
    implicit none

    ! imported parameter
    type(Tsilam_nl_item_ptr), intent(inout) :: nlItemPtr

    nlItemPtr%nli%chContent = ''
    nlItemPtr%nli%chName = ''

  end subroutine clean_itemPtr
  


  !******************************************************************

  subroutine clean_namelist(nl)
    !
    ! Cleaning without destroying memory
    !
    implicit none

    ! imported parameter
    type(Tsilam_namelist), intent(inout) :: nl

    ! Local variables
    integer :: i

    do i=1,nl%iLastFilledItem
      nl%items(i)%nli%chContent = ''
      nl%items(i)%nli%chName = ''
    end do

  end subroutine clean_namelist


  !******************************************************************

  subroutine clean_namelist_group(nlGrp)
    !
    ! Cleaning without destroying memory
    !
    implicit none

    ! imported parameter
    type(Tsilam_namelist_group), intent(inout) :: nlGrp

    ! Local variables
    integer :: i, iList

    do iList = 1, nlGrp%iLastFilledNL
      do i=1,nlGrp%lists(iList)%nl%iLastFilledItem
        nlGrp%lists(iList)%nl%items(i)%nli%chContent = ''
        nlGrp%lists(iList)%nl%items(i)%nli%chName = ''
      end do
    end do

  end subroutine clean_namelist_group


  !******************************************************************

  subroutine report_namelist (nl)
    !
    ! Just writes the report for the namelist
    !
    implicit none

    type(Tsilam_namelist), intent(in) :: nl

    call write_namelist(run_log_funit, nl,.true., .true.) ! header and public

  end subroutine report_namelist


  !******************************************************************

  logical function fu_namelist_defined(nl) 
    !
    ! Checks whether the namelist is reasonable
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nl

    if(.not.associated(nl))then
      fu_namelist_defined = .false.
    else
      fu_namelist_defined = nl%defined == silja_true
    endif

  end function fu_namelist_defined


  !******************************************************************

  logical function fu_namelist_empty(nl) 
    !
    ! Checks whether the namelist is reasonable
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nl

    if(.not.associated(nl))then
      fu_namelist_empty = .true.
    else
      fu_namelist_empty = nl%iLastFilledItem <= 0
    endif

  end function fu_namelist_empty


  !*******************************************************************

  integer function fu_nbr_of_items_of_namelist(nl)
    !
    ! Returns the number of filled items in the namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nl

    if(.not.associated(nl))then
      fu_nbr_of_items_of_namelist = 0
    else
      fu_nbr_of_items_of_namelist = nl%iLastFilledItem
    endif

  end function fu_nbr_of_items_of_namelist



!  type silam_namelist_item
!    private
!    character(len=clen) :: chName
!    character(len=fnlen) :: chContent
!  end type silam_namelist_item
!
!  type silam_nl_item_ptr
!    private
!    type(silam_namelist_item), pointer :: nli
!  end type silam_nl_item_ptr
!
!  !
!  ! Namelist itself
!  !
!  type Tsilam_namelist
!    private
!    type(silam_nl_item_ptr), dimension(:), pointer :: items
!    integer :: iLastFilledItem, iTotalSize
!    type(silja_logical) :: defined
!  end type Tsilam_namelist




END MODULE silam_namelist
