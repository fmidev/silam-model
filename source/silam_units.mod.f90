module silam_units
  !
  ! This module handles the units and unit conversion of SILAM. 
  ! The concept is: evey physical quantity stored as a SILAM field with field_id
  ! must have a unit, which follows the field transformations.
  !
  use toolbox

  implicit none
  
  ! Public routines of this module
  !
  public fu_conversion_factor
  public fu_unit_type
  public select_basic_unit
  !public fu_factor_to_basic_unit
  public fu_SI_unit
  public fu_set_named_value
  public set_named_value_and_unit
  public setNamedValue
  !
  ! Private routines of this module
  !
  private fu_conversion_factor_general
  private fu_factor_to_basic_unit_gen
  !
  ! Interface 
  !
  interface fu_conversion_factor
    module procedure fu_conversion_factor_general
  end interface

  !interface fu_factor_to_basic_unit
  !  module procedure fu_factor_to_basic_unit_gen
  !end interface

  
  ! Types of units used in the unit conversion routine
  !
  integer, public, parameter :: amount_unit = 10001 ! Mass, radioactivity
  integer, public, parameter :: time_unit = 10002  ! Just time
  integer, public, parameter :: size_unit = 10003  ! Just size
  integer, public, parameter :: fraction_unit = 10004  ! Just fraction of something, e.g. grid area
  integer, public, parameter :: temperature_unit = 10005  ! Temperature
  integer, public, parameter :: energy_unit = 10006  ! solar radiation, energy
  integer, public, parameter :: power_unit = 10007  ! a derivative of energy


  
  CONTAINS
  
  ! ******************************************************************

  real function fu_conversion_factor_general(chUnitFrom, chUnitTo)
    !
    ! This function makes the conversion from one unit to another.
    ! Idea: conversion between the units is performed in two
    ! steps. First, unitFrom is converted to some basic SI unit,
    ! second, the UnitTo is obtained from that basic SI unit.  Here we
    ! treat only conversions that does not need material.
    !
    implicit none

    ! Imported parameters with intent IN
    !
    character(len=*), intent(in) :: chUnitFrom, chUnitTo

    ! Local variables
    real :: tmpFactor
    character(len=clen) :: chSubUnit1, chSubUnit2, chSubUnit3, chSubUnit4, chBasicUnit
    integer :: unitType, iSlashFrom, iSlashTo

    if(chUnitFrom == "" .or. chUnitTo == "")then
      call set_error('Empty unit given','fu_conversion_factor_general')
      return
    endif
    tmpFactor = 1.
    iSlashFrom = index(chUnitFrom,"/")
    iSlashTo = index(chUnitTo,"/")
    if(iSlashFrom == 0)then ! [unit1] -> [unit3]
      if(iSlashTo == 0)then
        chSubUnit1 = chUnitFrom
        chSubUnit3 = chUnitTo
        chSubUnit2 =''  ! No denominators
        chSubUnit4 =''
      else
        call set_error('Incompatible units:' + chUnitFrom + '<->' + chUnitTo, &
                     & 'fu_conversion_factor_general')
        return
      endif
    else   ! [unit_1]/[unit_2] -> [unit_3]/[unit_4]
      if(iSlashTo == 0)then
        call set_error('Incompatible units' + chUnitFrom + '<->' + chUnitTo, &
                     & 'fu_conversion_factor_general')
        return
      else
        chSubUnit1 = chUnitFrom(1:iSlashFrom-1)
        chSubUnit3 = chUnitTo(1:iSlashTo-1)
        chSubUnit2 = chUnitFrom(iSlashFrom+1 : )
        chSubUnit4 = chUnitTo(iSlashTo+1 : )
      endif
    endif

    !------------------------------------------------------
    ! 
    ! Unit nominators: select basic unit and transform each given unit to basic one
    ! Reason for that is: sometimes properly selected basic unit will make the task 
    ! much simpler and, e.g., independent from the availability of material
    !
    call select_basic_unit(chSubUnit1, chSubUnit3, chBasicUnit)
    
    tmpFactor = tmpFactor * fu_factor_to_basic_unit_gen(chSubUnit1, chBasicUnit)    
    unitType = fu_unit_type(chSubUnit1)

    if(unitType /= fu_unit_type(chSubUnit3))then
      call set_error('Incompatible units:' + chUnitFrom + '_and_' + chUnitTo, &
                   & 'fu_conversion_factor_general')
      return
    endif

    tmpFactor = tmpFactor / fu_factor_to_basic_unit_gen(chSubUnit3, chBasicUnit)


    !------------------------------------------------------
    !
    ! Unit denominators:
    !
    if(iSlashFrom /= 0)then

      call select_basic_unit(chSubUnit2, chSubUnit4, chBasicUnit)

      tmpFactor = tmpFactor / fu_factor_to_basic_unit_gen(chSubUnit2, chBasicUnit)
      unitType = fu_unit_type(chSubUnit2)

      if(unitType /= fu_unit_type(chSubUnit4))then
        call set_error('Incompatible units:' + chUnitFrom + '_and_' + chUnitTo, &
                     & 'fu_conversion_factor_general')
        return
      endif

      tmpFactor = tmpFactor * fu_factor_to_basic_unit_gen(chSubUnit4, chBasicUnit)

    endif  ! Denominators exist

    if(error)then
      fu_conversion_factor_general = -1.
    else
      fu_conversion_factor_general = tmpFactor
    endif

  end function fu_conversion_factor_general


  !******************************************************************

  recursive integer function fu_unit_type(chUnit) result(unit_type)
    !
    ! Determines the type of the unit - amount or time (so far)
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chUnit

    !Local variables
    integer :: i, status, iLen
    character(len=1) :: chLast
    
    iLen  = len_trim(chUnit)
    unit_type = int_missing

    if(iLen == 0)then
      call set_error('Empty unit given','fu_unit_type')
      return

    elseif(iLen == 1)then
      chLast = chUnit
      if(chLast == 'm')then
        unit_type = size_unit
      elseif(chLast == 'g' .or. chLast == 't')then
        unit_type = amount_unit
      elseif(chLast == '%')then
        unit_type = fraction_unit
      elseif(chLast == 'W')then
        unit_type = power_unit
      elseif(chLast == 'J')then
        unit_type = energy_unit
      elseif(chLast == 'K' .or. chLast == 'C' .or. chLast == 'F')then
        unit_type = temperature_unit
      else
        call set_error('Unknown 1-char unit:' + chUnit, 'fu_unit_type')
        return
      endif

    elseif(iLen == 2)then
      chLast = chUnit(iLen:iLen) ! Get the last letter
      if(fu_if_digit(chLast))then
        unit_type = fu_unit_type(chUnit(1:len_trim(chUnit)-1))  ! something like m2 ?
      elseif(index(chUnit,'kg') == len_trim(chUnit)-1)then
        unit_type = amount_unit
      elseif(index(chUnit,'Bq') == len_trim(chUnit)-1)then
        unit_type = amount_unit
      elseif(trim(chUnit) == 'hr' .or. trim(chUnit) == 'yr')then
        unit_type = time_unit
      else
        unit_type = fu_unit_type(chUnit(2:2))  ! no ideas on 2-char unit, may be, "km"?
        if(unit_type == int_missing) &
                       & call set_error('Unknown 2-char unit:' + chUnit,'fu_unit_type')
      endif

    elseif(iLen == 3)then
      if(index(chUnit,'ton') == len_trim(chUnit)-2)then
        unit_type = amount_unit
      elseif(trim(chUnit) == 'sec' .or. trim(chUnit) == 'min' .or. &
         & trim(chUnit) == 'day' .or. trim(chUnit) == 'mon')then
        unit_type = time_unit
      else
        unit_type = fu_unit_type(chUnit(2:3))  ! no ideas on 3-char unit, may be, "kBq"?
        if(unit_type == int_missing) & 
                       & call set_error('Unknown 3-digit unit:' + chUnit,'fu_unit_type')
      endif

    elseif(iLen == 4)then
      if(index(chUnit,'mole') == len_trim(chUnit)-3)then
        unit_type = amount_unit
      elseif(trim(chUnit) == 'hour' .or. trim(chUnit) == 'year')then
        unit_type = time_unit
      else
        unit_type = fu_unit_type(chUnit(2:4))  ! no ideas on 4-char unit, may be, "kton"?
        if(unit_type == int_missing) &
                       & call set_error('Unknown 4-digit unit:' + chUnit,'fu_unit_type')
      endif

    elseif(iLen > 4)then

      if(trim(chUnit) == 'number')then
        unit_type = amount_unit

      elseif(trim(chUnit) == 'month')then
        unit_type = time_unit

      elseif(trim(chUnit) == 'fraction')then
        unit_type = fraction_unit
      else
        unit_type = fu_unit_type(chUnit(2:len_trim(chUnit)))  ! no ideas on x-char unit, may be, "k<whatever>"?
        if(unit_type == int_missing) &
                       & call set_error('Unknown >4 digit unit:' + chUnit,'fu_unit_type')
      endif

    else
      call set_error('Unknown unit:' + chUnit,'fu_unit_type')
    endif

  end function fu_unit_type


  !****************************************************************************

  subroutine select_basic_unit(chUnit1, chUnit2, chBasicUnit)
    !
    ! Takes care of choosing an optimal basic unit for the two ones given
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chUnit1, chUnit2
    character(len=*), intent(out) :: chBasicUnit

    chBasicUnit = ''
    if(fu_unit_type(chUnit1) /= fu_unit_type(chUnit2))then
      call set_error('Incompatible unit types','select_basic_unit')
      return
    endif

    select case(fu_unit_type(chUnit1))
      case(time_unit)
        chBasicUnit = 'sec'

      case(size_Unit)
        chBasicUnit = 'm'

      case(fraction_unit)
        chBasicUnit = 'fraction'

      case(temperature_unit)
        chBasicUnit = 'K'

      case(energy_unit)
        chBasicUnit = 'J'

      case(power_unit)
        chBasicUnit = 'W'

      case(amount_unit)
        !
        ! May be, we transfer kg to g or alike? Then basic unit is not mole
        ! We may have a few basic units: kg, number, Bq, and mole as default

        ! one unit number -> number
        ! one unit mole -> mole
        ! one unit Bq -> Bq
        ! else: both must be either ton or g, basic unit kg

        if (fu_if_ends_with(chUnit1, 'number') .or. fu_if_ends_with(chUnit2, 'number')) then
          chBasicUnit = 'number'
        else if (fu_if_ends_with(chUnit1, 'mole') .or. fu_if_ends_with(chUnit2, 'mole')) then
          chBasicUnit = 'mole'
        else if (fu_if_ends_with(chUnit1, 'Bq') .or. fu_if_ends_with(chUnit2, 'Bq')) then
          chBasicUnit = 'Bq'
        else
          if (.not. (fu_if_ends_with(chUnit1, 'g') .or. fu_if_ends_with(chUnit1, 'ton').or. fu_if_ends_with(chUnit1, 't')) .or. &
            & .not. (fu_if_ends_with(chUnit2, 'g') .or. fu_if_ends_with(chUnit2, 'ton').or. fu_if_ends_with(chUnit2, 't'))) then
            call set_error('Unable to select basic unit: ' + chUnit1 + ' -> ' + chUnit2, &
                         & 'select_basic_unit')
            return
          end if
          !chBasicUnit = 'kg' changed basic unit of mass to g to handle conversion of mkg to kg
          chBasicUnit = 'g'
        end if

      case default
        call set_error('Unknown unit type:' + chUnit1,'select_basic_unit')
        return
    end select ! Type of Unit1

  end subroutine select_basic_unit


  !****************************************************************************

  recursive real function fu_factor_to_basic_unit_gen(chUnit, chBasicUnit) &
                                         & result(factor_to_basic_unit)
    !
    ! Finds the conversion factor from some unit to the appropriate basic unit
    !
    implicit none
    character(len=*), intent(in) :: chUnit, chBasicUnit

    ! Local variables
    integer :: iPower, status
    real :: fTmp
    character(len=unitNmLen) :: chScale
    !
    ! The most-trivial case: unit and basic unit are the same
    !
    if(trim(chUnit) == trim(chBasicUnit))then
      factor_to_basic_unit = 1.0
      return
   endif    
    !
    ! Check first if the unit is in some power - 2, 3 or 4 are allowed. 
    !
    read(unit=chUnit(len_trim(chUnit):len_trim(chUnit)),fmt=*,iostat=status) iPower
    if(status == 0 .and. iPower > 1 .and. iPower < 5)then

      factor_to_basic_unit = fu_factor_to_basic_unit_gen(chUnit(1:len_trim(chUnit)-1), &
                                                       & chBasicUnit)
      !
      ! Take care of the power
      !
      fTmp = factor_to_basic_unit
      do status = 1, iPower-1
        factor_to_basic_unit = factor_to_basic_unit * fTmp
      end do
      return
    endif
    !
    ! Now we should have a basic unit. Be careful !!!!!!
    ! Every unit can be, in principle, trasnported to default basic unit (e.g. mole)
    ! In some cases, the "easier" basic units are used - like kg.
    !
!call msg(chUnit + ',' + chBasicUnit, index(trim(chUnit),trim(chBasicUnit),.true.))
    if(index(trim(chUnit),trim(chBasicUnit),.true.) > 0 .and. &
     & index(trim(chUnit),trim(chBasicUnit),.true.) == len_trim(chUnit)-len_trim(chBasicUnit)+1)then
      !
      ! Generic case: basic unit is a part of the unit. Deal with the case and return
      !
      chScale = chUnit(1:index(trim(chUnit),trim(chBasicUnit),.true.)-1)
      !call msg('chUnit, chScale, chBasicUnit: ' + chUnit + ',' + chScale + ',' + chBasicUnit)
      factor_to_basic_unit = 1.0

      if(trim(chScale) == 'p')then        ! The standard meaning of these scalings
        factor_to_basic_unit = 1.0e-12
      elseif(trim(chScale) == 'n')then
        factor_to_basic_unit = 1.0e-9
      elseif(trim(chScale) == 'mk')then
        factor_to_basic_unit = 1.0e-6
      elseif(trim(chScale) == 'u')then
        factor_to_basic_unit = 1.0e-6
      elseif(trim(chScale) == 'm')then
        factor_to_basic_unit = 1.0e-3
      elseif(trim(chScale) == 'c')then
        factor_to_basic_unit = 1.0e-2
      elseif(trim(chScale) == 'd')then
        factor_to_basic_unit = 1.0e-1
      elseif(trim(chScale) == 'k')then
        factor_to_basic_unit = 1.0e+3
      elseif(trim(chScale) == 'M')then
        factor_to_basic_unit = 1.0e+6
      elseif(trim(chScale) == 'G')then
        factor_to_basic_unit = 1.0e+9
      elseif(trim(chScale) == 'T')then
        factor_to_basic_unit = 1.0e+12
      else
        call set_error('Did not understand the unit scale:' + chUnit,'fu_factor_to_basic_unit_gen')
        return
      endif

      ! if(trim(chBasicUnit) == 'kg') then
      !    factor_to_basic_unit = factor_to_basic_unit * 1e-3
      ! end if
      
    else
      !
      ! More complicated case when the unit is a separate word, different from basic unit
      ! There are several cases where this can happen
      !
      select case(fu_unit_type(chBasicUnit))

        case(amount_unit)
         
          ! if(trim(chBasicUnit) == 'kg')then  ! kg has connection to ton
          !   !if(chUnit(len_trim(chUnit):len_trim(chUnit)) == 'g')then
          !   if (fu_if_ends_with(chUnit, 'g')) then
          !     factor_to_basic_unit = fu_factor_to_basic_unit_gen(chUnit, 'g') * 1.0e-3
          !   elseif(trim(chUnit) == 'ton' .or. trim(chUnit) == 't')then
          !     factor_to_basic_unit = 1.0e+3
          !   elseif(trim(chUnit) == 'kton' .or. trim(chUnit) == 'kt')then
          !     factor_to_basic_unit = 1.0e+6
          !   elseif(trim(chUnit) == 'Mton' .or. trim(chUnit) == 'Mt')then
          !     factor_to_basic_unit = 1.0e+9
          !   elseif(trim(chUnit) == 'Gton' .or. trim(chUnit) == 'Gt')then
          !     factor_to_basic_unit = 1.0e+12
          !   else
          !     call set_error('Unknown amount_unit:' + chUnit,'')
          !     return
          !   endif
          ! else
          !   call set_error('Unknown combination of amount unit:' + chUnit + ', and basic unit:' + chBasicUnit, &
          !                & 'fu_factor_to_basic_unit_gen')
          !   return
          ! endif
         
          if(trim(chBasicUnit) == 'g')then  ! kg has connection to ton
            !if(chUnit(len_trim(chUnit):len_trim(chUnit)) == 'g')then
            if (trim(chUnit) == 'kg') then
              factor_to_basic_unit = 1.0e+3
            elseif(trim(chUnit) == 'ton' .or. trim(chUnit) == 't')then
              factor_to_basic_unit = 1.0e+6
            elseif(trim(chUnit) == 'kton' .or. trim(chUnit) == 'kt')then
              factor_to_basic_unit = 1.0e+9
            elseif(trim(chUnit) == 'Mton' .or. trim(chUnit) == 'Mt')then
              factor_to_basic_unit = 1.0e+12
            elseif(trim(chUnit) == 'Gton' .or. trim(chUnit) == 'Gt')then
              factor_to_basic_unit = 1.0e+15
            else
              call set_error('Unknown amount_unit:' + chUnit,'')
              return
            endif
          else
            call set_error('Unknown combination of amount unit:' + chUnit + ', and basic unit:' + chBasicUnit, &
                         & 'fu_factor_to_basic_unit_gen')
            return
          endif
        
        case(time_unit)
          if(trim(chUnit) == 'sec')then
            factor_to_basic_unit = 1.
          elseif(trim(chUnit) == 'min')then
              factor_to_basic_unit = 60.
          elseif(trim(chUnit) == 'hr' .or. trim(chUnit) == 'hour')then
            factor_to_basic_unit = 3600.
          elseif(trim(chUnit) == 'day' .or. trim(chUnit) == 'd')then
            factor_to_basic_unit = 86400.
          elseif(trim(chUnit) == 'mon' .or. trim(chUnit) == 'month')then
            factor_to_basic_unit = 2635200.  ! 30.5 days is assumed in a month
          elseif(trim(chUnit) == 'yr' .or. trim(chUnit) == 'year')then
            factor_to_basic_unit = 31557600. ! 365.25 days per year is assumed
          else
            call set_error('Unknown time unit:' + chUnit,'fu_factor_to_basic_unit_gen')
            return
          endif

        case(fraction_unit)
          if(trim(chUnit) =='%')then
            factor_to_basic_unit = 0.01
          elseif(trim(chUnit) == 'fraction')then
            factor_to_basic_unit = 1.0
          else
            call set_error('Unknown time unit:' + chUnit,'fu_factor_to_basic_unit_gen')
          endif

        case default
          call set_error('Unknown type of the basic unit:' + chBasicUnit,'fu_factor_to_basic_unit_gen')
          return
      end select  ! unit type of basic unit

    endif  ! whether chBasicUnit is part of chUnit

  end function fu_factor_to_basic_unit_gen


  !***********************************************************************

  function fu_SI_unit(chUnit) result(chSI_unit)
    !
    ! Returns the SI unit corresponding to the given derived one
    !
    implicit none

    ! Return value
    character(len=10) :: chSI_unit

    ! Imported parameters
    character(len=*), intent(in) :: chUnit

    ! Local variable
    character(len=1) :: chPower

    if(len_trim(chUnit) == 0)then
      call set_error('Empty unit given','fu_SI_unit')
      return
    endif
    
    chPower = chUnit(len_trim(chUnit):len_trim(chUnit))
    if(.not. fu_if_digit(chPower)) chPower = ''
    
    select case(fu_unit_type(chUnit))
      case(time_unit)
        chSI_unit = 'sec'

      case(size_Unit)
        chSI_unit = 'm'

      case(fraction_unit)
        chSI_unit = 'fraction'

      case(temperature_unit)
        chSI_unit = 'K'

      case(amount_unit)
        !
        ! May be, we transfer kg to g or alike? Then basic unit is not mole
        ! We may have a few basic units: kg, number, Bq, and mole as default
        !
        select case(chUnit)
          case('mkmole','mmole','mole','kmole','Mmole','Gmole')
            chSI_unit = 'mole'

          case('number')
            chSI_unit = 'number'

          case('Bq','kBq','MBq')
            chSI_unit = 'Bq'

          case('pg','ng','mkg', 'ug', 'mg','g','kg','Mg','Gg','ton','kton','Mton','Gton','t','kt','Mt','Gt')
            chSI_unit = 'kg'

          case default
            call set_error(fu_connect_strings('Unknown unit:',chUnit),'fu_SI_unit')
        end select  ! value of Unit1

      case default
        call set_error(fu_connect_strings('Unknown unit:',chUnit),'fu_SI_unit')
    end select ! Type of Unit1

    chSI_unit = chSI_unit + chPower

  end function fu_SI_unit


  !***********************************************************************

  real function fu_set_named_value(chInputLine, ifSilent)
    !
    ! Get a line, such as "10 kg" or "1 mkm" and returns the given vaue in SI unit. 
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chInputLine  ! number and unit
    logical, intent(in), optional :: ifSilent

    ! Local variables
    real :: fIn, factor
    integer :: iStat
    character(len=clen) :: chGivenUnit
    !
    ! Just read the line and, find out the conversion factor and return the result
    !
    fu_set_named_value = real_missing
    
    if(len_trim(chInputLine) == 0)then  ! the most-frequent error situation
      if(present(ifSilent))then
        if(ifSilent)return
      endif
      call set_error('Empty line is given','fu_set_named_value')
      return
    endif

    read(unit=chInputLine,fmt=*,iostat=iStat)fIn   ! get the value
    if(iStat /= 0)then
      if(present(ifSilent))then
        if(ifSilent)return
      endif
      call set_error('Failed the line:' + chInputLine,'fu_set_named_value')
      return
    endif

    chGivenUnit = chInputLine(index(adjustl(chInputLine),' ')+1:)
    chGivenUnit = adjustl(chGivenUnit)

    if(index(chGivenUnit,'/') == 0)then
      factor = fu_conversion_factor(chGivenUnit, fu_SI_unit(chGivenUnit))
    else
      factor = fu_conversion_factor(chGivenUnit(1:index(chGivenUnit,'/')-1), &
                                  & fu_SI_unit(chGivenUnit(1:index(chGivenUnit,'/')-1))) / &
             & fu_conversion_factor(chGivenUnit(index(chGivenUnit,'/')+1:), &
                                  & fu_SI_unit(chGivenUnit(index(chGivenUnit,'/')+1:)))
    endif
    if(error .or. (factor .eps. real_missing)) then
         call msg("chGivenUnit : '" // chGivenUnit // "'") 
         call msg("fIn",fIn)
         call set_error("Failed with line '"//trim(chInputLine)//"'","fu_set_named_value")
         return
    endif

    fu_set_named_value = fIn * factor

  end function fu_set_named_value


  !***********************************************************************
  
  subroutine set_named_value_and_unit(chInputLine, fValue, chUnit)
    !
    ! Similar to the fu_set_named_value but returns both the value and the basic unit
    !
    implicit none
    
    ! Improted parameters
    character(len=*), intent(in) :: chInputLine
    character(len=*), intent(out) :: chUnit
    real, intent(out) :: fValue
    
    ! Local variables
    real :: factor
    integer :: iStat
    
    fValue = real_missing
    chUnit = ''

    read(unit=chInputLine,fmt=*,iostat=iStat)fValue
    if(iStat /= 0)then
      call set_error(fu_connect_strings('Failed the line:',chInputLine),'set_named_value_and_unit')
      return
    endif

    chUnit = chInputLine(index(chInputLine,' ')+1:)
    chUnit = adjustl(chUnit)

    if(index(chUnit,'/') == 0)then
      factor = fu_conversion_factor(chUnit, fu_SI_unit(chUnit))
      chUnit = fu_SI_unit(chUnit)
    else
      factor = fu_conversion_factor(chUnit(1:index(chUnit,'/')-1), &
                                  & fu_SI_unit(chUnit(1:index(chUnit,'/')-1))) / &
             & fu_conversion_factor(chUnit(index(chUnit,'/')+1:), &
                                  & fu_SI_unit(chUnit(index(chUnit,'/')+1:)))
      chUnit = fu_SI_unit(chUnit(1:index(chUnit,'/')-1)) + '/' + fu_SI_unit(chUnit(index(chUnit,'/')+1:))
    endif
    if(error .or. (factor .eps. real_missing))return

    fValue = fValue * factor

  end subroutine set_named_value_and_unit
 

  !*******************************************************************

  subroutine setNamedValue(chValueStr, chUnitTo, fValue)
    !
    ! The routine makes a simple unit conversion for the value presented as a string
    ! with its own unit. Contrary to the corresponding routine in module materials, 
    ! here we apply only decimal sub-unit conversions.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chValueStr, chUnitTo
    real, intent(out) :: fValue

    ! Local variables
    integer :: iStatus
    character(len=20) :: chUnitFrom

    !
    ! First, let's get the value
    !
    read(unit=chValueStr,fmt=*,iostat=iStatus)fValue
    if(iStatus /= 0)then
      call set_error(fu_connect_strings('Failed to read the value from the string:', &
                                      & chValueStr),'setNamedValue')
      fValue = real_missing
      return
    endif
    !
    ! Now, let's get the unit and make conversion
    !
    chUnitFrom = chValueStr(index(adjustl(chValueStr),' ')+1 :)
    chUnitFrom = adjustl(chUnitFrom)

    fValue = fValue * fu_conversion_factor(chUnitFrom, chUnitTo)
    if(error)then
      fValue = real_missing
      return
    endif

  end subroutine setNamedValue

end module silam_units
