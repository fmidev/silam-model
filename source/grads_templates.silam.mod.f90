!----------------------------------------------------------------------------
!
! This module contains all necessary objects and procedures for operating with
! GrADS-type templates in strings (actually, file names). 
! Two-way operations are possible: decode template string to the template
! structure with further possibility to generate an appropriate file name; and 
! an inverse task (tricky, as usual) - to compute the template string from
! a set of file names
!
! Language: FORTRAN-90
!
! Author: M.Sofiev
!
module grads_templates

use silam_times

implicit none

public MakeTemplate ! Main subroutine for computing template string
public decode_template_string ! decodes template string to template strucure

public sort_times_ascend ! Inside time_file_set

public fu_type
public fu_start_pos
public fu_start_pos_ar
public fu_n_items
public fu_n_times
public fu_item
public fu_collection
public set_collection
public fu_anal_time
public fu_valid_time
public fu_forec_len
public fu_anal_time_vary
public fu_forec_len_vary
public fu_template_string  ! e.g. %y2
public fu_pure_grads_template ! Our template class is more rich than pure GrADS...
public fu_same_template_start_time
public fu_template_timestep   !! 

public fu_FNm ! From tf_set

public add_time_file_to_tf_set
public set_what_vary
public get_what_vary

public fu_silam_template_item ! Either %case or %source
public fu_silam_source_template_item ! %source as the output

public fu_forecast_length_templ_str
public fu_valid_time_hour_templ_str
public fu_valid_time_day_templ_str
public fu_valid_time_month_templ_str
public fu_valid_time_year_templ_str
public fu_anal_time_hour_templ_str
public fu_anal_time_day_templ_str
public fu_anal_time_month_templ_str
public fu_anal_time_year_templ_str

public defined
public report

private string_from_template
private check_string
private fu_master_items
private fu_acronym
private fu_item_maxlength
private compress_template_structure
private fu_FNm_from_template ! Create a file name having template and all times
private fu_template_string_for_item
private fu_template_string_for_template

private fu_item_type_of_template_item
private fu_start_pos_of_template_item
private fu_start_pos_ar_of_templ_item
private fu_n_items_of_template
private fu_item_of_template
private fu_template_collection
private set_template_collection

private fu_analysis_time_of_tf_set
private fu_valid_time_of_tf_set
private fu_forec_len_of_tf_set
private fu_anal_time_of_tf_set_vary
private fu_forec_len_of_tf_set_vary
private fu_FNm_of_tf_set
private fu_n_times_of_tf_set
!private exchange_times_of_tf_set
! Added by Jukkis 08.04.2008

private fu_compare_templates_eq  ! Interface for operator ==
private fu_compare_template_items_eq  ! Interface for operator ==
private grads_template_defined
private report_grads_template

interface fu_type
  module procedure fu_item_type_of_template_item
end interface
interface fu_start_pos
  module procedure fu_start_pos_of_template_item
end interface
interface fu_start_pos_ar
  module procedure fu_start_pos_ar_of_templ_item
end interface
interface fu_n_items
  module procedure fu_n_items_of_template
  module procedure fu_n_times_of_tf_set
end interface
interface fu_n_times
  module procedure fu_n_times_of_tf_set
end interface
interface fu_item
  module procedure fu_item_of_template
end interface
interface fu_collection
  module procedure fu_template_collection
end interface
interface set_collection
  module procedure set_template_collection
end interface
interface fu_anal_time
  module procedure fu_analysis_time_of_tf_set
end interface
interface fu_valid_time
  module procedure fu_valid_time_of_tf_set
end interface
interface fu_forec_len
  module procedure fu_forec_len_of_tf_set
end interface
interface fu_anal_time_vary
  module procedure fu_anal_time_of_tf_set_vary
end interface
interface fu_forec_len_vary
  module procedure fu_forec_len_of_tf_set_vary
end interface
interface fu_FNm
  module procedure fu_FNm_of_tf_set
  module procedure fu_FNm_from_template
end interface
interface fu_template_string
  module procedure fu_template_string_for_item
  module procedure fu_template_string_for_template
end interface

interface operator (==)
  module procedure fu_compare_templates_eq
  module procedure fu_compare_template_items_eq
end interface

interface defined
  module procedure grads_template_defined
end interface

interface report
  module procedure report_grads_template
end interface

integer, private,parameter :: max_pos_of_template = 20
integer, private,parameter :: max_nbr_of_templates = 40
integer, private,parameter :: max_nbr_of_times= 1000


! The main template-defining structure
!
type template_item
  private
  integer :: item_type, min_len, max_len ! Type, minimum and maximum length of the string
  integer, dimension(max_pos_of_template) :: start_pos=0 ! There may be several similar items
end type template_item

type (template_item), public, parameter :: template_item_missing = &
        &template_item(int_missing, int_missing, int_missing, 0)

type grads_template
  private
  integer :: n_items=0
  type(template_item), dimension(max_nbr_of_templates) :: items
  character(len=fnlen) :: ini_string = '' ! Initial string with template acronyms
  character(len=fnlen) :: templ_collection='' ! Final string of the templates
end type grads_template

type (grads_template), public, parameter ::  grads_template_missing = &
                  &  grads_template(0, template_item_missing, '', '')

!
! Input for the template-making routines: a set of times and file names
!
type time_file_set
  private
  type(silja_time), dimension(max_nbr_of_times):: analysis_times ! Reference(initialization) times
  type(silja_interval), dimension(max_nbr_of_times)::forec_len ! forecast length
  character(len=fnlen),dimension(max_nbr_of_times):: FNms
  integer :: n_times=0
  logical :: analysis_time_vary, forecast_length_vary
end type time_file_set

type ( time_file_set), public, parameter :: time_file_set_missing = &
              & time_file_set(time_missing, interval_missing, '', 0, .false.,.false.)


! Different types of the template items. Not all of them are allowed - 
! stuff like Th1 imply 1-2 characters length, which is forbidden. Length
! of the file names must be constant.
!
! Template parameters are defined via its type, min and max possible length
! AFTER substitution and GrADS string denoting the template
!
type template_param
  private
  integer::tType, tLenMin,tLenMax ! type and length AFTER substitution min/max
  character(len=7) :: tString ! GrADS format of a particular template
end type template_param

! All possible template types are below. 
! 
! ATTENTION. A trick.
!
! GRIB has analysis time as an independent variable, which is below treated
! via Ta* templates. But GrADS does not have such template. As a result, file
! with analysis time in the name and non-zero forecast length can never be
! catched by template options. The way out is: ignore GrADS forecast and assume
! that analysis time is valid time with corresponding templates. It leads to 
! time shift for the forecast length, but if it is fixed for the given file set
! it is not too bad (however, must not be forgotten). Algorithm below is:
! - find all possible templates, including Ta*.
! - if analysis time varies AND Ta* templates found (which means that file names
!   contain varying parts described by Ta* templates only) - create a time shift
!   for the forecast_length, which would make valid_time = analysis_time and
!   enable to use valid_time templates.
!
! A TRICK. Last two templates are made to support SILAM case name and source name
!          However, they must never appear in other work. 

integer,private,parameter :: Ty2 =          1001, Iy2 = 1
integer,private,parameter :: Ty4 =          1002, Iy4 = 2
integer,private,parameter :: Tm1 =          1003, Im1 = 3
integer,private,parameter :: Tm2 =          1004, Im2 = 4
integer,private,parameter :: Tmc =          1005, Imc = 5
integer,private,parameter :: Tmc_capital =  1006, Imc_capital = 6
integer,private,parameter :: Td1 =          1007, Id1 = 7
integer,private,parameter :: Td2 =          1008, Id2 = 8
integer,private,parameter :: Th1 =          1009, Ih1 = 9
integer,private,parameter :: Th2 =          1010, Ih2 = 10
integer,private,parameter :: Th3 =          1011, Ih3 = 11
integer,private,parameter :: Tf2 =          1012, If2 = 12
integer,private,parameter :: Tf3 =          1013, If3 = 13 ! Forecast time hour
integer,private,parameter :: Tn2 =          1014, In2 = 14  ! Forecast time minute
integer,private,parameter :: Tiy2 =         1015, Iiy2 = 15
integer,private,parameter :: Tiy4 =         1016, Iiy4 = 16
integer,private,parameter :: Tim1 =         1017, Iim1 = 17
integer,private,parameter :: Tim2 =         1018, Iim2 = 18 !!!month of ini time
integer,private,parameter :: Tin2 =         1019, Iin2 = 19 !!!minuts of ini time
integer,private,parameter :: Timc =         1020, Iimc = 20
integer,private,parameter :: Timc_capital = 1021, Iimc_capital = 21
integer,private,parameter :: Tid1 =         1022, Iid1 = 22
integer,private,parameter :: Tid2 =         1023, Iid2 = 23
integer,private,parameter :: Tih1 =         1024, Iih1 = 24
integer,private,parameter :: Tih2 =         1025, Iih2 = 25
integer,private,parameter :: Tih3 =         1026, Iih3 = 26
integer,private,parameter :: Tay2 =         1027, Iay2 = 27
integer,private,parameter :: Tay4 =         1028, Iay4 = 28
integer,private,parameter :: Tam1 =         1029, Iam1 = 29 
integer,private,parameter :: Tam2 =         1030, Iam2 = 30  !!analysis month
integer,private,parameter :: Tan2 =         1031, Ian2 = 31
integer,private,parameter :: Tamc =         1032, Iamc = 32  !!
integer,private,parameter :: Tamc_capital = 1033, Iamc_capital = 33
integer,private,parameter :: Tad1 =         1034, Iad1 = 34
integer,private,parameter :: Tad2 =         1035, Iad2 = 35
integer,private,parameter :: Tah1 =         1036, Iah1 = 36
integer,private,parameter :: Tah2 =         1037, Iah2 = 37
integer,private,parameter :: Tah3 =         1038, Iah3 = 38
integer,private,parameter :: Tcase =        1039, Icase = 39
integer,private,parameter :: Tsource =      1040, Isource = 40
integer,private,parameter :: Tclock =       1041, Iclock = 41 !!Wallclock
integer,private,parameter :: Tjd3 =         1042, Ijd3 = 42     ! Julian days
integer,private,parameter :: Tijd3 =        1043, Iijd3 = 43
integer,private,parameter :: Tajd3 =        1044, Iajd3 = 44
integer,private,parameter :: Ttask =        1055, Itask = 45
integer,private,parameter :: Tfs2 =         1056, Ifs2 = 46 ! Forecast time second Only two chars (minutes and hours present)
integer,private,parameter :: Tfm2 =         1057, Ifm2 = 47 ! Forecast time minute Only two chars (hour present)
integer,private,parameter :: Tfh2 =         1058, Ifh2 = 48 ! Forecast time hour (day persent) Only two chars
integer,private,parameter :: Tfd2 =         1059, Ifd2 = 49 ! Forecast time day  at least two chars

integer, private, parameter ::  nbr_of_templ_params = 49

type(template_param),dimension(nbr_of_templ_params),private,parameter :: &
  & template_values = &
       !               Template length??  
    & (/template_param (Ty2, 2,2, '%y2'), &
      & template_param (Ty4, 4,4, '%y4'), &
      & template_param (Tm1, 1,2, '%m1'), &
      & template_param (Tm2, 2,2, '%m2'), &
      & template_param (Tmc, 3,3, '%mc'), &
      & template_param (Tmc_capital, 3,3, '%MC'), &
      & template_param (Td1, 1,2, '%d1'), &
      & template_param (Td2, 2,2, '%d2'), &
      & template_param (Th1, 1,2, '%h1'), &
      & template_param (Th2, 2,2, '%h2'), &
      & template_param (Th3, 3,3, '%h3'), &
      & template_param (Tf2, 2,3, '%f2'), &
      & template_param (Tf3, 3,3, '%f3'), &
      & template_param (Tn2, 2,2, '%n2'), &
      & template_param (Tiy2, 2,2, '%iy2'), &
      & template_param (Tiy4, 4,4, '%iy4'), &
      & template_param (Tim1, 1,2, '%im1'), &
      & template_param (Tim2, 2,2, '%im2'), &
      & template_param (Tin2, 2,2, '%in2'), &
      & template_param (Timc, 3,3, '%imc'), &
      & template_param (Timc_capital, 3,3, '%iMC'), &
      & template_param (Tid1, 1,2, '%id1'), &
      & template_param (Tid2, 2,2, '%id2'), &
      & template_param (Tih1, 1,2, '%ih1'), &
      & template_param (Tih2, 2,2, '%ih2'), &
      & template_param (Tih3, 3,3, '%ih3'), &
      & template_param (Tay2, 2,2, '%ay2'), &
      & template_param (Tay4, 4,4, '%ay4'), &
      & template_param (Tam1, 1,2, '%am1'), &
      & template_param (Tam2, 2,2, '%am2'), &
      & template_param (Tan2, 2,2, '%an2'), &
      & template_param (Tamc, 3,3, '%amc'), &
      & template_param (Tamc_capital, 3,3, '%aMC'), &
      & template_param (Tad1, 1,2, '%ad1'), &
      & template_param (Tad2, 2,2, '%ad2'), &
      & template_param (Tah1, 1,2, '%ah1'), &
      & template_param (Tah2, 2,2, '%ah2'), &
      & template_param (Tah3, 3,3, '%ah3'), &
      & template_param (Tcase, 1,clen, '%case'), &
      & template_param (Tsource, 1,clen, '%source'), &
      & template_param (Tclock, 19,19, '%clock'), &  ! String made from wall clock time
      & template_param (Tjd3, 3,3, '%jd3'), &
      & template_param (Tijd3, 3,3, '%ijd3'), &
      & template_param (Tajd3, 3,3, '%ajd3'), &
      & template_param (Ttask, 3,3, '%task'), &
      & template_param (Tfs2, 2,2, '%fs2'), &
      & template_param (Tfm2, 2,2, '%fm2'), &
      & template_param (Tfh2, 2,2, '%fh2'), &
      & template_param (Tfd2, 2,2, '%fd2') &
      /)

type(template_item), parameter, public :: item_missing = template_item(&
          & int_missing, int_missing, int_missing, max_pos_of_template*int_missing)

type(grads_template),parameter, public :: template_missing = &
           & grads_template(0,item_missing,'','')

CONTAINS


!******************************************************************************

subroutine MakeTemplate(times_files, template, ifAnalysisTimeAllowed)
  !
  ! Main subroutine, which drives filling-in the template structure
  ! for a given arrays of times and corresponding to file names
  ! Algorithm:
  ! After a few checkings - all templates one-by-one are checked for
  ! existence in the file names. So stupid...
  ! A trick - template item can be in several places - e.g., names of 
  ! directories, etc. So, there is an array of starting points for each 
  ! template item.
  ! For each item type the first time and file name determine all
  ! possible locations of the template, while all others may only filter
  ! out some of them.
  !
  ! Author: M.Sofiev
  !
  implicit none

  ! Imported parameters
  logical, intent(in) :: ifAnalysisTimeAllowed ! GrADS does not have this time....
  type(time_file_set),intent(inout),target::times_files ! stands for times_files
  type(grads_template),intent(out)::template

  ! Local variables
  integer :: i,t,yr,mon,day,hr,min, iItem, iMaster,iSlave,iMPos, iSPos
  real :: sec
  integer, dimension(max_pos_of_template) :: indexTmp
  integer, dimension(5) :: tParts ! Year-mon-day-hour-minute
  integer, dimension(max_nbr_of_templates,2) :: positions ! Positions of the template items
  type(time_file_set),pointer::t_f ! stands for times_files
  character(len=10) templ_str, templ_str_chk ! String from time substituted into template
  logical :: ifVarying

  t_f => times_files
  if(t_f%n_times < 2)return

  ! Check the length of the file names and nullify the template structure
  !
  iItem = 1
  do while (len_trim(t_f%FNms(iItem)) < 2)
    iItem = iItem + 1
    if(iItem > t_f%n_times)then
      call set_error('All file names are empty','MakeTemplate')
      return
    end if
  end do
  do t=1,t_f%n_times
    t_f%FNms(t) = trim(adjustl(t_f%FNms(t)))
    if(len_trim(t_f%FNms(t)) > 1)then
      if(len_trim(t_f%FNms(iItem)) /= len_trim(t_f%FNms(t)))then
        call set_error('Varying length of the file names','MakeTemplate')
        return
      end if
    end if
  end do

  template%n_items=0
  do i=1,size(template%items)
    template%items(i)%item_type = int_missing
    template%items(i)%start_pos = int_missing
    template%items(i)%min_len = 999
    template%items(i)%max_len = 0
  end do

  ! Scan all supporting template items one-by-one in a cycle. Those, which
  ! are found are stored into a template structure with their suspected 
  ! positions.
  ! Steps are: get time, then create a template string as if this time is
  ! substituted to a particular template item and then check where this item
  ! can be found in the particular string
  ! Obs. There may an added analysis time with empty file name - see 
  ! arrange_times for details. So, we have to take care of possibly missing
  ! files allowing times with no files at all. GRIB allows it, so no problem
  !
  do iItem = 1, SIZE(template_values)

    if(fu_silam_template_item(template_values(iItem)%tType)) cycle

    if(.not.ifAnalysisTimeAllowed)then
      if(template_values(iItem)%tType == Tay2 .or. &
       & template_values(iItem)%tType == Tay4  .or. &
       & template_values(iItem)%tType == Tam1  .or. &
       & template_values(iItem)%tType == Tam2  .or. &
       & template_values(iItem)%tType == Tan2  .or. &
       & template_values(iItem)%tType == Tamc  .or. &
       & template_values(iItem)%tType == Tamc_capital  .or. &
       & template_values(iItem)%tType == Tad1  .or. &
       & template_values(iItem)%tType == Tad2  .or. &
       & template_values(iItem)%tType == Tajd3  .or. &
       & template_values(iItem)%tType == Tah1  .or. &
       & template_values(iItem)%tType == Tah2  .or. &
       & template_values(iItem)%tType == Tah3) cycle
    endif

    indexTmp = int_missing ! Temporary index of template string positions

    do i=1,t_f%n_times ! Cycle through all times and related file names
      call string_from_template(t_f%analysis_times(i),t_f%analysis_times(i), &
                              & t_f%forec_len(i), &
                              & template_values(iItem)%tType,templ_str)
      if(error)then
        call unset_error('MakeTemplate')
        exit
      end if

      if(len_trim(t_f%FNms(i)) > 1)then ! Skip empty file names
        call check_string(templ_str,t_f%FNms(i), indexTmp)
        if(ALL(indexTmp == int_missing)) exit ! Nothing for one time => void item
      end if
    end do

    if(ANY(indexTmp /= int_missing))then ! something found for this template item
!      print *,'Template',template%n_items+1,template_values(iItem)%tString
      template%n_items = template%n_items + 1
      if(template%n_items  > max_nbr_of_templates)then
        call set_error('Too many templates found','MakeTemplate')
      end if
      template%items(template%n_items)%item_type = template_values(iItem)%tType
      template%items(template%n_items)%min_len = min(len_trim(templ_str), &
                                       & template%items(template%n_items)%min_len)
      template%items(template%n_items)%max_len = max(len_trim(templ_str), &
                                       & template%items(template%n_items)%min_len)
      t=1
      do i=1,max_pos_of_template ! Store positions of the template item
        if(indexTmp(i) /= int_missing)then
          template%items(template%n_items)%start_pos(t)=indexTmp(i)
          t=t+1
        end if
      end do
    end if

  end do ! cycle through template items

  if(template%n_items == 0) then
    call set_error('Failed to find any template','MakeTemplate')
    return
  end if

  ! There are a few items in a hierarchy - e.g., Ty2 depends on Ty4. Depends 
  ! means: if Ty4 exists than Ty2 will be found too.
  ! This has to be checked. However, the elimination should be done
  ! only for those positions of slave template items, which overlap with the 
  ! master ones.
  !
  do iSlave=1,template%n_items
    !
    ! Check if there is a master item for this "slave" one
    !
    if(all(fu_master_items(template%items(iSlave)%item_type) == int_missing)) cycle
    !
    ! Check if the master item is in the list
    !
    do iMaster=1,template%n_items 
      if(any(template%items(iMaster)%item_type == fu_master_items( &
                               & template%items(iSlave)%item_type)))then
        !
        ! Cycling over master positions and delete all overlapping slaves
        ! i.e. slave item starting at or after master starts but before
        ! the master ends
        !
        do iMPos = 1,size(template%items(iMaster)%start_pos)
          if(template%items(iMaster)%start_pos(iMPos) == int_missing) cycle
          do iSPos = 1,size(template%items(iSlave)%start_pos)
            if(template%items(iSlave)%start_pos(iSPos) == int_missing) cycle
            if(template%items(iSlave)%start_pos(iSPos) >= & 
             & template%items(iMaster)%start_pos(iMPos) .and. &
             & template%items(iSlave)%start_pos(iSPos) <= &
             & template%items(iMaster)%start_pos(iMPos)+ &
                    & fu_item_maxlength(template%items(iMaster)%item_type))then
              template%items(iSlave)%start_pos(iSPos) = int_missing
            end if
          end do  ! cycle through slave item
        end do  ! cycle through master item
        !
        ! If there are no valid positions of the slave item - delete it
        !
        if(all(template%items(iSlave)%start_pos(:) == int_missing))then
          template%items(iSlave)%item_type = int_missing !Delete item
          exit
        end if
      end if
    end do
  end do
  !
  ! Now there are holes in the template list - thanks for deleted items
  ! Let's compress the structure. Below algorithm mixes-up the templates, but
  ! their order is anyway not important because does not reflect actual 
  ! positions in the file names.
  !
  call compress_template_structure(template)


!  print *, 'After slave-master cleaning: ',template%n_items,'  templates:'
!  do iItem = 1, template%n_items
!    print *,'Template ',fu_acronym(template%items(iItem)%item_type)
!  end do

  ! Previously, we have forbidden varying ;ength of the file names. So, all
  ! templates with varying string length must be excluded. Note - if some
  ! template in principle can be varying but in current case it does not - 
  ! it can stay
  !
  do iItem = 1, template%n_items

    i = 1 ! Skip empty file names
    do while (len_trim(t_f%FNms(i)) < 2)
      i = i + 1
    end do

    call string_from_template(t_f%analysis_times(i),t_f%analysis_times(i), &
                            & t_f%forec_len(i), &
                            & template%items(iItem)%item_type,templ_str_chk)
    if(error)return

    do t=i+1,t_f%n_times
      call string_from_template(t_f%analysis_times(t),t_f%analysis_times(t), &
                              & t_f%forec_len(t), &
                              & template%items(iItem)%item_type,templ_str)
      if(error)return
      if(len_trim(templ_str) /= len_trim(templ_str_chk)) then ! Varying string
        print *, 'Killing ', template%items(iItem)%item_type
        template%items(iItem)%item_type = int_missing
        exit
      end if
    end do

  end do
  
  call compress_template_structure(template)

  !
  ! We do not need to substitute all parts of the file name with templates. 
  ! For example, if some template happens to coinside with some part of the 
  ! name. There is actually no check if such coinsidence is actual (reflects 
  ! time dimension) unless this part is varying. So, if some template
  ! is not varying for given set of times - it should be removed.
  !
  do iItem = 1, template%n_items

    i = 1 ! Skip empty file names
    do while (len_trim(t_f%FNms(i)) < 2)
      i = i + 1
    end do

    call string_from_template(t_f%analysis_times(i),t_f%analysis_times(i), &
                            & t_f%forec_len(i), &
                            & template%items(iItem)%item_type,templ_str_chk)
    if(error)return

	ifVarying = .false.

    do t=i+1,t_f%n_times
      call string_from_template(t_f%analysis_times(t),t_f%analysis_times(t), &
                              & t_f%forec_len(t), &
                              & template%items(iItem)%item_type,templ_str)
      if(error)return
      if(trim(templ_str) /= trim(templ_str_chk)) then
        ifVarying = .true.
        exit 
      end if
    end do

    if(.not.ifVarying) template%items(iItem)%item_type = int_missing
  end do

  call compress_template_structure(template)

  !
  ! There may be cases when one template exactly coindsides with the other.
  ! It would mean that its start and end positions are exactly the same.
  ! This is an ambiguity, which can not be resolved - we have to just 
  ! delete one of them arbitrarily.
  !
  ifVarying = .false.
  do iMaster = 1, template%n_items
    do iMPos = 1, size(template%items(iMaster)%start_pos)
      if(template%items(iMaster)%start_pos(iMPos) == int_missing) cycle
      do iSlave = iMaster+1, template%n_items
        do iSPos = 1, size(template%items(iSlave)%start_pos)
          if(template%items(iSlave)%start_pos(iSPos) == int_missing) cycle
          if(template%items(iMaster)%start_pos(iMPos) == &
                          & template%items(iSlave)%start_pos(iSPos) .and. &
           & template%items(iMaster)%start_pos(iMPos) + template%items(iMaster)%max_len == &
           & template%items(iSlave)%start_pos(iSPos) + template%items(iSlave)%max_len)then
            call msg_warning('Coinsiding templates','MakeTemplate')
            print *,'Templates ',trim(fu_acronym(template%items(iMaster)%item_type)),', ', &
                               & trim(fu_acronym(template%items(iSlave)%item_type)), &
                               & ' seem to coinside. The last one is deleted'
            ifVarying = .true.
            template%items(iSlave)%start_pos(iSPos) = int_missing
            if(all(template%items(iSlave)%start_pos(:) == int_missing)) &
                            & template%items(iSlave)%item_type = int_missing !Delete item
          end if
        end do  ! iSPos
      end do  ! iSlave
    end do  ! iMPos
  end do  ! iMaster

  if(ifVarying)then ! If the template set was altered due to coinsiding items
    print *
    print *,'ATTENTION !'
    print *,'ATTENTION !'
    print *,'There were templates removed due to exact coinsidence with the others.'
    print *,'The choice what to delete was arbitrary! GrADS will work, but meaning'
    print *,'of the template string may be different from the reality.'
    print *
  end if

  call compress_template_structure(template)

  if(error)return

!  print *, 'After varying filter: ',template%n_items,'  templates:'
  print '(A,I3,A)', ' Finally, ',template%n_items,' template(s):'
  do iItem = 1, template%n_items
    print *,'Template ',fu_acronym(template%items(iItem)%item_type)
  end do

  !
  ! The last check - there must be no overlap of template items.
  !
  do iMaster=1,template%n_items
    do iSlave=iMaster+1,template%n_items 
      !
      ! Cycling over "master" positions and check all "slaves".
      ! Three types of overlapping:
      ! - master's begin is inside slave
      ! - master's end is insied slave
      ! - complete slave is inside master (<=> slave's start is inside master)
      !
      do iMPos = 1,size(template%items(iMaster)%start_pos)
        if(template%items(iMaster)%start_pos(iMPos) == int_missing) cycle
        do iSPos = 1,size(template%items(iSlave)%start_pos)
          if(template%items(iSlave)%start_pos(iSPos) == int_missing) cycle

          ! Option 1
          if(template%items(iMaster)%start_pos(iMPos) >= &       
              & template%items(iSlave)%start_pos(iSPos) .and. &
              & template%items(iMaster)%start_pos(iMPos) <= &
              & template%items(iSlave)%start_pos(iSPos)+template%items(iSlave)%max_len-1)then
            call set_error('Overlapping templates. Sorry','MakeTemplate')
            return
          end if

          ! Option 2
          if(template%items(iMaster)%start_pos(iMPos)+template%items(iMaster)%max_len-1 >= &       
              & template%items(iSlave)%start_pos(iSPos) .and. &
              & template%items(iMaster)%start_pos(iMPos)+template%items(iMaster)%max_len-1 <= &
              & template%items(iSlave)%start_pos(iSPos)+template%items(iSlave)%max_len-1)then
            call set_error('Overlapping templates. Sorry','MakeTemplate')
            return
          end if

          ! Option 3
          if(template%items(iMaster)%start_pos(iMPos) <= &       
              & template%items(iSlave)%start_pos(iSPos) .and. &
              & template%items(iMaster)%start_pos(iMPos)+template%items(iMaster)%max_len-1 >= &
              & template%items(iSlave)%start_pos(iSPos))then
            call set_error('Overlapping templates. Sorry','MakeTemplate')
            return
          end if

        end do  ! cycle through slave item
      end do  ! cycle through master item
    end do
  end do

  !
  ! If final list of templates includes any of %Ta... - it means that analysis time varies
  ! inside the file set and it is reflected in the names. Then - set analysis_time_vary,
  ! which later will be read by write_ctl to determine the time shift
  !
  call get_what_vary(template, times_files)

  !
  ! Copy all template items positions to the positions array and sort them acsending
  !
  positions = int_missing
  i=1
  do t=1,template%n_items
    do iSPos=1,size(template%items(t)%start_pos)
      if(template%items(t)%start_pos(iSPos) == int_missing)cycle
      positions(i,1) = template%items(t)%start_pos(iSPos)
      positions(i,2) = t
      i=i+1
      if (i > size(positions,dim=1))then
        call set_error('Too many template positions','MakeTemplate')
        return
      end if
    end do 
  end do

  iMPos = i-1 ! Total number of active template positions 
  
  ifVarying = .true. ! If order of positions is varying - sort them
  do while (ifVarying)
    ifVarying = .false.
    do i=1,iMPos-1
      if(positions(i,1)==int_missing .or. positions(i+1,1)==int_missing)cycle
      if(positions(i,1) > positions(i+1,1))then
        t=positions(i,1)
        positions(i,1)=positions(i+1,1)
        positions(i+1,1)=t
        t=positions(i,2)
        positions(i,2)=positions(i+1,2)
        positions(i+1,2)=t
        ifVarying = .true.
      end if
    end do
  end do

  !
  ! Final step - create a templ_collection string from the first-time file name
  ! and the set of found templates for this time.
  !
  ! First - prepare indices and fill-in permanent part
  !
  t=1
  do while (len_trim(t_f%FNms(t)) < 2)
    t=t+1
  end do
  positions(iMPos+1,1) = len_trim(t_f%FNms(t))+1
  template%templ_collection = t_f%FNms(t)(1:positions(1,1)-1)
  !
  ! Now - a mixture of templates and, if any, permanent characters.
  ! To do it - run through all positions
  !
  iSPos = positions(1,1) ! iSPos is a flying position in a string here
  do i=1,iMPos
    call string_from_template(t_f%analysis_times(t),t_f%analysis_times(t), &
                            & t_f%forec_len(t), &
                            & template%items(positions(i,2))%item_type,templ_str)
    if(error)return
    template%templ_collection = fu_connect_strings(template%templ_collection, &
                        & fu_acronym(template%items(positions(i,2))%item_type))
    iSPos = iSPos + len_trim(templ_str)
    if(iSPos > len_trim(t_f%FNms(t)))exit
    template%templ_collection = fu_connect_strings(template%templ_collection, &
                                          & t_f%FNms(t)(iSPos:positions(i+1,1)-1))
    iSPos = positions(i+1,1)
  end do
  return

end subroutine MakeTemplate



!*********************************************************************************

subroutine decode_template_string(string_in, gr_template)
  !
  ! Decodes the template string to the grads_template structure.
  ! Method - just scans the string looking for % with further scan of all existing
  ! template values, which one fits. For each found template a new item is
  ! created with only one starting position, which corresponds to start
  ! of the template code in the templ_collection string.
  !
  ! JV addition: the string is first checked for environmental variables written 
  ! as ${SOMETHING}, which are expanded.
  !
  ! Author: M.Sofiev

  implicit none

  ! Imported parameters 
  character(len = *),intent(in)::string_in
  type(grads_template),intent(out):: gr_template

  ! Local variables
  character(len = fnlen) :: string
  integer :: i, iIndex, iPos, strLen

  string = fu_expand_environment(string_in)
  
  ! Preparation to the string analysis. In particular, if there is no templates
  ! the string will just be stored into gr_template.

  gr_template = template_missing
  gr_template%ini_string = string
  gr_template%n_items = 0
  gr_template%templ_collection = trim(adjustl(string))
  strLen = len_trim(gr_template%templ_collection)

  iPos = index(string(1:strLen),'%')
  iIndex = iPos

  do while (iIndex > 0)

    ! Check all types of templates - which one fits
    ! For versions of SILAM before v.3.2 - switch the comments of the below two lines
    !
    do i=1,nbr_of_templ_params
!      if(index(fu_str_u_case(string(iPos:iPos+len_trim(template_values(i)%tString))), &
!             & trim(fu_str_u_case(template_values(i)%tString))) /= 0)then ! found
      if(index((string(iPos:iPos+len_trim(template_values(i)%tString))), &
             & trim((template_values(i)%tString))) /= 0)then ! found
        gr_template%n_items = gr_template%n_items +1
        if(gr_template%n_items >= max_nbr_of_templates)then
          call set_error('too many templates','decode_template_string')
          return
        end if
        gr_template%items(gr_template%n_items)%item_type = template_values(i)%tType
        gr_template%items(gr_template%n_items)%min_len = template_values(i)%tLenMin
        gr_template%items(gr_template%n_items)%max_len = template_values(i)%tLenMax
        gr_template%items(gr_template%n_items)%start_pos(1) = iPos
        iIndex = index(string(iPos + 1:strLen),'%')
        iPos = iPos + iIndex
        exit
      end if
    end do

    if(i > nbr_of_templ_params)then ! Unknown template code in string
      call set_error(fu_connect_strings('Unknown template code in string:', &
                                      & string(iPos:strLen)), &
                   & 'decode_template_string')
      return
    end if
  end do

  gr_template%items(gr_template%n_items+1)%start_pos(1) = len_trim(gr_template%templ_collection)+1

end subroutine decode_template_string



!*********************************************************************************

subroutine add_time_file_to_tf_set(tf_set,analysis_time, forecast_len, FNm)
  !
  ! Adds a new time and correcponding file to the time_file_set
  !
  ! M.Sofiev
  !
  implicit none

  ! Imported variables
  type(time_file_set), intent(inout)::tf_set
  type(silja_time), intent(in) :: analysis_time
  type(silja_interval), intent(in) :: forecast_len
  character(len=*), intent(in) :: FNm

  if(tf_set%n_times >= 0) then ! Stupid, but there may be -1 as initialisation...
    tf_set%n_times = tf_set%n_times+1
  else
    tf_set%n_times = 1
  end if

  if(tf_set%n_times > size(tf_set%analysis_times))then
    call set_error('Too many time moments','add_time_file_to_tf_set')
    return
  end if

  tf_set%analysis_times(tf_set%n_times) = analysis_time
  tf_set%forec_len(tf_set%n_times) = forecast_len
  tf_set%FNMs(tf_set%n_times) = FNm

end subroutine add_time_file_to_tf_set



!*********************************************************************************

subroutine sort_times_ascend(tf_set)
  !
  ! Sorts the times in the set in ascending order. Very straightforward method
  !
  ! M.Sofiev, FMI
  !
  implicit none
  
  !Imported parameters
  type(time_file_set), intent(inout) :: tf_set

  ! Local variables
  logical :: OK
  integer :: t
  type(silja_time):: tTmp
  type(silja_interval):: inTmp
  character(len=fnlen):: chTmp

  OK =.false.
  do while(.not.OK)
    OK = .true.
    do t=1,fu_n_times(tf_set)-1
      if(fu_Valid_Time(tf_set,t+1) < fu_Valid_Time(tf_set,t))then
        tTmp = tf_set%analysis_times(t)
        tf_set%analysis_times(t) = tf_set%analysis_times(t+1)
        tf_set%analysis_times(t+1) = tTmp

        inTmp = tf_set%forec_len(t)
        tf_set%forec_len(t) = tf_set%forec_len(t+1)
        tf_set%forec_len(t+1) = inTmp

        chTmp = tf_set%FNms(t)
        tf_set%FNms(t) = tf_set%FNms(t+1)
        tf_set%FNms(t+1) = chTmp

        OK = .false.
      end if
    end do
  end do
end subroutine sort_times_ascend


!********************************************************************************

function fu_pure_grads_template(Template, chCaseNm, chSrcNm)
  !
  ! Our template list is more rich than that of GrADS. E.g., it has %source, %time,
  ! %a* templates, etc. But none of them is allowed in the ctl file. So, we may 
  ! want to create the pure GrADS template or at least check that the given one
  ! belongs to that class.
  ! Method - check that time dimension does not contain illegal stuff and explore
  ! those template items, which can be explored
  ! 
  implicit none
  !
  ! Return value of the function
  character(len=fnlen) :: fu_pure_grads_template

  ! Imported variables
  type(grads_template), intent(in) :: Template
  character(len=*), intent(in), optional :: chCaseNm, chSrcNm

  ! Local variables
  integer :: iTmp, iEndPos
  character(len = clen) :: strTempl
  character(len = fnlen) :: strTmp

  fu_pure_grads_template = ''

  ! Set the beginning of the file name up to the first template (if any)
  !
  if(template%n_items == 0)then
    fu_pure_grads_template = trim(template%templ_collection)
    return
  else 
    strTmp = template%templ_collection(1:template%items(1)%start_pos(1)-1)
  end if

  ! Scan all templates and stuff between them one-by-one
  !
  do iTmp = 1,template%n_items

    ! Always - some stupidity check
    !
    select case(template%items(iTmp)%item_type)
      !
      ! Allowed tempalte items: 
      !
      case (Ty2,Ty4,Tm1,Tm2,Tmc,Tmc_capital,Td1,Td2,Th1,Th2,Th3, Tn2, &
          & Tiy2,Tiy4,Tim1,Tim2,Timc,Timc_capital,Tid1,Tid2,Tih1,Tih2,Tih3, Tin2, &
          & Tf2, Tf3, Tfs2, Tfm2, Tfh2, Tfd2)

         strTempl = fu_template_string(template%items(iTmp)%item_type)
      !
      ! Forbidden tempalte items:
      !
      case (Tay2,Tay4,Tam1,Tam2,Tamc,Tamc_capital,Tad1,Tad2,Tah1,Tah2,Tah3, Tan2, Tjd3, Tajd3, Tijd3)
        call set_error('Analysis time and/or julian day is not allowed in the output template', &
                     & 'fu_pure_grads_template')
        fu_pure_grads_template = ''
        return
      !
      ! Items to be explored:
      !
      case(Tcase)
        if(present(chCaseNm))then
          strTempl = chCaseNm
        else
          call set_error('Case name requested but absent','fu_pure_grads_template')
          return
        endif

      case(Tsource)
        if(present(chSrcNm))then
          strTempl = chSrcNm
        else
          call set_error('Source name requested but absent','fu_pure_grads_template')
          return
        endif

      case(Tclock) ! We shall need just wall clock
        call string_from_template(time_missing, time_missing, interval_missing, &
                                & template%items(iTmp)%item_type, &
                                & strTempl)

      case(Ttask)
        write(unit=strTempl,fmt='(i03)') smpi_global_rank

      case default
        call set_error(fu_connect_strings('unknown template type:', &
                                        & fu_acronym(template%items(iTmp)%item_type)), &
                     & 'fu_pure_grads_template')
        return
    end select

    if(error)return

    ! Connect the item to the main srting
    !
    strTmp = fu_connect_strings(strTmp, strTempl)

    ! Are there fixed characters between the current and the next template items ?
    ! If any - connect.
    ! Note, n_item+1 has after-the-end start position
    !
    iEndPos = template%items(iTmp)%start_pos(1) + &
            & len_trim(fu_acronym(template%items(iTmp)%item_type))

    if(iEndPos < template%items(iTmp+1)%start_pos(1))then
      strTmp=fu_connect_strings(strTmp, &
                              & template%templ_collection(iEndPos: &
                                         & template%items(iTmp+1)%start_pos(1)-1))
    end if
  end do

  fu_pure_grads_template = strTmp

!  !
!  ! Copy the chTemplateStr to the output fu_pure_grads_template string exploring
!  ! the %case, %source and %clock items, if any
!  !
!  iTmp = index(chTemplateStr, '%')
!  jTmp = 1
!  fu_pure_grads_template = chTemplateStr(1:iTmp-1)
!
!  do while (iTmp > 0)
!
!    ! Add non-template part of the string
!    !
!    fu_pure_grads_template = fu_connect_strings(fu_pure_grads_template, &
!                                              & chTemplateStr(iTmp : iTmp+jTmp-1))
!    !
!    ! If template item - check and explore if needed
!    !
!    if(chTemplateStr(iTmp:iTmp+5) == '%case')then
!      if(.not.present(chCaseNm))then
!        call set_error('Case name is required but abssent','fu_pure_grads_template')
!        return
!      endif
!      iTmp = iTmp+5
!      fu_pure_grads_template = fu_connect_strings(fu_pure_grads_template, chCaseNm)
!
!    elseif(chTemplateStr(iTmp:iTmp+7) == '%source')then
!      if(.not.present(chSrcNm))then
!        call set_error('Source name is required but abssent','fu_pure_grads_template')
!        return
!      endif
!      iTmp = iTmp+7
!      fu_pure_grads_template = fu_connect_strings(fu_pure_grads_template, chSrcNm)
!
!    elseif(chTemplateStr(iTmp:iTmp+6) == '%clock')then
!      iTmp = iTmp+6
!      fu_pure_grads_template = fu_connect_strings(fu_pure_grads_template, &
!                                                & fu_time_fname_string_utc(fu_wallclock()))
!    else
!      ! Allowed GrADS template item. Do nothing...
!    endif
!
!    jTmp = index(chTemplateStr(iTmp : len_trim(chTemplateStr)),'%') ! Next template position
!
!  enddo

end function fu_pure_grads_template


!****************************************************************************

function fu_same_template_start_time(grTemplate, valid_time, anal_time)result(start_time)
  !
  ! Finds the earliest time for the given valid_time and analysis_time until the 
  ! template stays the same
  !
  implicit none

  ! Returned variable
  type(silja_time) :: start_time

  ! Imported parameters
  type(grads_template), intent(in) :: grTemplate
  type(silja_time), intent(in) :: valid_time, anal_time

  ! Local variables
  integer :: iTemplPeriod, iItem
  !
  ! Just find the fastest-varying time parameter in this template. This and all higher-rank
  ! parameters must stay the same, the lower-rank parameters go to zero
  !
  if(grTemplate%n_items == 0)then
    start_time = anal_time  ! non-varying template means that analysis time is the first one
    return
  else
    iTemplPeriod = -1
    do iItem = 1, grTemplate%n_items

      select case(grTemplate%items(iItem)%item_type)

        case(Tiy2, Tiy4, Tay2, Tay4)
          iTemplPeriod = max(iTemplPeriod, 1)  ! with regard to valid_time - no dependance
        case(Tim1, Tim2, Timc, Timc_capital, Tam1, Tam2, Tamc, Tamc_capital)
          iTemplPeriod = max(iTemplPeriod, 2)
        case(Tin2, Tan2)
          iTemplPeriod = max(iTemplPeriod, 3)
        case(Tid1, Tid2, Tijd3, Tad1, Tad2, Tajd3)
          iTemplPeriod = max(iTemplPeriod, 4)
        case(Tih1, Tih2, Tih3, Tah1, Tah2, Tah3)
          iTemplPeriod = max(iTemplPeriod, 5)

        case(Ty2, Ty4)
          iTemplPeriod = max(iTemplPeriod, 11)
        case(Tm1, Tm2, Tmc, Tmc_capital)
          iTemplPeriod = max(iTemplPeriod, 12)
        case(Td1, Td2, Tjd3, Tfd2)
          iTemplPeriod = max(iTemplPeriod, 13)
        case(Th1, Th2, Th3, Tf2, Tf3, Tfh2)
          iTemplPeriod = max(iTemplPeriod, 14)
        case(Tn2, Tfm2)
          iTemplPeriod = max(iTemplPeriod, 15)
        case(Tfs2)
          iTemplPeriod = max(iTemplPeriod, 16)
        case default
          call set_error('Unknown GrADS template','fu_same_template_start_time')
          return
      end select
    end do
    !
    ! Having the fastest-varying parameter, just get the time
    !
    select case(iTemplPeriod)
      case(1)
        start_time = fu_set_time_utc(fu_year(anal_time), 1, 1, 0, 0, 0.)
      case(2)
        start_time = fu_set_time_utc(fu_year(anal_time), fu_mon(anal_time), &
                               & 1, 0, 0, 0.)
      case(3)
        start_time = fu_set_time_utc(fu_year(anal_time), fu_mon(anal_time), &
                               & fu_day(anal_time), 0, 0, 0.)
      case(4)
        start_time = fu_set_time_utc(fu_year(anal_time), fu_mon(anal_time), &
                               & fu_day(anal_time), fu_hour(anal_time), 0, 0.)
      case(5)
        start_time = fu_set_time_utc(fu_year(anal_time), fu_mon(anal_time), &
                               & fu_day(anal_time), fu_hour(anal_time), fu_min(anal_time), &
                               & 0.)
      case(11)
        start_time = fu_set_time_utc(fu_year(valid_time), 1, 1, 0, 0, 0.)
      case(12)
        start_time = fu_set_time_utc(fu_year(valid_time), fu_mon(valid_time), &
                               & 1, 0, 0, 0.)
      case(13)
        start_time = fu_set_time_utc(fu_year(valid_time), fu_mon(valid_time), &
                               & fu_day(valid_time), 0, 0, 0.)
      case(14)
        start_time = fu_set_time_utc(fu_year(valid_time), fu_mon(valid_time), &
                               & fu_day(valid_time), fu_hour(valid_time), 0, 0.)
      case(15)
        start_time = fu_set_time_utc(fu_year(valid_time), fu_mon(valid_time), &
                               & fu_day(valid_time), fu_hour(valid_time), fu_min(valid_time), &
                               & 0.)
      case(16)
        start_time = fu_set_time_utc(fu_year(valid_time), fu_mon(valid_time), &
                               & fu_day(valid_time), fu_hour(valid_time), fu_min(valid_time), &
                               & real(int(fu_sec(valid_time))))
      case default
        call set_error('Unknown template validity length','fu_same_template_start_time')
    end select

  endif  ! if template is real

end function fu_same_template_start_time



!*************************************************************************
!
!   PRIVATE stuff
!
!*************************************************************************

!*************************************************************************

subroutine string_from_template(ini_time, anal_time, forecast_len, &
                              & template_type, &
                              & str, &
                              & chCase, chSource, chSector)
!
! Creates a string from given analysis time and forecast length 
! for particular template item.
! Rules are: 
!      analysis time determines %i* templates
!      forecast length determines %f* templates
!      together analysis time and forecast length determine templates
!           like %y*, %m*, %d*, %h*, %n* as they correspond to VALID TIME.
! 
! Substitutes time to t template item
!
! Author: M.Sofiev
!
  implicit none

  ! Result value
  character(len=*),intent(out) :: str

  ! Imported variables
  type(silja_time),intent(in):: ini_time, anal_time
  type(silja_interval), intent(in)::forecast_len
  integer, intent(in) :: template_type
  character(len=*), intent(in), optional :: chCase, chSource, chSector

  ! Local declarations
  integer :: year, month, day, hour, min
  real :: sec
  character(len=3),dimension(12):: mon_str = (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)
  character(len=3),dimension(12):: mon_str_capital = (/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)

  str = ''

  select case (template_type)
    case(Ty2)         !!!!!!!!!!!!!!!!!!!!  Valid-time templates ------*****
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') mod(year,100)
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Ty4)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i4)') year

    case(Tm1)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      if(month <10)then
        write(unit=str,fmt='(i1)') month
      else
        write(unit=str,fmt='(i2)') month
      end if

    case(Tm2)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') month
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tmc)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(a3)') mon_str(month)

    case(Tmc_capital)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(a3)') mon_str_capital(month)

    case(Td1)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      if(day <10)then
        write(unit=str,fmt='(i1)') day
      else
        write(unit=str,fmt='(i2)') day
      end if

    case(Td2)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') day
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tjd3)
      write(unit=str,fmt='(i3)') fu_julian_date(anal_time+forecast_len)
      if(str(1:1) == ' ')str(1:1) = '0'
      if(str(2:2) == ' ')str(2:2) = '0'

    case(Th1)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      if(hour <10)then
        write(unit=str,fmt='(i1)') hour
      elseif(hour < 100)then
        write(unit=str,fmt='(i2)') hour
      else
        write(unit=str,fmt='(i3)') hour
      end if

    case(Th2)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      if(hour < 100)then
        write(unit=str,fmt='(i2)') hour
        if(str(1:1) == ' ') str(1:1) = '0'
      else
        write(unit=str,fmt='(i3)') hour
      end if

    case(Th3)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i3)') hour
      if(str(1:1) == ' ')str(1:1) = '0'
      if(str(2:2) == ' ')str(2:2) = '0'

    case(Tn2)
      call get_time_components_utc(anal_time+forecast_len,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') min
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tf2)       !!!!!!!!!!!!!!!!!!!!!!!! Forecast length-------******
      if(fu_hours(forecast_len) >= 1000)then
        call set_error('Forecast length is too long','string_from_template')
        call report(forecast_len)
        return
      endif
      if(fu_hours(forecast_len) >= 100)then
        write(unit=str,fmt='(i3.3)') nint(fu_sec8(forecast_len)) / 3600
      else
        write(unit=str,fmt='(i3.2)') nint(fu_sec8(forecast_len)) / 3600
      endif

    case(Tf3)
      if(fu_hours(forecast_len) >= 1000)then
        call set_error('Forecast length is too long','string_from_template')
        call report(forecast_len)
        return
      endif
      write(unit=str,fmt='(i3.3)')  nint(fu_sec8(forecast_len)) / 3600

    case(Tfd2)
      write(unit=str,fmt='(i2.2)') nint(fu_sec8(forecast_len)) / (24*3600)

    case(Tfh2)
      write(unit=str,fmt='(i2.2)') mod( nint(fu_sec8(forecast_len)) / 3600, 24)  

    case(Tfm2)
      write(unit=str,fmt='(i2.2)') mod( nint(fu_sec8(forecast_len)) / 60, 60)

    case(Tfs2)
      write(unit=str,fmt='(i2.2)') mod( nint(fu_sec8(forecast_len)), 60)

    case(Tiy2)       !!!!!!!!!!!!!!!!!!!!! Initial times-------******
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') mod(year,100)
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tiy4)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i4)') year

    case(Tim1)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      if(month <10)then
        write(unit=str,fmt='(i1)') month
      else
        write(unit=str,fmt='(i2)') month
      end if

    case(Tim2)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') month
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Timc)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(a3)') mon_str(month)

    case(Timc_capital)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(a3)') mon_str_capital(month)

    case(Tid1)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      if(day <10)then
        write(unit=str,fmt='(i1)') day
      else
        write(unit=str,fmt='(i2)') day
      end if

    case(Tid2)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') day
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tijd3)
      write(unit=str,fmt='(i3)') fu_julian_date(ini_time)
      if(str(1:1) == ' ')str(1:1) = '0'
      if(str(2:2) == ' ')str(2:2) = '0'

    case(Tih1)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      if(hour <10)then
        write(unit=str,fmt='(i1)') hour
      elseif(hour < 100)then
        write(unit=str,fmt='(i2)') hour
      else
        write(unit=str,fmt='(i3)') hour
      end if

    case(Tih2)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      if(hour < 100)then
        write(unit=str,fmt='(i2)') hour
        if(str(1:1) == ' ') str(1:1) = '0'
      else
        write(unit=str,fmt='(i3)') hour
      end if

    case(Tih3)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i3)') hour
      if(str(1:1) == ' ')str(1:1) = '0'
      if(str(2:2) == ' ')str(2:2) = '0'

    case(Tin2)
      call get_time_components_utc(ini_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') min
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tay2)       !!!!!!!!!!!!!!!!!!!!! Analysis times-------******
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') mod(year,100)
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tay4)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i4)') year

    case(Tam1)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      if(month <10)then
        write(unit=str,fmt='(i1)') month
      else
        write(unit=str,fmt='(i2)') month
      end if

    case(Tam2)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') month
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tamc)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(a3)') mon_str(month)

    case(Tamc_capital)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(a3)') mon_str_capital(month)

    case(Tad1)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      if(day <10)then
        write(unit=str,fmt='(i1)') day
      else
        write(unit=str,fmt='(i2)') day
      end if

    case(Tad2)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') day
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tajd3)
      write(unit=str,fmt='(i3)') fu_julian_date(anal_time)
      if(str(1:1) == ' ')str(1:1) = '0'
      if(str(2:2) == ' ')str(2:2) = '0'

    case(Tah1)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      if(hour <10)then
        write(unit=str,fmt='(i1)') hour
      elseif(hour < 100)then
        write(unit=str,fmt='(i2)') hour
      else
        write(unit=str,fmt='(i3)') hour
      end if

    case(Tah2)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      if(hour < 100)then
        write(unit=str,fmt='(i2)') hour
        if(str(1:1) == ' ') str(1:1) = '0'
      else
        write(unit=str,fmt='(i3)') hour
      end if

    case(Tah3)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i3)') hour
      if(str(1:1) == ' ')str(1:1) = '0'
      if(str(2:2) == ' ')str(2:2) = '0'

    case(Tan2)
      call get_time_components_utc(anal_time,year,month,day,hour,min,sec)
      write(unit=str,fmt='(i2)') min
      if(str(1:1) == ' ')str(1:1) = '0'

    case(Tcase)
      if(present(chCase))then
        write(unit=str,fmt='(A)') trim(chCase)
      else
        call set_error('Case name is needed but absent','string_from_template')
        return
      endif

    case(Tsource)
      if(present(chSource))then
        write(unit=str,fmt='(A)') trim(chSource)
      else
        call set_error('Source name is needed but absent','string_from_template')
        return
      endif

!    case(Tsector)
!      if(present(chSector))then
!        write(unit=str,fmt='(A)') trim(chSector)
!      else
!        call set_error('Sector name is needed but absent','string_from_template')
!        return
!      endif

    case(TClock)
      str = fu_time_fname_string_utc(fu_wallclock())

    case(Ttask)
      write(unit=str,fmt='(i03)') smpi_global_rank

    case default
      call set_error('Unknown template type','string_from_template')
      return
  end select

  if(index(str,'*') /= 0)then  ! Too large number for example
    call report(ini_time)
    call report(anal_time)
    call report(forecast_len)
    print *, 'Template type: ', template_type
    call set_error(fu_connect_strings('Strange time leads to * in final string:',str), &
                 & 'string_from_template')
    str = ''
  endif

end subroutine string_from_template


!************************************************************************

function fu_acronym(ItemType)
  !
  ! Returns GrADS acronym of the template
  !
  implicit none

  ! Return value and imported parameter
  character(len=10) :: fu_acronym
  integer, intent(in) :: ItemType

  ! Local declaration
  integer :: i

  if(ItemType == int_missing)then
    fu_acronym = ''
    return
  end if

  do i=1,size(template_values)
    if(ItemType == template_values(i)%tType)then
      fu_acronym = template_values(i)%tString
      return
    end if
  end do

  fu_acronym = ''
  call set_error('Unknown template type','fu_acronym')

end function fu_acronym


!*************************************************************************

subroutine check_string(templ_str, FNm, indexAr)
!
! Checks all occasions of templ_str in FNm. Works differently depending
! on IndexTmp. If it is empty - it is filled in with all template locations
! if it is not empty - the locations stored in this IndexTmp are checked
! for existence in current FNm.
!
! Author: M.Sofiev
!
  implicit none

  ! Imported parameters
  character(len=*),intent(in) :: templ_str, FNm
  integer, dimension(:) :: indexAr

  ! Local variables
  integer :: iBeg, iEnd, iInd

  ! If templ_str does not exist in FNm string - nullify indexAr and quit
  !
!  print *, 'index(FNm,trim(templ_str))=','*',FNm,'*',trim(templ_str),'*',index(FNm,trim(templ_str))
  if(index(FNm,trim(templ_str)) == 0)then
    indexAr = int_missing
    return
  end if

  !
  ! Templ_str found in FNm string - now all depends on indexAr
  !
  iBeg = 1 ! Position in FNm
  iEnd = len_trim(FNm) ! Position in FNm
  if(ALL(indexAr == int_missing))then
    !
    ! index is empty - fill it in with all occasions
    !
    iInd = 1
    do while (index(FNm(iBeg:iEnd), trim(templ_str)) /= 0)
      indexAr(iInd) = index(FNm(iBeg:iEnd), trim(templ_str)) +iBeg-1
      iBeg = indexAr(iInd) + 1
      iInd = iInd + 1
    end do
  else
    !
    ! indexAr contains something - check it 
    !
    iBeg = 1
    do while(ANY(indexAr(iBeg:SIZE(indexAr)) /= int_missing))
      if(indexAr(iBeg) /= int_missing)then
!        print *, indexAr(iBeg), indexAr(iBeg)+len_trim(templ_str), '  ', FNm(indexAr(iBeg):indexAr(iBeg)+len_trim(templ_str)-1)
        if(index(FNm(indexAr(iBeg):indexAr(iBeg)+len_trim(templ_str)-1),&
               & trim(templ_str)) == 0) indexAr(iBeg) = int_missing
      end if
      iBeg = iBeg + 1
    end do
  end if

end subroutine check_string




!*************************************************************************

function fu_master_items (itemType)
  ! 
  ! Checks if given itemType has a "master" type, which is stronger. For
  ! example, Ty2 has a master type Ty4. Consequence - if master item is
  ! found than the slave one can be deleted from the list if their positions
  ! overlap.
  !
  ! Author: M.Sofiev
  !
  implicit none

  ! Return value is an array of suspected master items
  integer, dimension(5) :: fu_master_items

  ! Imported parameters
  integer, intent(in):: itemType

  fu_master_items = int_missing

  select case (itemType)
    case(Ty2)  !---------------- valid time
      fu_master_items(1) = Ty4
    case(Ty4)
    case(Tm1)
      fu_master_items(1) = Tm2
    case(Tm2)
    case(Tmc)
    case(Td1)
      fu_master_items(1) = Td2
    case(Td2)
    case(Tjd3)
    case(Th1)
      fu_master_items(1) = Th2
      fu_master_items(2) = Th3
    case(Th2)
      fu_master_items(1) = Th3
    case(Th3)
    case(Tn2)

    case(Tf2)  !------------------- forecast
      fu_master_items(1) = Tf3
    case(Tf3)

    case(Tiy2) !------------------- initial time
               ! can always be substituted by analysis time
      fu_master_items(1) = Tiy4
      fu_master_items(2) = Tay2
      fu_master_items(3) = Tay4
    case(Tiy4)
      fu_master_items(1) = Tay4
    case(Tim1)
      fu_master_items(1) = Tim2
      fu_master_items(2) = Tam1
      fu_master_items(3) = Tam2
    case(Tim2)
      fu_master_items(1) = Tam2
    case(Tin2)
      fu_master_items(1) = Tan2
    case(Timc)
      fu_master_items(1) = Tamc
    case(Tid1)
      fu_master_items(1) = Tid2
      fu_master_items(2) = Tad1
      fu_master_items(3) = Tad2
    case(Tid2)
      fu_master_items(2) = Tad2
    case(Tijd3)
      fu_master_items(2) = Tajd3
    case(Tih1)
       fu_master_items(1) = Tih2
       fu_master_items(2) = Tih3
       fu_master_items(3) = Tah1
       fu_master_items(4) = Tah2
       fu_master_items(5) = Tah3
    case(Tih2)
      fu_master_items(1) = Tih3
      fu_master_items(2) = Tah2
      fu_master_items(3) = Tah3
    case(Tih3)
      fu_master_items(1) = Tah3

    case(Tay2)  !------------------- analysis time
      fu_master_items(1) = Tay4
    case(Tay4)
    case(Tam1)
      fu_master_items(1) = Tam2
    case(Tam2)
    case(Tan2)
    case(Tamc)
    case(Tad1)
      fu_master_items(1) = Tad2
    case(Tad2)
    case(Tajd3)
    case(Tah1)
       fu_master_items(1) = Tah2
       fu_master_items(2) = Tah3
    case(Tah2)
      fu_master_items(1) = Tah3
    case(Tah3)
    case(Tcase, Tsource)
!    case(Tcase, Tsource, Tsector)
    case default
      call set_error('Unknown template type','fu_master_item')
      fu_master_items = int_missing
  end select
end function fu_master_items



!*************************************************************************

integer function fu_item_maxlength (itemType)
! 
! For a given item type returns its length after the substitution of a 
! values. Just searches the template_values array...
!
! Author: M.Sofiev
!
implicit none
  ! Imported parameters
  integer, intent(in):: itemType
  ! Local declarations
  integer :: i

  do i=1,size(template_values)
    if(template_values(i)%tType == itemType)then
      fu_item_maxLength = template_values(i)%TLenMax
      return
    end if
  end do
  call set_error('Unknown template type','fu_item_maxlength')
  fu_item_maxlength = int_missing

end function fu_item_maxlength



!***************************************************************************

subroutine compress_template_structure (template)
  !
  ! Removes holes from the list of templates and recalculates the 
  ! number of templates in the list
  !
  implicit none

  ! Imported parameter
  type(grads_template), intent(inout) :: template

  ! Local variables
  integer :: iItem

  template%n_items = size(template%items) ! Set the maximum first

  do while(template%items(template%n_items)%item_type == int_missing)
    template%n_items = template%n_items - 1
    if(template%n_items < 1)then
      call set_error('No template items found','compress_template_structure')
      return
    end if
  end do
  
  iItem = 1
  do while (iItem < template%n_items)
    if(template%items(iItem)%item_type == int_missing)then
      template%items(iItem) = template%items(template%n_items) !Delete item
      template%items(template%n_items)%item_type = int_missing
      template%n_items = template%n_items - 1
      do while(template%items(template%n_items)%item_type == int_missing)
        template%n_items = template%n_items - 1
      end do
    end if
    iItem = iItem +1
  end do

end subroutine compress_template_structure


!**************************************************************************

function fu_FNm_from_template(gr_templ, &
                            & ini_time, anal_time, forec_len, &
                            & chCase, chSource, chSector) result (strFNm)
  !
  ! Makes up the file name from the given template and times / intervals.
  ! Note that NONE of them is optional, so if you are sure that particular
  ! time can never be used - just put time_missing or interval_missing
  !
  ! If there is a template item requiring a certain time, but corresponding 
  ! time is missing - an error will be generated. Note: %i* templates require
  ! initial ini_time, %a* templates require analysis anal_time, %y/m/d/h/n* 
  ! require valid_time, and %f* require forecast forec_len
  !
  ! Author: Mikhail Sofiev
  !
  implicit none

  ! Result of the function - just file name string
  character(len = fnlen) :: strFNm

  ! Imported parameters
  type(grads_template), intent(in) :: gr_templ
  type(silja_time), intent(in) :: ini_time, anal_time
  type(silja_interval), intent(in) :: forec_len
  character(len=*), intent(in), optional :: chCase, chSource, chSector

  ! Local variables
  integer :: i, iEndPos
  character(len = fnlen) :: strTmp
  character(len = clen) :: strTempl

  ! Set the beginning of the file name up to the first template (if any)
  !
  if(gr_templ%n_items == 0)then
    strFNm = trim(gr_templ%templ_collection)
    return
  else 
    strTmp = gr_templ%templ_collection(1:gr_templ%items(1)%start_pos(1)-1)
!    gr_templ%items(gr_templ%n_items+1)%start_pos(1)=len_trim(gr_templ%templ_collection)+1
  end if

  ! Scan all templates and stuff between them one-by-one
  !
  do i = 1,gr_templ%n_items

    ! Always - some stupidity check
    !
    select case(gr_templ%items(i)%item_type)
      case (Ty2,Ty4,Tm1,Tm2,Tmc,Tmc_capital,Td1,Td2,Tjd3,Th1,Th2,Th3, Tn2)
        if(.not.(defined(anal_time).and.defined(forec_len)))then
          call set_error('can not compute valid time','fu_FNm_from_template')
          return
        end if
        call string_from_template(ini_time, anal_time, forec_len, &
                                & gr_templ%items(i)%item_type, &
                                & strTempl)

      case (Tiy2,Tiy4,Tim1,Tim2,Timc,Timc_capital,Tid1,Tid2,Tijd3,Tih1,Tih2,Tih3, Tin2)
        if(.not.defined(ini_time))then
          call set_error('undefined initial time','fu_FNm_from_template')
          return
        end if
        call string_from_template(ini_time, anal_time, forec_len, &
                                & gr_templ%items(i)%item_type, &
                                & strTempl)

      case (Tay2,Tay4,Tam1,Tam2,Tamc,Tamc_capital,Tad1,Tad2,Tajd3,Tah1,Tah2,Tah3, Tan2)
        if(.not.defined(anal_time))then
          call set_error('undefined analysis time','fu_FNm_from_template')
          return
        end if
        call string_from_template(ini_time, anal_time, forec_len, &
                                & gr_templ%items(i)%item_type, &
                                & strTempl)

      case (Tf2, Tf3, Tfs2, Tfm2, Tfd2, Tfh2)
        if(.not.defined(forec_len))then
          call set_error('undefined forecast length','fu_FNm_from_template')
          return
        end if
        call string_from_template(ini_time, anal_time, forec_len, &
                                & gr_templ%items(i)%item_type, &
                                & strTempl)

      case(Tcase)
        if(present(chCase))then
          call string_from_template(ini_time, anal_time, forec_len, &
                                  & gr_templ%items(i)%item_type, &
                                  & strTempl, chCase)
        else
          call set_error('Case name requested but absent','fu_FNm_from_template')
          return
        endif

      case(Tsource)
        if(present(chSource))then
          call string_from_template(ini_time, anal_time, forec_len, &
                                  & gr_templ%items(i)%item_type, &
                                  & strTempl,'',chSource)
        else
          call set_error('Source name requested but absent','fu_FNm_from_template')
          return
        endif

!      case(Tsector)
!        if(present(chSector))then
!          call string_from_template(ini_time, anal_time, forec_len, &
!                                  & gr_templ%items(i)%item_type, &
!                                  & strTempl,'','', chSector)
!        else
!          call set_error('Sector name requested but absent','fu_FNm_from_template')
!          return
!        endif

      case(Tclock) ! We shall need just wall clock
        call string_from_template(time_missing, time_missing, interval_missing, &
                                & gr_templ%items(i)%item_type, &
                                & strTempl)

      case(Ttask)
        call string_from_template(ini_time, anal_time, forec_len, &
            & gr_templ%items(i)%item_type, strTempl)

      case default
        call set_error(fu_connect_strings('unknown template type:', &
                                        & fu_acronym(gr_templ%items(i)%item_type)), &
                     & 'fu_FNm_from_template')
        return
    end select

    if(error)return

    ! Connect the item to the main srting
    !
    strTmp = strTmp + strTempl

    ! Are there fixed characters between the current and the next template items ?
    ! If any - connect.
    ! Note, n_item+1 has after-the-end start position
    !
    iEndPos = gr_templ%items(i)%start_pos(1) + len_trim(fu_acronym(gr_templ%items(i)%item_type))

    if(iEndPos < gr_templ%items(i+1)%start_pos(1))then
      strTmp = strTmp + gr_templ%templ_collection(iEndPos: gr_templ%items(i+1)%start_pos(1)-1)
    end if

  end do

  strFNm = strTmp

end function fu_FNm_from_template


!**************************************************************************

subroutine set_what_vary(tf_set, ifAnalysisTime, ifForecastLen)
  !
  ! Sets analysis_time_vary and forecast_len_vary.
  !
  implicit none
  logical, intent(in) :: ifAnalysisTime, ifForecastLen
  type(time_file_set),intent(inout) :: tf_set

  tf_set%Analysis_time_vary = ifAnalysisTime
  tf_set%forecast_length_vary = ifForecastLen
end subroutine set_what_vary



!**************************************************************************

subroutine get_what_vary(template, tf_set)
  !
  ! Sets tf_set%analysis_time_vary and tf_set%forecast_len_vary form the list
  ! of templates.
  !
  implicit none
  !
  ! Imported parameter
  type(grads_template),intent(in) :: template
  type(time_file_set), intent(inout) :: tf_set

  ! Local variables
  integer :: i

  do i=1,template%n_items
    select case(template%items(i)%item_type)
    case (Ty2, Ty4, Tm1, Tm2, Tmc, Td1, Td2, Tjd3, Th1, Th2, Th3, Tn2)
      ! Do nothing - valid_time will, of course, vary
    case (Tf2, Tf3, Tfs2, Tfm2, Tfh2, Tfd2)
      tf_set%forecast_length_vary = .true.
      tf_set%analysis_time_vary = .false.
      return
    case (Tay2, Tay4, Tam1, Tam2, Tan2, Tamc, Tad1, Tad2, Tajd3, Tah1, Tah2, Tah3)
      tf_set%forecast_length_vary = .false.
      tf_set%analysis_time_vary = .true.
      return
    case default
      call set_error('Unknown template type','get_what_vary')
      return
    end select
  end do
  
end subroutine get_what_vary

!**************************************************************************

type (silja_interval) function fu_template_timestep(template) 
  !
  ! Returns maximum interval that will not jump over a file in a given template
  ! Purpose -- proper step to check all possible files matching the template
  !
  implicit none
  !
  ! Imported parameter
  type(grads_template),intent(in) :: template

  ! Local variables
  integer :: i, itmpl
  type (silja_interval) :: tTmp

  fu_template_timestep = very_long_interval
  do i=1,template%n_items
    itmpl = template%items(i)%item_type
    select case(itmpl)
        case(Ty2, Ty4,  Tiy2, Tiy4, Tay2, Tay4)
          tTmp = one_day*365
        case(Tm1, Tm2, Tmc, Tmc_capital, Tam1, Tam2, Tamc, Tamc_capital, Tim1, Tim2, Timc, Timc_capital)
           tTmp = one_day * 28
        case(Td1, Td2, Tjd3, Tad1, Tad2, Tajd3, Tid1, Tid2,Tijd3, Tfd2)
           tTmp = one_day
        case(Th1, Th2, Th3, Tah2, Tah3,  Tf2, Tf3, Tfh2)
          tTmp = one_hour
        case(Tfm2, Tn2)
          tTmp = one_minute
        case(Tfs2)
          tTmp = one_second
        case default
          call set_error('Unknown GrADS template type: '//trim(fu_str(itmpl)),'fu_template_timestep')
          return
    end select
    if (fu_template_timestep > tTmp ) fu_template_timestep = tTmp
  end do
  
end function fu_template_timestep

!************************************************************************

logical function fu_compare_templates_eq(t1, t2)
  !
  ! Compares two templates for their identity. Reloaded operator ==
  !
  implicit none

  ! Imported variables
  type(grads_template), intent(in) :: t1, t2

  ! Local variables
  integer :: i

  fu_compare_templates_eq = .false.

  if(t1%n_items /= t2%n_items) return

  if(t1%ini_string /= t2%ini_string) return

  if(t1%templ_collection /= t2%templ_collection) return

  do i=1, t1%n_items
    if(.not. (t1%items(i) == t2%items(i)))return
  end do

  fu_compare_templates_eq = .true.

end function fu_compare_templates_eq


!******************************************************************

logical function fu_silam_template_item(item)
  !
  ! Returns true is the item is either %case or %source, which are both
  ! serve SILAM internal naming convention - case name and emission source
  ! name correspondingly. Here these items are not served because they need
  ! SILAM environment around
  !
  implicit none

  integer, intent(in) :: item

  fu_silam_template_item = (item == Tcase) .or. (item == TSource) .or. &
                         & (item == TClock)
!                         & (item == TSector) .or. (item == TClock)
end function fu_silam_template_item


!******************************************************************

function fu_silam_source_template_item()result(chTemplate)
  !
  ! Returns the string title of the SILAM source template: "%source"
  ! In fact, just incapsulation
  !
  implicit none
  character(len=7) :: chTemplate

  chTemplate = template_values(Isource)%tString

end function fu_silam_source_template_item


!!******************************************************************
!
!function fu_silam_sector_template_item()result(chTemplate)
!  !
!  ! Returns the string title of the SILAM source template: "%source"
!  ! In fact, just incapsulation
!  !
!  implicit none
!  character(len=7) :: chTemplate
!
!  chTemplate = template_values(Isector)%tString
!
!end function fu_silam_sector_template_item


!******************************************************************

logical function fu_compare_template_items_eq(item1, item2)
  !
  ! Reloaded operator == for template items
  ! 
  implicit none

  ! Imported parameters
  type(template_item), intent(in) :: item1, item2

  ! Local variables
  integer :: i

  fu_compare_template_items_eq = .false.

  if(item1%item_type /= item2%item_type) return
  if(item1%min_len /= item2%min_len) return
  if(item1%max_len /= item2%max_len) return

  do i = 1, size(item1%start_pos)
    if(item1%start_pos(i) /= item2%start_pos(i)) return
  end do

  fu_compare_template_items_eq = .true.

end function fu_compare_template_items_eq




!***************************************************************************
!***************************************************************************
!
!  Encapsulation stuff
!
!***************************************************************************
!***************************************************************************

! type template_item
integer function fu_item_type_of_template_item(item) ! Returns item type 
  implicit none
  type(template_item),intent(in)::item
  fu_item_type_of_template_item = item%item_type
end function fu_item_type_of_template_item

integer function fu_start_pos_of_template_item(item,iPos) ! Returns particular position
  implicit none
  type(template_item),intent(in)::item
  integer, intent(in)::iPos
  if(iPos > size(item%start_pos))then
    call set_error('Position index is too big','fu_start_pos_of_template_item')
  else
    fu_start_pos_of_template_item = item%start_pos(iPos)
  end if
end function fu_start_pos_of_template_item

function fu_start_pos_ar_of_templ_item(item) result(posAr) ! returns position array
  implicit none
  integer,dimension(max_pos_of_template):: posAr
  type(template_item),intent(in)::item
  posAr = item%start_pos
end function fu_start_pos_ar_of_templ_item


!***************************************************************************

! type grads_template

integer function fu_n_items_of_template(template) ! Returns number of items
  implicit none
  type(grads_template),intent(in)::template
  fu_n_items_of_template = template%n_items
end function fu_n_items_of_template

function fu_item_of_template(template,iItem)result(item) ! Returns item 
  implicit none
  type(template_item)::item
  type(grads_template),intent(in)::template
  integer, intent(in):: iItem
  if(iItem > size(template%items))then
    call set_error('Item index is too big','fu_item_of_template')
  else
    item = template%items(iItem)
  end if
end function fu_item_of_template

function fu_template_collection(template) ! Returns the final template collection
  implicit none
  character(len=fnlen)::fu_template_collection
  type(grads_template),intent(in)::template
  fu_template_collection = template%templ_collection
end function fu_template_collection

subroutine set_template_collection(template, chColl) ! Sets the template collection by force
  implicit none
  type(grads_template),intent(out)::template
  character(len=*), intent(in) :: chColl
  template%templ_collection = adjustl(chColl)
end subroutine set_template_collection


!************************************************************************
!
! type template_param
!
function fu_template_string_for_item(templateType)
  !
  ! Returns the strings tString of the template  - e.g. "%y2"
  !
  implicit none
  character(len=7) :: fu_template_string_for_item ! Return value
  integer, intent(in) :: templateType

  integer :: i

  do i=1, nbr_of_templ_params
    if(template_values(i)%tType == templateType)then
      fu_template_string_for_item = template_values(i)%tString
      return
    endif
  end do

  call set_error('Unknown template type','fu_template_string_for_item')
  print *,'Template type: ', templateType
  fu_template_string_for_item = ''
  
end function fu_template_string_for_item

!
! Below encapsulation functions return characteristic parts of the 
! tempalte strings, so that it can be quickly checked that 
! particular type of the template exists in the string
!
function fu_forecast_length_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = "%f"
end function fu_forecast_length_templ_str

function fu_valid_time_hour_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = template_values(Ih2)%tString
end function fu_valid_time_hour_templ_str

function fu_valid_time_day_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = template_values(Id2)%tString
end function fu_valid_time_day_templ_str

function fu_valid_time_month_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = '%m'
end function fu_valid_time_month_templ_str

function fu_valid_time_year_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = '%y'
end function fu_valid_time_year_templ_str

function fu_anal_time_hour_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = template_values(Iah2)%tString
end function fu_anal_time_hour_templ_str

function fu_anal_time_day_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = template_values(Iad2)%tString
end function fu_anal_time_day_templ_str

function fu_anal_time_month_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = '%am'
end function fu_anal_time_month_templ_str

function fu_anal_time_year_templ_str() result(str)
  implicit none
  character(len=7) :: str
  str = '%ay'
end function fu_anal_time_year_templ_str


!**************************************************************************

function fu_template_string_for_template(template) result(str)
  !
  ! A function returning a full string with template acronyms just what was
  ! send as input for decode_template_string subroutine
  !
  implicit none
  character(len=fnlen) :: str
  type(grads_template), intent(in) :: template

  str = template%ini_string

end function fu_template_string_for_template

!***************************************************************************

logical function grads_template_defined(templ)
  implicit none
  type(grads_template), intent(in) :: templ
  grads_template_defined = .not.(templ == grads_template_missing .or. &
                               & (templ%n_items == 0 .and. templ%ini_string == ''))
end function grads_template_defined

!***************************************************************************

subroutine report_grads_template(templ)
  implicit none
  type(grads_template), intent(in) :: templ
  call msg('Grads template report -------------------------')
  call msg('Number of items:', templ%n_items)
  call msg('Initial string:' + templ%ini_string)
  call msg('Final collection now:' + templ%templ_collection)
  call msg('End of grads template report ------------------')
end subroutine report_grads_template


!***************************************************************************
!
! Input for the template-making routines: a set of times and file names
!
function fu_analysis_time_of_tf_set(tf,indTime) result(time)
  implicit none
  type(silja_time)::time
  type(time_file_set),intent(in)::tf
  integer, intent(in)::indTime
  if(indTime > size(tf%analysis_times))then
    call set_error('Time index is too big','fu_ref_time_of_tf_set')
    time = time_missing
  else
    time = tf%analysis_times(indTime)
  end if
end function fu_analysis_time_of_tf_set

function fu_valid_time_of_tf_set(tf,indTime) result(time)
  implicit none
  type(silja_time)::time
  type(time_file_set),intent(in)::tf
  integer, intent(in)::indTime
  if(indTime > size(tf%analysis_times))then
    call set_error('Time index is too big','fu_valid_time_of_tf_set')
    time = time_missing
  else
    time = tf%analysis_times(indTime)+tf%forec_len(indTime)
  end if
end function fu_valid_time_of_tf_set

function fu_forec_len_of_tf_set(tf,indTime) result(length)
  implicit none
  type(silja_interval)::length
  type(time_file_set),intent(in)::tf
  integer, intent(in)::indTime
  if(indTime > size(tf%analysis_times))then
    call set_error('Time index is too big','fu_forec_len_of_tf_set')
    length = interval_missing
  else
    length = tf%forec_len(indTime)
  end if
end function fu_forec_len_of_tf_set

logical function fu_anal_time_of_tf_set_vary(tf)
  implicit none
  type(time_file_set),intent(in)::tf
  fu_anal_time_of_tf_set_vary = tf%analysis_time_vary
end function fu_anal_time_of_tf_set_vary

logical function fu_forec_len_of_tf_set_vary(tf)
  implicit none
  type(time_file_set),intent(in)::tf
  fu_forec_len_of_tf_set_vary = tf%forecast_length_vary
end function fu_forec_len_of_tf_set_vary

function fu_FNm_of_tf_set(tf,iNm) result(name)
  implicit none
  character(len=fnlen):: name
  type(time_file_set),intent(in)::tf
  integer, intent(in)::iNm
  if(iNm > size(tf%analysis_times))then
    call set_error('File name index is too big:' + fu_str(iNm),'fu_FNm_of_tf_set')
    name = ''
  else
    name = tf%FNms(iNm)
  end if
end function fu_FNm_of_tf_set

integer function fu_n_times_of_tf_set(tf)
  implicit none
  type(time_file_set),intent(in)::tf
  fu_n_times_of_tf_set = tf%n_times
end function fu_n_times_of_tf_set


!********************************************************************************


  
end module grads_templates
