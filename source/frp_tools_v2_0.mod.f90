MODULE frp_tools
  !
  ! This module contains the  tools for handling modis frp product
  ! Three processing stages are served.
  ! 1. Parsing the MODIS HDF files for FRP and text files for TA. Four sets of data
  !    are obtained: FRP Aqua/Terra and TA Aqua/Terra
  ! 2. The four datasets of stage one are merged together. If needed, TA is converted to FRP
  !    and used to patch holes in the FRP time series. Aqua and Terra are written in one file
  ! 3. The merged dataset of stage 2 is analysed for spatial amd temporal overlaps and 
  !    complemetarity and processed into SILAM fire source v.1
  !
  ! Language: ANSI FORTRAN 90
  !
  ! Author: Mikhail Sofiev FMI, mikhail.sofiev@fmi.fi
  !

  !USE DFLIB
  !USE DFPORT

  use grids_geo
  use grads_templates
  use field_identifications
  use grads_io

  !    USE IFPORT
  !    USE PORTLIB
  !    use grids_geo
  !    use grads_templates

  implicit none

  public make_silam_fire_src_v1
  public mask_wrong_fires
  public find_wrong_fires
  public parse_modis_ta

  type TSatellitedata
    character(len=clen) :: chInstrument, chProduct
    integer :: nVals
    real, dimension(:), pointer :: lon, lat, val, dx, dy, T21, T31, T21mean, T31mean, DTmean
    integer, dimension(:), pointer :: hhmm
  end type TSatellitedata

   
   character(len=*), parameter :: chLabelStr = '#         index  -date-                              lon       lat      dx    dy unit     FRP         T4    T4b      T11    T11b      TA     MCE  areaFr'
                                           ! fire =       1 2009    1    1    0   10  0.0  151.2433  -24.9140   1.004 1.002 km    40.874 MW 346.752 309.340 316.984 306.880  10.104   0.500   0.010
  CONTAINS



  !********************************************************************************************

 subroutine parse_modis_ta(chInFile, chOutFile)
     !
     ! Reads the MODIS FRP and TA for one satellite, selects the better dataset, writes it to
     ! uFRPintermediate with conversion, if needed
     !
     implicit none

     ! Imported parameters
     character(len=*), intent(in) :: chInFile, chOutFile

     ! Local variables
     integer :: iInput, iTmp, jTmp, status, uFRPraw, uIn, iPos, imm, ihh, idd, iyyyy, hhmm, iconf
     character(len=20) :: chInstrument, chQuantity, sat, chInstrTmp, chDateTmp
     character(len=fnlen) :: chInputTemplate, chOutFNm, chLine
     type(grads_template) :: gr_template
     type(TSatellitedata) :: pArFRP
     logical :: eof
     real :: lat, lon, val, dx, dy
     character(len=*), parameter :: sub_name="parse_modis_ta"


     ! Input file
     uIn = fu_next_free_unit()
     open(uIn,file=chInFile,status='old',action='read',iostat = status)
     if(status /= 0)then
       call set_error('Failed to open the TA file: '//trim(chInFile), sub_name)
       return
     endif

     !
     ! The output file, which will collect FRP from the TA files
     !
     uFRPraw = fu_next_free_unit()
     open(uFRPraw, file=chOutFile,iostat = status)
     if(status /= 0)then
       call set_error('Failed to open the FRP raw file for output: '//trim(chOutFile), sub_name)
       return
     endif
     write(uFRPraw, '(A)')'#'
     write(uFRPraw, '(2A)')'# Infile:  ', trim(chInFile)
     write(uFRPraw, '(2A)')'# OutFile: ', trim(chOutFile)
     write(uFRPraw, '(4A)')'# Created with is4fires ',sub_name, '  from ', revision_str

     !
     ! Cycle through the From file copying it to To file, converting stuff
     ! Three different formats (att least) supported
     !

     do while (.True.)
       READ(unit = uIn, fmt = '(A)', iostat = status) chLine
       if (status /= 0) exit

       if(index(chLine,'latitude') > 0)cycle
       if (len(trim(chLine)) == 0) cycle

             do iPos = 1, len_trim(chLine)
               if(chLine(iPos:iPos) == ',')chLine(iPos:iPos) = ' '
               if(chLine(iPos:iPos) == '/')chLine(iPos:iPos) = ' '
             end do

             iPos = index(chLine,':')
             if(iPos > 0)then
               chLine(iPos:iPos) = ' '
!latitude,longitude,brightness,scan,track,acq_date,acq_time,satellite,confidence,version,bright_t31,frp
!-1.405,29.201,301.7,2.1,1.4,2016-03-06,00:00,A,44,6.0NRT,277.2,22.2
!or
!latitude,longitude,brightness,scan,track,acq_date,acq_time,satellite,confidence,version,bright_t31,frp
!-29.169,26.208,306.8,3,1.7,2012-03-11, 00:05,A,49,5.0       ,283.7,77
               read(unit=chLine, fmt= *, iostat=status) &
                                          & lat, lon, val, dx, dy, chDateTmp, ihh, imm, sat, iconf
               hhmm = ihh * 100 + imm
               chDateTmp(5:5) = " " !Convert YYYY-mm-dd to  YYYY mm dd
               chDateTmp(8:8) = " "
               
             else

!-6.435,-35.043,305.7,4.1,1.9,12/04/2007,0015,T,35

               read(unit=chLine, fmt=*, iostat=status) lat, lon, val, dx, dy, &
                                                   & imm, idd, iyyyy, hhmm, sat, iconf
               WRITE(chDateTmp, fmt = '(I4.4,X,I2.2,X,I2.2)') iyyyy,imm,idd
             endif
             if(status /= 0)then
               call msg_warning('Failed to read the line:' + chLine,'parse_MODIS_files')
               continue
             endif
             if(hhmm < 0 .or. hhmm > 2400)then
               call set_error('Wrong line:'+chLine,'convert_modis_ta_2_frp')
               call msg('********** TA Problematic file: ' + chInFile)
               call unset_error('convert_modis_ta_2_frp')
               cycle
             endif

             !
           ! Convert to FRP and write down
             !
             val = ((((1.040451202E-006 * val - 0.001403323036) * &
                                 & val + 0.7336164398) * val - 174.051266) * val + 15664.59653)

!            # longitude, latitude, FRP_value, T21, T31, T21mean, T31mean, DeltaTmean, dx, dy, year, mon, day, hhmm, instrument
!             print *,"Ops", jTmp

             if(fu_str_u_case(trim(sat)) == 'T')then
               !jTmp = jTmp
               write(uFRPraw,fmt='(3(F10.3,1x),A, 2(F6.2,1x),1x,A,I5,1x,A,I5)', iostat=status) &
               !  print '(3(F10.3,1x),A, 2(F6.2,1x),1x, A,I5,1x,A,I5)',  &
                                             & lon, lat, val, '-1 -1 -1 -1 -1 ', dx, dy, &
                                             & trim(chDateTmp), hhmm, &
                                             & 'MODIS_TERRA_TA_as_FRP', 100
             else
               write(uFRPraw,fmt='(3(F10.3,1x),A, 2(F6.2,1x),1x,A,I5,1x,A,I5)', iostat=status) &
               !  print *,  &
                                             & lon, lat, val, '-1 -1 -1 -1 -1 ', dx, dy, &
                                             & trim(chDateTmp), hhmm, &
                                             & ' MODIS_AQUA_TA_as_FRP', 100
                
             endif
              jTmp = jTmp + 1
          end do  ! through file

          call msg('Lines in TA file:'+chInFile+':',jTmp)

          close(uIn)
          close(uFRPraw)


    end subroutine parse_modis_ta 



  !********************************************************************************************


  subroutine make_silam_fire_src_v1(arFRP_FNm, outFNm)
    !
    ! Merges several observation flows into a single dataset
    ! Aggregates space- and time- wise related frp points and writes down the SILAM fire source file
    !
    implicit none

    character(len=*), dimension(:), intent(in) :: arFRP_FNm
    character(len=*), intent(in) :: outFNm

    !local variables

    real, dimension(:), allocatable :: frp, lat, lon, dx, dy, T21, T31, T21mean, T31mean, DTmean
    integer , dimension(:),allocatable :: mark, confidence
    character(len=clen), dimension(:), pointer :: label
    integer :: i,j,iStatus, npnts,nzap, uFRPintermediate, uOut, indFire, npnts_max, iPoint, &
             & iSkip, nTimeSeries, nFiresInTimeSeries, iTmp, iTmp1, iTmp2, iTmp3, iTmp4, nzapTA, &
             & day, mon, year, hhmm
    real, dimension(:,:), pointer, save :: mfrp,mlat,mlon,mT21, mT31, mT21mean, mT31mean, mDTmean, mdx,mdy
    integer, dimension(:,:), pointer, save :: mmark
    character(len=fnlen) :: line
    type(silja_time), dimension(:,:), pointer, save :: timeFire
    type(silja_time), dimension(:), pointer, save :: timeTmp
    logical :: eof, ifSplit, done_with_frp
    logical, save :: ifFirst = .true.
    real :: FRP_grand_total, FRP_eaten, FRP_SILAM_file, fTmp1, fTmp2, fTmp3, fTmp4, fTmp_eliminated

    allocate(frp(worksize), lat(worksize), lon(worksize), dx(worksize), dy(worksize), &
       & T21(worksize), T31(worksize), T21mean(worksize), T31mean(worksize), DTmean(worksize), &
       & mark(worksize), confidence(worksize))
    !
    ! read from files.
    ! Two groups of files: FRP and TA-as-FRP
    ! The FRP files are taken as they are and merged, TA-as-FRP are only used in points where FRP is 
    ! not available
    !
    npnts=0
    nTimeSeries = 0
    nFiresInTimeSeries = 0

    if(error)return

    mark = 0
    FRP_grand_total = 0.0
    FRP_eaten = 0.0
    FRP_SILAM_file = 0.0
    frp = 0.0; lon = 0.0; lat = 0.0; dx = 0.0; dy = 0.0
    nzap = 0
    !
    ! Count the grand total number of lines that might be needed for the exercise
    ! FRP files
    !
    do iTmp1 = 1, size(arFRP_FNm)
      inquire(file=arFRP_FNm(iTmp1), exist=eof)
      if(eof)then
        nzap = nzap + fu_number_of_lines_in_file(arFRP_FNm(iTmp1))
      else
        call msg_warning('Could not find the file:' + arFRP_FNm(iTmp1),'make_silam_fire_src_v1')
      endif
      if(error)return
    end do ! FRP files

    !
    ! Having the max imaginable number of lines, allocate the label array - the rest are work arrays
    !
    allocate(label(nzap), timeTmp(nzap), stat=iStatus)
    if(fu_fails(iStatus == 0,'Failed to allocate intermediate arrays','make_silam_fire_src_v1'))return
    timeTmp = time_missing

    uOut = fu_next_free_unit()
    open(uOut,file=outFNm,iostat=iStatus)
    if(fu_fails(iStatus == 0,'Cannot open the output file','make_silam_fire_src_v1'))return

    !
    ! Read the FRP files and collect the info into one long input array
    !
    uFRPintermediate = fu_next_free_unit()
    i = 1
    do iTmp1 = 1, size(arFRP_FNm)
      call msg('Reading FRP input:' + arFRP_FNm(iTmp1))
      open(uFRPintermediate, file = arFRP_FNm(iTmp1), action='read', iostat=iStatus)
      if(iStatus /= 0)then
        call msg_warning('Cannot open the FRP input file:' + arFRP_FNm(iTmp1), &
                       & 'make_silam_fire_src_v1')
        cycle
      endif
      !
      !--------------------------------------------------------------------------
      !
      ! READ the fires: FRP
      !
      eof = .false.
      done_with_frp = .false.
      do while (.not. eof)
        call next_line_from_input_file(uFRPintermediate, line, eof)
        if(eof)exit

        read(unit=line, fmt=*, iostat=iStatus) lon(i),lat(i),frp(i), &
                                             & T21(i), T31(i), &
                                             & T21mean(i), T31mean(i), DTmean(i), &
                                             & dx(i), dy(i), &
                                             & year, mon, day, hhmm, label(i), confidence(i)
        if(fu_fails(iStatus==0,'Failed intermediate file reading:' + line + &
                             & ', input file:' + arFRP_FNm(iTmp1),'make_silam_fire_src_v1'))return
        iTmp=len(trim(label(i)))
        if (label(i)(iTmp-8:iTmp) == "TA_as_FRP") done_with_frp = .True.
        if (done_with_frp .and. label(i)(iTmp-4:iTmp) == "A_FRP") then
          call set_error("True FRP must be processed first", 'make_silam_fire_src_v1')
          return
        endif

        if(frp(i) < 0.000001)then
call msg('Eliminated zero FRP in:' + line)
          mark(i) = 0
          cycle
        endif
        timeTmp(i) = fu_set_time_utc(year, mon, day, int(hhmm / 100), mod(hhmm, 100), 0.)
        if(error .or. .not. defined(timeTmp(i)))then
          mark(i) = 0
          if(error)call unset_error('After problematic line:' + line)
          cycle
        endif
        mark(i) = 1
        FRP_grand_total = FRP_grand_total + frp(i)
        i = i + 1
      end do  ! single file
      close(uFRPintermediate)
      call msg('Total FRP after: '//trim(arFRP_FNm(iTmp1)), FRP_grand_total)

    end do  ! FRP input files
    !
    ! The TA-based representation of fires is taken only if there is no FRP stuff.
    ! So, we read the TA files forllowing the same template but add the fire into the list
    ! only if there is no one in FRP covering that time and place. Boring exercise but 
    ! nothing to do...
    !
do iTmp3 = 1, nzap
if(mark(iTmp3) < 1)cycle
if(frp(iTmp3) < 0.001)call set_error('Low FRP -1' + fu_str(iTmp3),'after FRP reading')
enddo

    frp(i) = 0.0  ! for the case if last line was rejected but not emptied

    !call msg('Grand total FRP:', FRP_grand_total)


    !
    ! Two steps of the merging.
    ! 1. The pixels inside a single frame can be overlapping, thus seeing the same fire 
    !    more than once. In this case, the max-FRP pixel is selected, with care taken to avoid 
    !    a chain effect when every other pixel sequentially eliminates the previous one
    ! 2. In different frames (or frames from different instruments) the overlapping pixels
    !    mean several-time-moment observations of the same fire. In this case, the diurnal 
    !    variation is accounted for and the daily-mean FRPs are averaged.
    !
    ! Step 1. Same-frame repetitive fire images. Take max, preclude chain effect
    !
    if(.not. ifFirst)then
      if(nzap > size(mfrp,2))then
        deallocate(mfrp,mlon,mlat,mT21, mT31, mT21mean, mT31mean, mDTmean, timeFire, mmark)
        ifFirst = .true. 
      endif
    endif
    if(ifFirst)then
      allocate(mfrp(50,nzap),mlon(50,nzap),mlat(50,nzap),timeFire(50,nzap),mdx(50,nzap),mdy(50,nzap), &
             & mT21(50,nzap), mT31(50,nzap), mT21mean(50,nzap), mT31mean(50,nzap), mDTmean(50,nzap), &
             & mmark(50,nzap))
      ifFirst = .false.
    endif
    timeFire = time_missing
    mfrp = 0.
    mlon = 0.
    mlat = 0.
    mT21 = 0.
    mT31 = 0.
    mT21mean = 0.
    mT31mean = 0.
    mDTmean = 0.
    mdx = 0.
    mdy = 0.
    mmark = 0
    iSkip = 0

    !-------------------------------------------------------------------------------
    !
    ! Step 1. Remove the duplicated observations from the same frame
    !
    call remove_duplicates(timeTmp, frp, lon, lat, T21, T31, T21mean, T31mean, DTmean, dx, dy, mark, nzap, FRP_eaten)
    if(error)return

fTmp1 = 0.0
fTmp2 = 0.0
iTmp1 = 0
iTmp2 = 0
do i = 1, nzap
  if(mark(i) > 0)then
    fTmp1 = fTmp1 + frp(i)
    iTmp1 = iTmp1 + 1
  else
    fTmp2 = fTmp2 + frp(i)
    iTmp2 = iTmp2 + 1
  endif
enddo
call msg('After removing duplicates, valid fires and FRP:',iTmp1, fTmp1)
call msg('... and eliminated fires and FRP:',iTmp2, fTmp2)
call msg('Total:',iTmp1+iTmp2, fTmp1+fTmp2)
fTmp_eliminated = fTmp2
if(abs(fTmp1+fTmp2-FRP_grand_total)/(FRP_grand_total +1e-5) > 0.001)then
  call set_error('High error','main')
  stop
endif
if(abs(fTmp1+FRP_eaten-FRP_grand_total)/(FRP_grand_total + 1e-5) > 0.001)then
  call set_error('High error 2','main')
  stop
endif

    !---------------------------------------------------------------------------------
    !
    ! Step 2. Scan the fires one by one, if found several ones for the same fire with different 
    ! time stamps, write to that very fire index, otherwise take next
    !
    indFire = 0
    call make_time_series(timeTmp, frp, lon, lat, T21, T31, T21mean, T31mean, DTmean, dx, dy, mark, nzap, indFire)

    nTimeSeries = indFire  ! store the actual number of independent fires

fTmp1 = 0.0
fTmp2 = 0.0
fTmp3 = 0.0
fTmp4 = 0.0
iTmp1 = 0
iTmp2 = 0
iTmp3 = 0
iTmp4 = 0
do i = 1, nzap
  if(mark(i) > 0)then
    fTmp1 = fTmp1 + frp(i)  ! valid raw FRP
    iTmp1 = iTmp1 + 1
  else
    fTmp2 = fTmp2 + frp(i)
    iTmp2 = iTmp2 + 1
  endif
enddo
do i = 1, nTimeSeries
  do j = 1, size(mmark,1)
    if(mmark(j,i) > 0)then
      fTmp3 = fTmp3 + mfrp(j,i)
      iTmp3 = iTmp3 + 1
    else
      fTmp4 = fTmp4 + mfrp(j,i)
      iTmp4 = iTmp4 + 1
    endif
  end do
end do
call msg('After making time series, valid raw fires and FRP:',iTmp1, fTmp1)
call msg('... and eliminated fires and FRP:',iTmp2, fTmp2)
call msg('... and valid time series and FRP:',iTmp3, fTmp3)
call msg('... and void time series and FRP:',iTmp4, fTmp4)
call msg('Total N time series, FRP with eliminated duplicates:',iTmp1+iTmp2+iTmp3+iTmp4, fTmp3 + fTmp_eliminated)
if(abs(fTmp3+fTmp_eliminated-FRP_grand_total)/(FRP_grand_total + 1e-5) > 0.001)then
  call set_error('High error 2','main')
  stop
endif


    !--------------------------------------------------------------------------------------
    !
    ! Step 3. Analyse each time fire line to identify the bunches of spots with the same time.
    !         These are the high-resolution pixels collected by a single (few) large pixels
    !         from another overpass.
    ! Options are:
    ! - kill the large pixel and re-analyse the time line in order to split the fires. 
    !   +: more detailed land-use, -:loss of temporal information
    ! - split the large pixel by "projecting" it to small ones proportionaly to their FRPs
    !   +: no info loss, -: artificial trick, voluntary definition of pixels to split
    ! - sum-up small pixels. +: plume-rise for close plumes, -: loss of land-use
    !
    ! Below call realises the option 2.
    !
    indFire = 1
    do while(indFire <= nTimeSeries)

!if(indFire > 3384)then
! call msg('1:',indFire)
! do i = 1, size(mlon,1)
!  if(.not. defined(timeFire(i,indFire)))cycle
!  call msg('i,lon,lat:'+ fu_str(i),mlon(i,indFire),mlat(i,indFire))
!  call msg('frp,mark',mfrp(i,indFire),mmark(i,indFire))
!  call msg('dx,dy',mdx(i,indFire),mdy(i,indFire))
! enddo
!endif

      call analyse_time_series(timeFire(:,indFire), mfrp(:,indFire), &
                             & mlon(:,indFire), mlat(:,indFire), &
                             & mT21(:,indFire), mT31(:,indFire), &
                             & mT21mean(:,indFire), mT31mean(:,indFire), &
                             & mDTmean(:,indFire), &
                             & mdx(:,indFire), mdy(:,indFire), mmark(:,indFire), ifSplit, nFiresInTimeSeries)

!fTmp1 = 0.0
!fTmp2 = 0.0
!fTmp3 = 0.0
!fTmp4 = 0.0
!iTmp1 = 0
!iTmp2 = 0
!iTmp3 = 0
!iTmp4 = 0
!do i = 1, nzap
!  if(mark(i) > 0)then
!    fTmp1 = fTmp1 + frp(i)
!    iTmp1 = iTmp1 + 1
!  else
!    fTmp2 = fTmp2 + frp(i)
!    iTmp2 = iTmp2 + 1
!  endif
!enddo
!do i = 1, nTimeSeries
!  do j = 1, size(mmark,1)
!    if(mmark(j,i) > 0)then
!      fTmp3 = fTmp3 + mfrp(j,i)
!      iTmp3 = iTmp3 + 1
!    else
!      fTmp4 = fTmp4 + mfrp(j,i)
!      iTmp4 = iTmp4 + 1
!    endif
!  end do
!end do
!if(abs(fTmp3+fTmp_eliminated-FRP_grand_total)/FRP_grand_total > 0.001)then
!  call msg('')
!  call msg('After ANALYSING time series, valid raw fires and FRP:',iTmp1, fTmp1)
!  call msg('... and eliminated fires and FRP:',iTmp2, fTmp2)
!  call msg('... and valid time series and FRP:',iTmp3, fTmp3)
!  call msg('... and void time series and FRP:',iTmp4, fTmp4)
!  call msg('Total N time series and FRP with eliminated:',iTmp1+iTmp2+iTmp3+iTmp4, fTmp3+fTmp_eliminated)
!  call msg('Problem at indFire:', indFire)
!  call set_error('High error 3','main')
!  stop
!endif


      if(ifSplit)then
call msg('Splitting time series:', indFire, nTimeSeries)


        call remove_duplicates(timeFire(:,indFire), mfrp(:,indFire), &
                             & mlon(:,indFire), mlat(:,indFire), &
                             & mT21(:,indFire), mT31(:,indFire), &
                             & mT21mean(:,indFire), mT31mean(:,indFire), mDTmean(:,indFire), &
                             & mdx(:,indFire), mdy(:,indFire), &
                             & mmark(:,indFire), nFiresInTimeSeries, FRP_eaten)
if(fTmp1 > 1.0)then
  call msg('indFire: duplicates are removed:',fTmp1)
endif


!fTmp1 = 0.0
!fTmp2 = 0.0
!fTmp3 = 0.0
!fTmp4 = 0.0
!iTmp1 = 0
!iTmp2 = 0
!iTmp3 = 0
!iTmp4 = 0
!do i = 1, nzap
!  if(mark(i) > 0)then
!    fTmp1 = fTmp1 + frp(i)
!    iTmp1 = iTmp1 + 1
!  else
!    fTmp2 = fTmp2 + frp(i)
!    iTmp2 = iTmp2 + 1
!  endif
!enddo
!do i = 1, nTimeSeries
!  do j = 1, size(mmark,1)
!    if(mmark(j,i) > 0)then
!      fTmp3 = fTmp3 + mfrp(j,i)
!      iTmp3 = iTmp3 + 1
!    else
!      fTmp4 = fTmp4 + mfrp(j,i)
!      iTmp4 = iTmp4 + 1
!    endif
!  end do
!end do
!if(abs(fTmp3+fTmp_eliminated-FRP_grand_total)/FRP_grand_total > 0.001)then
!  call msg('')
!  call msg('After REMOVING DUPLICATES-2, valid raw fires and FRP:',iTmp1, fTmp1)
!  call msg('... and eliminated fires and FRP:',iTmp2, fTmp2)
!  call msg('... and valid time series and FRP:',iTmp3, fTmp3)
!  call msg('... and void time series and FRP:',iTmp4, fTmp4)
!  call msg('Total N time series and FRP with eliminated:',iTmp1+iTmp2+iTmp3+iTmp4, fTmp3+fTmp_eliminated)
!  call msg('Problem at indFire:', indFire)
!  call set_error('High error 4','main')
!  stop
!endif


        call make_time_series(timeFire(:,indFire), mfrp(:,indFire), &
                            & mlon(:,indFire), mlat(:,indFire), &
                            & mT21(:,indFire), mT31(:,indFire), mT21mean(:,indFire), &
                            & mT31mean(:,indFire), mDTmean(:,indFire), &
                            & mdx(:,indFire), mdy(:,indFire), &
                            & mmark(:,indFire), nFiresInTimeSeries, nTimeSeries)    ! the last valid series

!fTmp1 = 0.0
!fTmp2 = 0.0
!fTmp3 = 0.0
!fTmp4 = 0.0
!iTmp1 = 0
!iTmp2 = 0
!iTmp3 = 0
!iTmp4 = 0
!do i = 1, nzap
!  if(mark(i) > 0)then
!    fTmp1 = fTmp1 + frp(i)
!    iTmp1 = iTmp1 + 1
!  else
!    fTmp2 = fTmp2 + frp(i)
!    iTmp2 = iTmp2 + 1
!  endif
!enddo
!do i = 1, nTimeSeries
!  do j = 1, size(mmark,1)
!    if(mmark(j,i) > 0)then
!      fTmp3 = fTmp3 + mfrp(j,i)
!      iTmp3 = iTmp3 + 1
!    else
!      fTmp4 = fTmp4 + mfrp(j,i)
!      iTmp4 = iTmp4 + 1
!    endif
!  end do
!end do
!if(abs(fTmp3+fTmp_eliminated-FRP_grand_total)/FRP_grand_total > 0.001)then
!  call msg('')
!  call msg('After MAKING TIME SERIES-2, valid raw fires and FRP:',iTmp1, fTmp1)
!  call msg('... and eliminated fires and FRP:',iTmp2, fTmp2)
!  call msg('... and valid time series and FRP:',iTmp3, fTmp3)
!  call msg('... and void time series and FRP:',iTmp4, fTmp4)
!  call msg('Total N time series and FRP with eliminated:',iTmp1+iTmp2+iTmp3+iTmp4, fTmp3+fTmp_eliminated)
!  call msg('Problem at indFire:', indFire)
!  call set_error('High error 4','main')
!  stop
!endif
      endif

      indFire = indFire + 1
    end do  !
    !
    ! Get the max number of obs in time series and the actual number of time series
    !
    npnts_max = 1
    i = 0
    line = ''
    do indFire = 1, nTimeSeries
      npnts = 0
      iStatus = 0
      do j = 1, size(timeFire,1)
        if(mmark(j,indFire) > 0)then
          npnts = npnts + 1
          iStatus = 1
          if(len_trim(line) == 0) line = 'fire_' + &
                   & fu_time_fname_string_utc(fu_set_time_utc(fu_year(timeFire(j,indFire)), &
                                                            & fu_mon(timeFire(j,indFire)), &
                                                            & fu_day(timeFire(j,indFire)), &
                                                            & 0, 0, 0.0))
        endif
      end do
      i = i + iStatus
      if(npnts_max < npnts) npnts_max = npnts
    enddo

call msg('Number of time series and max length of time series:', i, npnts_max)


    write(uOut, '(A)')'FRP_DATASET_V1'
    write(uOut, '(2A)')'# ', trim(line)
    write(uOut, '(2A)')'# FRP_DATASET_V1 Created with is4fires from ', revision_str
    do iTmp1 = 1, size(arFRP_FNm)
      write(uOut, '(A,X,I3,A,A)')'# Input file ', iTmp1,":", trim(arFRP_FNm(iTmp1))
    enddo  

    write(uOut, '(A,I8)')'number_of_fires = ', i
    write(uOut, '(A,I8)')'max_number_of_same_fire_observations = ', npnts_max
!    write(uOut, '(A)')'fire_metadata_file = d:\!model\silam_v5_2\ini\fire_metadata.ini'

    write(uOut, '(A)') ""
    write(uOut, '(A)') chLabelStr

    !
    ! COunt the 
    !
    j = 1
    do indFire = 1, nTimeSeries
      iStatus = 0
      do i = 1, npnts_max
        if(mmark(i,indFire) > 0)then
          iStatus = 1
          write(uOut, &
              & fmt='(A,I7,1x,5(I4,1x),F4.1,1x,2(F9.4,1x),2(F7.3,1x),A,1x,F9.3,1x,A,1x,7(F7.3,1x))') &
                                         & 'fire = ', &
                                         & j, &  !indFire, &
                                         & fu_year(timeFire(i,indFire)), &
                                         & fu_mon(timeFire(i,indFire)), &
                                         & fu_day(timeFire(i,indFire)), &
                                         & fu_hour(timeFire(i,indFire)), &
                                         & fu_min(timeFire(i,indFire)), &
                                         & fu_sec(timeFire(i,indFire)), &
                                         & mlon(i,indFire), &
                                         & mlat(i,indFire), &
                                         & mdx(i,indFire), mdy(i,indFire), 'km', &
                                         & mfrp(i,indFire), 'MW', &
                                         & mT21(i,indFire), mT31(i,indFire), mT21mean(i,indFire), &
                                         & mT31mean(i,indFire), mDTmean(i,indFire), &
!                                         & 300., 300., 300., 300., 300., &
                                         & 0.5, 0.01
          FRP_SILAM_file = FRP_SILAM_file + mfrp(i,indFire)
        endif
      end do  ! npnts_max
      j = j + iStatus
    end do  ! nzap

    write(uOut, '(A)')'END_FRP_DATASET_V1'
    close(uOut)

    call msg('Day processed, FRP original and in SILAM file:', FRP_grand_total, FRP_SILAM_file)
    call msg('FRP in duplicates and total absolute error:', FRP_eaten, &
                                                   & FRP_grand_total - FRP_SILAM_file - FRP_eaten)
    if(abs(FRP_grand_total - FRP_SILAM_file - FRP_eaten) / (FRP_grand_total +1e-5) > 0.001)then
      call set_error('Absolute error > 0.1% of total FRP','main')
    endif


    deallocate(label, timeTmp)

    CONTAINS

    !=============================================================================

    logical function fu_if_cover(lon,lat,dx,dy, lon_new,lat_new)
      implicit none
      real, intent(in) :: lon,lat,dx,dy, lon_new,lat_new
      real :: dist_2

      dist_2 = ((fu_dx_deg_to_m(lon-lon_new,lat)/dx)**2 + &
              & (fu_dy_deg_to_m(lat-lat_new)/dy)**2) * 1e-6
      fu_if_cover = dist_2 < 0.0625  !0.25
      
    end function fu_if_cover


    !=============================================================================
    
    real function fu_fraction_covered(lon_cover, lat_cover, dx_cover, dy_cover, &
                                    & lon_new, lat_new, dx_new, dy_new)result(fract)
      !
      ! Roughly estimates a fraction of the _new pixel covered by the _cover pixel
      !
      implicit none

      real, intent(in) :: lon_cover, lat_cover, dx_cover, dy_cover, &
                        & lon_new, lat_new, dx_new, dy_new

      real :: dist_2, dx, dy

      dx = 0.5 * (dx_cover + dx_new)
      dy = 0.5 * (dy_cover + dy_new)
      dist_2 = ((fu_dx_deg_to_m(lon_cover-lon_new,lat_cover)/dx)**2 + &
            & (fu_dy_deg_to_m(lat_cover-lat_new)/dy)**2) * 1e-6

      fract = max(0.,min(1., dx_cover*dy_cover / (dx_new*dy_new) * (1.-dist_2)))

    end function fu_fraction_covered


    !=============================================================================

    subroutine remove_duplicates(timeFire, frp, lon, lat, T21, T31, T21mean, T31mean, DTmean, dx, dy, mark, nzap, FRP_removed)
      implicit none

      ! Imported parameters
      type(silja_time), dimension(:), intent(in) :: timeFire
      real, dimension(:), intent(inout) :: frp, lon, lat, T21, T31, T21mean, T31mean, DTmean, dx, dy
      integer, dimension(:), intent(inout) :: mark
      integer, intent(in) :: nzap
      real, intent(inout) :: FRP_removed

      ! Local variables
      integer :: i, j
      real :: fTmp

do i=1,nzap
if(mark(i) < 1) cycle  ! already excluded fire
if(isnan(frp(i)))call msg('nan-1',i)
if(frp(i) < 0.001)call msg('small frp 1',i,frp(i))
end do


!fTmp = 0.0
      do i=1,nzap
        if(mark(i) < 1) cycle  ! already excluded fire

!call msg('Next:'+fu_str(hhmm(i)),lon(i),lat(i))
        !
        ! Check if there are same-location same-frame fires below
        !
        do j = 1, nzap

!call msg('To check:'+fu_str(hhmm(j)),frp(j),mark(j))

          if(mark(j) < 1 .or. i==j)cycle ! already eaten or the same

          if(timeFire(i) == timeFire(j))then ! year(i)==year(j) .and. mon(i)==mon(j) .and. day(i)==day(j) .and. hhmm(i)==hhmm(j))then
            !
            ! Check for overlap
            !
            if(fu_if_cover(lon(i),lat(i),dx(i),dy(i), lon(j),lat(j))) then
              !
              ! Duplicate. Very strong one will always eat much smaller one. Also, if the small
              ! has not eaten more than one, it wil be eaten too. Otherwise, chain-resolving needed.
              !
if(frp(i)+frp(j) < 0.0001)then
call msg('Very small frps:',frp(i),frp(j))
endif
!call msg('Neighbours:',lon(j),lat(j))
              if(frp(i) > frp(j))then
                if(mark(j) < 3 .or. frp(i) > 5.*frp(j))then
                  mark(j) = 0
                  mark(i) = mark(i) + 1
!fTmp = fTMp + frp(j)
                  fTmp = fu_fraction_covered(lon(i),lat(i),dx(i),dy(i), lon(j),lat(j),dx(j),dy(j))
                  frp(i) = frp(i) + frp(j) * (1.-fTmp)
                  frp(j) = frp(j) * fTmp
                  FRP_removed = FRP_removed + frp(j)
!call msg('Neighbour has been eaten-1. FRP, tot:',frp(j),fTmp)
                else
                  call msg('Chain resolving -1. FRPs:',frp(i),frp(j))
                  if(mark(i) > 1 .or. frp(i) > 2.*frp(j))then  ! comparable weight. Average locations!
                    lon(j) = (lon(j)*frp(j) + lon(i)*frp(i)) / (frp(j) + frp(i))
                    lat(j) = (lat(j)*frp(j) + lat(i)*frp(i)) / (frp(j) + frp(i))
                    dx(j) = (dx(j)*frp(j) + dx(i)*frp(i)) / (frp(j) + frp(i))
                    dy(j) = (dy(j)*frp(j) + dy(i)*frp(i)) / (frp(j) + frp(i))
                  endif
                  fTmp = frp(j)    ! swap FRPs
                  frp(j) = frp(i)
                  frp(i) = fTmp
                  mark(j) = mark(j) + mark(i)
                  mark(i) = 0         ! whether the location is averaged or not, kill i-th pixel
                  fTmp = fu_fraction_covered(lon(j),lat(j),dx(j),dy(j), lon(i),lat(i),dx(i),dy(i))
                  frp(j) = frp(j) + frp(i) * (1.-fTmp)
                  frp(i) = frp(i) * fTmp
                  FRP_removed = FRP_removed + frp(i)
                endif
              else
                if(mark(i) < 2 .or. frp(j) > 5.*frp(i))then
                  mark(i) = 0           ! mark(i)>1 => fire has "eaten" some neighbour(s). Keep it!
                  mark(j) = mark(j) + 1 ! j-th pixel has eaten the i-th one
                  fTmp = fu_fraction_covered(lon(j),lat(j),dx(j),dy(j), lon(i),lat(i),dx(i),dy(i))
                  frp(j) = frp(j) + frp(i) * (1.-fTmp)
                  frp(i) = frp(i) * fTmp
                  FRP_removed = FRP_removed + frp(i)
!call msg('Neighbour has been eaten-2. FRPs, tot:',frp(i),fTmp)
!fTmp = fTMp + frp(i)
                  exit
                else
                  call msg('Chain resolving-2. FRPs:',frp(i),frp(j))
                  if(mark(j) > 1 .or. frp(j) > 2.*frp(i))then  ! comparable weight. Average locations!
                    lon(i) = (lon(j)*frp(j) + lon(i)*frp(i)) / (frp(j) + frp(i))
                    lat(i) = (lat(j)*frp(j) + lat(i)*frp(i)) / (frp(j) + frp(i))
                    dx(i) = (dx(j)*frp(j) + dx(i)*frp(i)) / (frp(j) + frp(i))
                    dy(i) = (dy(j)*frp(j) + dy(i)*frp(i)) / (frp(j) + frp(i))
                  endif
                  fTmp = frp(i)
                  frp(i) = frp(j)
                  frp(j) = fTmp
                  mark(i) = mark(i) + mark(j)
                  mark(j) = 0         ! whether the location is averaged or not, kill i-th pixel
                  fTmp = fu_fraction_covered(lon(i),lat(i),dx(i),dy(i), lon(j),lat(j),dx(j),dy(j))
                  frp(i) = frp(i) + frp(j) * (1.-fTmp)
                  frp(j) = frp(j) * fTmp
                  FRP_removed = FRP_removed + frp(j)
                endif
              endif  ! frp(i) <> frp(j)
            endif  ! duplicate
          end if  ! same frame
        end do  ! j: 1:nzap
      end do   ! i: 1:nzap


do i=1,nzap
if(mark(i) < 1) cycle  ! already excluded fire
if(isnan(frp(i)))call msg('nan-2',i)
if(frp(i) < 0.001)call msg('small frp 2',i,frp(i))
end do

    end subroutine remove_duplicates



    !=============================================================================
    
    subroutine make_time_series(timeTmp, frp, lon, lat, T21, T31, T21mean, T31mean, DTmean, dx, dy, mark, nzap, nTimeSeries)
      !
      ! Analyses the list of fires and stores them in a bunch of time series.
      ! Time series is a sequence of observations of one fire performed during one day.
      !
      implicit none

      ! Imported parameters
      type(silja_time), dimension(:), intent(inout) :: timeTmp
      real, dimension(:), intent(inout) :: frp, lon, lat, T21, T31, T21mean, T31mean, DTmean, dx, dy
      integer, dimension(:), intent(inout) :: mark
      integer, intent(in) :: nzap
      integer, intent(inout) :: nTimeSeries

      ! Local declarations
      integer :: npnts, indFire, i, j, iPoint

      indFire = nTimeSeries

do i=1,nzap
if(mark(i) < 1) cycle  ! already excluded fire
if(isnan(frp(i)))call msg('nan-3',i)
if(frp(i) < 0.001)call msg('small frp 3',i,frp(i))
end do


      do i=1,nzap
        if(mark(i) < 1) cycle
        npnts=1
        indFire = indFire + 1
        timeFire(npnts,indFire) = timeTmp(i) 
        mfrp(npnts,indFire)=frp(i)  ! / arDailyVariation(fu_hour_local_time(timeFire) + 1)
        mlon(npnts,indFire)=lon(i)
        mlat(npnts,indFire)=lat(i)
        mT21(npnts,indFire)=T21(i)
        mT31(npnts,indFire)=T31(i)
        mT21mean(npnts,indFire)=T21mean(i)
        mT31mean(npnts,indFire)=T31mean(i)
        mDTmean(npnts,indFire)=DTmean(i)
        mdx(npnts,indFire)=dx(i)
        mdy(npnts,indFire)=dy(i)
        mmark(npnts,indFire) = 1
        mark(i) = 0

!call msg('Next2:'+fu_str(mhhmm),mlon(npnts,indFire),mlat(npnts,indFire))

        !
        ! Check if there are same-location fires in different time-frames
        !
        do j=1,nzap
          if(mark(j) < 1)cycle
          if(i==j)cycle

          do iPoint = 1, npnts

            if(fu_if_cover(mlon(iPoint,indFire),mlat(iPoint,indFire), &
                         & mdx(iPoint,indFire), mdy(iPoint,indFire), lon(j),lat(j))) then

              if(timeTmp(j) == timeFire(iPoint,indFire)) then
call msg('N:',npnts+1)
call msg('Fire itself:'+fu_str(timeTmp(i)) + ', index:' + fu_str(i),lon(i),lat(i))
call msg('Checking nearby2:'+fu_str(timeTmp(j)) + ', index:' + fu_str(j),lon(j),lat(j))
call msg('Relative distance2:',((fu_dx_deg_to_m(lon(j)-lon(i),lat(i))/dx(i))**2 + &
                                        & (fu_dy_deg_to_m(lat(j)-lat(i))/dy(i))**2)*1e-6)
call msg('Longitude and dist,[km]:',mlon(1,indFire), 0.001*fu_dx_deg_to_m(lon(j)-mlon(1,indFire),mlat(1,indFire)))
call msg('Latitude and dist,[km]:',mlat(1,indFire), 0.001*fu_dy_deg_to_m(lat(j)-mlat(1,indFire)))
call msg('dx,dy,[km]:',dx(i),dy(i))
call msg('FRPs:',frp(i),frp(j))
call msg('Place fits, times:' + fu_str(timeTmp(j)) + &
                              & fu_str(timeFire(iPoint,indFire)))
                call set_error('problem with timing:'+fu_str(timeTmp(j)),'make_time_series')
                return
              endif

              npnts=npnts+1

              timeFire(npnts,indFire) = timeTmp(j)
              mfrp(npnts,indFire) = frp(j)
              mlat(npnts,indFire)=lat(j)
              mlon(npnts,indFire)=lon(j)
              mT21(npnts,indFire)=T21(j)
              mT31(npnts,indFire)=T31(j)
              mT21mean(npnts,indFire)=T21mean(j)
              mT31mean(npnts,indFire)=T31mean(j)
              mDTmean(npnts,indFire)=DTmean(j)
              mdx(npnts,indFire)=dx(j)
              mdy(npnts,indFire)=dy(j)
              mmark(npnts,indFire) = 1
              mark(j)=0                  ! mark as used
              exit
            end if
          enddo  ! npnts
        end do  ! j=1:nzap

      end do   ! i=1:nzap

      nTimeSeries = indFire

      if(iSkip > 0)then
        call msg('Skipped large cells, initial and resulting nbr of lines:' + &
                           & fu_str(iSkip),nzap,indFire)
      else
        call msg('Number of lines:',nzap,indFire)
      endif

do i=1,nzap
if(mark(i) < 1) cycle  ! already excluded fire
if(isnan(frp(i)))call msg('nan-4',i)
if(frp(i) < 0.001)call msg('small frp 4',i,frp(i))
end do

    end subroutine make_time_series

    
    !=============================================================================
    
    subroutine analyse_time_series(pTimeFire, pFrp, pLon, pLat, &
                                 & pT21, pT31, pT21mean, pT31mean, pDTmean, pdx, pdy, mark, ifSplit, nP)
      !
      ! Checks if several small pixels in one time slot are "collected" inside one large
      ! pixel in another time slot. 
      !
      ! Options are:
      ! - kill the large pixel and re-analyse the time line in order to split the fires. 
      !   +: more detailed land-use, -:loss of temporal information
      ! - split the large pixel by "projecting" it to small ones proportionaly to their FRPs
      !   +: no info loss, -: artificial trick, voluntary definition of pixels to split
      ! - sum-up small pixels. +: plume-rise for close plumes, -: loss of land-use
      !
      ! So far, splitting seems to be the best option. Some artificial noise is added but
      ! no information lost.
      !
      implicit none

      ! Imported parameters
      type(silja_time), dimension(:), intent(inout) :: pTimeFire
      real, dimension(:), intent(inout) :: pFrp, pLat, pLon, pT21, pT31, pT21mean, pT31mean, pDTmean, pdx, pdy
      integer, dimension(:), intent(inout) :: mark
      logical, intent(out) :: ifSplit
      integer, intent(out) :: nP

      ! Local variables
      integer :: iFire, j, iArMax, nSTF, nTimes, iLargestSet, iPixelToSplit, maxPixelsCovered
      type(silja_time), dimension(50) :: times
      integer, dimension(50) :: nSameTimeFires, nCoveredFires
      integer, dimension(50,50) :: pSameTimeFire, pCoveredFire
      real :: frpTot

      !
      ! Get the size of time line and pixel size variability
      !
      nP = 1
      nSTF = 1  ! same-time fires
      ifSplit = .false.
      nSameTimeFires = 0
      pSameTimeFire = 0
      pCoveredFire = 0
      nCoveredFires = 0
      nTimes = 0
      do nP = 1, size(pFrp)
        if(mark(nP) == 0)exit
      end do ! iFire
      nP = nP - 1
      !
      ! Is this time series a collection around one big pixel?
      !
      do iFire = 1, nP
        !
        ! Get the time, check whether it is duplicated
        ! Also, get the pixels that are covered by this one
        !
        do j = 1, nP
          if(pTimeFire(j) == pTimeFire(iFire))then         ! same time fires
            nSameTimeFires(iFire) = nSameTimeFires(iFire) + 1
            pSameTimeFire(iFire,nSameTimeFires(iFire)) = j
          endif
        enddo
      end do
      !
      ! Take the largest set of same-time fires 
      !
      iLargestSet = 1
      do iFire = 1, nP
        if(nSameTimeFires(iFire) > nSameTimeFires(iLargestSet)) iLargestSet = iFire
      end do
      !
      ! Finally, we can decide whether this series is worth splitting
      !
      if(nSameTimeFires(iLargestSet) < 2) return
      !
      ! There is such set. Find the pixels that cover them
      !
      do iFire = 1, nSameTimeFires(iLargestSet)          ! scanning the largest set
        !
        ! Who covers them? Apart from themselves, of course.
        !
        do j = 1, nP
          if(j == pSameTimeFire(iLargestSet,iFire))cycle

          if(fu_if_cover(pLon(j),pLat(j),pdx(j),pdy(j), &
                       & pLon(pSameTimeFire(iLargestSet,iFire)), &
                       & pLat(pSameTimeFire(iLargestSet,iFire)))) then
            nCoveredFires(j) = nCoveredFires(j) + 1
            pCoveredFire(j,nCoveredFires(j)) = pSameTimeFire(iLargestSet,iFire)
          endif
        enddo
      end do  ! same-time fires of the largest set
      !
      ! Find the pixel that covers max fraction of this set. Note that pCoveredFires refer
      ! only to this subset.
      !
      iPixelToSplit = 0
      maxPixelsCovered = 0
      do j = 1, nP
        if(nCoveredFires(j) > maxPixelsCovered)then
          maxPixelsCovered = nCoveredFires(j)
          iPixelToSplit = j
        endif
      end do
      if(maxPixelsCovered < 1)then
        call msg('Same-time fires are not covered by anyone. Split the set!')
        ifSplit = .true.
        return
      endif

      ifSplit = .true.
      !
      ! Last check for weirds: the pixel to split must be larger than all those to be 
      ! projected onto it. If this is not true, average all the covered pixels
      ! and be done with it.
      !
      do iFire = 1, nCoveredFires(iPixelToSplit)
        if(pdx(pCoveredFire(iPixelToSplit,iFire)) > pdx(iPixelToSplit) .or. &
         & pdx(pCoveredFire(iPixelToSplit,iFire)) > pdx(iPixelToSplit))then
          !
          ! Eliminate all covered pixels: sum them up
          !
          frpTot = 0.
          do j = 1, nCoveredFires(iPixelToSplit)
            frpTot = frpTot + pFrp(pCoveredFire(iPixelToSplit,j))
          end do
call msg('Summing-up N same-time pixels, adding to the end:',nCoveredFires(iPixelToSplit))
          pTimeFire(nP + 1) = pTimeFire(pCoveredFire(iPixelToSplit,1))
          pFrp(nP + 1) = frpTot
          pLat(nP + 1) = pLat(iPixelToSplit)
          pLon(nP + 1) = pLon(iPixelToSplit)
          pT21(nP + 1) = pT21(iPixelToSplit)
          pT31(nP + 1) = pT31(iPixelToSplit)
          pT21mean(nP + 1) = pT21mean(iPixelToSplit)
          pT31mean(nP + 1) = pT31mean(iPixelToSplit)
          pDTmean(nP + 1) = pDTmean(iPixelToSplit)
          pdx(nP + 1) = pdx(iPixelToSplit)
          pdy(nP + 1) = pdy(iPixelToSplit)
          mark(nP +1) = 1
          do j = 1, nCoveredFires(iPixelToSplit)
            mark(pCoveredFire(iPixelToSplit,j)) = 0
          enddo
          nP = nP + 1
          return
        endif
      end do  ! n covered fires

      !
      ! Finally, we have the pixel that will be split and the pixels that are from the same time
      ! and covered by this one. They will be projected to it.
      !
      frpTot = 0.
      do iFire = 1, nCoveredFires(iPixelToSplit)
        frpTot = frpTot + pFrp(pCoveredFire(iPixelToSplit,iFire))
      end do
      !
      ! Split the pixel
      !
call msg('Splitting the pixel, adding to the end:',iPixelToSplit, nP)
      do iFire = 1, nCoveredFires(iPixelToSplit)
        pTimeFire(nP + iFire) = pTimeFire(iPixelToSplit)
        pFrp(nP + iFire) = pFrp(iPixelToSplit) / frpTot * pFrp(pCoveredFire(iPixelToSplit,iFire))
        pLat(nP + iFire) = pLat(pCoveredFire(iPixelToSplit,iFire))
        pLon(nP + iFire) = pLon(pCoveredFire(iPixelToSplit,iFire))
        pT21(nP + iFire) = pT21(pCoveredFire(iPixelToSplit,iFire))
        pT31(nP + iFire) = pT31(pCoveredFire(iPixelToSplit,iFire))
        pT21mean(nP + iFire) = pT21mean(pCoveredFire(iPixelToSplit,iFire))
        pT31mean(nP + iFire) = pT31mean(pCoveredFire(iPixelToSplit,iFire))
        pDTmean(nP + iFire) = pDTmean(pCoveredFire(iPixelToSplit,iFire))
        pdx(nP + iFire) = pdx(pCoveredFire(iPixelToSplit,iFire))
        pdy(nP + iFire) = pdy(pCoveredFire(iPixelToSplit,iFire))
        mark(nP + iFire) = 1
!call msg('New pixel:'+fu_time_to_string(pTimeFire(nP+iFire))+', FRP:' + &
!                                        & fu_real2str(pFrp(nP+iFire)), pLon(nP+iFire), pLat(nP+iFire))
      enddo
call msg('Eliminating pixel:'+fu_str(pTimeFire(iPixelToSplit))+', FRP:' + &
                          & fu_str(pFrp(iPixelToSplit)), pLon(iPixelToSplit),pLat(iPixelToSplit))
      mark(iPixelToSplit) = 0

      nP = nP + nCoveredFires(iPixelToSplit)

    end subroutine analyse_time_series


  end subroutine make_silam_fire_src_v1


  !*******************************************************************************
  
  subroutine mask_wrong_fires(SrcNames, chOutDir, chMaskFNm, nYearsMax)
    implicit none
    ! Imported parameters
    character(len=*), dimension(:), intent(in) :: SrcNames
    character(len=*), intent(in) ::  chOutDir, chMaskFNm
    integer , intent(in) :: nYearsMax
    ! Local variables
    type(silja_field_id) :: id
    type(silja_grid) :: MaskGrid
    integer :: iUnit, igf, nx, ny, i1d, iUnitSrc, iUnitDest,iFile, status
    real :: x,y, lon, lat
    real, dimension(:), pointer :: mask_array
    character(len = fnlen) :: chLine, chSrcFNm, chTmp
    character(len = clen) :: chLabelStr = '#  index  -date-  lon  lat   dx dy unit  FRP   TA  T4  T4b T11 T11b  MCE  areaFr'
    
    call init_grads_io(10) 
    !
    ! Get the mask field
    !
    igf =  fu_open_gradsfile_i(chMaskFNm)
    MaskGrid = fu_silamGrid_of_grads(igf)
    call grid_dimensions(MaskGrid,nx,ny)
    mask_array => fu_work_array(nx*ny)
    mask_array = 0.0

    !call read_field_from_grads_indices(igf, indVar_, indLev_, indTime_, grid_data, fill_value_)
    call read_field_from_grads_indices(igf, 1, 1, 1, mask_array, 0.)
    call close_gradsfile_i(igf)
    !
    ! The rest is simple: read the fire files through copying the lines if they do not fall into 
    ! forbidden cells and commenting them out if they do.
    !
    do iFile=1,size(SrcNames)
      call msg("processing source of total ",  iFile, size(SrcNames) )
      chSrcFNm = SrcNames(iFile)
      if(len_trim(chSrcFNm) < 2)cycle
      iUnitSrc = fu_next_free_unit()
      call msg('Starting:' + chSrcFNm)
      open(iUnitSrc,file=chSrcFNm,status = 'old')
      iUnitDest = fu_next_free_unit()
      call msg(chOutDir + dir_slash + chSrcFNm(index(chSrcFNm,dir_slash,.true.)+1:len_trim(chSrcFNm)))
      open(iUnitDest,file = chOutDir + dir_slash + chSrcFNm(index(chSrcFNm,dir_slash,.true.)+1:len_trim(chSrcFNm)))
      write(iUnitDest,'(A)')'# mask_wrong_fires observed more than '//fu_str(nYearsMax)//' years'
      do while (.True.)
        READ(unit = iUnitSrc, fmt = '(A)', iostat = status) chLine
        if (status /= 0) exit
         !        call next_line_from_input_file(iUnitSrc,chLine,eof2)
        if(index(chLine,'fire = ') > 0) then
          read(unit=chLine,fmt=*)chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,lon,lat
          call project_point_to_grid(lon, lat, MaskGrid, x, y)
          i1d = nint(x) + (nint(y) - 1) * nx
          if(mask_array(i1d) > nYearsMax)then
            call msg('Non-fire: ' + chLine)
            write(iUnitDest,'(A,A)')'##$$##$$##  ',trim(chLIne)
          else
            write(iUnitDest,'(A)')trim(chLine)
          endif
        else
          write(iUnitDest,'(A)')trim(chLine)
        endif  ! fire line
      end do
      close(iUnitSrc)
    end do
    
  end subroutine mask_wrong_fires
  
  
  !*****************************************************************************************************

  subroutine find_wrong_fires(FireSrcTemplIn, timeStart, timeEnd, chFrequentFireMaskFNm, &
                            & iAnnualCountThresh) !, nYearsAboveThreshold)
    implicit none
    ! Imported parameters
    type(grads_template), intent(in) :: FireSrcTemplIn
    type(silja_time), intent(in) :: timeStart, timeEnd
    character(len=*), intent(in) :: chFrequentFireMaskFNm
    integer, intent(in) :: iAnnualCountThresh !, nYearsAboveThreshold
  
    ! Local variables
    
    type(silja_field_id) :: id
    type(silja_time) :: now
    type(silja_grid) :: gridFire
    integer :: i1d, year, nYears, indYear, startYear, endYear, ix, iy
    real :: x,y, lon, lat, frp
    logical :: eof2, ifFirst, ifFireLine
    real, dimension(:,:), pointer :: count_array, frp_array
    integer, dimension(:), pointer :: arUnit
    integer*1, dimension(:), pointer :: cntTmp
    real, dimension(:), pointer :: frpTmp, f_cntTmp
    character(len = fnlen) :: chLine, chSrcFNm, chTmp, chOutDir
    character(len = clen) :: chLabelStr = '#         index  -date-                              lon       lat      dx    dy unit     FRP         T4    T4b      T11    T11b      TA     MCE  areaFr'
                                           ! fire =       1 2009    1    1    0   10  0.0  151.2433  -24.9140   1.004 1.002 km    40.874 MW 346.752 309.340 316.984 306.880  10.104   0.500   0.010
    real, parameter :: xStart = -179.985, dx = 0.03, yStart = -89.985, dy = 0.03
    integer, parameter :: nx = 12000, ny = 6000
  
    !
    ! Make the fire mask field: create the grid and reserve space for it
    !
    gridFire = fu_set_grid('Fire grid', lonlat, pole_geographical, &
                         & xStart, yStart, & !-179.985, -89.985, &      !sw_corner_modlon, sw_corner_modlat, &
                         & nx, ny, dx, dy) !12000, 6000, 0.03, 0.03)  !nx, ny, dx, dy)
    
    chOutDir = chFrequentFireMaskFNm(1:index(chFrequentFireMaskFNm,dir_slash,.true.)-1)
    call create_directory_tree(chOutDir)
    
    startYear = fu_year(timeStart)
    endYear = fu_year(timeEnd)
    nYears = endYear - startYear + 1
    
    allocate(count_array(nYears, nx*ny), frp_array(nYears, nx*ny), &
        & cntTmp(nx*ny), f_cntTmp(nx*ny), frpTmp(nx*ny), stat=i1d)
    count_array = 0.
    frp_array = 0.
    if(fu_fails(i1d==0,'Failed mask allocation','find_wrong_fires'))return
    !
    ! The rest is simple: read the fire files projecting the points to the grid and increasing the counters
    !
    arUnit => fu_work_int_array()
    if(error)return

    now = timeStart
    do while(now <= timeEnd)
      chSrcFNm = fu_FNm(FireSrcTemplIn, now, now, zero_interval)
      call msg(chSrcFNm)
      !
      ! Read this source file
      !
      arUnit(101) = fu_next_free_unit()
      open(arUnit(101),file=chSrcFNm,status = 'old',iostat=ix)
      if(fu_fails(ix == 0,'Absent file:'+chSrcFNm,'find_wrong_fires'))return
      eof2 = .false.
      cntTmp = 0
      frpTmp = 0.
      do while(.not.eof2)
        call next_line_from_input_file(arUnit(101),chLine,eof2)
        if(index(chLine,'fire = ') > 0) then
          read(unit=chLine,fmt=*)chTmp,chTmp,chTmp,year,chTmp,chTmp,chTmp,chTmp,chTmp,lon,lat, chTMp, chTmp, chTmp, frp
          if(year < startYear .or. year > endYear)exit
          indYear = year - startYear + 1
          call project_point_to_grid(lon, lat, gridFire, x, y)
          i1d = nint(x) + (nint(y) - 1) * nx
          cntTmp(i1d)  = cntTmp(i1d) + 1
          frpTmp(i1d) = frpTmp(i1d) + frp
        endif  ! fire line
      end do
      do i1d = 1, nx*ny
        if(cntTmp(i1d) > 0)then
          count_array(indYear,i1d)  = count_array(indYear,i1d) + 1.0     ! no more than 1 fire per day
          frp_array(indYear,i1d) = frp_array(indYear,i1d) + frpTmp(i1d) / real(cntTmp(i1d))  ! mean daily frp
        endif
      end do
      close(arUnit(101))
      now = now + one_day
    end do    ! fire files list <=> considered period
    close(arUnit(100))
    !
    ! Fires are counted. Write them down to grads file
    !
    call msg('Writing grads...')
    open(arUnit(1),file=chFrequentFireMaskFNm + '_count.grads',access='direct',recl=nx*ny*4_8,form='unformatted')
    i1d = nx*ny
    do indYear = 1, nYears
      write(arUnit(1),rec=indYear)count_array(indYEar,1:i1d)
    end do
    close(arUnit(1))
    open(arUnit(1),file=chFrequentFireMaskFNm + '_frp.grads',access='direct',recl=nx*ny*4_8,form='unformatted')
    do indYear = 1, nYears
      write(arUnit(1),rec=indYear)frp_array(indYEar,1:i1d)
    end do
    close(arUnit(1))
    
    call write_ctl(chFrequentFireMaskFNm + '_count.grads', 'fire_count', "Yearly fire count", nYears)
    call write_ctl(chFrequentFireMaskFNm + '_frp.grads', 'frp', "Yearly total FRP", nYears)

    !
    ! Process what has been read and recorded
    ! Stage 1. Get the cells with > threshold fires for each year
    !
    do indYear = 1, nYears
      call msg('Count stuff above threshold for year:',indYear)
      open(arUnit(1),file=chFrequentFireMaskFNm + '_' + fu_str(indYear + startYear - 1))
      write(arUnit(1),*)'Fires with counts above threshod:', iAnnualCountThresh, 'lon, lat, count, total frp_MW'
      do iy = 1, ny
        do ix = 1, nx
          i1d = ix + (iy-1)*nx
          if(nint(count_array(indYear,i1d)) <= iAnnualCountThresh)cycle
          write(arUnit(1),'(2(F9.4,2x),i5,1x,f9.2)') fu_lon_geographical_from_grid(real(ix), real(iy), gridFire), &
                       & fu_lat_geographical_from_grid(real(ix), real(iy), gridFire), &
                       & nint(count_array(indYear,i1d)), frp_array(indYear,i1d)
        end do
      enddo
      close (arUnit(1))
    enddo  ! years
    !
    ! If we have several years in line, process them all and list the points with >1 year high counts
    !
    if(nYears == 1)return
  
    call msg('Count cells across years')
    cntTmp = 0
    frpTmp = 0.
    do i1d = 1, nx*ny
      do indYear = 1, nYears
        if(count_array(indYear,i1d) > iAnnualCountThresh)then
          cntTmp(i1d) = cntTmp(i1d) + 1
          frpTmp(i1d) = frpTmp(i1d) + frp_array(indYear,i1d)
        endif
      end do
    end do  ! i1d

    ! Write year-count array
    f_cntTmp(1:nx*ny) = real(cntTmp(1:nx*ny))
    open(arUnit(1),file=chFrequentFireMaskFNm + '_yrs_above_thr.grads',access='direct',recl=nx*ny*4_8,form='unformatted')
    write(arUnit(1),rec=1) f_cntTmp(1:nx*ny)
    close(arUnit(1))
    call write_ctl(chFrequentFireMaskFNm + '_yrs_above_thr.grads', 'yrcnt', &
          & "Years with nFires above "//trim(fu_str(iAnnualCountThresh)),1)

    !
    ! Now simply write down a few files with counts etc
    !
    do indYear = 1, nYears
      call msg('Writing down cumulatives for N years:',indYear)
  
      open(arUnit(1),file=chFrequentFireMaskFNm + '_seen_in_' + fu_str(indYear) + '_years_or_more.txt')
      write(arUnit(1),'(A,i3,A,i5,A)')'Cells FRP with fires seen >=', iAnnualCountThresh, ' times in at least', indYear, ' years. Lon,lat,frp'
      arUnit(2) = fu_next_free_unit()
      open(arUnit(2),file=chFrequentFireMaskFNm + '_seen_in_' + fu_str(indYear) + '_years_or_more.kml')
      write(arUnit(2),'(A)')'<kml xmlns="http://earth.google.com/kml/2.2">'
      write(arUnit(2),'(A)')'<Document>'
      arUnit(3) = fu_next_free_unit()
      open(arUnit(3),file=chFrequentFireMaskFNm + '_seen_in_' + fu_str(indYear) + '_years.txt')
      write(arUnit(3),'(A,i3,A,i5,A)')'Cells FRP with fires seen >=', iAnnualCountThresh, ' times in', indYear, ' years. Lon,lat,frp'
      arUnit(4) = fu_next_free_unit()
      open(arUnit(4),file=chFrequentFireMaskFNm + '_seen_in_' + fu_str(indYear) + '_years.kml')
      write(arUnit(4),'(A)')'<kml xmlns="http://earth.google.com/kml/2.2">'
      write(arUnit(4),'(A)')'  <Document>'
      
      do iy = 1, ny
        do ix = 1, nx
          i1d = ix + (iy-1)*nx
          !
          ! Many fires seen in at least N years
          !
          if(cntTmp(i1d) >= indYear)then
            write(arUnit(1),'(2(F9.4,2x),f10.3)') fu_lon_geographical_from_grid(real(ix), real(iy), gridFire), &
                                                & fu_lat_geographical_from_grid(real(ix), real(iy), gridFire), &
                                                & frpTmp(i1d)
            write(arUnit(2),'(A)')'    <Placemark>'
            write(arUnit(2),'(A)')'      <name>'
            write(arUnit(2),'(A)') '         _' + fu_str(int(cntTmp(i1d))) + 'yrs_' + fu_str(int(frpTmp(i1d))) + 'MWday'
            write(arUnit(2),'(A)')'      </name>'
            write(arUnit(2),'(A)')'      <Point>'
            write(arUnit(2),'(A)')'        <coordinates>'
            write(arUnit(2),'(10x,2(F9.4,2x))') fu_lon_geographical_from_grid(real(ix), real(iy), gridFire), &
                                           & fu_lat_geographical_from_grid(real(ix), real(iy), gridFire)
            write(arUnit(2),'(A)')'        </coordinates>'
            write(arUnit(2),'(A)')'      </Point>'
            write(arUnit(2),'(A,I3,F10.3,A)')'      <description />'
            write(arUnit(2),'(A)')'    </Placemark> '
            !
            ! Many fires seen in exactly N years
            !
            if(cntTmp(i1d) == indYear)then
              write(arUnit(3),'(2(F9.4,2x),f10.3)') fu_lon_geographical_from_grid(real(ix), real(iy), gridFire), &
                                               & fu_lat_geographical_from_grid(real(ix), real(iy), gridFire), &
                                               & frpTmp(i1d)
              write(arUnit(4),'(A)')'    <Placemark>'
              write(arUnit(4),'(A)')'      <name>'
              write(arUnit(4),'(A)') fu_str(int(frpTmp(i1d))) + 'MWday'
              write(arUnit(4),'(A)')      '</name>'
              write(arUnit(4),'(A)')'      <Point>'
              write(arUnit(4),'(A)')'        <coordinates>'
              write(arUnit(4),'(10x,2(F9.4,2x))') fu_lon_geographical_from_grid(real(ix), real(iy), gridFire), &
                                             & fu_lat_geographical_from_grid(real(ix), real(iy), gridFire)
              write(arUnit(4),'(A)')'        </coordinates>'
              write(arUnit(4),'(A)')'      </Point>'
              write(arUnit(4),'(A,F10.3,A)')'      <description />'
              write(arUnit(4),'(A)')'    </Placemark> '
            endif
          endif
        end do
      end do
      write(arUnit(2),'(A)')'  </Document>'
      write(arUnit(2),'(A)')'</kml> '
      write(arUnit(4),'(A)')'  </Document>'
      write(arUnit(4),'(A)')'</kml> '
      close(arUnit(1))
      close(arUnit(2))
      close(arUnit(3))
      close(arUnit(4))
    end do
    !
    ! The rest is simple: read the fire files one by one copying the lines if they do not fall into 
    ! forbidden cells and commenting them out if they do.
    ! To save the pain, we shall write many output directories: for >1, >2, >3, etc >12 years counts,
    ! then choose the right ones.
    !
    now = timeStart
    do while(now <= timeEnd)
      chSrcFNm = fu_FNm(FireSrcTemplIn, now, now, zero_interval)
      if(len_trim(chSrcFNm) < 2)cycle
      arUnit(101) = fu_next_free_unit()
      call msg('Starting:' + chSrcFNm)
      open(arUnit(101),file=chSrcFNm,status = 'old', action='read')
      do indYear = 1, nYears
        arUnit(indYear) = fu_next_free_unit()
        call create_directory_tree(chOutDir + '_exclude_' + fu_str(indYear) + '_yrs_or_more')
        open(arUnit(indYear),file = chOutDir + '_exclude_' + fu_str(indYear) + '_yrs_or_more' + dir_slash + chSrcFNm(index(chSrcFNm,dir_slash,.true.)+1:len_trim(chSrcFNm)), action='write')
      end do
      call msg(chOutDir + dir_slash + chSrcFNm(index(chSrcFNm,dir_slash,.true.)+1:len_trim(chSrcFNm)))
      eof2 = .false.
      ifFIrst = .true.
      do while(.not.eof2)
        call next_line_from_input_file(arUnit(101),chLine,eof2)
        ifFireLine = index(chLine,'fire = ') > 0
        if(ifFireLine) then
          read(unit=chLine,fmt=*)chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,chTmp,lon,lat
          call project_point_to_grid(lon, lat, gridFire, x, y)
          i1d = nint(x) + (nint(y) - 1) * nx
        endif
        do indYear = 1, nYears
          if(ifFirst) write(arUnit(indYear),'(A)')trim(chLabelStr)
          if(ifFireLine)then             ! not automatically but after checking...
            if(cntTmp(i1d) >= indYear)then
              if(cntTmp(i1d) == indYear)call msg('Non-fire:' + fu_str(indYear) + ':' + chLine)
              write(arUnit(indYear),'(A,A)')'##$$##$$##  ',trim(chLIne)
            else
              write(arUnit(indYear),'(A)')trim(chLine)
            endif
          else
              write(arUnit(indYear),'(A)')trim(chLine)
          endif  ! ifFireLine
        end do  ! nYears or more
        ifFirst = .false.
      end do  ! eof2
      close(101)
      do indYear = 1, nYears
        close(arUnit(indYear))
      enddo
      now = now + one_day
    end do  !  source file lines
    
  contains
                            
    subroutine write_ctl(datFNm, chVar, chLongVar, nT)
      implicit none
      character(len=*), intent(in) :: datFNm, chVar, chLongVar
      integer, intent(in) :: nT
      
      integer :: iUnit,iTmp
      

      iUnit = fu_next_free_unit()

      ! First character of the basename
      iTmp=index(datFNm, dir_slash, BACK=.True.)
      if (iTmp < 1) then
        iTmp = 1 
      else
        iTmp = iTmp + 1
      endif

      open(iUnit, file=datFNm + '.ctl')
      WRITE(iUnit,'(A)') 'DSET  ^'//trim(datFNM(iTmp:))
      WRITE(iUnit,'(A)') 'TITLE SILAM GrADS output'
      WRITE(iUnit,'(A)') 'OPTIONS LITTLE_ENDIAN'
      WRITE(iUnit,'(A)') 'UNDEF  0'
      WRITE(iUnit,'(A5,I5,A8,2(F15.7,1x))') 'XDEF ',nx,' LINEAR ', xStart, dx
      WRITE(iUnit,'(A5,I5,A8,2(F15.7,1x))') 'YDEF ',ny,' LINEAR ', yStart, dy
      WRITE(iUnit,'(A)') 'ZDEF 1 LEVELS 0'
      WRITE(iUnit,'(A,I3,3A)') 'TDEF ', nT,' LINEAR ', fu_time_to_grads_string(timeStart), ' 1yr'
      WRITE(iUnit,'(A)')'VARS 1'
      WRITE(iUnit,'(A,A,A)') trim(chVar), ' 0 99 99 0 ', trim(chLongVar)
      WRITE(iUnit,'(A7)') 'ENDVARS'
      close(iUnit)

      open(iUnit, file=datFNm + '.super_ctl')
      WRITE(iUnit,'(A)') 'LIST = general'
      WRITE(iUnit,'(A,A)') 'ctl_file_name = ^'//trim(datFNM(iTmp:))//'.ctl'
      WRITE(iUnit,'(A)') ' grid_title = global fire grid'
      WRITE(iUnit,'(A)') ' grid_type = LON_LAT'
      WRITE(iUnit,'(A)') ' lon_s_pole =   0.0'
      WRITE(iUnit,'(A)') ' lat_s_pole =   -90.0'
      WRITE(iUnit,'(A)') ' lon_start =    '//fu_str(xStart)
      WRITE(iUnit,'(A)') ' lat_start =    '//fu_str(yStart)
      WRITE(iUnit,'(A,I6)') ' nx = ', nx
      WRITE(iUnit,'(A,I6)') ' ny = ', ny
      WRITE(iUnit,'(A,F15.7)') ' dx = ', dx
      WRITE(iUnit,'(A,F15.7)') ' dy = ', dy
      WRITE(iUnit,'(A)') ' resol_flag = 128'
      WRITE(iUnit,'(A)') ' ifReduced = 0'
      WRITE(iUnit,'(A)') ' earth_flag = 0'
      WRITE(iUnit,'(A)') ' wind_component = 0 '
      WRITE(iUnit,'(A)') ' reduced_nbr_str = 0'
      WRITE(iUnit,'(A)') ' lat_pole_stretch = 0. '
      WRITE(iUnit,'(A)') ' lon_pole_stretch = 0.'
      WRITE(iUnit,'(A)') ' number_of_levels = 1'
      WRITE(iUnit,'(A)') ' vertical_method = SURFACE_LEVEL'
      WRITE(iUnit,'(A)') ' level_type = SURFACE_LEVEL'
      WRITE(iUnit,'(A)') ' time_label_position  =  END_of_period  '
      WRITE(iUnit,'(A)') ' data_time_features = static_data'
      WRITE(iUnit,'(A)') 'END_LIST = general'
      WRITE(iUnit,'(A)') 'LIST = '//trim(chVar)
      WRITE(iUnit,'(A)') '   quantity_short_name = '//trim(chVar)
      WRITE(iUnit,'(A)') 'END_LIST = '//trim(chVar)
      close(iUnit)
    end subroutine write_ctl
    
  end subroutine find_wrong_fires


END MODULE frp_tools

