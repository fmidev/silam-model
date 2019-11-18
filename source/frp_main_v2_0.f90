!****************************************************************************
!
!  PROGRAM: FRP  processor
!
!  PURPOSE:  reads FRP, makes freqient_fire mask, applies it
!
!****************************************************************************

program is4fires20

  use frp_tools
  use grids_geo
  use grads_templates

  implicit none

  character(len=*), parameter :: usage = "progname make_src infile [infile [...]] outfile.fs1"//NEW_LINE('A')//&
                &" or "//NEW_LINE('A')//&
                &"progname make_mask  FireSrcTempl_all_fires yymmddhh-start yymmddhh-end  chFrequentFireMaskFNm"//NEW_LINE('A')//&
                &" or "//NEW_LINE('A')//&
                &"progname parse_modis_ta infile outfile"//NEW_LINE('A')//&
                &" or "//NEW_LINE('A')//&
                &"progname mask_fires maskfile infile.fs1 [infile [...]] outdir"


  character(len=clen) :: chCaseNm
  character(len=fnlen) :: chHDFParsingExeFNm, chEmisCoefFileBin, &
                         & chIniFNm, chCocktailName, chReleaseRateUnit, &
                         & chEcotype, chAverFuelLoad, chOutDir, chIntermediateTemplate, &
                         & chSilamEmissionTemplate, chGradsOutputTemplate
  character(len=fnlen), dimension(:), allocatable :: chInputFNames, arFRP_FNm
  real :: scaling_factor, fMinFileSizeRelative, fMaxFileSizeRelative
  integer :: iTmp, Nargs,  nFireDaysThreshold !, nYearsAboveThreshold

  type(silja_time) :: timeStart, timeEnd
!  type(grads_template) :: logTempl
!  type(grads_template), dimension(:), allocatable :: FRPtempl
  type(grads_template) :: FireSrcTempl_all_fires
!  logical :: ifParseHDF, ifMakeSILAMFireSource,  ifMaskWrongFires, ifUseOldMask
  !logical :: ifParseHDF, ifProcessFRP, ifConvert_2_emis, ifAnalyse_emis_files, &
  !         & ifRemoveSmallFiles, ifRemoveLargeFiles, ifFillGaps,  ifMaskWrongFires
!  real, dimension(:), pointer :: daily_totals


  !----------------------------------------
  !
  !  Open the global error and warning recording file
  !  This very file will be used for error messages, warnings and, later,
  !  for other log functions
  !  However, when the output directory becomes known - it will be transferred to it
  !  and its number will be changed
  !


  run_log_funit = fu_next_free_unit()
  open(run_log_funit, file="/dev/null", action="write") !No need to keep a separate log


  Nargs = command_argument_count()

  

  if (Nargs < 3 ) then !--------- No arguments
    CALL get_command_argument(0, chIniFNm)
    print *, chIniFNm, " called with too little arguments"
    print *, usage
    call exit(-1)
  END IF

  !Read the arguments
  allocate(chInputFNames(Nargs))!, stat=iStat)
  do iTmp=1,Nargs
    call get_command_argument(iTmp, chInputFNames(iTmp))
    !    print *, "arg", iTmp, ":", trim(chInputFNames(iTmp))
  enddo


  select case (trim(chInputFNames(1)))
    case ("make_src")
                      !Input files for one day              !Output file

      call make_silam_fire_src_v1(chInputFNames(2:Nargs-1), chInputFNames(Nargs))

    case ("make_mask")

      !  frp20  make_mask  FireSrcTempl_all_fires yymmddhh-start yymmddhh-end  chFrequentFireMaskFNm

      call decode_template_string(chInputFNames(2), FireSrcTempl_all_fires)   ! all-fire source template
      timeStart = fu_hirlam_string_to_time(chInputFNames(3))
      timeEnd   = fu_hirlam_string_to_time(chInputFNames(4))
      nFireDaysThreshold = 50
      if (error) call exit(-1)

      call find_wrong_fires( FireSrcTempl_all_fires, &    ! all-fire source template
                          & timeStart, timeEnd, &
                          & chInputFNames(5), &  ! chFrequentFireMaskFNm
                          & nFireDaysThreshold)  !, nYearsAboveThreshold)

!!!!!      call find_wrong_fires(chFireSourceOutDir, FireSrcTempl_all_fires, &    ! all-fire source template
!!!!!                          & timeStart, timeEnd, &
!!!!!                          & chFrequentFireMaskFNm, &
!!!!!                          & nFireDaysThreshold)  !, nYearsAboveThreshold)

    case ("mask_fires")
      call mask_wrong_fires(chInputFNames(3:Nargs-1), & ! Fire sources unmasked
                          & chInputFNames(Nargs), &     ! Output dir
                          & chInputFNames(2), 6)           ! chFrequentFireMaskFNm

    case ("parse_modis_ta")
      if (Nargs /= 3) then
        CALL get_command_argument(0, chIniFNm)
        print *, trim(chIniFNm), "parse_modis_ta called with strange arguments"
        print *, "Usage" 
        print *, usage
        call exit(-1)
      endif
        
      call parse_modis_ta(chInputFNames(2), & ! Fire sources unmasked
                          & chInputFNames(3))

    case default
      CALL get_command_argument(0, chIniFNm)
      print *, trim(chIniFNm), " called with strange arguments"
      print *, usage
      call exit(-1)
  end select
  if(error) call exit(-1)



  



  ! call create_directory_tree(chFireSourceOutDir)
  stop
  

  


end program is4fires20

