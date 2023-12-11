MODULE dispersion_server
!
! This module includes several services for the dispersion model:
!
! 1. Interpolation of maps of cocktails (interpolation of plain fields is
!    covered in field_buffer)
! 2. Initialization routines for maps of cocktails
!
! The position of this server is specific: it knows ALL basic classes
! except for emission soruce terms and IO - where specialized servers exist
! Therefore, all the stuff connecting chemistry and meteorology can be 
! placed here
!
! Language: ANSI FORTRAN-90
!
! Units: SI
!
! Author: Mikhail Sofiev, mikhail.sofiev@fmi.fi
!

use diagnostic_variables
use lagrange_particles

implicit none

public write_mass_map_links_to_files
public write_lpSet_links_to_files
public interp_3d_mass_map_to_fld3d
public interp_2d_mass_map_to_fld
public add_field_to_mass_map_interp
public set_LMML_missing
public set_MML_missing

public update_mass_map_from_namelist
public update_mass_map_from_file
public fu_if_eulerian_present
public fu_if_lagrangian_present
public fu_if_hybrid_present

public merge_mass_map_2_mass_map
public merge_lagrPartSet_2_mass_map
public merge_mass_map_2_stack
public merge_vector_data_to_stack
public start_new_output_period_massmap
public start_new_output_period_lpSet

private get_scaling
private find_mass_map_4_input_id
private write_mass_map_to_files
private prepare_mass_map_writing
private convert_fld_to_mass_map_unit
private merge_inst_and_aver_2_instant
private merge_cumul_2_instant
!private merge_inst_2_cumulative   No instant fields to merge somewhwew
private merge_interval_2_cumulative  
private merge_cumul_2_cumulative
private merge_cumul_MMdata

!  !
!  ! Indices in the momentum map field - for momentums along x, y, z
!  !
!  integer, public, parameter :: indXC = 1
!  integer, public, parameter :: indYC = 2
!  integer, public, parameter :: indZC = 3

! Unit conversions for initializing mass maps:
! 
integer, private, parameter :: unity = 1001, &
                             & volume_ratio = 1002, &
                             & per_area = 1003, per_area_inv = 1004, &
                             & per_volume = 1005, per_volume_inv = 1006, &
                             & vmr_2_cnc = 1007, cnc_2_vmr = 1008, &
                             & vmr_2_mass = 1009, mass_2_vmr = 1010
!
! Variables needed for the data assimilation module.
!
integer, public, save :: DA_NUM_SOURCES = int_missing
integer, public, save :: DA_POSITIVE = int_missing, DA_NEGATIVE = int_missing, DA_ZERO = int_missing
! if da_zero == da_positive set da_zero_coef = 0 to avoid double counting, if
! da_negative /= da_positive, set da_negt_coef = -1
real, public, save :: DA_ZERO_COEF_GRAD = int_missing, &
                    & DA_NEGT_COEF_OBS = int_missing, &
                    & DA_NEGT_COEF_GRAD = int_missing

  ! type of dynamics environment the run can deal with
  !
  integer, parameter, public :: lagrangian_flag = 2201
  integer, parameter, public :: eulerian_flag = 2202
  integer, parameter, public :: hybrid_flag = 2203


  !
  ! For the needs of fast IO, intermediate mass maps can be needed for averaging purposes
  ! Here is the structure that links them
  !
  integer, public, parameter :: nMassMapLinks = 10
  integer, public, parameter :: nLagrPart2MassMapLinks = 10

  type TMassMapLink
    integer, dimension(nMassMapLinks) :: iVerticalTreatment, iAveragingType
    type(silja_interval), dimension(nMassMapLinks) :: AveragingPeriod
    type(silja_time), dimension(nLagrPart2MassMapLinks) :: IntegrationStart
    type(chemical_adaptor), dimension(nMassMapLinks) :: adaptor
    type(TMass_Map_ptr), dimension(nMassMapLinks) :: pMMIn
    type(TMass_Map), dimension(nMassMapLinks) :: pMMOut
  end type TMassMapLink
  public TMassMapLink

  !
  ! Similar to MassMap to MassMap, we have LagrangeParticles to MassMap links
  ! Note that internal MassMap is not needed: lpSet does not hold cumulative quantities
  !
  type TLagrangeToMassMapLink
    integer, dimension(nLagrPart2MassMapLinks) :: iVerticalTreatment, iAveragingType
    type(silja_interval), dimension(nLagrPart2MassMapLinks) :: AveragingPeriod
    type(silja_time), dimension(nLagrPart2MassMapLinks) :: IntegrationStart
    type(chemical_adaptor), dimension(nLagrPart2MassMapLinks) :: adaptor
    type(Tlagrange_particles_set_ptr), dimension(nLagrPart2MassMapLinks) :: pLpSetIn
    integer, dimension(nLagrPart2MassMapLinks) :: iLPSetType   ! main, short-lived or aerosol
    type(TMass_Map), dimension(nLagrPart2MassMapLinks) :: pMMOut
  end type TLagrangeToMassMapLink
  public TLagrangeToMassMapLink


CONTAINS

  subroutine set_MML_missing(MML)
    type(TMassMapLink), intent(inout) :: MML
    MML%iVerticalTreatment = int_missing
    MML%iAveragingType = int_missing
    MML%AveragingPeriod = interval_missing
    MML%IntegrationStart = time_missing
    MML%adaptor = chemical_adaptor_missing()
    MML%pMMIn = Mass_Map_ptr_missing
    MML%pMMOut = Mass_Map_missing
  end subroutine set_MML_missing
  
  
  ! **************************************************************************
  
  subroutine set_LMML_missing(LMML)
    type(TLagrangeToMassMapLink), intent(inout) :: LMML
    
    LMML%iVerticalTreatment = int_missing
    LMML%iAveragingType = int_missing
    LMML%AveragingPeriod = interval_missing
    LMML%IntegrationStart = time_missing
    LMML%adaptor = chemical_adaptor_missing()
    LMML%pLpSetIn = lagrange_particles_set_ptr_missing
    LMML%iLPSetType = int_missing
    LMML%pMMOut = Mass_Map_missing
  end subroutine set_LMML_missing
  
  
  !***************************************************************************

  subroutine write_mass_map_links_to_files(MassMapLinks, iSource, &
                                         & gridOut, vertOut, data_buf, now, output_step, &
                                         & iGrads, iNetcdf, iGrib, iGribType, ifRegular,trimfactor)
    implicit none
    
    ! Imported parameters
    type(TMassMapLink), intent(in) :: MassMapLinks
    type(silja_grid), intent(in) :: gridOut
    type(silam_vertical), intent(in) :: vertOut
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: output_step
    type(Tfield_buffer), intent(in) :: data_buf
    integer, intent(in) :: iSource, iGrads, iNetcdf, iGrib, iGribType
    logical, intent(in) :: ifRegular
    real, intent(in) :: trimfactor

    ! local variables
    integer :: iLink

    do iLink = 1, nMassMapLinks
      if(MassMapLinks%iVerticalTreatment(iLink) == int_missing)exit  ! all done
!call msg('Writing mass map link:', iLink)   
      call write_mass_map_to_files(MassMapLinks%pMMOut(iLink), &
                                 & MassMapLinks%iAveragingType(iLink), &
                                 & now - MassMapLinks%IntegrationStart(iLink), &
                                 & iSource, &
                                 & gridOut, vertOut, data_buf, now, output_step, &
                                 & iGrads, iNetcdf, iGrib, iGribType, ifRegular, trimfactor)
      if(error)return
    enddo

  end subroutine write_mass_map_links_to_files

  !***************************************************************************

  subroutine write_lpSet_links_to_files(lpSetLinks, iSource, &
                                      & gridOut, vertOut, data_buf, now, output_step, &
                                      & iGrads, iNetcdf, iGrib, iGribType, ifRegular, trimfactor)
    implicit none
    
    ! Imported parameters
    type(TLagrangeToMassMapLink), intent(in) :: lpSetLinks
    type(silja_grid), intent(in) :: gridOut
    type(silam_vertical), intent(in) :: vertOut
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: output_step
    type(Tfield_buffer), intent(in) :: data_buf
    integer, intent(in) :: iSource, iGrads, iNetcdf, iGrib, iGribType
    logical, intent(in) :: ifRegular
    real, intent(in) :: trimfactor

    ! local variables
    integer :: iLink

    do iLink = 1, nLagrPart2MassMapLinks
      if(lpSetLinks%iVerticalTreatment(iLink) == int_missing)exit  ! all done
call msg('Writing lpSet link:', iLink)
      call write_mass_map_to_files(lpSetLinks%pMMOut(iLink), &
                                 & lpSetLinks%iAveragingType(iLink), &
                                 & now - lpSetLinks%IntegrationStart(iLink), &
                                 & iSource, &
                                 & gridOut, vertOut, data_buf, now, output_step, &
                                 & iGrads, iNetcdf, iGrib, iGribType, ifRegular,trimfactor)
      if(error)return
    enddo

  end subroutine write_lpSet_links_to_files


  !***************************************************************************

  subroutine write_mass_map_to_files(pMM, iAverType, integration_period, &
                                    & iSrc, &
                                    & gridOut, vertOut, data_buf, now, output_step, &
                                    & iGrads, iNetcdf, iGrib, iGribType, ifRegular, trimfactor)
    !
    ! Wwrites down all the output mass maps in the given collection of mass map links
    ! to the output files: grads, netcdf and grib, whcih are pointed out by their indices
    ! Note that the mass maps are in their own grid, whereas the output files should be in 
    ! the given grid and vertical. Therefore, 3D interpolation is possible
    !
    implicit none
    
    ! Imported parameters
    type(TMass_map), intent(in) :: pMM
    type(silja_interval), intent(in) :: integration_period
    type(silja_grid), intent(in) :: gridOut
    type(silam_vertical), intent(in) :: vertOut
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: output_step
    type(Tfield_buffer), intent(in) :: data_buf
    integer, intent(in) :: iAverType, iSrc, iGrads, iNetcdf, iGrib, iGribType
    logical, intent(in) :: ifRegular
    real, intent(in) :: trimfactor

    ! Local variables
    integer :: iz, ix, iy,  iSpecies, i1d, iFieldKind, fs
   !    integer :: iSrc No loop over sources
    logical :: ifVertInterp, ifHorizInterp
    type(field_3d_data_ptr) :: Field3dTmp
!    type(TMass_map), pointer :: pMM
    real, dimension(:,:,:,:,:), pointer :: pMMDat
    type(THorizInterpStruct), pointer :: pHorizInterpStruct, pHorizInterpFromMeteo
    type(TVertInterpStruct), pointer :: pVertInterpStruct, pVertInterpFromMeteo

      pMMdat => pMM%Arm
      call grid_dimensions(gridOut, ix,iy)
      fs = ix*iy !! gridOut!
      ! 
      ! May be, we do not need interpolation?
      !
      if(fu_NbrOfLevels(pMM%vertTemplate) > 1)then
        ifVertInterp = .not. fu_cmp_verts_eq(pMM%vertTemplate, vertOut)
      else
        ifVertInterp = .false.
      endif
      ifHorizInterp = .not. pMM%gridTemplate == gridOut
      select case(iAverType)
        case(iMeanLastHrs)
          iFieldKind = averaged_flag
        case(iAverage)
          iFieldKind = averaged_flag
        case(iCumulative)
          iFieldKind = accumulated_flag
        case(iInstant)
          iFieldKind = forecast_flag
        case(iAsIs)
          if(fu_accumulated_quantity(pMM%quantity))then
            iFieldKind = accumulated_flag
          else
            iFieldKind = forecast_flag
          endif
        case default
          call msg("Requested averagng type", iAverType)
          call set_error('Unsupported averaging type','write_mass_map_to_files')
          return
      end select
      !
      ! Depending on the interpolation necessity, try to find shortcuts
      !
      if(ifVertInterp)then
        !
        ! Have to interpolate the vertical. Most probably, horizontal too, so just brute-force
        !
        pHorizInterpStruct => fu_horiz_interp_struct(pMM%gridTemplate, gridOut, linear, .true., &
                                                   & iOutside = setMissVal)
        pVertInterpStruct => fu_vertical_interp_struct(pMM%vertTemplate, vertOut, pMM%gridTemplate, &
                                                     & linear, output_step, 'outMassMap_to_file')
        !
        ! Note that the grid and vertical are most-probably not the same as meteo, so have to 
        ! use supplementary structures
        !
        pHorizInterpFromMeteo => fu_horiz_interp_struct(meteo_grid, pMM%gridTemplate, linear, .true., &
                                                      & iOutside = setMissVal)
        pVertInterpFromMeteo => fu_vertical_interp_struct(meteo_vertical, pMM%vertTemplate, meteo_grid, &
                                                        & linear, output_step, 'meteo_to_output')
        if(error)return
!        !
!        ! Refine coefficients if vertical intrpolation is needed. If they have
!        ! already been refined, the function will just return
!        !
! Handled in diagnostic_quantities for the whole pool     
!        call refine_interp_vert_coefs_v2(pVertInterpStruct, & ! Structure to fill
!                                       & data_buf, &       ! Meteo data buffer, pointers are ready
!                                       & now, &               ! Current time
!                                       & pHorizInterpFromMeteo, pVertInterpFromMeteo) ! grid /= meteo_grid:
!        if(error)return

        ! Prepare the 3D cube to receive reprojected data
        !
        fs = pMM%nx*pMM%ny
        do iz = 1, fu_nbrOfLevels(vertOut)
          field3dTmp%p2d(iz)%ptr => fu_work_array(fs)
          field3dTmp%p2d(iz)%ptr(1:fs) = 0.
          allocate(field3dTmp%p2d(iz)%idPtr)
        end do
        !
        ! Now, scan the whole mass map doing the job
        !
        do iSpecies = 1, pMM%nSpecies
            !
            ! Prepare field ids for all levels
            !
            do iz = 1, fu_nbrOfLevels(vertOut)
              field3dTmp%p2d(iz)%idPtr = fu_set_field_id(fmi_silam_src, &
                                             & pMM%quantity, &
                                             & now - integration_period, & !ini-time
                                             & integration_period, &    ! forecast length
                                             & gridOut, &
                                             & fu_level(vertOut,iz), &
                                             & integration_period, &  ! accum length
                                             & zero_interval, &  ! validity length
                                             & iFieldKind, &     ! 
                                             & pMM%species(iSpecies))
              IF (error) RETURN
            end do  ! iz 1:n3d
            !
            ! ... and do the 3d reprojection
            !
            call interp_3d_mass_map_to_fld3d(pMM, &  ! a set of concentr maps
                                           & iSrc, &   ! source ID of conc map to be interpolated
                                           & pMM%species(iSpecies), &
                                           & ifHorizInterp, &        ! if horizontal interp needed
                                           & .true., &               ! if vertical interp needed
                                           & pHorizInterpStruct, &   ! horizontal interp coefs
                                           & pVertInterpStruct, &    ! vertical interp coefs
                                           & do_nothing_flag, &      ! iVerticalTreatment: all done
                                           & .false., &              ! ifPerVolume: all is done
                                           & field3dTmp)       ! Place for the interpolated fields
            if(error)return
            !
            ! ... and write the stuff down to the file, with cleaning afterwards
            !
            if (trimfactor > 0) then
               DO iz = 1, fu_nbrOfLevels(vertOut)
                     call trim_precision(field3dTmp%p2d(iz)%ptr(1:fs), trimfactor, real_missing)
               END DO  ! iz
            endif
            DO iz = 1, fu_nbrOfLevels(vertOut)
              if(iGrads /= int_missing) CALL write_next_field_to_gradsfile(iGrads, &
                                                                        & field3dTmp%p2d(iz)%idPtr, &
                                                                        & field3dTmp%p2d(iz)%ptr, &
                                                                        & if_regular_output_times=ifRegular)
              if(iNetcdf /= int_missing) CALL write_next_field_to_netcdf_file(iNetcdf, &
                                                                        & field3dTmp%p2d(iz)%idPtr, &
                                                                        & field3dTmp%p2d(iz)%ptr)
              if(iGrib /= int_missing) CALL write_next_field_to_gribfile(iGribType, iGrib, &
                                                                        & field3dTmp%p2d(iz)%idPtr, &
                                                                        & field3dTmp%p2d(iz)%ptr)
              if(error)return
              field3dTmp%p2d(iz)%ptr(1:fs) = 0.
            END DO  ! iz
        end do  ! iSpecies
        !
        ! OK, this mass map is handled. Release the memory: next mass map may have different n3d
        !
        do iz = 1, fu_nbrOfLevels(vertOut)
          call free_work_array(field3dTmp%p2d(iz)%ptr)
          deallocate(field3dTmp%p2d(iz)%idPtr)
        end do

      else
        !
        ! No vertical interpolation needed. Go input layer by input layer.
        ! Note that input may be 2D or integrated columns - no interpolation in that case either
        !
        if(ifHorizInterp)then
          !
          ! Still, horizontal interpolation needed
          !
          pHorizInterpStruct => fu_horiz_interp_struct(pMM%gridTemplate, gridOut, linear, .true., &
                                                     & iOutside = setMissVal) 

          ! Prepare the 2D cube to receive reprojected data
          !
          field3dTmp%p2d(1)%ptr => fu_work_array(fs)
          field3dTmp%p2d(1)%ptr(1:fs) = 0.
          allocate(field3dTmp%p2d(1)%idPtr)
          !
          ! Now, scan the whole mass map doing the job
          !
          do iSpecies = 1, pMM%nSpecies
              !
              ! Prepare field ids for all levels.
              !
              do iz = 1, pMM%n3d
                field3dTmp%p2d(1)%idPtr = fu_set_field_id(fmi_silam_src, &
                                             & pMM%quantity, &
                                             & now - integration_period, & !ini-time
                                             & integration_period, &    ! forecast length
                                             & gridOut, &
                                             & fu_level(pMM%vertTemplate,iz), &
                                             & integration_period, &  ! accum length
                                             & zero_interval, &  ! validity length
                                             & iFieldKind, &     ! 
                                             & pMM%species(iSpecies))
                IF (error) RETURN
                !
                ! ... and do the 2d reprojection
                !
                call interp_2d_mass_map_to_fld(pMM, & ! a set of 2d maps
                                             & iSrc, & ! source ID of the map to pick
                                             & iz, &    ! level in ptrMassMap to pick
                                             & pMM%species(iSpecies), &   ! species to pick
                                             & .true., &  ! if horizontal interp needed
                                             & pHorizInterpStruct, & ! horizontal interp coefs
                                             & field3dTmp%p2d(1))
                if(error)return
                !
                ! ... and write the stuff down to the file
                !
                if (trimfactor > 0) then
                    call trim_precision(field3dTmp%p2d(1)%ptr(1:fs), trimfactor, real_missing)
                endif

                if(iGrads /= int_missing) CALL write_next_field_to_gradsfile(iGrads, &
                                                                        & field3dTmp%p2d(1)%idPtr, &
                                                                        & field3dTmp%p2d(1)%ptr, &
                                                                        & if_regular_output_times=ifRegular)
                if(iNetcdf /= int_missing) CALL write_next_field_to_netcdf_file(iNetcdf, &
                                                                        & field3dTmp%p2d(1)%idPtr, &
                                                                        & field3dTmp%p2d(1)%ptr)
                if(iGrib /= int_missing) CALL write_next_field_to_gribfile(iGribType, iGrib, &
                                                                        & field3dTmp%p2d(1)%idPtr, &
                                                                        & field3dTmp%p2d(1)%ptr)
                if(error)return
                field3dTmp%p2d(1)%ptr(1:fs) = 0.
              END DO  ! iz
          end do  ! iSpecies

          call free_work_array(field3dTmp%p2d(1)%ptr)
          deallocate(field3dTmp%p2d(1)%idPtr)

        else
          !
          ! No transformation at all. Just reorder the indices
          !
          ! Prepare the 2D cube to receive reprojected data
          !
          field3dTmp%p2d(1)%ptr => fu_work_array(fs)
          allocate(field3dTmp%p2d(1)%idPtr)
          !
          ! Now, scan the whole mass map doing the job
          !
          do iSpecies = 1, pMM%nSpecies
              !
              ! Prepare field ids for all levels
              !
              do iz = 1, pMM%n3d
                field3dTmp%p2d(1)%idPtr = fu_set_field_id(fmi_silam_src, &
                                             & pMM%quantity, &
                                             & now - integration_period, & !ini-time
                                             & integration_period, &    ! forecast length
                                             & gridOut, &
                                             & fu_level(pMM%vertTemplate,iz), &
                                             & integration_period, &  ! accum length
                                             & zero_interval, &  ! validity length
                                             & iFieldKind, &     ! 
                                             & pMM%species(iSpecies))
!call msg('')
!call report(field3dTmp%p2d(1)%idPtr)
                IF (error) RETURN
                !
                ! ... and do the mass-map to field rearrangement
                !
                do iy = 1, pMM%ny
                  do ix = 1, pMM%nx
                    i1d = ix + (iy-1) * pMM%nx
                    field3dTmp%p2d(1)%ptr(i1d) = pMMdat(iSpecies,iSrc,iz,ix,iy)
                  end do
                end do  
                !
                ! ... and write the stuff down to the file
                !
                if (trimfactor > 0) then
                    call trim_precision(field3dTmp%p2d(1)%ptr(1:fs), trimfactor, real_missing)
                endif

                if(iGrads /= int_missing) CALL write_next_field_to_gradsfile(iGrads, &
                                                                        & field3dTmp%p2d(1)%idPtr, &
                                                                        & field3dTmp%p2d(1)%ptr, &
                                                                        & if_regular_output_times=ifRegular)
                if(iNetcdf /= int_missing) CALL write_next_field_to_netcdf_file(iNetcdf, &
                                                                        & field3dTmp%p2d(1)%idPtr, &
                                                                        & field3dTmp%p2d(1)%ptr)
                if(iGrib /= int_missing) CALL write_next_field_to_gribfile(iGribType, iGrib, &
                                                                        & field3dTmp%p2d(1)%idPtr, &
                                                                        & field3dTmp%p2d(1)%ptr)
                if(error)return
              END DO  ! iz
            end do  ! iSpecies
          
          call free_work_array(field3dTmp%p2d(1)%ptr)
          deallocate(field3dTmp%p2d(1)%idPtr)

        endif  ! if horizontal interpolation
      endif  ! if vertical interpolation

  end subroutine write_mass_map_to_files


  !***************************************************************************
  !
  !
  !  INTERPOLATION of maps of cocktails 
  !
  !
  !***************************************************************************
  ! 
  ! Here we utilize the shortcut given by interpolation coefficient structures,
  ! which allows, in principle, a considerable speedup because for each output
  ! grid cell it provides weighting coefficients for all affecting input cells.
  !
  !***************************************************************************



  subroutine interp_3d_mass_map_to_fld3d(MMin, &  ! a set of concentr maps
                                         & iSourceId, &  ! source ID of the map to pick
                                         & species, &
                                         & ifHorizInterp, & ! if horizontal interp needed
                                         & ifVertInterp, &  ! if vertical interp needed
                                         & interpCoefHoriz, &  ! horizontal interp coefs
                                         & interpCoefVert, &   ! vertical interp coefs
                                         & iColumnTreatment, & ! if interate/average into 2D fld
                                         & ifPerVolume, &  ! if mass, to be converted to density
                                         & field3d)    ! Place for interpolated fields
    !
    ! Takes the given 3D map of cocktail, picks the given source, substance and size class 
    ! (and gas/aerosol phase) and interpolated the obtained map to the 3D field structure.
    ! All the interpolation rules are in the interpCoefHoriz and ifHorizInterp, as well as in
    ! interpCoefVert and ifVertInterp. Checking is limited but still exists
    !
    implicit none

    type(Tmass_map), intent(in) :: MMin    ! a set of 2d maps
    integer, intent(in) :: iSourceId          ! source ID of the map to pick
    type(silam_species), intent(in) :: species   ! species to pick
    logical, intent(in) :: ifHorizInterp, ifVertInterp    ! if horizontal and vertical interp needed
    type(THorizInterpStruct), pointer :: interpCoefHoriz  ! horizontal interp coefs
    type(TVertInterpStruct), pointer :: interpCoefVert    ! vertical interp coefs
    integer, intent(in) :: iColumnTreatment
    logical, intent(in) :: ifPerVolume  ! if mass that will be converted to density
    type(field_3d_data_ptr), intent(out) :: field3d ! Place for interpolated field

    ! Local variables
    integer :: ixTo, iyTo, iCoef, nxTo, nyTo, nxFrom, nyFrom, iLevTo, nLevsTo, iLevFrom, nLevsFrom, &
             & nSpeciesTmp, iSpecies, ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
             & ixStartTo, iyStartTo, ixEndTo, iyEndTo, &
             & iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo
    real, dimension(max_levels) :: fTmpFrom
    real :: fTmp, fWave
    type(THorizInterpCells) :: ptrCoefHoriz
    type(TVertInterpOneCell) :: ptrCoefVert
    integer, dimension(max_species) :: spcTmp
    character(len=substNmLen) :: chSubstNm ! substance to pick
    type(TAerosol_mode) :: aerMode  ! size mode to pick
    logical :: ifSumSubst, ifSumAerMode, ifSumWave, ifValid

    !
    ! A bit of checking...
    !
    if(ifHorizInterp)then
      if(.not. (MMin%gridTemplate == fu_gridFrom(interpCoefHoriz)))then
        call set_error('Cocktail map grid and interpolation grid_from do not correspond', &
                     & 'interp_3d_mass_map_to_fld3d')
        return
      endif
    endif
!    if(.not.associated(fields))then
!      call set_error('Fields pointer is not associated','interp_mass_map_to_fld3d')
!      return
!    endif

    if(iColumnTreatment == integrate_column_flag &
            &.or. iColumnTreatment == lowest_level_flag)then
      nLevsFrom = MMin%n3D
      nLevsTo = 1
    elseif(ifVertInterp)then
      nLevsFrom = fu_nbrOfLevels(fu_vertFrom(interpCoefVert))
      nLevsTo = fu_nbrOfLevels(fu_vertTo(interpCoefVert))
      call get_vertical_range_vert_interp(interpCoefVert, &
                                        & iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo)
    else
      nLevsFrom = MMin%n3D
      nLevsTo = nLevsFrom
    endif

!    if(size(fields) < nLevsTo)then
!      call set_error('Too small size of fields array','interp_mass_map_to_fld3d')
!      return
!    endif

    !
    ! Depending on what is given as real indices and what is int_missing, some summing may be needed
    ! Handled by an array of species indices to be picked from the mass map
    !
    nSpeciesTmp = 0
    ifSumSubst = len_trim(fu_substance_name(species)) == 0
    ifSumAerMode = .not. defined(fu_mode(species))
    ifSumWave = (fu_optical_wave_length(species) .eps. real_missing)

    if(ifSumSubst .and. ifSumAerMode .and. ifSumWave)then
      !
      ! Sum over the whole mass map
      !
      nSpeciesTmp = MMin%nSpecies
      do iSpecies = 1, nSpeciesTmp
        spcTmp(iSpecies) = iSpecies
      end do

    elseif(ifSumAerMode .and. ifSumWave)then
      !
      ! Sum all for the given substance: modes and wave lengths
      !
      chSubstNm = fu_substance_name(species)
      do iSpecies = 1, MMin%nSpecies
        if(fu_substance_name(MMin%species(iSpecies)) == chSubstNm)then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do

    elseif(ifSumAerMode)then
      !
      ! Sum over all modes for the given substance and wave lenght
      !
      chSubstNm = fu_substance_name(species)
      fWave = fu_optical_wave_length(species)
      do iSpecies = 1, MMin%nSpecies
        if(fu_substance_name(MMin%species(iSpecies)) == chSubstNm .and. &
         & (fWave .eps. fu_optical_wave_length(MMin%species(iSpecies))))then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do

    elseif(ifSumWave)then
      !
      ! For the given subsatnce and mode sum all wave lengths
      !
      chSubstNm = fu_substance_name(species)
      aerMode = fu_mode(species)
      do iSpecies = 1, MMin%nSpecies
        if(fu_substance_name(MMin%species(iSpecies)) == chSubstNm .and. &
         & (aerMode == fu_mode(MMin%species(iSpecies))))then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do

    else
      !
      ! Exact species is requested
      !
      do iSpecies = 1, MMin%nSpecies
        if(MMin%species(iSpecies) == species)then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do
    endif    !  what to sum-up

    !----------------------------------------------------------------------------
    !
    ! The main cycle over the vertical_To
    !
    if(ifHorizInterp)then
      call grid_dimensions(fu_gridTo(interpCoefHoriz), nxTo, nyTo)
      call grid_dimensions(fu_gridFrom(interpCoefHoriz), nxFrom, nyFrom)
    else
      nxFrom = MMin%nx
      nyFrom = MMin%ny
      nxTo = nxFrom
      nyTo = nyFrom
    endif
    fTmpFrom(1:nLevsFrom) = 0.

    if(ifHorizInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      call get_area_limits(interpCoefHoriz, ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
                                          & ixStartTo, iyStartTo, ixEndTo, iyEndTo)
      call get_coefs(interpCoefHoriz, ptrCoefHoriz)
      if(error)return
      !
      ! Means of vertical treatment: collaps to single layer, interpolate, do-nothing
      !
      if(iColumnTreatment == integrate_column_flag)then
        !
        ! The whole vertical is collapsed into a single level, interated or averaged over
        !
        if(.not. ifPerVolume)then
          do iLevFrom = 1, nLevsFrom
            fTmpFrom(iLevFrom) = fu_layer_thickness_m(fu_level(MMin%vertTemplate,iLevFrom))
          end do
          fTmpFrom(nLevsFrom+1) = sum(fTmpFrom(1:nLevsFrom))
        endif

        do iyTo = iyStartTo, iyEndTo
          do ixTo = ixStartTo, ixEndTo
            if(ptrCoefHoriz%indX(1,ixTo,iyTo) == 0)cycle
            !
            ! Prepare the fTmpFrom - the integral over vertical column that corresponds to the 
            ! current cell of gridTo: do the horizontal interpolation
            !
            fTmp = 0.0

            do iLevFrom = iLevStartFrom, iLevEndFrom
              do iSpecies = 1, nSpeciesTmp           ! Summing up or picking the single species
                if(ifPerVolume)then
                  fTmp = fTmp + ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                                               & MMin%arM(spcTmp(iSpecies), &
                                                              & iSourceId, &
                                                              & iLevFrom, &
                                                              & ptrCoefHoriz%indX(iCoef,ixTo,iyTo), &
                                                              & ptrCoefHoriz%indY(iCoef,ixTo,iyTo))
                else
                  fTmp = fTmp + fTmpFrom(iLevFrom) * ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                                               & MMin%arM(spcTmp(iSpecies), &
                                                              & iSourceId, &
                                                              & iLevFrom, &
                                                              & ptrCoefHoriz%indX(iCoef,ixTo,iyTo), &
                                                              & ptrCoefHoriz%indY(iCoef,ixTo,iyTo))
                endif
              end do   ! species
            end do ! Summing-up the vertFrom column but for gridTo
            !
            ! Put the value to right element
            !
!            if(iColumnTreatment == integrate_column_flag)then
              field3d%p2d(1)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp   ! The integral...
!            else
!              fields(1)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp / fTmpFrom(nLevsFrom+1) ! ...or mean over vertical
!            endif

          end do   ! ixTo
        end do   ! iyTo

      elseif(ifVertInterp)then
        !
        ! Full-blown vertical interpolation.
        ! Cycle over the horizontal grid _To
        !
        do iyTo = iyStartTo, iyEndTo
          do ixTo = ixStartTo, ixEndTo
            if(ptrCoefHoriz%indX(1,ixTo,iyTo) == 0)cycle
            !
            ! Prepare the fTmpFrom - the vertical column in the From coordinates
            ! that corresponds to the current cell of gridTo: do the horizontal interpolation
            ! Note: there can be plenty of vertical elements to skip in the From vertical.
            ! We skip them by limiting the variation of the iLevFrom index. Note: interpolation
            ! structure is made so that indices are made in ascending order, so the lowest index
            ! in the vertFrom that will be needed is pointed by the first index at the lowest 
            ! level of the interpolation structure, etc.
            !
            fTmpFrom(1:nLevsFrom) = 0.0

            do iLevFrom = iLevStartFrom, iLevEndFrom
              do iCoef = 1, fu_nCoefs(interpCoefHoriz)
                if(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) == 0)exit ! end of useful indices
                do iSpecies = 1, nSpeciesTmp           ! Summing up or picking the single species
                  fTmpFrom(iLevFrom) = fTmpFrom(iLevFrom) + ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                                               & MMin%arM(spcTmp(iSpecies), &
                                                              & iSourceId, &
                                                              & iLevFrom, &
                                                              & ptrCoefHoriz%indX(iCoef,ixTo,iyTo), &
                                                              & ptrCoefHoriz%indY(iCoef,ixTo,iyTo))
                end do   ! species
              end do   ! coefs
            end do ! Filling-up the vertFrom column but for gridTo

            !
            ! The horizontal interpolation made the vertical column in vertFrom
            ! Now it can be interpolated to vertTo
            !
            do iLevTo = iLevStartTo, iLevEndTo
              !
              ! Fill-in the fTmpFrom(:) column in the vertical_From coordinates
              !
              call get_coefs(interpCoefVert,ixTo,iyTo,iLevTo,ptrCoefVert)
              fTmp = 0.
              do iCoef = 1, fu_nCoefs(interpCoefVert)
                if(ptrCoefVert%indLev(iCoef) == 0)exit ! end of useful indices
                fTmp = fTmp + ptrCoefVert%weight(iCoef) * fTmpFrom(ptrCoefVert%indLev(iCoef))
              end do
              field3d%p2d(iLevTo)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp  ! Put the value to right element
            end do  ! cycle over the vertical To
          end do   ! ixTo
        end do   ! iyTo

      else
        !
        ! No vertical interpolation or integration/averaging is needed
        !
        call get_coefs(interpCoefHoriz, ptrCoefHoriz)
        if(error)return

        do iyTo = 1, nyTo
          do ixTo = 1, nxTo
            if(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) == 0)cycle
            !
            ! Vertical interpolation is void - the horizontal step is the only one
            !
            do iLevFrom = iLevStartFrom, iLevEndFrom
              fTmp = 0.
              do iCoef = 1, fu_nCoefs(interpCoefHoriz)
                if(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) == 0)exit
                do iSpecies = 1, nSpeciesTmp
                  fTmp = fTmp + ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                              & MMin%arM(spcTmp(iSpecies), &
                                             & iSourceId, &
                                             & iLevFrom, &
                                             & ptrCoefHoriz%indX(iCoef,ixTo,iyTo), &
                                             & ptrCoefHoriz%indY(iCoef,ixTo,iyTo))
                end do
              end do
              if(fTmp > 0.001)then
                call msg('ix,iy:',ixTo,real(iyTo))
              endif
              field3d%p2d(iLevFrom)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp
            end do ! Filling-up the vertFrom column but for gridTo

          end do   ! ixTo
        end do   ! iyTo
      endif ! if vertical interpolation

      !
      ! Apply missing-value mask 
      if (associated( ptrCoefHoriz%ifValid)) then ! mask out missings
        do iyTo = 1, nyTo
          do ixTo = 1, nxTo
               if (.not. ptrCoefHoriz%ifValid(ixTo,iyTo)) then
                   do iLevFrom = iLevStartFrom, iLevEndFrom
                     field3d%p2d(iLevFrom)%ptr(ixTo + nxTo*(iyTo-1)) = real_missing
                   enddo
               endif
          end do   ! ixTo
        end do   ! iyTo
        
      endif
    else  !-----------------------------------------------------------------------------
      !
      ! Horizontal interpolation is void, just copy the values, possibly,
      ! with vertical interpolation or column integration
      !
      if(iColumnTreatment == integrate_column_flag)then
        !
        ! The whole vertical is collapsed into a single level, interated or averaged over
        !
        if(.not. ifPerVolume)then
          do iLevFrom = 1, nLevsFrom
            fTmpFrom(iLevFrom) = fu_layer_thickness_m(fu_level(MMin%vertTemplate,iLevFrom))
          end do
          fTmpFrom(nLevsFrom+1) = sum(fTmpFrom(1:nLevsFrom))
        endif

        do iyTo = 1, nyTo
          do ixTo = 1, nxTo
            !
            ! Prepare the fTmp - the integral over vertical column. Note that fTmpFrom now has 
            ! the layer thicknesses
            !
            fTmp = 0.0

            do iLevFrom = 1, nLevsFrom
              do iSpecies = 1, nSpeciesTmp           ! Summing up or picking the single species
                if(ifPerVolume)then
                  fTmp = fTmp + MMin%arM(spcTmp(iSpecies), iSourceId, iLevFrom, ixTo,iyTo)
                else
                  fTmp = fTmp + fTmpFrom(iLevFrom) * &
                              & MMin%arM(spcTmp(iSpecies), iSourceId, iLevFrom, ixTo,iyTo)
                endif
!if(fTmp > 1)then
!  print *
!endif
              end do   ! species
            end do ! Summing-up the vertFrom column but for gridTo
            !
            ! Put the value to right element. Note that we have summed-up the mass. 
            !
!            if(iColumnTreatment == integrate_column_flag)then
              field3d%p2d(1)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp   ! The integral...
!            else
!              fields(1)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp / fTmpFrom(nLevsFrom+1)  ! ...or mean over vertical
!            endif

          end do   ! ixTo
        end do   ! iyTo

      elseif(ifVertInterp)then
        !
        ! Full-blown vertical interpolation (still no horizontal one)
        ! Cycle over the horizontal grid _To
        !
        do iyTo = 1, nyTo
          do ixTo = 1, nxTo
            !
            ! The horizontal interpolation is void, so ixFrom==ixTo, iyFrom==iyTo
            ! It should only be interpolated to vertTo
            !
            do iLevTo = iLevStartTo, iLevEndTo
              !
              ! Fill-in the fTmpFrom(:) column in the vertical_From coordinates
              !
              call get_coefs(interpCoefVert,ixTo,iyTo,iLevTo,ptrCoefVert)
              if(ptrCoefVert%indLev(1) == 0)cycle
              fTmp = 0.

              do iCoef = 1, fu_nCoefs(interpCoefVert)
                if(ptrCoefVert%indLev(iCoef) == 0)exit
                do iSpecies = 1, nSpeciesTmp
                   fTmp = fTmp + ptrCoefVert%weight(iCoef) * MMin%arM(spcTmp(iSpecies), &
                                                                           & iSourceId, &
                                                                           & ptrCoefVert%indLev(iCoef), &
                                                                           & ixTo, &
                                                                           & iyTo)
                end do
              end do
              field3d%p2d(iLevTo)%ptr(ixTo + nxTo*(iyTo-1)) = fTmp  ! Put the value to right element
            end do  ! cycle over the vertical To

          end do   ! ixTo
        end do   ! iyTo

      else
        !
        ! No vertical, no horizontal interpolation is needed - just copy the requested field
        !
        call start_count('cocktail_map_to_fld_no_interp')
        do iyTo = 1, nyTo
          do ixTo = 1, nxTo
            do iLevTo = 1, nLevsTo
              field3d%p2d(iLevTo)%ptr(ixTo + nxTo*(iyTo-1)) = 0.
              do iSpecies = 1, nSpeciesTmp
                field3d%p2d(iLevTo)%ptr(ixTo+nxTo*(iyTo-1)) = &
                                            & field3d%p2d(iLevTo)%ptr(ixTo+nxTo*(iyTo-1)) + &
                                            & MMin%arM(spcTmp(iSpecies), &
                                                                      & iSourceId,iLevTo,ixTo,iyTo)
              end do
            end do
          end do   ! ixTo
        end do   ! iyTo
        call stop_count('cocktail_map_to_fld_no_interp')
                
      endif ! if vertical interpolation 

    endif ! if horizontal interpolation


  end subroutine interp_3d_mass_map_to_fld3d



  !**********************************************************************************
  
  subroutine interp_2d_mass_map_to_fld(MassMapIn, & ! a set of 2d maps
                                     & iSourceId, & ! source ID of the map to pick
                                     & iLevel, &    ! level in MassMapIn to pick
                                     & species, &   ! species to pick
                                     & ifInterp, &  ! if horizontal interp needed
                                     & interpCoefHoriz, & ! horizontal interp coefs
                                     & field)       ! Place for interpolated field
    !
    ! Takes the given 2D map of cocktail, picks the given source, substance and size class 
    ! (and gas/aerosol phase) and interpolated the obtained map to the field.
    ! All the interpolation rules are in the interpCoefHoriz and ifInterp. Checking is 
    ! limited but still exists
    !
    implicit none

    type(Tmass_map), intent(in) :: MassMapIn    ! a set of 2d maps
    integer, intent(in) :: iSourceId, iLevel  ! source ID and level of the map to pick
    type(silam_species), intent(in) :: species
    logical, intent(in) :: ifInterp  ! if horizontal interp needed
    type(THorizInterpStruct), intent(in) :: interpCoefHoriz ! horizontal interp coefs
    type(field_data_ptr), intent(inout) :: field         ! Place for interpolated field

    ! Local variables
    integer :: ixTo, iyTo, iCoef, nxTo, nyTo, nxFrom, nyFrom, nSpeciesTmp, iSpecies, &
             & ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, ixStartTo, iyStartTo, ixEndTo, iyEndTo
    real :: fTmp, fWave
    type(THorizInterpCells) :: ptrCoefHoriz
    integer, dimension(max_species) :: spcTmp
    character(len=substNmLen) :: chSubstNm ! substance index to be picked
    type(TAerosol_mode) :: aerMode      ! index of size mode to pick
    logical :: ifSumSubst, ifSumAerMode, ifSumWave, ifValid

!    mapper => MassMapIn%mapper
    !
    ! A bit of checking...
    !
    if(ifInterp)then
      if(.not. (MassMapIn%gridTemplate == fu_gridFrom(interpCoefHoriz)))then
        call set_error('Cocktail map grid and interpolation grid_from do not correspond', &
                     & 'interp_2d_mass_map_to_fld')
        return
      endif
      if(.not. (fu_grid(field%idPtr) == fu_gridTo(interpCoefHoriz)))then
        call set_error('field grid and interpolation grid_to do not correspond', &
                     & 'interp_2d_mass_map_to_fld')
        return
      endif
    endif
    !
    ! Depending on what is given as real indices and what is int_missing, some summing may be needed
    ! Handled by an array of species indices to be picked from the mass map
    !
    nSpeciesTmp = 0
    ifSumSubst = len_trim(fu_substance_name(species)) == 0
    ifSumAerMode = .not. defined(fu_mode(species))
    ifSumWave = (fu_optical_wave_length(species) .eps. real_missing)

    if(ifSumSubst .and. ifSumAerMode .and. ifSumWave)then
      !
      ! Sum over the whole mass map
      !
      nSpeciesTmp = MassMapIn%nSpecies
      do iSpecies = 1, nSpeciesTmp
        spcTmp(iSpecies) = iSpecies
      end do

    elseif(ifSumAerMode .and. ifSumWave)then
      !
      ! Sum all for the given substance: modes and wave lengths
      !
      chSubstNm = fu_substance_name(species)
      do iSpecies = 1, MassMapIn%nSpecies
        if(fu_substance_name(MassMapIn%species(iSpecies)) == chSubstNm)then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do

    elseif(ifSumAerMode)then
      !
      ! Sum over all modes for the given substance and wave lenght
      !
      chSubstNm = fu_substance_name(species)
      fWave = fu_optical_wave_length(species)
      do iSpecies = 1, MassMapIn%nSpecies
        if(fu_substance_name(MassMapIn%species(iSpecies)) == chSubstNm .and. &
         & (fWave .eps. fu_optical_wave_length(MassMapIn%species(iSpecies))))then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do

    elseif(ifSumWave)then
      !
      ! For the given subsatnce and mode sum all wave lengths
      !
      chSubstNm = fu_substance_name(species)
      aerMode = fu_mode(species)
      do iSpecies = 1, MassMapIn%nSpecies
        if(fu_substance_name(MassMapIn%species(iSpecies)) == chSubstNm .and. &
         & (aerMode == fu_mode(MassMapIn%species(iSpecies))))then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do

    else
      !
      ! Exact species is requested
      !
      do iSpecies = 1, MassMapIn%nSpecies
        if(MassMapIn%species(iSpecies) == species)then
          nSpeciesTmp = nSpeciesTmp + 1
          spcTmp(nSpeciesTmp) = iSpecies
        end if
      end do
    endif

    !----------------------------------------------------------------------------
    !
    ! The main cycle over the output grid. Careful! In case of no interpoaltion needed
    ! the interpolation structure is nullified, so do NOT use it.
    !
    if(ifInterp)then
      !
      ! Interpolation needed, structure is defined
      !
      call get_area_limits(interpCoefHoriz, ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
                                          & ixStartTo, iyStartTo, ixEndTo, iyEndTo)
      call grid_dimensions(fu_gridTo(interpCoefHoriz), nxTo, nyTo)

      call get_coefs(interpCoefHoriz, ptrCoefHoriz)

      do iyTo = iyStartTo, iyEndTo
        do ixTo = ixStartTo, ixEndTo
          ifValid=.true.
          if (associated( ptrCoefHoriz%ifValid)) ifValid = ptrCoefHoriz%ifValid(ixTo,iyTo)
          if (ifValid)then
             fTmp = 0.  ! Get a single value in the output grid 
             do iCoef = 1, fu_nCoefs(interpCoefHoriz)
               if(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) == 0)exit
               do iSpecies = 1, nSpeciesTmp
                 fTmp = fTmp + ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                                     & MassMapIn%arM(spcTmp(iSpecies), &
                                                    & iSourceId, &
                                                    & iLevel, &
                                                    & ptrCoefHoriz%indX(iCoef,ixTo,iyTo), &
                                                    & ptrCoefHoriz%indY(iCoef,ixTo,iyTo))
               end do
             end do
             field%ptr(ixTo + nxTo*(iyTo-1)) = fTmp  ! Put the value to right element of the field
          else
              field%ptr(ixTo + nxTo*(iyTo-1)) = real_missing
          endif
        end do
      end do
      
    else
      !
      ! No interpolation needed, just sum-up the arrays. Dimensions can be taken from the 
      ! cocktail map array because grids are the same
      !
      do iyTo = 1, MassMapIn%ny
        do ixTo = 1, MassMapIn%nx
          field%ptr(ixTo + MassMapIn%nx*(iyTo-1)) = 0.
          do iSpecies = 1, nSpeciesTmp
            field%ptr(ixTo+MassMapIn%nx*(iyTo-1)) = field%ptr(ixTo+MassMapIn%nx*(iyTo-1)) + &
                                     & MassMapIn%arM(spcTmp(iSpecies), iSourceId, iLevel,ixTo,iyTo)
          end do
        end do
      end do

    endif  ! if the interpolation is needed

  end subroutine interp_2d_mass_map_to_fld


  !*************************************************************************************************
  
  subroutine add_field_to_mass_map_interp(field_3d, indSpeciesOut, fMassFractionDescr, &
                                        & nSpeciesDescr, fMassAdded, &
                                        & fRateScaling, &
                                        & ifHorizInterp, ifVertInterp, interpCoefHoriz, interpCoefVert,&
                                        & iSourceId, &
                                        & pMassMap, pMomentumMapX, pMomentumMapY, pMomentumMapZ)
    !
    ! Adds a field to mass map with possible distribution of the species to the map
    ! Also copies momentum if corresponding mass map is available
    !
    implicit none
    
    ! Imported parameters
    type(silja_3d_field), pointer :: field_3d             ! Input data. NO CHECKING
    integer, dimension(:), intent(in) :: indSpeciesOut    ! indices of the species in the mass map (1:nSpeciesDescr)
    real, dimension(:), intent(in) :: fMassFractionDescr     ! what goes into each species
    real(kind=8),  dimension(:), intent(inout) :: fMassAdded  ! Actually transferred mass (massmap species indices)
    integer, intent(in) :: nSpeciesDescr                     ! Dimension of species and fMassFraction 
    real, intent(in) :: fRateScaling                      ! Scaling factor for the added mass 
    logical, intent(in) :: ifHorizInterp, ifVertInterp    ! if horizontal and vertical interp needed
    type(THorizInterpStruct), intent(in) :: interpCoefHoriz  ! horizontal interp coefs
    type(TVertInterpStruct), intent(in) :: interpCoefVert    ! vertical interp coefs
    integer, intent(in) :: iSourceId                      ! source ID of the map to pick
    type(Tmass_map), intent(inout), target :: pMassMap    ! place for output mass
    type(Tmass_map),  intent(inout), target, optional :: pMomentumMapX, pMomentumMapY, pMomentumMapZ  ! output momentum

    ! Local parameters
    real, dimension(3*max_levels) :: fTmpFrom
    integer :: nxFrom, nyFrom, nxTo, nyTo, nLevsFrom, nLevsTo, ixTo, iyTo, ixFrom, iyFrom, &
             & iLevFrom, iLevTo, i1dFrom, iSpecies, iCoef, ixStartFrom, iyStartFrom, &
             & ixEndFrom, iyEndFrom, ixStartTo, iyStartTo, ixEndTo, iyEndTo, &
             & iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo
    type(THorizInterpCells) :: ptrCoefHoriz
    type(TVertInterpOneCell) :: ptrCoefVert
    type(silja_rp_1d), dimension(max_levels) :: pDataFrom
    logical :: ifMoments
    real :: fTmp, fTmp2, fTmp3, fTmp4
    real, pointer :: fPtr

!integer :: iSpLoc, iSrcLoc, ixLoc,iyLoc,izLoc
!real, dimension(:), pointer :: fArHist
!integer, dimension(:), pointer, save :: iArHist

    !
    ! A bit of preparation
    !
    if(present(pMomentumMapX))then
      if(.not. (present(pMomentumMapY) .and. present(pMomentumMapZ)))then
        call set_error('MomentumX present, other(s) not','add_field_to_mass_map_interp')
        return
      endif
      ifMoments = .true.
    else
      ifMoments = .false.
    endif

    if(ifVertInterp)then
      nLevsFrom = fu_nbrOfLevels(fu_vertFrom(interpCoefVert))
      nLevsTo = fu_nbrOfLevels(fu_vertTo(interpCoefVert))
      call get_vertical_range_vert_interp(interpCoefVert, &
                                        & iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo)
    else
      nLevsFrom = pMassMap%n3D
      nLevsTo = nLevsFrom
    endif

    if(ifHorizInterp)then
      if (associated(interpCoefHoriz%ifValid)) then
            call set_error("Handling missing values not implemented",&
                        &'add_field_to_mass_map_interp')
            return
      endif
        
      call grid_dimensions(fu_gridTo(interpCoefHoriz), nxTo, nyTo)
      call grid_dimensions(fu_gridFrom(interpCoefHoriz), nxFrom, nyFrom)
      call get_area_limits(interpCoefHoriz, ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
                                          & ixStartTo, iyStartTo, ixEndTo, iyEndTo)
    else
      nxFrom = pMassMap%nx
      nyFrom = pMassMap%ny
      nxTo = nxFrom
      nyTo = nyFrom
    endif

    do iLevFrom = 1, nLevsFrom
      pDataFrom(iLevFrom)%pp => fu_grid_data_from_3d(field_3d,iLevFrom)
    enddo

    !----------------------------------------------------------------------------
    !
    ! The main cycle
    !
    if(ifHorizInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      call get_coefs(interpCoefHoriz, ptrCoefHoriz)
      if(error)return

      if(ifVertInterp)then
        !
        ! Full-blown vertical interpolation.
        ! Cycle over the horizontal grid _To
        !
        do iyTo = iyStartTo, iyEndTo
          do ixTo = ixStartTo, ixEndTo
            if(ptrCoefHoriz%indX(1,ixTo,iyTo) == 0)cycle ! if void cell
            !
            ! Prepare the fTmpFrom - the vertical column in the From coordinates
            ! that corresponds to the current cell of gridTo: do the horizontal interpolation
            ! Note: there can be plenty of vertical elements to skip in the From vertical.
            ! We skip them by limiting the variation of the iLevFrom index. Note: interpolation
            ! structure is made so that indices are in ascending order, so the lowest index
            ! in the vertFrom that will be needed is pointed by the first index at the lowest 
            ! level of the interpolation structure, etc.
            !
            do iSpecies = 1, nSpeciesDescr           ! Distribute the values
              !
              ! Fill-in vertical column in the (ixTo,iyTo) of gridTo
              !
              fTmpFrom(1:nLevsFrom*3) = 0.
              do iLevFrom = iLevStartFrom, iLevEndFrom
                do iCoef = 1, fu_nCoefs(interpCoefHoriz)
                  i1dFrom = ptrCoefHoriz%indX(iCoef,ixTo,iyTo) + &
                                                & (ptrCoefHoriz%indY(iCoef,ixTo,iyTo)-1) * nxFrom
                  if(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) == 0)exit
                  fTmpFrom(iLevFrom) = fTmpFrom(iLevFrom) + ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                                                          & pDataFrom(iLevFrom)%pp(i1dFrom)
!if(ixTo == 41 .and. iyTo == 46)then
!  call msg('Summing iCoef, ixFrom, iyFrom, mass, Px, Py: (' + fu_str(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) ) + ',' + &
!                                                          & fu_str(ptrCoefHoriz%indY(iCoef,ixTo,iyTo)) + ')', &
!    & (/ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * pDataFrom(iLevFrom)%pp(i1dFrom), &
!      & ptrCoefHoriz%weight_X(iCoef,ixTo,iyTo) * pDataFrom(iLevFrom)%pp(i1dFrom), &
!      & ptrCoefHoriz%weight_Y(iCoef,ixTo,iyTo) * pDataFrom(iLevFrom)%pp(i1dFrom)/))
!endif
                  !
                  ! Moments: x and y are treated with their weights, z is intact
                  if(ifMoments)then
                    fTmpFrom(iLevFrom+nLevsFrom) = fTmpFrom(iLevFrom+nLevsFrom) + &     
                            & ptrCoefHoriz%weight_X(iCoef,ixTo,iyTo) * pDataFrom(iLevFrom)%pp(i1dFrom)
                    fTmpFrom(iLevFrom+2*nLevsFrom) = fTmpFrom(iLevFrom+2*nLevsFrom) + &     
                            & ptrCoefHoriz%weight_Y(iCoef,ixTo,iyTo) * pDataFrom(iLevFrom)%pp(i1dFrom)
                  endif
                end do   ! coefs
              end do ! Filling-up the vertFrom column but for gridTo
              !
              ! The horizontal interpolation created the vertical column in vertFrom
              ! Now it can be interpolated to vertTo
              ! Horizontal moments are interpolated same way as masses. Vertical one takes
              ! separate stand.
              !
              do iLevTo = iLevStartTo, iLevEndTo
                !
                ! Fill-in the fTmpFrom(:) column in the vertical_From coordinates
                !
                call get_coefs(interpCoefVert,ixTo,iyTo,iLevTo,ptrCoefVert)
                if(ptrCoefVert%indLev(1) == 0)cycle
                fTmp = 0.
                fTmp2 = 0.
                fTmp3 = 0.
                fTmp4 = 0.
                do iCoef = 1, fu_nCoefs(interpCoefVert)
                  if(ptrCoefVert%indLev(iCoef) == 0)exit
                  fTmp = fTmp + ptrCoefVert%weight(iCoef) * fTmpFrom(ptrCoefVert%indLev(iCoef))
                  !
                  ! Moments: horizontal are summed-up just as mass, vertical is treated with weight_Z
                  if(ifMoments)then
                    fTmp2 = fTmp2 + ptrCoefVert%weight(iCoef) * &
                                                  & fTmpFrom(ptrCoefVert%indLev(iCoef)+nLevsFrom)
                    fTmp3 = fTmp3 + ptrCoefVert%weight(iCoef) * &
                                                  & fTmpFrom(ptrCoefVert%indLev(iCoef)+2*nLevsFrom)
                    fTmp4 = fTmp4 + ptrCoefVert%weight_Z(iCoef) * fTmpFrom(ptrCoefVert%indLev(iCoef))
                  endif
                end do

                fTmp =  fTmp * fMassFractionDescr(iSpecies) * fRateScaling
                fMassAdded(indSpeciesOut(iSpecies)) = fMassAdded(indSpeciesOut(iSpecies)) + fTmp 
                fPtr =>  pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo)
                fPtr = fPtr + fTmp

                pMassMap%ifGridValid(ilevTo,iSourceId) = .true.
                pMassMap%ifColumnValid(iSOurceId,ixTo,iyTo) = .true.
                if(ifMoments)then
                   fTmp = fMassFractionDescr(iSpecies) * fRateScaling

                  fPtr => pMomentumMapX%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo) 
                  fPtr = fPtr + fTmp2  * fTmp

                  fPtr => pMomentumMapY%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo)
                  fPtr = fPtr + fTmp3 *ftmp

                  fPtr => pMomentumMapZ%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo)
                  fPtr = fPtr + fTmp4 *fTmp

                  pMomentumMapX%ifGridValid(ilevTo,iSourceId) = .true.
                  pMomentumMapX%ifColumnValid(iSourceId,ixTo,iyTo) = .true.
                  pMomentumMapY%ifGridValid(ilevTo,iSourceId) = .true.
                  pMomentumMapY%ifColumnValid(iSourceId,ixTo,iyTo) = .true.
                  pMomentumMapZ%ifGridValid(ilevTo,iSourceId) = .true.
                  pMomentumMapZ%ifColumnValid(iSourceId,ixTo,iyTo) = .true.

call check_mass_moment_point(pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo), &
                      & pMomentumMapX%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo) , &
                      & 'X-monent wrong after full reprojection','add_field_to_mass_map_interp')
call check_mass_moment_point(pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo), &
                      & pMomentumMapY%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo) , &
                      & 'Y-monent wrong after full reprojection','add_field_to_mass_map_interp')
call check_mass_moment_point(pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo), &
                      & pMomentumMapZ%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo) , &
                      & 'Z-monent wrong afterifull reprojection','add_field_to_mass_map_interp')
!if(pMassMap%arM(iSpecies, iSourceId, iLevTo, ixTo,iyTo) > 0)then
!ixLoc=ixTo; iyLoc=iyTo; izLoc=iLevTo
!print *,'*'
!endif
                endif  ! ifMoments
              end do  ! cycle over the vertical To
            end do   ! species
          end do   ! ixTo
        end do   ! iyTo

!open(50,file='histogram.txt',access='append')
!write(50,fmt='(110(i5,1x))')(iArHist(ixTo),ixTo=1,105)
!close(50)
      else
        !
        ! No vertical interpolation is needed
        !
        do iyTo = iyStartTo, iyEndTo
          do ixTo = ixStartTo, ixEndTo
            if(ptrCoefHoriz%indX(1,ixTo,iyTo) == 0)cycle
            !
            ! Vertical interpolation is void - the horizontal step is the only one
            !
            do iSpecies = 1, nSpeciesDescr
              do iLevFrom = 1, nLevsFrom
                fTmp = 0.
                fTmp2 = 0.
                fTmp3 = 0.
                do iCoef = 1, fu_nCoefs(interpCoefHoriz)
                  if(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) == 0)exit
                  fTmp = fTmp + ptrCoefHoriz%weight(iCoef,ixTo,iyTo) * &
                              & pDataFrom(iLevFrom)%pp(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) + &
                                                   & (ptrCoefHoriz%indY(iCoef,ixTo,iyTo)-1) * &
                                                   & nxFrom)
                  if(ifMoments)then
                    fTmp2 = fTmp2 + ptrCoefHoriz%weight_X(iCoef,ixTo,iyTo) * &
                                & pDataFrom(iLevFrom)%pp(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) + &
                                                     & (ptrCoefHoriz%indY(iCoef,ixTo,iyTo)-1) * &
                                                     & nxFrom)
                    fTmp3 = fTmp3 + ptrCoefHoriz%weight_Y(iCoef,ixTo,iyTo) * &
                                & pDataFrom(iLevFrom)%pp(ptrCoefHoriz%indX(iCoef,ixTo,iyTo) + &
                                                     & (ptrCoefHoriz%indY(iCoef,ixTo,iyTo)-1) * &
                                                     & nxFrom)
                  endif
                end do

                fTmp =  fTmp * fMassFractionDescr(iSpecies) * fRateScaling
                fMassAdded(indSpeciesOut(iSpecies)) = fMassAdded(indSpeciesOut(iSpecies)) + fTmp 
                fPtr =>  pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo)
                fPtr = fPtr + fTmp

                pMassMap%ifGridValid(ilevFrom,iSourceId) = .true.
                pMassMap%ifColumnValid(iSOurceId,ixTo,iyTo) = .true.
                if(ifMoments)then
                  pMomentumMapX%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo) = &
                     & pMomentumMapX%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo) + &
                     & fTmp2 * fMassFractionDescr(iSpecies) * fRateScaling
                  pMomentumMapY%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo) = &
                     & pMomentumMapY%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo) + &
                     & fTmp3 * fMassFractionDescr(iSpecies) * fRateScaling
                  pMomentumMapX%ifGridValid(ilevFrom,iSourceId) = .true.
                  pMomentumMapX%ifColumnValid(iSourceId,ixTo,iyTo) = .true.
                  pMomentumMapY%ifGridValid(ilevFrom,iSourceId) = .true.
                  pMomentumMapY%ifColumnValid(iSOurceId,ixTo,iyTo) = .true.
call check_mass_moment_point(pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo), &
                      & pMomentumMapX%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo) , &
                      & 'X-monent wrong after horizontal-only reprojection','add_field_to_mass_map_interp')
call check_mass_moment_point(pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo), &
                      & pMomentumMapY%arM(indSpeciesOut(iSpecies), iSourceId, iLevFrom, ixTo,iyTo) , &
                      & 'Y-monent wrong after horizontal-only reprojection','add_field_to_mass_map_interp')
                endif  ! ifMoment
              end do ! LevFrom
            end do  ! species
          end do   ! ixTo
        end do   ! iyTo
      endif ! if vertical interpolation 

    else
      !
      ! Horizontal interpolation is void, just copy the values, possibly,
      ! with vertical interpolation
      !
      if(ifVertInterp)then
        !
        ! Full-blown vertical interpolation.
        ! Cycle over the horizontal grid _To
        !
        do iyTo = 1, nyTo
          do ixTo = 1, nxTo
            !
            ! The horizontal interpolation is void, so ixFrom==ixTo, iyFrom==iyTo
            ! It should only be interpolated to vertTo
            !
!if(pDataFrom(1)%pp(ixTo + nxTo * (iyTo-1)) > 1e-10)then
!  call msg('pDataFrom(1)%pp(ixTo + nxTo * (iyTo-1))',pDataFrom(1)%pp(ixTo + nxTo * (iyTo-1)))
!endif
            i1dFrom = ixTo + nxTo * (iyTo-1)   ! same as To

            !!! Anything in this column?
            do iLevFrom = 1, nLevsFrom
              if (pDataFrom(iLevFrom)%pp(i1dFrom) /= 0.)  exit
            enddo
            if (iLevFrom > nLevsFrom) cycle !! Loop above finished -- some mass exists

            do iSpecies = 1, nSpeciesDescr
              do iLevTo = iLevStartTo, iLevEndTo
                !
                ! Fill-in the fTmpFrom(:) column in the vertical_From coordinates
                !
                call get_coefs(interpCoefVert,ixTo,iyTo,iLevTo,ptrCoefVert)
                if(ptrCoefVert%indLev(1) <= 0) cycle   ! may be, this is void cell?
                fTmp = 0.
                fTmp2 = 0.
                do iCoef = 1, fu_nCoefs(interpCoefVert)
                  if(ptrCoefVert%indLev(iCoef) == 0)exit
                  fTmp = fTmp + ptrCoefVert%weight(iCoef) * &
                                                    & pDataFrom(ptrCoefVert%indLev(iCoef))%pp(i1dFrom)
                  if(ifMoments)fTmp2 = fTmp2 + ptrCoefVert%weight_Z(iCoef) * &
                                                    & pDataFrom(ptrCoefVert%indLev(iCoef))%pp(i1dFrom)
                end do

                fTmp =  fTmp *  fMassFractionDescr(iSpecies) * fRateScaling
                fPtr => pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo)
                fPtr = fPtr + fTmp
                fMassAdded(indSpeciesOut(iSpecies)) = fMassAdded(indSpeciesOut(iSpecies)) + fTmp 

                pMassMap%ifGridValid(ilevTo,iSourceId) = .true.
                pMassMap%ifColumnValid(iSOurceId,ixTo,iyTo) = .true.

                if(ifMoments)then
                  fPtr => pMomentumMapZ%arM(indSpeciesOut(iSpecies),iSourceId,iLevTo,ixTo,iyTo)
                  fPtr = fPtr + fTmp2* fMassFractionDescr(iSpecies) * fRateScaling

                  pMomentumMapZ%ifGridValid(ilevTo,iSourceId) = .true.
                  pMomentumMapZ%ifColumnValid(iSourceId,ixTo,iyTo) = .true.
call check_mass_moment_point(pMassMap%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo), &
                      & pMomentumMapZ%arM(indSpeciesOut(iSpecies), iSourceId, iLevTo, ixTo,iyTo) , &
                      & 'Z-monent wrong after vertical-only reprojection','add_field_to_mass_map_interp')
                endif
              end do  ! iLevTo
            end do  ! iSpecies
          end do   ! ixTo
        end do   ! iyTo

      else
        !
        ! No interpolation is needed at all - just copy the requested field. Moments are intact.
        ! Note that From can still be a subset of To, they just both have to start from 1
        !
        do iyFrom = 1, nyFrom
          do ixFrom = 1, nxFrom
            do iLevFrom = 1, nLevsFrom
              do iSpecies = 1, nSpeciesDescr
                fTmp = pDataFrom(iLevFrom)%pp(ixFrom+nxFrom*(iyFrom-1)) * &
                     & fMassFractionDescr(iSpecies) * fRateScaling
                fPtr => pMassMap%arM(indSpeciesOut(iSpecies),iSourceId,iLevFrom,ixFrom,iyFrom)
                fPtr = fPtr + fTmp 
                fMassAdded(indSpeciesOut(iSpecies)) = fMassAdded(indSpeciesOut(iSpecies)) + fTmp 
                pMassMap%ifGridValid(iLevFrom,iSourceId) = .true.
                pMassMap%ifColumnValid(iSourceId,ixFrom,iyFrom) = .true.
              end do
            end do
          end do   ! ixFrom
        end do   ! iyFrom
        call stop_count('cocktail_map_to_fld_no_interp')
                
      endif ! if vertical interpolation 

    endif ! if horizontal interpolation

  end subroutine add_field_to_mass_map_interp


  !**************************************************************************

  subroutine check_mass_moment_point(mass, moment, chmsg, chplace)
    !
    ! More fancy check for wrong moments.....
    !
    implicit none
    real, intent(in) :: mass
    real, intent(inout) :: moment
    character(len=*) :: chmsg, chplace

    if( .not. (abs(moment) <= 0.49999*abs(mass) )) then !! 0 vs 0 is also fine
       if( abs(moment) < 0.5001*abs(mass) ) then
!         ! Bark and fix it
!          call msg(chmsg//"("//chplace//")"//" before fixing mass, cm", mass, moment/mass)
          moment = moment * 0.9997 
       else
          call msg(chmsg//"("//chplace//")"//" mass, cm", mass, moment/mass)
          call set_error(chmsg,chplace)
       endif
    endif

  end subroutine check_mass_moment_point


  !*************************************************************************************************
  !*************************************************************************************************
  !
  !  Initialization of maps of cocktails
  !
  !*************************************************************************************************
  !*************************************************************************************************

  subroutine update_mass_map_from_namelist(nlIn, chItemTitle, & ! namelist and item name
                                         & ptrMap, nMaps, &  ! array of mass maps and their number
                                         & valid_time, &     ! time to search
                                         & ifRandomise, &    ! if smooth the reprojection aliasing
                                         & nFields_updated)  ! how many were set
    !
    ! Handles the update of the mass map from a list of files given in a form of namelist. 
    ! Actually decodes the namelist and possible templates down to a single file and
    ! call the next routine in chain, which does more specific work. 
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(Tsilam_namelist), pointer :: nlIn
    character(len=*), intent(in) :: chItemTitle
    type(Tmass_map_ptr), dimension(:), intent(in) :: ptrMap
    type(silja_time), intent(in) :: valid_time
    integer, intent(in) :: nMaps
    logical, intent(in) :: ifRandomise
    integer, intent(out) :: nFields_updated

    ! Local variables
    type(silam_sp), dimension(:), save, pointer :: fnames
    type(grads_template) :: fname_template
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    type(silam_fformat) :: fform
    type(silam_sp) :: sp
    integer :: iFNm, iItem, nItems

    ! Check input
    !
    IF (.NOT.associated(nlIn)) THEN
      CALL set_error('Input namelist is not associated','update_mass_map_from_namelist')
      RETURN
    END IF
    nFields_updated = 0
    
    ! Get the items with given title and scan the content to determine the 
    ! format of the input file and its name as decoded from template
    nullify(ptrItems)
    call get_items(nlIn, chItemTitle, ptrItems, nItems)
    if(error .or. nItems < 1)then
      call set_error(fu_connect_strings('Could not extract namelist items:',chItemTitle), &
                   & 'update_mass_map_from_namelist')
      return
    endif
    sp%sp => fu_work_string()
    
    ! Scan the items and update the mass map values
    do iItem = 1, nItems
      
      sp%sp = fu_content(ptrItems(iItem))
      fname_template = template_missing

      fform = fu_input_file_format(sp%sp)
      if(error)return
      sp%sp=adjustl(sp%sp)
      sp%sp=adjustl(sp%sp(index(sp%sp,' ')+1:))
      if(fform%iFormat == test_field_value_flag)then
        
        CALL update_mass_map_from_file(sp%sp, fform, ptrMap, nMaps, &
                                     & valid_time, ifRandomise, nFields_updated)

      else
        
        call decode_template_string(sp%sp, fname_template)
        nullify(fnames)
        call FNm_from_single_template(fname_template, &
                                      & valid_time, fnames, &
                                      & ifAdd = .false., &
                                      & ifStrict = .false., &
                                      & ifAllowZeroFcLen = .true., &
                                      & ifWait = .false., & 
                                      & fc_step = interval_missing) 
        IF (error) return

        do iFNm = 1, size(fnames)
          if(fnames(iFNm)%sp=='')exit
          CALL update_mass_map_from_file(fnames(iFNm)%sp, fform, ptrMap, nMaps, valid_time, &
                                       & ifRandomise, nFields_updated)
        end do
      endif

   end do

   call free_work_array(sp%sp)

  end subroutine update_mass_map_from_namelist


  !****************************************************************************

  subroutine update_mass_map_from_file(chFName, FFormat, & ! File name and its format
                                     & ptrMap, nMaps, timeOfMap, ifRandomise, &
                                     & nFields_updated)

    ! Finds the necessary data and gets them to the given set of mass maps. Maps must be
    ! fully initialized. This routine just updates the fields
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chFName
    type(Tmass_map_ptr), dimension(:), intent(in) :: ptrMap
    type(silja_time), intent(in) :: timeOfMap
    type(silam_fformat) :: fformat
    integer, intent(in) :: nMaps
    logical, intent(in) :: ifRandomise
    integer, intent(inout) :: nFields_updated

    ! Local variables
    integer :: uIn, iMap, iVar,  iSubst, iMode, iWave, ix, iy, fields_found, fields_accepted, &
             & indTime, iTmp, ispecies, conversion_type, iLevIn, iLevMap, nTmp, &
             & n_out_of_range, npatched, n_failed
    real, dimension(:), pointer :: dataIn, dataMap, weightSum_mapVert
    real, dimension(:,:), pointer :: fractions, centres
    logical :: ifTmp, ifSameGrids, ifSameVerticals, ifProjectionOver, eof, ifUpsideDown
    type(silja_field_id) :: id, id_prev
    type(silja_grid), pointer :: pGrid
    type(Tmass_map), pointer :: pMap
    type(silam_vertical) :: vertIn, vertTmp
    type(silam_sp) :: spTmp
    type(silam_species) :: species
    type(silja_time), dimension(:), pointer ::  timeLst
    type(silja_field_id),dimension(:), pointer :: idList
    type(silam_vertical), dimension(:), pointer :: verticals
    character(len = *), parameter :: sub_name = 'update_mass_map_from_file'
    
    
    dataIn => fu_work_array()
    weightSum_mapVert => fu_work_array()
    ifUpsideDown = .False.
    
    IF (nMaps < 1) THEN
      CALL set_error('Mass map is not ready', sub_name)
      RETURN
    END IF

    call msg('Updating the mass map from:' + chFName)
    !
    fields_found = 0
    fields_accepted = 0

    select case(fformat%iformat)

      case(grib_file_flag) ! we haven't needed it this far ..
          call set_error('Initialization from grib is not implemented', sub_name)
          return

      case(ascii_file_flag)
          
        fractions => fu_work_array_2d()
        centres => fu_work_array_2d()
        if(error)return

        call open_ascii_file_i(chFName, uIn)
        eof = .false.
        do while(.not. eof)
          
          call read_next_field_from_ascii_file(uIn, eof, id, dataIn, fFillValue=0.)
          if(eof .or. error)exit

          fields_found = fields_found + 1
            
          call find_mass_map_4_input_id(ptrMap, nMaps, timeOfMap, id, iMap, &
                                      & iSubst, iMode, iWave, conversion_type)
          if(error)return
          if(iMap == int_missing)cycle  ! this variable variable is not needed

          call msg('Accepting the ASCII field:')
          call report(id)

          pMap => ptrMap(iMap)%ptrMassMap
          ispecies = fu_isp(isubst, imode, iwave, pMap%mapper)             

          if(fu_grid(id) == pMap%gridTemplate)then
            dataMap => dataIn
          else    ! Reproject the data to the map grid
            ! 
            if(error)return
            dataMap => fu_work_array(fu_number_of_gridpoints(pMap%gridTemplate))
            if(error)return

            call grid_data_horizontal_select(fu_grid(id), &
                                           & dataIn,&
                                           & ifTmp, &
                                           & pMap%gridTemplate,&
                                           & dataMap, &
                                           & ifRandomise, &
                                           & 5, & ! iAccuracu
                                           & fu_regridding_method(fu_quantity(id)), & !method
                                           & 0.0, & !fMissingValue
                                           & nearestPoint) !iOutside
            if(error)return

            do iTmp = 1, fu_number_of_gridpoints(pMap%gridTemplate)
              if(dataMap(iTmp) < 0)then
                call msg('Negative reprojected input at the  1D index:',iTmp,dataMap(iTmp))
                dataMap(iTmp) = 0.
              endif
            end do
          endif
          call set_vertical(fu_level(id), vertIn) 
          call convert_fld_to_mass_map_unit(dataMap, &              ! data array
                                          & conversion_type, &      ! quantity
                                          & pMap%gridTemplate, &    ! grid 
                                          & vertIn, 1)              ! vertical and level
          if(error)return
          !
          ! Get the right layer. 
          ! ATTENTION. Here we set the concentrations out of a single layer.
          ! If a fraction <1 there may be other fields with overlapping layers, so we have to 
          ! average along the vertical. In this case, the later consumed field will set the 
          ! layer. No other way: we do not know if these fields will or will not show up. 
          ! All what can be suggested, is that each layer that overlaps with the current one, 
          ! determines it completely. Then the last take will prevail
          !
          call overlap_fraction_lyr_in_vert(vertIn, pMap%vertTemplate, fractions, centres)
          if(error)return

          do iLevMap = 1, fu_NbrOfLevels(pMap%vertTemplate)
            if(fractions(1, iLevMap) < 1e-5)cycle      ! do not touch the layers you do not overlap with
            do iy = 1, pMap%ny
              do ix = 1, pMap%nx
                pMap%arM(ispecies, 1, iLevMap, ix, iy) = dataMap(ix+(iy-1)*pMap%nx) ! * fractions(iLev,1)
              end do   ! nx
            end do   ! ny
          end do
          
          if(.not. fu_grid(id) == pMap%gridTemplate) call free_work_array(dataMap)
              
          fields_accepted = fields_accepted + 1
          
        end do ! while .not. eof

        call free_work_array(fractions)
        call free_work_array(centres)
        call close_ascii_file_i(uIn)

        
      case(grads_file_flag)  

        uIn = fu_open_gradsfile_i(chFName)
        if(error)return

        ! Turn the vertical to thick layers for vertical interpolation
        if(fu_if_layer(fu_level(fu_silamVert_of_grads(uIn),1)))then
          vertIn = fu_silamVert_of_grads(uIn)
        else
          call levels_to_layers(fu_silamVert_of_grads(uIn), vertIn)
          call msg_warning("Converting levels_to_layers for "//trim(chFName), sub_name)
        endif
        
        ! Find the right timestep in grads
        indTime = fu_grads_time_index(uIn, timeOfMap, single_time, .false.)   ! direction instant field
        if(indTime == int_missing) then ! GrADS file does not cover the needed time moment
          call set_error('Needed time:' + fu_str(timeOfMap) + &
                         & ', is not covered by GrADS file:' + chFName, &
                         & sub_name)
          return
        endif


        !
        ! Scan the GrADS variables one-by-one trying to find proper mass map for each of them
        !
        do iVar = 1, fu_n_gvars(uIn)
          call get_grads_var_metadata(uIn,iVar,1,indTime,id)  !indices of gfile,gvar,glev,gtime;SILAM-id
          if(error)return
          fields_found = fields_found + 1
          
          call find_mass_map_4_input_id(ptrMap, nMaps, timeOfMap, id, iMap, &
                                      & iSubst, iMode, iWave, conversion_type)
          if(error)return
          if(iMap == int_missing)cycle  ! this grads variable is not needed
            
          pMap => ptrMap(iMap)%ptrMassMap
          ispecies = fu_isp(isubst, imode, iwave, pMap%mapper)            
          ifSameGrids = pMap%gridTemplate == fu_grid(id)                   ! same grid ?

          call msg("Updating "//trim(fu_quantity_short_string(pMap%quantity))// &
                                      & '_' // trim(fu_str(fu_species(id))))
          if(fu_multi_level_quantity(pMap%quantity))then ! go along the vertical
            !
            ! Do we have to reproject the vertical?
            !
            if(fu_cmp_verts_eq(pMap%vertTemplate, vertIn))then               ! vertical same?
              weightSum_mapVert(1:fu_NbrOfLevels(pMap%vertTemplate)) = 1.0
              ifSameVerticals = .true.
            else
              weightSum_mapVert(1:fu_NbrOfLevels(pMap%vertTemplate)) = 0.0   ! will fill them in
              ifSameVerticals = .false.
            endif

            do iLevIn = 1, fu_NbrOfLevels(vertIn)  ! cycle levels of the grads file
              
              ! read the data
!                call set_level(id, fu_level(fu_silamVert_of_grads(uIn),iLevIn)) ! actual grads vertical here!!
!                call read_field_from_grads_id(uIn, id, dataIn, fill_value_=0.)

              call read_field_from_grads_indices(uIn, iVar, iLevIn, indTime, dataIn, fill_value_=0.)



              !!
              !! Depite all atributes say that it is moment, actually SILAM dumps and restores centers of mass!!! 
              !! We can check them here
              !!
              iTmp  = fu_quantity(id)
              if (iTmp == advection_moment_X_flag) iTmp = advection_cm_X_flag
              if (iTmp == advection_moment_Y_flag) iTmp = advection_cm_Y_flag
              if (iTmp == advection_moment_Z_flag) iTmp = advection_cm_Z_flag

              call grid_dimensions(fu_grid(id), ix, iy)
              call check_quantity_range(iTmp, dataIn, ix*iy, ix, &
                                       & .false., .false., &  ! ifRequireValidity, ifSilent
                                       & n_out_of_range, npatched, n_failed)

              if (n_out_of_range > 0) &
                & call msg( fu_quantity_short_string(iTmp)+ &
                   &": off-range, patched, total",(/n_out_of_range, npatched, ix*iy/))

              if(error)return
              !
              ! Bring them to the same grid
              !
              if(ifSameGrids)then
                dataMap => dataIn  ! already right grid
!call msg('Sum of GrADS:',sum(dataMap(1:fu_number_of_gridpoints(pMap%gridTemplate))))
              else

                dataMap => fu_work_array(fu_number_of_gridpoints(pMap%gridTemplate))
                if(error)return

                call grid_data_horizontal_select(fu_grid(id), &
                                               & dataIn,&
                                               & ifTmp, &
                                               & pMap%gridTemplate,&
                                               & dataMap, &
                                               & ifRandomise, &
                                               & 5, &
                                               & fu_regridding_method(fu_quantity(id)), & !pMap%quantity), & ! linear, &
                                               & real_missing, & !0.0, &
                                               & nearestPoint)
                if(error)return

              endif
#ifdef DEBUG_INITIAL
              call grid_dimensions(fu_grid(id), ix, iy)
              call check_quantity_range(iTmp, dataMap, ix*iy, ix, &
                                       & .false., .false., &  ! ifRequireValidity, ifSilent
                                       & n_out_of_range, npatched, n_failed)

              if (n_out_of_range > 0) &
                & call msg( fu_quantity_short_string(iTmp)+ &
                   &": after grid_data_horizontal_select off-range, patched, total",(/n_out_of_range, npatched, ix*iy/))

              if(error)return
#endif
              
              call convert_fld_to_mass_map_unit(dataMap, &           ! data array
                                              & conversion_type, &   ! what we actually do
!                                                & fu_grid(id), &       ! input grid from GrADS
                                              & pMap%gridTemplate, & ! gridOut of the mass map and the field
                                              & vertIn, iLevIn)      ! vertical and level
              if(error)return
#ifdef DEBUG_INITIAL
              call grid_dimensions(fu_grid(id), ix, iy)
              call check_quantity_range(iTmp, dataMap, ix*iy, ix, &
                                       & .false., .false., &  ! ifRequireValidity, ifSilent
                                       & n_out_of_range, npatched, n_failed)

              if (n_out_of_range > 0) &
                & call msg( fu_quantity_short_string(iTmp)+ &
                   &": after convert_fld_to_mass_map_unit off-range, patched, total",(/n_out_of_range, npatched, ix*iy/))

              if(error)return
#endif
              !
              ! Having grid reprojected and unit & area scaling applied, the final step is to project vertical
              !
              if(ifSameVerticals)then
                !
                ! Same vertical means that we just pick the same level and store the whole thing there
                !
                do iy = 1, pMap%ny
                  do ix = 1, pMap%nx
                    pMap%arM(ispecies, 1, iLevIn, ix, iy) = dataMap(ix+(iy-1)*pMap%nx)
                  end do
                end do
              else
                !
                ! Have to do the projection
                !
                call project_input_layer(pMap, vertIn, iLevIn, ifProjectionOver, &
                                       & iSpecies, weightSum_mapVert, dataMap)
                if(error)return
                if(ifProjectionOver)then
                  if(.not. ifSameGrids) call free_work_array(dataMap)
                  exit  ! no need to search for more input layers
                endif

              endif  ! ifSameVerticals

              if(.not. ifSameGrids) call free_work_array(dataMap)
              
            enddo ! lev in file

            if(.not. ifSameVerticals)then
              do iLevMap = 1, fu_NbrOfLevels(pMap%vertTemplate)
                if(weightSum_mapVert(iLevMap) > 1.05)then ! this should never happen
                  call set_error('Vertical reprojection fail 3', sub_name)
                  call msg('iLevMap, weightSum_mapVert(iLevMap)', iLevMap, weightSum_mapVert(iLevMap))
                  call report(fu_level(pMap%vertTemplate, iLevMap))
                  call report(vertIn)
                  return
                endif                
                if(weightSum_mapVert(iLevMap) < 0.95)then  ! dispersion level only partially overlaps with input vertical - allowed for upper levels
                  if(iLevMap == 1)then ! error if this happens for the lowest level
                    call set_error('Vertical reprojection fail 4, grads', sub_name)
                    call msg('iLevMap, weightSum_mapVert(iLevMap)', iLevMap, weightSum_mapVert(iLevMap))
                    call report(fu_level(pMap%vertTemplate, iLevMap))
                    call report(vertIn)
                  else ! upper levels not initialized 
                    call msg_warning('Level not fully initialised', sub_name)
                    call msg('iLevMap, weightSum_mapVert(iLevMap)', iLevMap, weightSum_mapVert(iLevMap))
                    call report(fu_level(pMap%vertTemplate, iLevMap))
                    call report(vertIn)
                  endif
                endif
              enddo
            endif  ! .not. ifSameVerticals
            
          else ! 2d

!              call read_field_from_grads_id(uIn, id, dataIn, fill_value_=0.)
            call read_field_from_grads_indices(uIn, iVar, 1, indTime, dataIn, fill_value_=0.)

            if(fu_grid(id) == pMap%gridTemplate)then
              dataMap => dataIn
            else    ! Reproject the data to the map grid

              dataMap => fu_work_array(fu_number_of_gridpoints(pMap%gridTemplate))
              if(error)return

              call grid_data_horizontal_select(fu_grid(id), &
                                             & dataIn,&
                                             & ifTmp, &
                                             & pMap%gridTemplate,&
                                             & dataMap, &
                                             & ifRandomise, &
                                             & 5, &
                                             & fu_regridding_method(fu_quantity(id)), & !pMap%quantity), & ! linear, &
                                             & real_missing, & !0.0, &
                                             & nearestPoint)
              if(error)return

              do iTmp = 1, fu_number_of_gridpoints(pMap%gridTemplate)
                if(dataMap(iTmp) < 0)then
                  call msg('Negative reprojected input2')
                  dataMap(iTmp) = 0.
                endif
              end do
            endif
#ifdef DEBUG                                                    
call msg('')
call msg("Quantity: "//fu_quantity_string(fu_quantity(id)))
call report(fu_species(id))
call msg('Ave of GRADS, local unit:',sum(dataMap(1:fu_number_of_gridpoints(pMap%gridTemplate)))/fu_number_of_gridpoints(pMap%gridTemplate))
#endif
              
            call convert_fld_to_mass_map_unit(dataMap, &                          ! data array
                                            & conversion_type, &                  ! quantity
                                            & pMap%gridTemplate, &                ! grid
                                            & fu_silamVert_of_grads(uIn), iLevIn) ! vertical and level
            if(error)return

#ifdef DEBUG                                                    
call msg('Ave of GRADS, converted:',sum(dataMap(1:fu_number_of_gridpoints(pMap%gridTemplate)))/fu_number_of_gridpoints(pMap%gridTemplate))
#endif

            do iy = 1, pMap%ny
              do ix = 1, pMap%nx
                  pMap%arM(ispecies, 1, 1, ix, iy) = dataMap(ix+(iy-1)*pMap%nx)
              end do   ! nx
            end do   ! ny
            
            if(.not. fu_grid(id) == pMap%gridTemplate) call free_work_array(dataMap)
          endif
          fields_accepted = fields_accepted + 1

        end do  ! grads variables

      case(netcdf_file_flag)  

        uIn = open_netcdf_file_i(chFName, fformat)
        if(error)return        
        
        ! Find the time
        call timelst_from_netcdf_file(uIn, timeLst, nTmp)
        indTime = int_missing
        do iTmp = 1, nTmp
          if(fu_between_times(timeLst(iTmp), timeOfMap + one_minute, timeOfMap - one_minute, .true.))then
            indTime = iTmp
            exit
          endif
        enddo  
        if(indTime == int_missing) then ! file does not cover the needed time
          call msg("Timelist:")  
          do iTmp = 1, nTmp
             call report(timeLst(iTmp))
          enddo
          call set_error('Needed time:' + fu_str(timeOfMap) + &
                         & ', is not covered by NETCDF file:' + chFName, &
                         & sub_name)
          return
        endif

        call id_list_from_netcdf_file(uIn, idList, fields_found)
        id_prev = field_id_missing
        do iVar = 1, fields_found
          id = idList(iVar)
          call set_level(id_prev, fu_level(id))
          if (id == id_prev) cycle
          id_prev = id
          
          call set_valid_time(id, timeLst(indTime))

          call find_mass_map_4_input_id(ptrMap, nMaps, timeOfMap, id, iMap, &
                                          & iSubst, iMode, iWave, conversion_type)
          if(error)return
          if(iMap == int_missing)cycle  ! this grads variable is not needed
          pMap => ptrMap(iMap)%ptrMassMap
          ispecies = fu_isp(isubst, imode, iwave, pMap%mapper)             
            
            
          if(fu_multi_level_quantity(pMap%quantity))then ! reproject vertical
            
            call get_netcdf_verticals(uIn, fu_quantity(id), fu_species(id), verticals, nTmp)
            if (nTmp /= 1)then
              call set_error('Strange number of verticals for a variable', sub_name)
              return
            endif
   
            ! Turn te vertical to thick layers for vertical interpolation
            if(fu_if_layer(fu_level(verticals(1),1)))then
              vertIn = verticals(1)
              ifUpsideDown = (fu_level(vertIn,1) > fu_level(vertIn, fu_NbrOfLevels(vertIn)))
            else
              ! netcdf verticals are not sorted, can be upside down 
              ! need to temporarily sort it here to make the layers
              vertTmp = verticals(1)  
              call arrange_levels_in_vertical(vertTmp, ifUpsideDown)
              call levels_to_layers(vertTmp, vertIn)
              if (ifUpsideDown) then ! if sorting changed anything, flip it back
                  vertTmp = vertIn
                  do iLevIn = 1, fu_NbrOfLevels(vertIn)
                      call set_level(vertIn, iLevIn, fu_level(vertTmp, fu_NbrOfLevels(vertIn)-iLevIn+1))
                  enddo
              endif   
            endif
            !
            ! Do we have to reproject the vertical?
            !
            vertTmp = pMap%vertTemplate
            if(fu_cmp_verts_eq(vertTmp, vertIn))then               ! vertical same?
              weightSum_mapVert(1:fu_NbrOfLevels(pMap%vertTemplate)) = 1.0
              ifSameVerticals = .true.
            else
              weightSum_mapVert(1:fu_NbrOfLevels(pMap%vertTemplate)) = 0.0   ! will fill them in
              ifSameVerticals = .false.
            endif

            do iLevIn = 1, fu_NbrOfLevels(vertIn)
                
            ! read the data
              call set_level(id, fu_level(verticals(1),iLevIn)) ! actual netcdf vertical here
              call read_field_from_netcdf_file(uIn, id, dataIn, fill_value=0.)
              if(error)return

              if(fu_grid(id) == pMap%gridTemplate)then
                dataMap => dataIn
                ifSameGrids = .true.
              else    ! Reproject the data to the map grid
                ifSameGrids = .false.
                !
                if(error)return
                dataMap => fu_work_array(fu_number_of_gridpoints(pMap%gridTemplate))
                if(error)return

                call grid_data_horizontal_select(fu_grid(id), &
                                               & dataIn,&
                                               & ifTmp, &
                                               & pMap%gridTemplate,&
                                               & dataMap, &
                                               & ifRandomise, &
                                               & 5, &
                                               & fu_regridding_method(fu_quantity(id)), &  !pMap%quantity), & ! linear, &
                                               & real_missing, & !0.0, &
                                               & nearestPoint)
                if(error)return

                do iTmp = 1, fu_number_of_gridpoints(pMap%gridTemplate)
                  if(dataMap(iTmp) < 0)then
                    call msg('Negative reprojected input3')
                    dataMap(iTmp) = 0.
                  endif
                end do
              endif  ! ifSameGrids

#ifdef DEBUG                                                    
call msg('')
call msg("Quantity: "//fu_quantity_string(fu_quantity(id)))
call report(fu_species(id))
call msg('Ave of NetCDF, local unit:',sum(dataMap(1:fu_number_of_gridpoints(pMap%gridTemplate)))/fu_number_of_gridpoints(pMap%gridTemplate))
#endif
                
              call convert_fld_to_mass_map_unit(dataMap, &                          ! data array
                                              & conversion_type, &                  ! quantity
                                              & pMap%gridTemplate, &                ! grid
                                              & vertIn, iLevIn)               ! vertical and level
              if(error)return

#ifdef DEBUG                                                    
call msg('Ave of NetCDF, converted:',sum(dataMap(1:fu_number_of_gridpoints(pMap%gridTemplate)))/fu_number_of_gridpoints(pMap%gridTemplate))
#endif

              if(ifSameVerticals)then
                !
                ! Same vertical means that we just pick the same level and store the whole thing there
                !
                do iy = 1, pMap%ny
                  do ix = 1, pMap%nx
                    pMap%arM(ispecies, 1, iLevIn, ix, iy) = dataMap(ix+(iy-1)*pMap%nx)
                  end do
                end do
              else
                !
                ! Have to do the projection
                !
                call project_input_layer(pMap, vertIn, iLevIn, ifProjectionOver, &
                                       & iSpecies, weightSum_mapVert, dataMap)
                if(error)return
                if(ifProjectionOver .and. .not. ifUpsideDown)then
                  if(.not. ifSameGrids) call free_work_array(dataMap)
                  exit  ! no need to search for more input layers
                endif
              endif   ! ifSameVertical

              if(.not. fu_grid(id) == pMap%gridTemplate) call free_work_array(dataMap)
                
            enddo ! lev in file

            if(.not. ifSameVerticals)then
              do iLevMap = 1, fu_NbrOfLevels(pMap%vertTemplate)
                if(weightSum_mapVert(iLevMap) > 1.05)then ! this should never happen
                  call set_error('Vertical reprojection fail 3', sub_name)
                  call msg('iLevMap, weightSum_mapVert(iLevMap)', iLevMap, weightSum_mapVert(iLevMap))
                  call report(fu_level(pMap%vertTemplate, iLevMap))
                  call report(vertIn)
                  return
                endif                
                if(weightSum_mapVert(iLevMap) < 0.95)then  ! dispersion level only partially overlaps with input vertical - allowed for upper levels
                  if(iLevMap == 1)then ! error if this happens for the lowest level
                    call set_error('Vertical reprojection fail 4', sub_name)
                    call msg('iLevMap, weightSum_mapVert(iLevMap)', iLevMap, weightSum_mapVert(iLevMap))
                    call report(fu_level(pMap%vertTemplate, iLevMap))
                    call report(vertIn)
                  else ! upper levels not initialized 
                    call msg_warning('Level not fully initialised', sub_name)
                    call msg('iLevMap, weightSum_mapVert(iLevMap)', iLevMap, weightSum_mapVert(iLevMap))
                    call report(fu_level(pMap%vertTemplate, iLevMap))
                    call report(vertIn)
                  endif
                endif
              enddo
            endif  ! .not. ifSameVerticals
              
          else ! 2d
            call read_field_from_netcdf_file(uIn, id, dataIn, fill_value=0.)
            if(fu_grid(id) == pMap%gridTemplate)then
              dataMap => dataIn
            else    ! Reproject the data to the map grid
              !
              if(error)return
              dataMap => fu_work_array(fu_number_of_gridpoints(pMap%gridTemplate))
              if(error)return

              call grid_data_horizontal_select(fu_grid(id), &
                                             & dataIn,&
                                             & ifTmp, &
                                             & pMap%gridTemplate,&
                                             & dataMap, &
                                             & ifRandomise, &
                                             & 5, &
                                             & fu_regridding_method(fu_quantity(id)), & !pMap%quantity), & ! linear, &
                                             & real_missing, & !0.0, &
                                             & nearestPoint)
              if(error)return

              do iTmp = 1, fu_number_of_gridpoints(pMap%gridTemplate)
                if(dataMap(iTmp) < 0)then
                  call msg('Negative reprojected input')
                  dataMap(iTmp) = 0.
                endif
              end do
            endif
                
            call convert_fld_to_mass_map_unit(dataMap, &                          ! data array
                                            & conversion_type, &                  ! quantity
                                            & pMap%gridTemplate, &                ! grid
                                            & fu_silamVert_of_grads(uIn), iLevIn) ! vertical and level
            if(error)return

            do iy = 1, pMap%ny
              do ix = 1, pMap%nx
                  pMap%arM(ispecies, 1, 1, ix, iy) = dataMap(ix+(iy-1)*pMap%nx)
              end do   ! nx
            end do   ! ny
              
            if(.not. fu_grid(id) == pMap%gridTemplate) call free_work_array(dataMap)
          endif

          fields_accepted = fields_accepted + 1
        end do  !  variables
        
      case(test_field_value_flag)   
      
        call msg('Making test field (update_mass_map_from_file):' + chFName)

        spTmp%sp => fu_work_string()
        
        ! Find the right map - need the quantity
        read(unit=chFName, iostat=uIn, fmt=*) spTmp%sp
        if(fu_fails(uIn == 0,'Cannot read quantity name from:'+chFName, sub_name))return
        call decode_id_params_from_io_str(spTmp%sp, &     ! string to decode
                                        & .false., &      ! ifMultiLevel
                                        & iVar, &         ! decoded quantity
                                        & species, &      ! decoded species
                                        & .true.)         ! scream if fail
        call free_work_array(spTmp%sp)
        if(fu_fails(iVar /= int_missing,'Strange quantity name in:' + chFName, sub_name))return

        id = fu_set_field_id_simple(fmi_silam_src, iVar, timeOfMap, surface_level, species)

        call find_mass_map_4_input_id(ptrMap, nMaps, timeOfMap, id, iMap, &
                                    & iSubst, iMode, iWave, conversion_type)
        if(error)return
       !        if(fu_fails(iMap /= int_missing, 'Unused test field:'+chFName,sub_name))return
        if (iMap /= int_missing) then

            pMap => ptrMap(iMap)%ptrMassMap
            ispecies = fu_isp(isubst, imode, iwave, pMap%mapper) 
            pGrid => pMap%gridTemplate

            ! All levels of the variable will be initiliased with this value
            do iLevMap = 1, fu_NbrOfLevels(pMap%vertTemplate)
              call make_test_field(chFName,  id,  dataIn,  timeOfMap,  zero_interval, pGrid)
              call convert_fld_to_mass_map_unit(dataIn,  conversion_type,  pMap%gridTemplate, &
                                              & pMap%vertTemplate, iLevMap) 
              do iy = 1, pMap%ny
                do ix = 1, pMap%nx
                   pMap%arM(ispecies,1,iLevMap,ix,iy) = dataIn(ix+(iy-1)*pMap%nx)
                end do   ! nx
              end do   ! ny
            end do   ! levels
            fields_accepted = fields_accepted + 1
        else
                call msg_warning('Unused test field:'+chFName,sub_name)
        endif
        
      case default   
        call msg('Not supported file type for initialization',fformat%iformat)
        call set_error('Not supported file type',sub_name)

    end select


    select case(fformat%iformat)
      case(grads_file_flag)
        call close_gradsfile_i(uIn)
      case(netcdf_file_flag)
        call close_netcdf_file(uIn)
    end select

    call msg('Variables found: ',fields_found)
!    if(fields_found > 0 .and. fields_accepted == 0)then
!      call set_error('Zero updated fields: update failed',sub_name)
!    else
      call msg('Variables updated: ',fields_accepted)
!    endif

    nFields_updated = nFields_updated + fields_accepted
    
    call free_work_array(dataIn)
    call free_work_array(weightSum_mapVert)

    
    CONTAINS
    
    subroutine project_input_layer(pMap, vertIn, iLevIn, ifProjectionOver, iSpecies, &
                                 & weightSum_mapVert, dataMap)
      !
      ! Projects a single input layer of vertIn pointed by iLevIn to the receiving vertical of 
      ! the mass map pMap
      !
      implicit none

      ! Imported parameters
      type(TMass_map), intent(inout) :: pMap
      type(silam_vertical), intent(in) :: vertIn
      integer, intent(in) :: iLevIn, iSpecies
      real, dimension(:), pointer :: weightSum_mapVert, dataMap
      logical, intent(out) :: ifProjectionOver

      ! Local variables
      integer :: ix, iy, iLevMap
      real :: weightSum_inLev, fWeight

      ifProjectionOver = .false.
      weightSum_inLev = 0.0
      do iLevMap = 1, fu_NbrOfLevels(pMap%vertTemplate)
        if(pMap%quantity == mass_in_air_flag)then  ! mass conservation
          fWeight = fu_vert_overlap_fraction(fu_level(vertIn, iLevIn), fu_level(pMap%vertTemplate, iLevMap))
        else ! averaging
          fWeight = fu_vert_overlap_fraction(fu_level(pMap%vertTemplate, iLevMap), fu_level(vertIn, iLevIn))
        endif
        if(error)return
        if(fWeight < 1.e-10)cycle
        do iy = 1, pMap%ny
          do ix = 1, pMap%nx
            if(weightSum_mapVert(iLevMap) .eps. 0.0)then
              pMap%arM(ispecies, 1, iLevMap, ix, iy) = dataMap(ix+(iy-1)*pMap%nx) * fWeight
            else
              pMap%arM(ispecies, 1, iLevMap, ix, iy) = pMap%arM(ispecies, 1, iLevMap, ix, iy) + &
                                                                    & dataMap(ix+(iy-1)*pMap%nx) * fWeight
            endif
          end do   ! nx
        end do   ! ny
                  
        weightSum_mapVert(iLevMap) = weightSum_mapVert(iLevMap) +  &
                                & fu_vert_overlap_fraction(fu_level(pMap%vertTemplate, iLevMap), &
                                                         & fu_level(vertIn, iLevIn))
        weightSum_inLev = weightSum_inLev + &
                                & fu_vert_overlap_fraction(fu_level(vertIn, iLevIn), &
                                                         & fu_level(pMap%vertTemplate, iLevMap))
      enddo !lev in mass map
      !
      ! Did we manage to project this whole level into the new vertical?
      !
      if(weightSum_inLev > 1.05)then ! this should never happen
        call set_error('Vertical reprojection fail 1', sub_name)
        call msg('iLevIn, weightSum_inLev', iLevIn, weightSum_inLev)
        call report(fu_level(vertIn, iLevIn))
        call report(pMap%vertTemplate)
        return
      endif                
      if(weightSum_inLev < 0.95)then 
        ! input level only partially overlaps with dispersion vertical - allowed for upper levels
        if(fu_NbrOfLevels(vertIn) > 1 .and. iLevIn == 1 .and. (.not. ifUpsideDown))then ! .or. &
!                                           & (iLevIn == fu_NbrOfLevels(vertIn) .and. ifUpsideDown)))then 
          ! error if this happens for the lowest level, unless there is just one
          call set_error('Vertical reprojection fail 2', sub_name)
          call msg('iLevIn, weightSum_inLev', iLevIn, weightSum_inLev)
          call report(fu_level(vertIn, iLevIn))
          call report(pMap%vertTemplate)
        else ! above the dispersion vertival - no need to continue
          ifProjectionOver = .true.
        endif
      endif
    end subroutine project_input_layer

  end subroutine update_mass_map_from_file


  !**********************************************************************************

  subroutine find_mass_map_4_input_id(ptrMap,nMaps,timeOfMap,id,iMap,&
                                    & iSubst,iMode,iWave, conversion)
    !
    ! Check-find the quantity in the given maps that satisfy the given id
    !
    implicit none

    ! Imported parameters
    type(Tmass_map_ptr), dimension(:), intent(in) :: ptrMap
    integer, intent (in) :: nMaps
    type(silja_time), intent(in) :: timeOfMap
    type(silja_field_id), intent(in) :: id
    integer, intent(out) :: iMap, iSubst,iMode,iWave
!    real, intent(out) :: fLevelIndex
    integer, intent(out) :: conversion
    
    ! Local variables
    integer :: iMapTmp, iTmp
    type(Tmass_map), pointer :: pMap
    real :: wave1, wave2

    iMap = int_missing

    iMapTmp = 1
    do while (iMapTmp <= nMaps)
      ! can the mass map quantity be directly diagnosed from the given one?
      if(fu_ifDiagnosticQuantity(fu_quantity(id), ptrMap(iMapTmp)%ptrMassMap%quantity))exit
      iMapTmp  = iMapTmp  +1
    end do
    if(error)return
    if(iMapTmp  > nMaps) return
    pMap => ptrMap(iMapTmp)%ptrMassMap

    !
    ! Valid period of the ID must cover the time of the map
    !
    if(.not. fu_between_times(timeOfMap, fu_valid_time(id), &
                                       & fu_valid_time(id) + fu_validity_length(id), .true.))return
    if(error)return
    !
    ! Find the substance name and index
    !
    iSubst = int_missing
    
    if (fu_cocktail_name(id) /= '') then !Some stupidity check
       call set_error("Trying to init massmap  from cocktail field", "find_mass_map_4_input_id")
      return
    endif

    do iTmp = 1, fu_nSubst(pMap%mapper)
      if(fu_substance_name(id) == fu_substance_name(iTmp, pMap%mapper)) then
        iSubst = iTmp
        exit
      endif
    end do
    if(iSubst == int_missing)return
    !
    ! With aerosol - be careful. Substance may not have the specific aerosol size class but it
    ! can still be reported to output file - with zeroes, of course.
    ! Should the specific mode not be found, skip the field!
    ! However, if the substance has only one state (a specific mode or gas-only), it can be accepted
    ! by-default because in this case the specific state and mode may be skipped from GrADS name as 
    ! self-evident
    !
    iMode = int_missing
    do iTmp = 1, fu_nModes(iSubst, pMap%mapper)
      if (fu_species_match( fu_species(id), pMap%species(fu_iSp(iSubst, iTmp, pMap%mapper))) ) then
        iMode = iTmp
        exit
      endif
    end do

    if(iMode == int_missing)return

    iWave = int_missing
    ! Note: in the current implementation, different modes can have
    ! different number of wavelengths. 
    do iTmp = 1, fu_nWaves(iSubst,iMode,pMap%mapper)
      wave1 = fu_optical_wave_length(id)
      wave2 = fu_optical_wave_length(pMap%species(fu_iSp(iSubst, iMode, iTmp, pMap%mapper)))
      if ((wave1 .eps. real_missing) .and. (wave2 .eps. real_missing)) then
        iwave = itmp
        exit
      end if
      if (1e6*fu_optical_wave_length(id) .eps. &
        & 1e6*fu_optical_wave_length(pMap%species(fu_iSp(iSubst, iMode, iTmp, pMap%mapper))))then
        iWave = iTmp
        exit
      endif
    end do
    if(iWave == int_missing)return

    iMap = iMapTmp

    ! Finally, unit conversion if needed. We have some implicit assumptions:
    ! depositions are per area in output
    !
    if (pMap%quantity == drydep_flag .or. pMap%quantity == wetdep_flag) then
      conversion = per_area_inv
    else if (pMap%quantity == mass_in_air_flag) then
      if (fu_quantity(id) == mass_in_air_flag) then
        conversion = volume_ratio
      else if (fu_quantity(id) == concentration_flag) then
        conversion = per_volume_inv
      else if (fu_quantity(id) == volume_mixing_ratio_flag) then
        conversion = vmr_2_mass      
      end if
    else if (pMap%quantity == concentration_flag) then
      if (fu_quantity(id) == concentration_flag) then
        conversion = unity
      else if (fu_quantity(id) == mass_in_air_flag) then
        conversion = per_volume
      else if (fu_quantity(id) == volume_mixing_ratio_flag) then
        conversion = vmr_2_cnc      
      end if
    else if (pMap%quantity == fu_quantity(id) .and. &
           & any(pMap%quantity == (/advection_moment_X_flag, &
                              & advection_moment_Y_flag, advection_moment_Z_flag/))  ) then
           conversion = volume_ratio  !! Should be same as mass                              
    else
      call msg_warning('Assuming unit conversion between ' &
                     & // fu_quantity_short_string(pMap%quantity) // ' and ' &
                     & // fu_quantity_short_string(fu_quantity(id)), 'find_mass_map_4_input_id')
      conversion = unity
    end if

  end subroutine find_mass_map_4_input_id


  !************************************************************************

  subroutine convert_fld_to_mass_map_unit(dataInOut, conversion, & !gridIn, verticalIn, &
                                                             & gridOut, verticalOut, iLev)
    !
    ! Depending on the given quantity, makes the scaling from input/output parameter
    ! to the corresponding mass-map quantity.
    ! The whole mess is caused by the ambiguity in the SILAM quantity definitions:
    ! advection works with masses in grid cell while the quantity is called concentration
    ! To the output and chemistry it certainly goes as true concentration. So far, the 
    ! storage follows the advection definition, which makes concentration quantity to
    ! have units of mass per grid cell. This scaling has to be done here.
    !
    implicit none

    ! Imported parameters
    real, dimension(:), pointer :: dataInOut
    integer, intent(in) :: conversion, iLev
    type(silja_grid), intent(in) :: gridOut  !, gridIn
    type(silam_vertical), intent(in) :: VerticalOut !, verticalIn

    ! Local variables
    real, dimension(:), pointer :: xSizeIn => null(), ySizeIn => null(), &
                        & xSizeOut => null(), ySizeOut => null()
    real :: fTmp, zSizeIn, zSizeOut, border, dx_deg, dy_deg, corner_lat_N, t, p, x0, y0
    integer :: iTmp, nx, ny, ix, iy, nxIn, nyIn
    type(silja_level) :: level
    logical :: ifAllocate = .False. , if_geo_grid

    ifAllocate = .false.

!    if(conversion == volume_ratio)then
!        if(fu_gridtype(gridIn) == lonlat)then
!          ifAllocate = .true.
!          xSizeIn => fu_work_array()
!          ySizeIn => fu_work_array()
!          CALL lonlat_grid_parameters(gridIn,& ! Get parameters of the grid
!                             & fTmp, corner_lat_N, if_geo_grid, &
!                             & nx, ny, &
!                             & fTmp, fTmp, & 
!                             & dx_deg, dy_deg)
!          DO iy=1,ny
!            DO ix=1,nx
!              xSizeIn(ix+(iy-1)*nx)= fu_dx_deg_to_m(dx_deg, (iy-1)*dy_deg+corner_lat_N)
!            END DO
!          END DO
!          DO iy=1,ny   ! Fill-in temporary array with y-size
!            DO ix=1,nx
!              ySizeIn(ix+(iy-1)*nx) = fu_dy_deg_to_m(dy_deg)
!            END DO
!          END DO
!        else
!          ifAllocate = .false.
!          xSizeIn => fu_dx_fld_m(gridIn)
!          ySizeIn => fu_dy_fld_m(gridIn)
!        endif
!      level = fu_level(verticalIn,iLev)
!      if(fu_if_layer(level)) then
!        zSizeIn = fu_layer_thickness_m(level)
!      else
!        border = 0.
!        do iTmp = 1, iLev
!          border = border + 2.*(fu_level_height(fu_level(verticalOut,iTmp)) - border)
!        end do
!        zSizeIn = 2. * (border - fu_level_height(level))
!      endif  ! IfLayer
!    endif  ! volume ratio conversion

    if(conversion == per_area .or. conversion == per_area_inv .or. &
     & conversion == per_volume .or. conversion == per_volume_inv .or. &
     & conversion == vmr_2_mass .or. conversion == mass_2_vmr .or. &
     & conversion == volume_ratio)then
        if(fu_gridtype(gridOut) == lonlat)then
          ifAllocate = .true.
          xSizeOut => fu_work_array()
          ySizeOut => fu_work_array()
          CALL lonlat_grid_parameters(gridOut,& ! Get parameters of the grid
                             & fTmp, corner_lat_N, if_geo_grid, &
                             & nx, ny, &
                             & x0, y0, & 
                             & dx_deg, dy_deg)
          DO iy=1,ny
            DO ix=1,nx
              xSizeOut(ix+(iy-1)*nx)= fu_dx_deg_to_m(dx_deg, (iy-1)*dy_deg+corner_lat_N)
            END DO
          END DO
          DO iy=1,ny   ! Fill-in temporary array with y-size
            DO ix=1,nx
              ySizeOut(ix+(iy-1)*nx) = fu_dy_deg_to_m(dy_deg)
            END DO
          END DO
        else
          xSizeOut => fu_dx_fld_m(gridOut)
          ySizeOut => fu_dy_fld_m(gridOut)
        endif
    endif
    zSizeOut = 1.
    if(conversion == per_volume .or. conversion == per_volume_inv .or. &
     & conversion == vmr_2_mass .or. conversion == mass_2_vmr .or. &
     & conversion == volume_ratio) then
      level = fu_level(verticalOut,iLev)
      if(fu_if_layer(level)) then
        zSizeOut = fu_layer_thickness_m(level)
      else
        border = 0.
        do iTmp = 1, iLev
          border = border + 2.*(fu_level_height(fu_level(verticalOut,iTmp)) - border)
        end do
        zSizeOut = 2. * (border - fu_level_height(level))
      endif  ! IfLayer
    endif

    if (conversion == vmr_2_cnc .or. conversion == vmr_2_mass .or. &
      & conversion == cnc_2_vmr .or. conversion == mass_2_vmr)then
         
      call msg_warning('Standard atmosphere pressure and temperature are used for conversion', &
                     & 'convert_fld_to_mass_map_unit')   
      
      select case(fu_leveltype(fu_level(verticalOut,iLev)))
        case(constant_pressure)
          call us_standard_atmosphere(fu_height_for_press(fu_pr_level_pressure(&
                                                          & fu_level(verticalOut,iLev))), fTmp, p, t)
        case(layer_btw_2_pressure)
          call us_standard_atmosphere(fu_height_for_press(fu_layer_centre_value( &
                                                          & fu_level(verticalOut,iLev))), fTmp, p, t)

        case(hybrid)
          call us_standard_atmosphere(fu_height_for_press(fu_hybrid_level_pressure( &
                                          & fu_level(verticalOut,iLev), std_pressure_sl)), fTmp, p, t)

        case(layer_btw_2_hybrid)
          call us_standard_atmosphere((fu_height_for_press(fu_hybrid_level_pressure( &
                                            & fu_upper_boundary_of_layer(fu_level(verticalOut,iLev)), &
                                            & std_pressure_sl)) + & 
                                     & fu_height_for_press(fu_hybrid_level_pressure( &
                                            & fu_lower_boundary_of_layer(fu_level(verticalOut,iLev)), &
                                            & std_pressure_sl))) / 2., fTmp, p, t)

        case(constant_height)
          call us_standard_atmosphere(fu_level_height(fu_level(verticalOut,iLev)), fTmp, p, t)

        case(layer_btw_2_height)
          call us_standard_atmosphere(fu_layer_centre_value(fu_level(verticalOut,iLev)), fTmp, p, t)

        case(surface)
          t = 1. 
          p = 1. 
        case default
          call set_error('Cannot compute standard pressure and temperature for this level', &
                       & 'convert_fld_to_mass_map_unit')
          call report(fu_level(verticalOut,iLev))
          return
      end select
      t = t * 288.15
      p = p * std_pressure_sl
    endif
   
    select case(conversion)
      case(per_area_inv, per_volume_inv)
!call msg('dataIn Ave:', sum(dataInOut(1:fu_number_of_gridpoints(gridOut)))/fu_number_of_gridpoints(gridOut))
        DO iTmp = 1, fu_number_of_gridpoints(gridOut)
!if(mod(iTmp, 1000) == 0) call msg('iTmp, x,y,z size:' + fu_str(xSizeOut(iTmp)) + ',' + &
!                                    & fu_str(ySizeOut(iTmp)) + ',' + fu_str(zSizeOut))
          dataInOut(iTmp) = dataInOut(iTmp) * zSizeOut * xSizeOut(iTmp)*ySizeOut(iTmp)
        END DO
!call msg('dataInOut Ave:', sum(dataInOut(1:fu_number_of_gridpoints(gridOut)))/fu_number_of_gridpoints(gridOut))
       
      case(per_area, per_volume)     
        DO iTmp = 1, fu_number_of_gridpoints(gridOut)
          dataInOut(iTmp) = dataInOut(iTmp) / (zSizeOut * xSizeOut(iTmp)*ySizeOut(iTmp))
        END DO

      case(volume_ratio)   ! so far doing nothing
!!        call grid_dimensions(gridIn, nxIn, nyIn)
!!        DO iy = 1, ny  ! in output grid
!!          do ix = 1, nx  ! in output grid
!!            iTmp = ix + (iy-1) * nx
!!            call project_point_to_grid_xy(gridOut, real(ix), real(iy), gridIn, ixIn, iyIn)
!!            jTmp = ixIn + (iyIn-1)* nxIn
!!            dataInOut(iTmp) = dataInOut(iTmp) / (zSizeIn * xSizeIn(jTmp) * ySizeIn(jTmp)) * &
!!                                              & zSizeOut * xSizeOut(iTmp) * ySizeOut(iTmp)
!!          END DO
!!        END DO
!!        call free_work_array(xSizeIn)
!!        call free_work_array(ySizeIn)

      case(vmr_2_cnc, vmr_2_mass)         
!call msg('vmr2cnc dataIn Ave:', sum(dataInOut(1:fu_number_of_gridpoints(gridOut)))/fu_number_of_gridpoints(gridOut))
        if(conversion == vmr_2_mass)then
          do iTmp = 1, fu_number_of_gridpoints(gridOut) 
            fTmp = dataInOut(iTmp)
            dataInOut(iTmp) = dataInOut(iTmp) * p /( gas_constant_uni * t) * zSizeOut * xSizeOut(iTmp)*ySizeOut(iTmp)
!if(mod(iTmp, 1000) == 0) call msg('iTmp, x,y,z size, In, Out:' + fu_str(iTmp), (/xSizeOut(iTmp),ySizeOut(iTmp),zSizeOut,fTmp,dataInOut(iTmp)/))
          enddo
        else
          do iTmp = 1, fu_number_of_gridpoints(gridOut) 
            dataInOut(iTmp) = dataInOut(iTmp) * p /( gas_constant_uni * t)
          enddo
        endif
!call msg('vmr2cnc dataInOut Ave:', sum(dataInOut(1:fu_number_of_gridpoints(gridOut)))/fu_number_of_gridpoints(gridOut))
      case(cnc_2_vmr, mass_2_vmr)         
        if(conversion == mass_2_vmr)then
          do iTmp = 1, fu_number_of_gridpoints(gridOut) 
            dataInOut(iTmp) = dataInOut(iTmp) * gas_constant_uni * t / (p * zSizeOut * xSizeOut(iTmp)*ySizeOut(iTmp))
          enddo
        else
          do iTmp = 1, fu_number_of_gridpoints(gridOut) 
            dataInOut(iTmp) = dataInOut(iTmp) * gas_constant_uni * t / p 
          enddo
        endif
      case (unity)
        ! no scaling
        
      case default
        call set_error('Strange conversion request', 'convert_fld_to_mass_map_unit')
        return
    end select
    
    if(ifAllocate)then
      call free_work_array(xSizeOut)
      call free_work_array(ySizeOut)
    endif
    
  end subroutine convert_fld_to_mass_map_unit


  !**********************************************************************************
  
  logical function fu_if_eulerian_present(simulation_type)
    !
    ! Answers whether the simulation include Eulerian dynamics
    !
    implicit none
    integer, intent(in) :: simulation_type
    fu_if_eulerian_present = simulation_type == eulerian_flag .or. &
                           & simulation_type == hybrid_flag
  end function fu_if_eulerian_present

  !**********************************************************************************
  
  logical function fu_if_lagrangian_present(simulation_type)
    !
    ! Answers whether the simulation include Eulerian dynamics
    !
    implicit none
    integer, intent(in) :: simulation_type
    fu_if_lagrangian_present = simulation_type == lagrangian_flag .or. simulation_type == hybrid_flag
  end function fu_if_lagrangian_present

  !**********************************************************************************
  
  logical function fu_if_hybrid_present(simulation_type)
    !
    ! Answers whether the simulation include Eulerian dynamics
    !
    implicit none
    integer, intent(in) :: simulation_type
    fu_if_hybrid_present = simulation_type == hybrid_flag
  end function fu_if_hybrid_present


  !**********************************************************************************
  !**********************************************************************************
  !
  ! Functions for merging the data from MassMap and Lagrangian arDyn & arMass into stack
  ! Used in the output.
  ! Stupid stuff !!!
  ! Have to repeat three times the whole logic in order to cover all three classes of data
  !
  !**********************************************************************************
  !**********************************************************************************

  !********************************************************************************
  
  !*****************************************************************************
  
  subroutine get_scaling(quantity, iVerticalTreatment, &
        & ifPerDensity, ifPerArea, ifPerDz, ifTimestepIntegrated, &
                       & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize, &
                       & data_buf)
    !
    ! Checks whether scaling with air density and/or cell size is needed and returns the 
    ! appropriate quantities
    !
    implicit none
    
    ! Imported parameters
    type(Tfield_buffer), intent(in), optional :: data_buf
    integer, intent(in) :: quantity, iVerticalTreatment
    logical, intent(out) :: ifPerDensity, ifPerArea, ifPerDz, ifTimestepIntegrated
    type(field_3d_data_ptr), pointer :: dz_past_3d, dz_future_3d, rho_past_3d, rho_future_3d
    real, dimension(:), pointer :: xSize, ySize
    
    ! Local variables
    integer :: iTmp

    ifPerDensity = .false.

    ifTimestepIntegrated = quantity == emission_intensity_flag .or. quantity == emission_flux_flag

    select case(quantity)
      case(concentration_flag, particle_counter_flag, optical_density_flag)
        ifPerArea = .true.
        ifPerDz = .true.
      case(optical_column_depth_flag, drydep_flag, wetdep_flag, concentration_2m_flag)
        ifPerArea = .true.
        ifPerDz = .false.
      case(volume_mixing_ratio_flag)
        ifPerArea = .true.
        ifPerDensity = .true.
        ifPerDz = .true.
      case(advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag, emission_intensity_flag)
        ifPerArea = .false.
        ifPerDz = .false.
      case(emission_flux_flag)
        ifPerArea = .true.
        ifPerDz = .false.
      case default
        call msg('Unknown  quantity:',quantity)
        call set_error('Unknown 3D(?) quantity','get_scaling')
        return
    end select ! once again, the quantity

    select case(iVerticalTreatment)
      case(do_nothing_flag,lowest_level_flag)          ! all are 3D
         !Relax...
      case(integrate_column_flag)
        ifPerDz = .false.
      case default
        call set_error('Strange vertical treatment of variable:' + fu_str(iVerticalTreatment), &
                     & 'get_scaling')
        return
    end select

    if(present(data_buf))then
      if(ifPerDensity) then
        iTmp = fu_index(data_buf, air_density_flag)
        if(fu_fails(iTmp >= 1,'No density in the buffer','get_scaling'))return
        rho_past_3d => data_buf%p4d(iTmp)%past
        rho_future_3d => data_buf%p4d(iTmp)%past
      else
       nullify(rho_past_3d)
       nullify(rho_future_3d)
      endif

      if(ifPerArea)then
        iTmp = fu_index(data_buf, cell_size_x_flag)
        if(fu_fails(iTmp >= 1,'No x-size in the buffer','get_scaling'))return
        xSize => data_buf%p2d(iTmp)%present%ptr
        iTmp = fu_index(data_buf, cell_size_y_flag)
        if(fu_fails(iTmp >= 1,'No y-size in the buffer','get_scaling'))return
        ySize => data_buf%p2d(iTmp)%present%ptr

      else
        nullify(xSize)
        nullify(ySize)
        nullify(dz_past_3d)
        nullify(dz_future_3d)
      endif
      if(ifPerDz)then
        iTmp = fu_index(data_buf, cell_size_z_flag)
        if(fu_fails(iTmp >= 1,'No layer thickness in the buffer','get_scaling'))return
        dz_past_3d => data_buf%p4d(iTmp)%past
        dz_future_3d => data_buf%p4d(iTmp)%future
      else
        nullify(dz_past_3d)
        nullify(dz_future_3d)
      endif
    else   ! no data buffer, nullify all pointers
      nullify(rho_past_3d)
      nullify(rho_future_3d)
      nullify(xSize)
      nullify(ySize)
      nullify(dz_past_3d)
      nullify(dz_future_3d)
    endif

  end subroutine get_scaling


  !********************************************************************************

  subroutine merge_mass_map_2_mass_map(pMassMapLinks, &      ! Structure to merge
                                     & data_buf, writingTime, now, timestepIn, ifFirstTime, ifLastOutput)
    !
    ! Merges the mass map to the cumulating mass map. Uses all species addressed by the adaptor,
    ! all levels or their sum or the fiels one (depending on vertical treatment). Normalization 
    ! allowed is per-volume and per-air-denstiy
    ! Does the same thing as merge_vector_data_to_stack below but for mass_map.
    ! ATTENTION.
    ! This sub takes ONLY input mass maps with instaint fields and serves ONLY merging that
    ! requires averaging, i.e. some action that calls for intermediate mass map. All simpler
    ! (or more complicated) requests are handled directly by other massMap to stack subroutines.
    !
    implicit none

    ! Imported parameters
    type(TMassMapLink), intent(inout) :: pMassMapLinks     ! Masses to merge, metadata
    type(TField_buffer), intent(in) :: data_buf
    type(silja_time), intent(in) :: writingTime, now
    type(silja_interval), intent(in) :: timestepIn
    logical, intent(in) :: ifFirstTime, ifLastOutput
    !
    ! Local variables
    !
    integer :: iLink, iTmp
    real :: fTmp
    logical :: ifPerDensity, ifPerArea, ifPerDz, ifPerDzDynamic, &
               & ifTimestepIntegrated, ifStartAveragingPeriod
    type(field_3d_data_ptr), pointer :: rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d
    real, dimension(:), pointer :: xSize, ySize

    !
    ! Grand cycle through the links
    !
    do iLink = 1, nMassMapLinks

      if(pMassMapLinks%iVerticalTreatment(iLink) == int_missing)exit ! all done

      !
      ! Some quantities need to be scaled from mass map for the output. Since scaling
      ! can be time-dependent (air density or dz), have to do it every time step
      !
      call get_scaling(pMassMapLinks%pMMOut(iLink)%quantity, &
                     & pMassMapLinks%iVerticalTreatment(iLink), &
                     & ifPerDensity, ifPerArea, ifPerDz, ifTimestepIntegrated, &
                     & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize, &
                     & data_buf)
      !
      ! Can we do it once at the writing moment or have to waste time at every collection?
      ! Basically, simple: we can ask for meteo dependence of projection from the MMOut vertical
      ! to metric height.
      !
      if(ifPerDz)then
        ifPerDzDynamic = fu_if_level_meteo_dependent(fu_leveltype( &
                                           & pMassMapLinks%pMMOut(iLink)%vertTemplate))
      else
        ifPerDzDynamic = .false.
      endif

      !
      ! Input mass map is either instant or cumulative (since the start of the run).
      ! Output can be as-is, instant, averaged, cumulative, mean-over-some-time
      ! Handle this all!
      !
      if(fu_accumulated_quantity(pMassMapLinks%pMMIn(iLink)%ptrMassMap%quantity))then
        !
        ! Input quantity is cumulative from the start of the run. Any operation hammers down to
        ! subtracting its value at some moment from its value at the output time moment.
        !
        if(pMassMapLinks%iAveragingType(iLink) == iMeanLastHrs)then
          if(fu_abs(now - writingTime) > pMassMapLinks%AveragingPeriod(iLink))cycle
        endif

        select case(pMassMapLinks%iAveragingType(iLink))
          case(iAsIs, iCumulative)                      ! for cumulative quantity, take the last value
            if(.not. now == writingTime)cycle
            call merge_instant_MMdata(pMassMapLinks%pMMIn(iLink)%ptrMassMap, &   ! input, instant
                                    & pMassMapLinks%pMMOut(iLink), &  ! output, meanXXHours
                                    & one_second, &                              ! no integration
                                    & pMassMapLinks%iVerticalTreatment(iLink), & ! what to do with vertical
                                    & pMassMapLinks%adaptor(iLink), &            ! link input->output species
                                    & ifPerDensity, ifPerDzDynamic, &
                                    & data_buf%weight_past)
          case(iMeanLastHrs)                ! start new period or subtract stored value and scale with period
            if(now == writingTime)then          ! finish the period and write the stuff
              ifStartAveragingPeriod = .false.
            elseif(timestepIn * 0.5 - one_second > &    ! start new averaging period
                 & now - (writingTime - pMassMapLinks%AveragingPeriod(iLink)))then
              pMassMapLinks%IntegrationStart(iLink) = now
              ifStartAveragingPeriod = .true.
            else
              cycle   ! neither start nor end of integration
            endif
            call merge_cumul_MMdata(pMassMapLinks%pMMIn(iLink)%ptrMassMap, &     ! input, instant
                                  & pMassMapLinks%pMMOut(iLink), &    ! output, meanXXHours
                                  & now - pMassMapLinks%IntegrationStart(iLink), &   ! integration period
                                  & pMassMapLinks%iVerticalTreatment(iLink), &   ! what to do with vertical
                                  & pMassMapLinks%adaptor(iLink), &              ! link input->output species
                                  & ifPerDensity, ifPerDzDynamic, &
                                  & data_buf%weight_past, &
                                  & ifStartAveragingPeriod, &         ! store the value-to-be-subtracted
                                  & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize)
          case(iAverage, iInstant)             ! subtract stored value and scale with period
            if(.not. now == writingTime)cycle
            call merge_cumul_MMdata(pMassMapLinks%pMMIn(iLink)%ptrMassMap, &     ! input, instant
                                  & pMassMapLinks%pMMOut(iLink), &    ! output, meanXXHours
                                  & now - pMassMapLinks%IntegrationStart(iLink), &   ! integration period
                                  & pMassMapLinks%iVerticalTreatment(iLink), &   ! what to do with vertical
                                  & pMassMapLinks%adaptor(iLink), &              ! link input->output species
                                  & ifPerDensity, ifPerDzDynamic, &
                                  & data_buf%weight_past, &
                                  & .false., &        ! start of period is done by start_new_output_period
                                  & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize)
          case default
            call set_error('Only Average,cumulative,mean-XX-hours allowed, not:' + &
                         & fu_str(pMassMapLinks%iAveragingType(iLink)),'merge_mass_map_2_mass_map')
            return
        end select   ! iAveraginType

      else
        !
        ! Input quantity is sometimes instant. We have to sum it up over the averaging interval.
        ! Outside it - do nothing
        !
        if(ifFirstTime .and. pMassMapLinks%iAveragingType(iLink) /= iCumulative)then
          iTmp = iAsIs      ! for the first time, cannot do anything: only instant is possible
        else
          iTmp = pMassMapLinks%iAveragingType(iLink)
        endif
        
        fTmp = 1.
        if (ifTimestepIntegrated) then ! Emission comes as mass per cell per timestep
                                       ! Have to scale it right
           fTmp = abs(1./fu_sec(timestepIn))
        endif

        select case(iTmp)
          case(iAsIs)                        ! take the last value
            if(.not. now == writingTime)cycle
            pMassMapLinks%IntegrationStart(iLink) = now
            call merge_instant_MMdata(pMassMapLinks%pMMIn(iLink)%ptrMassMap, &   ! input, instant
                                    & pMassMapLinks%pMMOut(iLink), &  ! output, meanXXHours
                                    & one_second, &                              ! no integration
                                    & pMassMapLinks%iVerticalTreatment(iLink), & ! what to do with vertical
                                    & pMassMapLinks%adaptor(iLink), &            ! link input->output species
                                    & ifPerDensity, ifPerDzDynamic, &
                                    & data_buf%weight_past)
          case(iInstant)                        ! for instant quantity, take the last value
            if(.not. now == writingTime)cycle
            pMassMapLinks%IntegrationStart(iLink) = now
            call merge_instant_MMdata(pMassMapLinks%pMMIn(iLink)%ptrMassMap, &   ! input, instant
                                    & pMassMapLinks%pMMOut(iLink), &  ! output, meanXXHours
                                    & one_second * fTmp, &                      ! no integration
                                    & pMassMapLinks%iVerticalTreatment(iLink), & ! what to do with vertical
                                    & pMassMapLinks%adaptor(iLink), &            ! link input->output species
                                    & ifPerDensity, ifPerDzDynamic, &
                                    & data_buf%weight_past)
          case(iMeanLastHrs, iAverage, iCumulative)    ! sum-up over the averaging period
            if(pMassMapLinks%iAveragingType(iLink) == iMeanLastHrs)then
              if(fu_abs(now - writingTime) > pMassMapLinks%AveragingPeriod(iLink))cycle
            endif
            call merge_instant_MMdata(pMassMapLinks%pMMIn(iLink)%ptrMassMap, &     ! input, instant
                                    & pMassMapLinks%pMMOut(iLink), &    ! output, meanXXHours
                                    & fu_abs(timestepIn) * fTmp, &                        ! model time step
                                    & pMassMapLinks%iVerticalTreatment(iLink), &   ! what to do with vertical
                                    & pMassMapLinks%adaptor(iLink), &              ! link input->output species
                                    & ifPerDensity, ifPerDzDynamic, &
                                    & data_buf%weight_past)
 !           call msg("Merging with timestep:"+fu_interval_string(timestepIn))
          case default
            call set_error('Only Average,cumulative,mean-XX-hours allowed, not:' + &
                         & fu_str(pMassMapLinks%iAveragingType(iLink)),'merge_mass_map_2_mass_map')
            return
        end select   ! iAveraginType
      endif  ! if cumulative quantity

      !
      ! At the output time moment, we should do the static part of the normalization:
      ! dx, dy, and, possibly, dz. Also, unless the stuff is cumulative, divide with 
      ! averaging period
      !
      if(now == writingTime)then
        if(pMassMapLinks%iAveragingType(iLink) == iCumulative)then
          call prepare_mass_map_writing(pMassMapLinks%pMMOut(iLink), &   ! map to normalise
                                      & interval_missing, &                         ! nothing for cumulatives
                                      & ifPerArea, &                                ! if dx,dy scalng
                                      & ifPerDz .and. .not. ifPerDzDynamic, &       ! if static dz scaling
                                      & dz_past_3d, xSize, ySize)
        else
          call prepare_mass_map_writing(pMassMapLinks%pMMOut(iLink), &   ! map to normalise
                                      & fu_abs(now - pMassMapLinks%IntegrationStart(iLink)), & ! time period to apply
                                      & ifPerArea, &                                 ! if dx,dy scalng
                                      & ifPerDz .and. .not. ifPerDzDynamic, &        ! if static dz scaling
                                      & dz_past_3d, xSize, ySize)
!            call msg("Finalizing with time interval:"+fu_interval_string(fu_abs(now - pMassMapLinks%IntegrationStart(iLink))))
        endif
      endif

    end do  ! iLink

    contains
    
      !===============================================================================

      subroutine merge_instant_MMdata(pMMIn, pMMOut, timestep, iVertTreat, adaptor, &
                                    & ifPerDensity, ifPerDzDyn, weight_past)
        !
        ! Add-up an instant mass map to the cumulative/average one. Integration goes over
        ! seconds, if output is to be written at this time moment, normalization with the 
        ! total cumulated time takes place
        !
        implicit none

        ! Imported parameters
        type(TMass_map), intent(in) :: pMMIn
        type(TMass_map), intent(inout) ::  pMMOut
        type(chemical_adaptor), intent(in) :: adaptor
        type(silja_interval), intent(in) :: timestep
        integer, intent(in) :: iVertTreat
        logical, intent(in) :: ifPerDensity, ifPerDzDyn
        real, intent(in) :: weight_past

        ! Local variables
        integer :: ix, iy, iz, iSpIn, iSrc, i1d, nz
        real :: seconds, dz, cnc_air
      
        seconds = fu_sec(timestep)
        cnc_air = 1.0
        dz = 1.0

        select case(iVertTreat)

          case(do_nothing_flag, lowest_level_flag)  ! merge the whole 5D or 4D cube
            if(iVertTreat == do_nothing_flag)then
              nz = pMMIn%n3d
            else
              nz = 1
            endif
            do iy = 1, pMMIn%ny
              do ix = 1, pMMIn%nx
                i1d = ix+(iy-1)*nx_dispersion
                do iz = 1, nz
                  if(ifPerDzDyn) dz = dz_past_3d%p2d(iz)%ptr(i1d)*weight_past + &
                                    & dz_future_3d%p2d(iz)%ptr(i1d)*(1.-weight_past)
                  if(ifPerDensity) cnc_air = (weight_past * rho_past_3d%p2d(iz)%ptr(i1d) + &
                                              & (1.-weight_past) * rho_future_3d%p2d(iz)%ptr(i1d)) / &
                                           & molecular_weight_air                  
                  do iSrc = 1, pMMIn%nSrc
                    do iSpIn = 1, adaptor%nSpecies
                      if(adaptor%iSp(iSpIn) > 0) pMMOut%arM(iSpIn,iSrc,iz,ix,iy) = &
                                                  & pMMOut%arM(iSpIn,iSrc,iz,ix,iy) + &
                                                  & pMMIn%arM(adaptor%iSp(iSpIn),iSrc,iz,ix,iy) * &
                                                  & seconds / (dz * cnc_air)
                    end do      ! iSpIn
                  end do     ! iSrc
                end do    ! iz
              end do    ! ix
            end do   ! iy

          case(integrate_column_flag)   ! sum over the vertical, cannot be for VMR-quantity
            do iy = 1, pMMIn%ny
              do ix = 1, pMMIn%nx
                i1d = ix+(iy-1)*nx_dispersion
                do iz = 1, pMMIn%n3d
                  if(ifPerDzDyn) dz = dz_past_3d%p2d(iz)%ptr(i1d)*weight_past + &
                                    & dz_future_3d%p2d(iz)%ptr(i1d)*(1.-weight_past)
                  if(ifPerDensity) cnc_air = (weight_past * rho_past_3d%p2d(iz)%ptr(i1d) + &
                                     & (1.-weight_past) * rho_future_3d%p2d(iz)%ptr(i1d)) / &
                                     & molecular_weight_air                  
                  do iSrc = 1, pMMIn%nSrc
                    do iSpIn = 1, adaptor%nSpecies
                      if(adaptor%iSp(iSpIn) > 0) pMMOut%arM(iSpIn,iSrc,1,ix,iy) = &
                                                   & pMMOut%arM(iSpIn,iSrc,1,ix,iy) + &
                                                   & pMMIn%arM(adaptor%iSp(iSpIn),iSrc,iz,ix,iy) * &
                                                   & seconds / dz
                    end do  ! speciesIn
                  end do  ! iSrc
                end do  ! iz
              end do  ! ix
            end do  ! iy

          case default
            call set_error('Strange vertical treatment:' + fu_str(iVertTreat),'merge_instant_MMdata')
            return
        end select

      end subroutine merge_instant_MMdata

  end subroutine merge_mass_map_2_mass_map


  !********************************************************************************
  
  subroutine merge_lagrPartSet_2_mass_map(pLPSet2MassMapLinks, &      ! Structure to merge
                                        & data_buf, writingTime, now, timestepIn)
    !
    ! Merges the Lagrange Particle set to the cumulating mass map. Uses all species addressed 
    ! by the adaptor, all levels or their sum or the fiels one (depending on vertical treatment). 
    ! Normalization allowed is per-volume and per-air-denstiy
    ! Does the same thing as merge_vector_data_to_stack below but for lpSet and mass_map.
    ! ATTENTION.
    ! This sub serves ONLY merging that requires averaging, i.e. some action that calls for 
    ! intermediate mass map. All simpler (or more complicated) requests are handled directly 
    ! by other lpSet to stack subroutines.
    !
    implicit none

    ! Imported parameters
    type(TLagrangeToMassMapLink), intent(inout) :: pLPSet2MassMapLinks  ! Masses to merge, metadata
    type(TField_buffer), intent(in) :: data_buf
    type(silja_time), intent(in) :: writingTime, now
    type(silja_interval), intent(in) :: timestepIn
    !
    ! Local variables
    !
    integer :: iLink
    logical :: ifPerDensity, ifPerArea, ifPerDz, ifPerDzDynamic, ifTimestepIntegrated
    type(field_3d_data_ptr), pointer :: rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d
    real, dimension(:), pointer :: xSize, ySize
    type(THorizInterpStruct), pointer :: pHorizIS
    type(TVertInterpStruct), pointer :: pVertIS
    type(silam_vertical) :: vertOfLayers
    type(silja_grid), save :: gridMeteoHighRes
    integer, save :: iDimScaling
    logical, save :: ifFirst = .true.

    if(pLPSet2MassMapLinks%iVerticalTreatment(1) == int_missing)return  ! nothing to do

    if(ifFirst)then
      !
      ! Prepare the single-layer z-type vertical and high-resolution meteo-type grid. 
      ! LP indices in that grid will be used instead of actual coordinates - that way, 
      ! we shall be able to involve the interpolation structures
      !
      call set_grid_from_lores_templ(pLPSet2MassMapLinks%pLPSetIn(1)%ptrLpSet%gridTemplate, &  ! grid to follow
                                 & pLPSet2MassMapLinks%pMMOut(1)%gridTemplate, & ! resolution to beat
                                 & 2.0, &         ! factor of resolution improvement
                                 & gridMeteoHighRes, iDimScaling) ! output grid and scaling
      ifFirst = .false.
    endif

    !
    ! Grand cycle through the links
    !
    do iLink = 1, nLagrPart2MassMapLinks

      if(pLPSet2MassMapLinks%iVerticalTreatment(iLink) == int_missing)exit ! all done

      if(.not. pLPSet2MassMapLinks%pMMOut(iLink)%gridTemplate == &
           &   pLPSet2MassMapLinks%pMMOut(1)%gridTemplate)then
        call msg_warning('Stragely, output grids differ between the links:','merge_lagrPartSet_2_mass_map')
        call report(pLPSet2MassMapLinks%pMMOut(iLink)%gridTemplate)
        call report(pLPSet2MassMapLinks%pMMOut(1)%gridTemplate)
        call set_error('Stragely, output grids differ between the links:','merge_lagrPartSet_2_mass_map')
      endif
      !
      ! Some quantities need to be scaled from mass map for the output. Since scaling
      ! can be time-dependent (air density or dz), have to do it every time step
      !
      call get_scaling(pLPSet2MassMapLinks%pMMOut(iLink)%quantity, &
                     & pLPSet2MassMapLinks%iVerticalTreatment(iLink), &
                     & ifPerDensity, ifPerArea, ifPerDz, ifTimestepIntegrated, &
                     & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize, &
                     & data_buf)
      if(error)return
      !
      ! Get the interpolation structures. 
      ! ATTENTION. 
      ! We need here a tricky structure: the one, which will allow "forward" projecting
      ! from the given LP coordinates to the output grid and vertical. Existing structures 
      ! work "backward" providing the indices and values for the given receptor points.
      ! We shall use nearest-point interpolation and reverse the direction.
      !
      call make_vertical_of_layers(pLPSet2MassMapLinks%pLPSetIn(iLink)%ptrLpSet%verticalTemplate, &
                                 & vertOfLayers)
      pVertIS => fu_vertical_interp_struct(pLPSet2MassMapLinks%pMMOut(iLink)%vertTemplate,&
                                         & vertOfLayers, &
                                         & pLPSet2MassMapLinks%pLPSetIn(iLink)%ptrLpSet%gridTemplate, &
!                                         & summation, &
                                         & nearest_point, &
                                         & one_hour, 'output_to_meteo_nearest_point')
      if(error)return
      if(fu_fails(associated(pVertIS),'Failed to get vertical Lagrangian interpolation structure', &
                                    & 'merge_lagrPartSet_2_mass_map'))return

      pHorizIS => fu_horiz_interp_struct(pLPSet2MassMapLinks%pMMOut(iLink)%gridTemplate, &
                                       & gridMeteoHighRes, &
!                                       & summation, &
                                       & nearest_point, .true., &
                                       & 5)
      if(error)return
      if(fu_fails(associated(pHorizIS),'Failed to get horizontal Lagrangian interpolation structure', &
                                     & 'merge_lagrPartSet_2_mass_map'))return
      !
      ! Can we do it once at the writing moment or have to waste time at every collection?
      ! Basically, simple: we can ask for meteo dependence of projection from the MMOut vertical
      ! to metric height.
      !
      if(ifPerDz)then
        ifPerDzDynamic = fu_if_level_meteo_dependent(fu_leveltype( &
                                      & pLPSet2MassMapLinks%pMMOut(iLink)%vertTemplate))
      else
        ifPerDzDynamic = .false.
      endif
      !
      ! Input mass map is always instant, output can be as-is, instant, averaged, cumulative, mean-over-some-time
      !
      select case(pLPSet2MassMapLinks%iAveragingType(iLink))
        !FIXME make emission working properly 

        case(iAsIs, iInstant)
          if(.not. (now == writingTime))cycle ! not the output time
          pLPSet2MassMapLinks%IntegrationStart(iLink) = now
          call merge_lp2MMdat(pLPSet2MassMapLinks%pLpSetIn(iLink)%ptrLpSet, &  ! input, instant
                            & pLPSet2MassMapLinks%iLPSetType(iLink), &       ! input
                            & pLPSet2MassMapLinks%pMMOut(iLink), &  ! output, meanXXHours
                            & one_second, &                                    ! averaging period
                            & pLPSet2MassMapLinks%iVerticalTreatment(iLink), & ! what to do with vertical
                            & pLPSet2MassMapLinks%adaptor(iLink), &  ! between input and output species
                            & ifPerDensity, ifPerDzDynamic, &
                            & pHorizIS, pVertIS, iDimScaling, &
                            & data_buf%weight_past)

        case(iMeanLastHrs, iAverage, iCumulative)

          if(pLPSet2MassMapLinks%iAveragingType(iLink) == iMeanLastHrs)then
            if(fu_abs(now - writingTime) > pLPSet2MassMapLinks%AveragingPeriod(iLink))cycle
          endif

          call merge_lp2MMdat(pLPSet2MassMapLinks%pLpSetIn(iLink)%ptrLpSet, &    ! input, instant
                            & pLPSet2MassMapLinks%iLPSetType(iLink), &       ! input
                            & pLPSet2MassMapLinks%pMMOut(iLink), &   ! output, meanXXHours
                            & timestepIn, &                           ! integrate over time step
                            & pLPSet2MassMapLinks%iVerticalTreatment(iLink), &    ! what to do with vertical
                            & pLPSet2MassMapLinks%adaptor(iLink), &  ! between input and output species
                            & ifPerDensity, ifPerDzDynamic, &
                            & pHorizIS, pVertIS, iDimScaling, &
                            & data_buf%weight_past)

        case default
          call set_error('Unsupported averaging type:' + &
                       & fu_str(pLPSet2MassMapLinks%iAveragingType(iLink)), &
                       & 'merge_inst_mass_map_2_mass_map')
          return
      end select   ! iAveraginType

      !
      ! At the output time moment, we should do the static part of the normalization:
      ! dx, dy, and, possibly, dz. Also, unless the stuff is cumulative, divide with 
      ! averaging period
      !
      if(writingTime == now)then
        if(pLPSet2MassMapLinks%iAveragingType(iLink) == iCumulative)then
          call prepare_mass_map_writing(pLPSet2MassMapLinks%pMMOut(iLink), &   ! map to normalise
                                      & interval_missing, &                    ! nothing for cumulatives
                                      & ifPerArea, &                           ! if dx,dy scalng
                                      & ifPerDz .and. .not. ifPerDzDynamic, &  ! if static dz scaling
                                      & dz_past_3d, xSize, ySize)
        else
          call prepare_mass_map_writing(pLPSet2MassMapLinks%pMMOut(iLink), &   ! map to normalise
                                      & now - pLPSet2MassMapLinks%IntegrationStart(iLink), & ! time period 
                                      & ifPerArea, &                                 ! if dx,dy scalng
                                      & ifPerDz .and. .not. ifPerDzDynamic, &        ! if static dz scaling
                                      & dz_past_3d, xSize, ySize)
        endif
      endif
      if(error)return
    end do  ! iLink

    contains
    
      !===============================================================================

      subroutine merge_lp2MMdat(pLpSetIn, iLPSetType, pMMOut, timestep, iVertTreat, &
                              & adaptor, &
                              & ifPerDensity, ifPerDzDyn, &
                              & pHorizIS, pVertIS, iDimScaling, &
                              & weight_past)
        !
        ! Add-up an instant mass map to the cumulative/average one. Integration goes over
        ! seconds, if output is to be written at this time moment, normalization with the 
        ! total cumulated time takes place
        !
        implicit none

        ! Imported parameters
        type(Tlagrange_particles_set), pointer :: pLpSetIn
        integer, intent(in) :: iLPSetType
        type(TMass_map), intent(inout) ::  pMMOut
        type(chemical_adaptor), intent(in) :: adaptor
        type(silja_interval), intent(in) :: timestep
        integer, intent(in) :: iVertTreat, iDimScaling
        logical, intent(in) :: ifPerDensity, ifPerDzDyn
        type(THorizInterpStruct), pointer :: pHorizIS
        type(TVertInterpStruct), pointer :: pVertIS
        real, intent(in) :: weight_past

        ! Local variables
        integer :: iPart, ixHiRes, iyHiRes, ixMeteo, iyMeteo, izMeteo, ixOut, iyOut, izOut, &
                 & iSpIn, iSrc, i1dMeteo, i1dOut, nxOut, iHorCoef, iVertCoef, nHorCoefs, nVertCoefs, &
                 & nActualCoefs
        real :: seconds, dz, cnc_air, fScale
        type(TVertInterpCells) :: coefsVertIS
        type(THorizInterpCells) :: coefsHorizIS
        real, dimension(:,:), pointer :: pLPMassIn
      
        seconds = fu_sec(timestep)
        cnc_air = 1.0
        dz = 1.0
        nxOut = pMMOut%nx
        nHorCoefs = fu_nCoefs(pHorizIS)
        nVertCoefs = fu_nCoefs(pVertIS)
        select case(iLPSetType)
          case(1)
            pLPMassIn => pLPSetIn%lpMassTrn
          case(2)
            pLPMassIn => pLPSetIn%lpMassSL
          case(3)
            pLPMassIn => pLPSetIn%lpMassAer
          case default
            call set_error('Unknown LPSetType:' + fu_str(iLPSetType),'merge_lp2MMdat')
            return
        end select

        select case(iVertTreat)

          case(do_nothing_flag)  ! merge the whole 5D cube

            call get_coefs(pHorizIS,coefsHorizIS)
            call get_coefs(pVertIS, coefsVertIS)
            if(error)return

            do iPart = 1, pLpSetIn%nop
              !
              ! Indices of the partiles in the high-resolution grid.
              ! The high-resolution grid is:
              ! - along x axis, it has nx_meteo * xDimScaling cells
              ! - along y axis, it has ny_meteo * yDimScaling cells
              ! - along z axis the index has not changed: nz_meteo cells
              !
              if(pLpSetIn%lpStatus(iPart) == int_missing)cycle
              ixHiRes = min(nx_meteo*iDimScaling,max(1, &
                     & nint((pLpSetIn%lpDyn(lp_x,iPart) - 0.5) * iDimScaling + 0.5)))
              iyHiRes = min(ny_meteo*iDimScaling,max(1, &
                     & nint((pLpSetIn%lpDyn(lp_y,iPart) - 0.5) * iDimScaling + 0.5)))

              ixMeteo = min(nx_meteo,max(1,nint(pLpSetIn%lpDyn(lp_x,iPart))))
              iyMeteo = min(ny_meteo,max(1,nint(pLpSetIn%lpDyn(lp_y,iPart))))
              izMeteo = min(nz_meteo,max(1,nint(pLpSetIn%lpDyn(lp_z,iPart))))   ! relative, meteo grid
              i1dMeteo = ixMeteo + nx_meteo * (iyMeteo -1)

              !
              ! Find the right vertical layer for the particle. Note that we have
              ! only one level: the method is nearest-point,
              ! vertical structure is in meteo grid, horizontal one if from output grid to 
              ! high-resolution one.
              ! Watchout the trick: e.g. horizontal hires grid is gridTo but 
              ! the nearest_point method gives 1:1 connection between output and hires
              ! grids. Then, knowing hires cell we know the output cell that is the closest.
              ! The other option is to use summation, which gives several target cells: lp is split
              !
!call msg('Meteo hires (' + fu_str(ixHiRes) + ',' + fu_str(iyHiRes) + ',' + &
!                                                     & fu_str(izMeteo) + ') -> output (' + &
!                         & fu_str(coefsHorizIS%indX(1,ixHiRes,iyHiRes)) + ',' + &
!                         & fu_str(coefsHorizIS%indY(1,ixHiRes,iyHiRes)) + ',' + &
!                         & fu_str(coefsVertIS%indLev(1,ixMeteo,iyMeteo,izMeteo)) + ')')
!              ixOut = coefsHorizIS%indX(1,ixHiRes,iyHiRes)
!              iyOut = coefsHorizIS%indY(1,ixHiRes,iyHiRes)
!              izOut = coefsVertIS%indLev(1,ixMeteo,iyMeteo,izMeteo)
              !
              ! Here we accept some uncertainty. The interpolation structure is not
              ! symmetrical: weight, in case of summation, is the fraction of From-cell that is
              ! to be transferred to To-cell. Whereas we need the other way round (remember,
              ! the grids are already swapped to catch-up with inverse-projection principle).
              ! Solution so far is:
              ! - horizontal reprojection turned into nearest-point. For nearly-homogeneous grids OK
              ! - vertical reprojection assumes that the particle is split to as many equal pieces
              !   as there are grid cells to serve. Crude but not very frequent case. 
              !
!              do iHorCoef = 1, nHorCoefs
!                ixOut = coefsHorizIS%indX(iHorCoef,ixHiRes,iyHiRes)
!                if(ixOut == 0)exit
!                iyOut = coefsHorizIS%indY(iHorCoef,ixHiRes,iyHiRes)
                ixOut = coefsHorizIS%indX(1,ixHiRes,iyHiRes)
                iyOut = coefsHorizIS%indY(1,ixHiRes,iyHiRes)
                i1dOut = ixOut + nxOut * (iyOut - 1)
 
                nActualCoefs = count(coefsVertIS%indLev(1:nVertCoefs,ixMeteo,iyMeteo,izMeteo)>0)
                fScale = 1. / real(nActualCoefs)
 
                do iVertCoef = 1, nActualCoefs !nVertCoefs
                  izOut = coefsVertIS%indLev(iVertCoef,ixMeteo,iyMeteo,izMeteo)
!                  if(izOut == 0)exit
                  iSrc = mod(pLpSetIn%lpStatus(iPart),100)

                  if(ifPerDzDyn) dz = dz_past_3d%p2d(izOut)%ptr(i1dOut)*weight_past + &
                                    & dz_future_3d%p2d(izOut)%ptr(i1dOut)*(1.-weight_past)
                  if(ifPerDensity) &
                         & cnc_air = (weight_past * rho_past_3d%p2d(izMeteo)%ptr(i1dMeteo) + &
                                    & (1.-weight_past) * rho_future_3d%p2d(izMeteo)%ptr(i1dMeteo))  &
                                   &  / molecular_weight_air                  
                  do iSpIn = 1, adaptor%nSpecies
                    if(adaptor%iSp(iSpIn) > 0) &
                        & pMMOut%arM(iSpIn,iSrc,izOut,ixOut,iyOut) = &
                                         & pMMOut%arM(iSpIn,iSrc,izOut,ixOut,iyOut) + &
                                         & pLPMassIn(adaptor%iSp(iSpIn),iPart) * &
                                         & seconds / (dz * cnc_air) * fScale
!                                         & coefsHorizIS%weight(iHorCoef,ixHiRes,iyHiRes) * &
!                                         & coefsVertIS%weight(iVertCoef,ixMeteo,iyMeteo,izMeteo)
                  end do      ! iSpIn
                end do  ! iVertCoefs
!              end do  ! iHorCoef
            end do   ! iPart

          case(integrate_column_flag)   ! sum over the vertical, cannot be for VMR-quantity

            call get_coefs(pHorizIS,coefsHorizIS)
            call get_coefs(pVertIS, coefsVertIS)
            if(error)return

            do iPart = 1, pLpSetIn%nop
              !
              ! Indices of the partiles in the high-resolution grid.
              ! The high-resolution grid is:
              ! - along x axis, it has nx_meteo * xDimScaling cells
              ! - along y axis, it has ny_meteo * yDimScaling cells
              ! - along z axis the index has not changed: nz_meteo cells
              !
              ixHiRes = nint((pLpSetIn%lpDyn(lp_x,iPart) - 0.5) * iDimScaling + 0.5)
              iyHiRes = nint((pLpSetIn%lpDyn(lp_y,iPart) - 0.5) * iDimScaling + 0.5)

              ixMeteo = nint(pLpSetIn%lpDyn(lp_x,iPart))
              iyMeteo = nint(pLpSetIn%lpDyn(lp_y,iPart))
              !
              ! Find the right vertical layer for the particle. Note that we have
              ! only one level: the method is nearest-point,
              ! vertical structure is in meteo grid, horizontal one if from output grid to 
              ! high-resolution one.
              ! Watchout the trick: e.g. horizontal hires grid is gridTo but since we
              ! use nearest_point method, we get 1:1 connection between output and hires
              ! grids. Then, knowing hires cell we know the output cell that is the closest.
              !
! For summation
!              do iHorCoef = 1, nHorCoefs
!                ixOut = coefsHorizIS%indX(iHorCoef,ixHiRes,iyHiRes)
!                if(ixOut == 0)exit
!                iyOut = coefsHorizIS%indY(iHorCoef,ixHiRes,iyHiRes)
! For nearest-point
                ixOut = coefsHorizIS%indX(1,ixHiRes,iyHiRes)
                iyOut = coefsHorizIS%indY(1,ixHiRes,iyHiRes)
                iSrc = mod(pLpSetIn%lpStatus(iPart),100)

                if(ifPerDzDyn)then
                  izMeteo = nint(pLpSetIn%lpDyn(lp_z,iPart))
                  izOut = coefsVertIS%indLev(1,ixMeteo,iyMeteo,izMeteo)
                  i1dOut = ixOut + nxOut * (iyOut -1)
                  dz = dz_past_3d%p2d(izOut)%ptr(i1dOut)*weight_past + &
                     & dz_future_3d%p2d(izOut)%ptr(i1dOut)*(1.-weight_past)
                endif
                if(ifPerDensity) &
                         & cnc_air = (weight_past * rho_past_3d%p2d(izMeteo)%ptr(i1dMeteo) + &
                                    & (1.-weight_past) * rho_future_3d%p2d(izMeteo)%ptr(i1dMeteo))  &
                                   &  / molecular_weight_air                  
                do iSpIn = 1, adaptor%nSpecies
                  if(adaptor%iSp(iSpIn) > 0) &
                           & pMMOut%arM(iSpIn,iSrc,1,ixOut,iyOut) = &
                                           & pMMOut%arM(iSpIn,iSrc,1,ixOut,iyOut) + &
                                           & pLPMassIn(adaptor%iSp(iSpIn),iPart) * &
                                           & seconds / (dz * cnc_air) !* &
!                                           & coefsHorizIS%weight(iHorCoef,ixHiRes,iyHiRes)
                end do      ! iSpIn
!              end do
            end do   ! iPart

          case(lowest_level_flag)     ! pick the first level

            call get_coefs(pHorizIS,coefsHorizIS)
            call get_coefs(pVertIS, coefsVertIS)
            if(error)return

            do iPart = 1, pLpSetIn%nop
              !
              ! Indices of the partiles in the high-resolution grid.
              ! The high-resolution grid is:
              ! - along x axis, it has nx_meteo * xDimScaling cells
              ! - along y axis, it has ny_meteo * yDimScaling cells
              ! - along z axis the index has not changed: nz_meteo cells
              !

              ixMeteo = nint(pLpSetIn%lpDyn(lp_x,iPart))
              iyMeteo = nint(pLpSetIn%lpDyn(lp_y,iPart))
              izMeteo = nint(pLpSetIn%lpDyn(lp_z,iPart))
              i1dMeteo = ixMeteo + nx_meteo * (iyMeteo -1)
              !
              ! Find the right vertical layer for the particle. Note that we have
              ! only one level: the method is nearest-point,
              ! vertical structure is in meteo grid, horizontal one if from output grid to 
              ! high-resolution one.
              ! Watchout the trick: e.g. horizontal hires grid is gridTo but since we
              ! use nearest_point method, we get 1:1 connection between output and hires
              ! grids. Then, knowing hires cell we know the output cell that is the closest.
              !
              izOut = 1   ! lowest level... not pHorizIS%indLev(1,ixMeteo,iyMeteo,izMeteo)
              if(izOut > 1)cycle

              ixHiRes = nint((pLpSetIn%lpDyn(lp_x,iPart) - 0.5) * iDimScaling + 0.5)
              iyHiRes = nint((pLpSetIn%lpDyn(lp_y,iPart) - 0.5) * iDimScaling + 0.5)
! For summation
!              do iHorCoef = 1, nHorCoefs
!                ixOut = coefsHorizIS%indX(iHorCoef,ixHiRes,iyHiRes)
!                if(ixOut == 0)exit
!                iyOut = coefsHorizIS%indY(iHorCoef,ixHiRes,iyHiRes)
                ixOut = coefsHorizIS%indX(1,ixHiRes,iyHiRes)
                iyOut = coefsHorizIS%indY(1,ixHiRes,iyHiRes)
                iSrc = mod(pLpSetIn%lpStatus(iPart),100)

                nActualCoefs = count(coefsVertIS%indLev(1:nVertCoefs,ixMeteo,iyMeteo,izMeteo)>0)
                fScale = 1. / real(nActualCoefs)
 
                do iVertCoef = 1, nActualCoefs !nVertCoefs
                  izOut = coefsVertIS%indLev(iVertCoef,ixMeteo,iyMeteo,izMeteo)
                  if(izOut > 1)cycle
!                  if(izOut == 0)exit

                  if(ifPerDzDyn)then
                    i1dOut = ixOut + nxOut * (iyOut -1)
                    dz = dz_past_3d%p2d(1)%ptr(i1dOut)*weight_past + &
                       & dz_future_3d%p2d(1)%ptr(i1dOut)*(1.-weight_past)
                  endif
                  if(ifPerDensity) &
                           & cnc_air = (weight_past * rho_past_3d%p2d(izMeteo)%ptr(i1dMeteo) + &
                                      & (1.-weight_past) * rho_future_3d%p2d(izMeteo)%ptr(i1dMeteo))  &
                                     &  / molecular_weight_air                  
                  do iSpIn = 1, adaptor%nSpecies
                    if(adaptor%iSp(iSpIn) > 0) pMMOut%arM(iSpIn,iSrc,1,ixOut,iyOut) = &
                                             & pMMOut%arM(iSpIn,iSrc,1,ixOut,iyOut) + &
                                             & pLPMassIn(adaptor%iSp(iSpIn),iPart) * &
                                             & seconds / (dz * cnc_air) * fScale
!                                             & coefsHorizIS%weight(iHorCoef,ixHiRes,iyHiRes)
                  end do      ! iSpIn
                end do    ! iVertCoefs
!              end do      ! iHorCoef
            end do   ! iPart

          case default
            call set_error('Strange vertical treatment:' + fu_str(iVertTreat),'merge_lp2MMdat')
            return
        end select

      end subroutine merge_lp2MMdat

  end subroutine merge_lagrPartSet_2_mass_map  


  !************************************************************************************

  subroutine merge_cumul_MMdata(pMMIn, pMMOut, IntegrPeriod, iVertTreat, adaptor, &
                              & ifPerDensity, ifPerDzDyn, weight_past, ifStartAveragingPeriod, &
                              & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize)
    !
    ! Add-up a cumulative mass map to the cumulative/average one. Integration goes over
    ! seconds, if output is to be written at this time moment, normalization with the 
    ! total cumulated time takes place
    !
    implicit none

    ! Imported parameters
    type(TMass_map), intent(in) :: pMMIn
    type(TMass_map), intent(inout) ::  pMMOut
    type(chemical_adaptor), intent(in) :: adaptor
    type(silja_interval), intent(in) :: IntegrPeriod
    integer, intent(in) :: iVertTreat
    logical, intent(in) :: ifPerDensity, ifPerDzDyn, ifStartAveragingPeriod
    real, intent(in) :: weight_past
    type(field_3d_data_ptr), pointer :: rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d
    real, dimension(:), pointer :: xSize, ySize

    ! Local variables
    integer :: ix, iy, iz, iSpIn, iSrc, i1d, nz
    real :: seconds_integr, dz, cnc_air
      
    seconds_integr = fu_sec(IntegrPeriod)
    cnc_air = 1.0
    dz = 1.0

    select case(iVertTreat)

      case(do_nothing_flag, lowest_level_flag)  ! merge the whole 5D or 4D cube
        if(iVertTreat == do_nothing_flag)then
          nz = pMMIn%n3d
        else
          nz = 1
        endif
        do iy = 1, pMMIn%ny
          do ix = 1, pMMIn%nx
            i1d = ix+(iy-1)*nx_dispersion
            do iz = 1, nz
              if(ifPerDzDyn) dz = dz_past_3d%p2d(iz)%ptr(i1d)*weight_past + &
                                & dz_future_3d%p2d(iz)%ptr(i1d)*(1.-weight_past)
              if(ifPerDensity) cnc_air =  (weight_past * rho_past_3d%p2d(iz)%ptr(i1d) + &
                                          & (1.-weight_past) * rho_future_3d%p2d(iz)%ptr(i1d)) / &
                                       & molecular_weight_air                  
              do iSrc = 1, pMMIn%nSrc
                do iSpIn = 1, adaptor%nSpecies
                  if(adaptor%iSp(iSpIn) > 0)then
                    if(ifStartAveragingPeriod)then
                      ! Store the current stage of the data
                      pMMOut%arM(iSpIn,iSrc,iz,ix,iy) = &
                                      & pMMIn%arM(adaptor%iSp(iSpIn),iSrc,iz,ix,iy) / (dz * cnc_air)
                    else
                      ! Subtract the stored data from the new piece. Scaling with averaging period is later
                      pMMOut%arM(iSpIn,iSrc,iz,ix,iy) = &
                                  & (pMMIn%arM(adaptor%iSp(iSpIn),iSrc,iz,ix,iy) / (dz * cnc_air) - &
                                      & pMMOut%arM(iSpIn,iSrc,iz,ix,iy)) !/ (seconds_integr + 0.001)
                    endif  ! ifStartAveragingPeriod
                  endif  ! adaptor species OK
                end do      ! iSpIn
              end do     ! iSrc
            end do    ! iz
          end do    ! ix
        end do   ! iy

      case(integrate_column_flag)   ! sum over the vertical, cannot be for VMR-quantity
        do iy = 1, pMMIn%ny
          do ix = 1, pMMIn%nx
            i1d = ix+(iy-1)*nx_dispersion
            if(ifStartAveragingPeriod) pMMOut%arM(:,:,1,ix,iy) = 0.0  ! prepare the place
            do iz = 1, pMMIn%n3d
              if(ifPerDzDyn) dz = dz_past_3d%p2d(iz)%ptr(i1d)*weight_past + &
                                & dz_future_3d%p2d(iz)%ptr(i1d)*(1.-weight_past)
              if(ifPerDensity) cnc_air = (weight_past * rho_past_3d%p2d(iz)%ptr(i1d) + &
                                 & (1.-weight_past) * rho_future_3d%p2d(iz)%ptr(i1d))  / &
                                 & molecular_weight_air                  
              do iSrc = 1, pMMIn%nSrc
                do iSpIn = 1, adaptor%nSpecies
                  if(adaptor%iSp(iSpIn) > 0)then 
                    if(ifStartAveragingPeriod)then
                      ! Store the current stage of the data
                      pMMOut%arM(iSpIn,iSrc,1,ix,iy) = pMMOut%arM(iSpIn,iSrc,1,ix,iy) + &
                                              & pMMIn%arM(adaptor%iSp(iSpIn),iSrc,iz,ix,iy) / dz
                    else
                       ! Subtract the stored data from the new piece and scale with averaging period
                      pMMOut%arM(iSpIn,iSrc,1,ix,iy) = pMMOut%arM(iSpIn,iSrc,1,ix,iy) - &
                                            & pMMIn%arM(adaptor%iSp(iSpIn),iSrc,iz,ix,iy) / dz
                    endif  ! ifStartAveragingPeriod
                  endif   ! adaptor species OK
                end do  ! speciesIn
              end do  ! iSrc
            end do  ! iz
            if(.not. ifStartAveragingPeriod) pMMOut%arM(:,:,1,ix,iy) = &
                                         & -pMMOut%arM(:,:,1,ix,iy) !/ (seconds_integr + 0.0001)
          end do  ! ix
        end do  ! iy

      case default
        call set_error('Strange vertical treatment:' + fu_str(iVertTreat),'merge_instant_MMdata')
        return
    end select

  end subroutine merge_cumul_MMdata


  !*****************************************************************************

  subroutine prepare_mass_map_writing(pMassMap, &          ! map to normalise
                                    & IntegrationPeriod, & ! time period to apply
                                    & ifPerArea, &         ! if dx,dy scalng
                                    & ifPerDzStatic, &     ! if static dz scaling
                                    & dz_3d, xSize, ySize)
    !
    ! Scales the mass map with the grid cell size. Done just before writing.
    !
    implicit none
 
    ! imported parameters
    type(TMass_map), intent(inout) :: pMassMap
    type(silja_interval), intent(in) :: IntegrationPeriod
    logical, intent(in) :: ifPerArea, ifPerDzStatic
    type(field_3d_data_ptr), pointer :: dz_3d
    real, dimension(:), pointer :: xSize, ySize

    ! local variables
    integer :: ix, iy, iz, iSrc, iSp, i1d
    real :: fTimeScale, area, dz

    !
    ! First, get the time scale
    !
    area = 1.0
    dz = 1.0
    fTimeScale = 1.0
    if(defined(IntegrationPeriod))then
      if(fu_abs(IntegrationPeriod) > one_second) fTimeScale = 1. / fu_sec(IntegrationPeriod)
    else
      if(.not. (ifPerDzStatic .or. ifPerArea))return ! no scaling whatsoever
    endif
    !
    ! Now, do the scaling with cell area, layer thickness, and time scale - all only if needed
    !
    do iy = 1, pMassMap%ny
      do ix = 1, pMassMap%nx
        i1d = ix+(iy-1)*pMassMap%nx
        if(ifPerArea) area = xSize(i1d) * ySize(i1d)
        do iz = 1, pMassMap%n3d
          if(ifPerDzStatic) dz = dz_3d%p2d(iz)%ptr(i1d)
          do iSrc = 1, pMassMap%nSrc
            do iSp = 1, pMassMap%nSpecies
              pMassMap%arM(iSp,iSrc,iz,ix,iy) = &
                                    & pMassMap%arM(iSp,iSrc,iz,ix,iy) * fTimeScale / (area * dz)
            end do      ! iSpIn
          end do     ! iSrc
        end do    ! iz
      end do    ! ix
    end do   ! iy

  end subroutine prepare_mass_map_writing


  !********************************************************************************

  subroutine start_new_output_period_massmap(MassMapLinks, &
                                           & data_buf, &          ! data buffer, also has time
                                           & next_output_time, & ! next output time
                                           & now)
    !
    ! Resets the mass map links for the new output period. Mass maps can be instant or 
    ! cumulative from the beginning of the run, so two options are to be considered
    !
    implicit none
    
    type(TMassMapLink), intent(inout) :: MassMapLinks
    type(Tfield_buffer), pointer :: data_buf
    type(silja_time), intent(in) :: next_output_time, now

    ! Local vairables
    integer :: iLink
    logical :: ifPerDz, ifPerArea, ifPerDensity, ifPerDzDynamic, ifTimestepIntegrated
    type(field_3d_data_ptr), pointer :: dz_past_3d, dz_future_3d, rho_past_3d, rho_future_3d
    real, dimension(:), pointer :: xSize, ySize

    !
    ! Somewhat boring stuff: have to all-again get the scaling and pass it to the mising function.
    ! 
    do iLink = 1, nMassMapLinks
      if(MassMapLinks%iVerticalTreatment(iLink) == int_missing)exit
      if(fu_accumulated_quantity(MassMapLinks%pMMIn(iLink)%ptrMassMap%quantity))then
        !
        ! Cumulative input
        !
        select case(MassMapLinks%iAveragingType(iLink))

          case(iAsIs, iCumulative) 

          case(iInstant, iMeanLastHrs)

          case(iAverage)
            call get_scaling(MassMapLinks%pMMOut(iLink)%quantity, &
                           & MassMapLinks%iVerticalTreatment(iLink), &
                           & ifPerDensity, ifPerArea, ifPerDz, ifTimestepIntegrated, &
                           & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize, &
                           & data_buf)
            if(ifPerDz)then
              ifPerDzDynamic = fu_if_level_meteo_dependent(fu_leveltype( &
                                           & MassMapLinks%pMMOut(iLink)%vertTemplate))
            else
              ifPerDzDynamic = .false.
            endif

            MassMapLinks%IntegrationStart(iLink) = now
            call merge_cumul_MMdata(MassMapLinks%pMMIn(iLink)%ptrMassMap, &   ! input, instant
                                  & MassMapLinks%pMMOut(iLink), &  ! output, meanXXHours
                                  & zero_interval, &                           ! integration period
                                  & MassMapLinks%iVerticalTreatment(iLink), & ! what to do with vertical
                                  & MassMapLinks%adaptor(iLink), &            ! link input->output species
                                  & ifPerDensity, ifPerDzDynamic, &
                                  & data_buf%weight_past, &     ! weight_past
                                  & .true., &                     ! Start new period
                                  & rho_past_3d, rho_future_3d, dz_past_3d, dz_future_3d, xSize, ySize)
          case default
            call set_error('Unknown averaging type in massmap links','start_new_output_period_massmap')
            return
        end select
 
      else
        !
        ! Instant input
        !
        select case(MassMapLinks%iAveragingType(iLink))
          case(iAsIs, iInstant, iCumulative) 
          case(iMeanLastHrs)
          case(iAverage)
            MassMapLinks%IntegrationStart(iLink) = now
            MassMapLinks%pMMOut(iLink)%arM = 0.0
          case default
            call set_error('Unknown averaging type in massmap links','start_new_output_period_massmap')
            return
        end select
      endif  ! if accumulation quantity
    end do  ! nMassMapLinks

  end subroutine start_new_output_period_massmap


  !************************************************************************

  subroutine start_new_output_period_lpSet(Lagr2MMLinks, now)
    !
    ! Resetting lagrange particles link structures for the next output period
    !
    implicit none

    type(TLagrangeToMassMapLink), intent(inout) :: Lagr2MMLinks
    type(silja_time), intent(in) :: now

    integer :: i
    !
    ! lpSet is always instant, so resetting the time pointer to the current moment and, in some cases,
    ! mass to zero is enough
    !
    do i = 1, nLagrPart2MassMapLinks
      if(Lagr2MMLinks%iVerticalTreatment(i) == int_missing)exit

      select case(Lagr2MMLinks%iAveragingType(i))
        case(iCumulative) 
        case(iAsIs, iInstant)
          Lagr2MMLinks%IntegrationStart(i) = time_missing
          Lagr2MMLinks%pMMOut(i)%arM = 0.0
        case(iAverage, iMeanLastHrs)
          Lagr2MMLinks%IntegrationStart(i) = now
          Lagr2MMLinks%pMMOut(i)%arM = 0.0
        case default
          call set_error('Unknown averaging type in lpSet-massmap links','start_new_output_period')
          return
      end select
    end do  ! nMassMapLinks

  end subroutine start_new_output_period_lpSet


  !********************************************************************************

  subroutine merge_mass_map_2_stack(pMassMap, &             ! Masses to merge
                                  & iSrc, iSpecies, iz, &   ! what to pick from mass map
                                  & idIn, &                 ! Main parameters of the fields
                                  & ifPerVolume, ifPerDensity, ifPerDz, &  ! normalization
                                  & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                  & weight_past, &
                                  & stackPtr, targetId)     ! Receiving side
    !
    ! Merges the mass map to the output stack. The sub picks out a requested species and
    ! emission source. Normalization allowed is per-volume and per-air-denstiy
    ! Does the same thing as merge_vector_data_to_stack below but for mass_map
    !
    implicit none
    
    ! Imported parameters
    type(TMass_map), pointer :: pMassMap             ! Masses to merge
    integer, intent(in) :: iSrc, iSpecies, iz        ! what to pick from mass map
    type(silja_field_id), intent(in) :: idIn         ! Main parameters of the input fields
    type(silja_field_id), intent(inout) :: targetId   ! Main parameters of the target output fields
    logical, intent(in) :: ifPerVolume, ifPerDensity, ifPerDz ! normalization
    real, dimension(:), pointer :: xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_
    real, intent(in) :: weight_past
    type(silja_stack), pointer :: stackPtr
    !
    ! Local variables
    !
    type(silja_field_id), pointer :: idTmp, idOut, idOut_internal
    type(silja_field), pointer :: field
    real, dimension(:), pointer :: dataOut, dataOut_internal
    integer :: ix, iy, i, timeDir, fs
    logical :: found, found_internal
    real, dimension(:,:,:,:,:), pointer :: arDat_

    call set_error("Obsolete... See code...","merge_mass_map_2_stack")
    !
    ! Can't handle fields with finite validity interval.
    ! 
    !
    return

call report (idin)

    !
    ! First of all, find the receiving fields
    !
    call find_receiving_fields_4_merge(idIn, stackPtr, targetId, &
                                     & found, found_internal, timeDir, &
                                     & idOut, dataOut, idOut_internal, dataOut_internal)
    if(error .or. .not. found)return

    arDat_ => pMassMap%arM
    fs = pMassMap%nx * pMassMap%ny
    !
    ! When the run just starts, there is nothing in the stack, therefore the
    ! validity_length is set to missing. Then the first field that comes puts itself
    ! into the data and overwrites the validity_length.
    ! In principle, all times must be missing but then it is nearly impossible to 
    ! get around such field
    !
    if(.not. defined(fu_validity_length(idOut)))then
      call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                           & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                           & weight_past, iSrc, iSpecies, iz, dataOut, 1.0, 0.0)  ! replace
      call set_validity_length(idOut, fu_validity_length(idIn))
      return
    endif
    !
    ! If valid times are already OK - fields are the same.
    ! The only thing, which might eventually be out of sense is validity length
    !
    if(fu_valid_time(idOut) == fu_valid_time(targetId))then
      call set_validity_length(idOut, fu_validity_length(targetId))
      return
    endif

    !-------------------------------------------------------------------
    !
    ! Having both ids known, let's check what should be done with the data,
    ! depending on the averaging type of idIn, idOut and targetId
    !
    select case(fu_field_kind(targetId))

      case(forecast_flag) 
        !
        ! Target ID is instant
        !
        select case(fu_field_kind(idIn))

          case(forecast_flag, averaged_flag) ! instant + aver -> instant
            if(fu_valid_time(idIn) == fu_valid_time(targetId)) then ! Last time ?
              call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                   & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                   & weight_past, iSrc, iSpecies, iz, dataOut, 1.0, 0.0)  ! replace
              idOut = idIn
            elseif(fu_accumulated_id(idIn))then
              if(fu_between_times(fu_valid_time(targetId), &  ! time
                                & fu_valid_time(idIn), &      ! limit 1
                                & fu_valid_time(idIn) + fu_accumulation_length(idIn), &  ! limit 2
                                & .true.)) then                ! if include borders
                call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                     & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                     & weight_past, iSrc, iSpecies, iz, dataOut, 1.0, 0.0)  ! replace
                idOut = idIn
              endif
            endif

          case(accumulated_flag) ! cumulative -> instant

            if(fu_valid_time(idIn) == fu_valid_time(targetId)) then  ! Last time ?

              if(fu_accumulation_start_time(idIn) == fu_accumulation_start_time(idOut))then
                !
                ! Same start of the fields, their difference gives mean over last timestep
                !
                call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                     & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                     & weight_past, iSrc, iSpecies, iz, dataOut, &
                                     & 1./ abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut))), &
                                     & -1./ abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut))))

!                dataOut(1:fs) = (dataIn(1:fs) - dataOut(1:fs)) / &
!                              & abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut)))

                idOut = targetId

              elseif(fu_accumulation_start_time(idIn) == fu_valid_time(idOut))then  
                !
                ! Last accumulation covers the needs
                !
                call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                     & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                     & weight_past, iSrc, iSpecies, iz, dataOut, &
                                     & 1. / abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut))), &
                                     & 0.0)
!                dataOut(1:fs) = dataIn(1:fs) / abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut)))
                idOut = targetId
              else
                call report_error('Cannot make instant field from accumulated', .false.)
                dataOut(1:fs) = 0.
                call set_valid_time(idOut, fu_valid_time(targetId))
                return
              endif

            else    ! Intermediate collection - just store
              call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                   & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                   & weight_past, iSrc, iSpecies, iz, dataOut, 1., 0.) ! dataOut(1:fs) = dataIn(1:fs)
              idOut = idIn
            endif
            
          case default
            dataOut(1:fs) = 0.
            call set_valid_time(idOut, fu_valid_time(targetId))
            call report_error('Unsupported input field type, instant target', .false.)
            return
        end select

      case(accumulated_flag, averaged_flag) !, accumulated_period_valid_flag) ! Target ID
        !
        ! TargetID is cumulative or averaged.
        ! Cumulative and averaged fields have to be integrated in time
        ! from beginning of the current accumulation period till its end.
        ! Average field will at the end be divided by length of integration.
        !
        select case(fu_field_kind(idIn))

          case (forecast_flag)  
            if(fu_between_times(fu_valid_time(targetId), &  ! time
                              & fu_accumulation_start_time(targetId), &  ! limit 1
                              & fu_valid_time(idIn), .false.) .or. &  ! limit 2, include borders
             & fu_accumulation_start_time(targetId) /= fu_accumulation_start_time(idOut)) then
              dataOut(1:fs) = 0.
              call set_accumulation_length(idOut, zero_interval)
              call set_analysis_time(idOut,fu_accumulation_start_time(targetId))
              call set_valid_time(idOut,fu_accumulation_start_time(targetId))
            endif
            !
            ! The main action happens if and only if the idIn is in-between the target valid_time
            ! and target accumulation start time
            !
            if(fu_between_times(fu_valid_time(idIn), &  ! time
                              & fu_accumulation_start_time(targetId), &  ! limit 1
                              & fu_valid_time(targetId), .true.)) then   ! limit 2, true to include borders
              !
              ! Of course, it is not always correct integration - in reality, the new
              ! field is valid at valid_time, while below it is used from 
              ! last_integrated time till valid_time. 
              ! Ideally, there should be symmetrical integration period, but for that
              ! I would need 3 fields - past, current and future, which is too complicated.
              ! See also io_server, where this error is somewhat coped in routine write_output.
              !
              call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                   & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                   & weight_past, iSrc, iSpecies, iz, dataOut, &
                                   & abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut))), &
                                   & 1.) !      dataOut(1:fs) = dataOut(1:fs) + dataIn(1:fs) * fTmp

              call set_accumulation_length(idOut, fu_valid_time(idIn) - fu_accumulation_start_time(idOut))
              call set_valid_time(idOut, fu_valid_time(idIn))
              call set_validity_length(idOut, fu_validity_length(idIn))
            endif

          case (accumulated_flag) !, accumulated_period_valid_flag)

            if(found_internal)then 
              !
              ! There will be several cases of mutual disposition of accumulation
              ! periods, and they will be considered one-by-one
              !
              call merge_cumulative_arrays(arDat_, fu_sec(fu_accumulation_length(idIn)), &
                                         & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                         & weight_past, iSrc, iSpecies, iz, &
                                         & ifPerVolume, ifPerDensity, ifPerDz, &
                                         & dataOut)
            else
              !
              ! May be, the in-field is just what we expect ? If we want it to be 
              ! just dropped into the stack - the internal field may be not needed
              !
              if(idIn == targetId)then
                idOut = targetId
                call set_validity_length(idOut, fu_validity_length(idIn))
                call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                     & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                     & weight_past, iSrc, iSpecies, iz, dataOut, 1., 0.)  ! replace
                return
              else
                call report_error('Cannot find internal field for merging', .false.)
                return
              endif
            endif     ! .not.found_internal

          case (averaged_flag) !, accumulated_period_valid_flag)
            !
            !
            ! Integrated or averaged - take with care in order not to 
            ! overlap the accumulation periods. Use intermediate field to store
            ! the prior-to-target accumulation.
            !
            if(found_internal)then
              !
              ! Convert average field to cumulative one and do not forget the 
              ! internal field needed to handle 2 cumulatives
              !
              call merge_cumulative_arrays(arDat_, fu_sec(fu_accumulation_length(idIn)), &
                                         & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                         & weight_past, iSrc, iSpecies, iz, &
                                         & ifPerVolume, ifPerDensity, ifPerDz, &
                                         & dataOut)
            else
              !
              ! May be, the in-field is just what we expect ? If we want it to be 
              ! just dropped into the stack - the internal field may be not needed
              !
              if(idIn == targetId)then
                idOut = targetId
                call set_validity_length(idOut, fu_validity_length(idIn))
                call merge_data_arrays(arDat_, ifPerVolume, ifPerDensity, ifPerDz, &
                                     & xSize_, ySize_, dz_past_, dz_future_, rho_past_, rho_future_, &
                                     & weight_past, iSrc, iSpecies, iz , dataOut, 1., 0.) !replace
                return
              else
                call report_error('Cannot find internal field for merging', .false.)
                return
              endif
            endif     ! .not.found_internal

          case default ! idIn
            !
            ! Types like difference_flag, period_averaging_flag
            !
            call set_error('Non-supported averaging type of input field','merge_mass_map_2_stack')
            return
        end select

        if(fu_field_kind(targetId) == averaged_flag)then
          !
          ! For averaged field, if this is the last data collection before the
          ! output - divide the integrated stuff by the accumulation length
          !
          if(fu_valid_time(targetId) == fu_valid_time(idIn))then
            if(fu_accumulation_length(idOut) == zero_interval)then
              dataOut(1:fs) = 0.
            else
              dataOut(1:fs) = dataOut(1:fs) / abs(fu_sec(fu_accumulation_length(idOut)))
            endif
            idOut = targetId
            call set_level(idOut, fu_level(idIn))
            call set_validity_length(idOut, fu_validity_length(idIn))
          endif
        endif

      case default ! target ID
        dataOut(1:fs) = 0.
        call set_valid_time(idOut, fu_valid_time(targetId))
        call set_error('Non-supported averaging type of output field','merge_mass_map_2_stack')
        return
    end select
    
    contains
    
    !========================================================================
    
    subroutine merge_data_arrays(arDat, ifPerVolume, ifPerDensity, ifPerDz, &
                               & xSize, ySize, dz_past, dz_future, rho_past, rho_future, &
                               & weight_past, iSrc, iSpecies, iz, &
                               & datOut, factor_in, factor_out)
      !
      ! Replaces the output data vector with the stuff from mass map, possibly with 
      ! scaling with const factors both in and out arrays, with grid cell volume, and 
      ! with air density.
      ! outNew = outOld * factorOut + in * factorIn (/cellArea) (/cell-thickness) (/air_density)
      ! It serves for replacing data, for adding and subtracting them with arbitrary scaling
      !
      implicit none

      ! Imported parameters
      real, dimension(:,:,:,:,:), pointer :: arDat
      logical, intent(in) :: ifPerVolume, ifPerDensity, ifPerDz
      real, dimension(:), pointer :: datOut
      real, dimension(:), pointer :: xSize, ySize, dz_past, dz_future, rho_past, rho_future
      real, intent(in) :: factor_in, factor_out, weight_past
      integer, intent(in) :: iSrc, iSpecies, iz 

      ! Local variables
      integer :: i1d, ix, iy
      real :: cnc_air

      if(ifPerVolume .and. ifPerDensity)then    ! normalised to cell volume and air density
        if(ifPerDz)then
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              i1d = ix+(iy-1)*nx_dispersion
              cnc_air = (weight_past*rho_past(i1d) + (1.0-weight_past)*rho_future(i1d)) / &
                      & molecular_weight_air
              datOut(i1d) = datOut(i1d) * factor_out + &
                           & arDat(iSpecies,iSrc,iz,ix,iy) * factor_in / &
                           & (cnc_air * xSize(i1d) * ySize(i1d) * &
                            & (dz_past(i1d)*weight_past + dz_future(i1d)*(1.-weight_past)))
            end do
          end do
        else    ! no dz
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              i1d = ix+(iy-1)*nx_dispersion
              cnc_air = (weight_past*rho_past(i1d) + (1.0-weight_past)*rho_future(i1d)) / &
                      & molecular_weight_air
              datOut(i1d) = datOut(i1d) * factor_out + &
                           & arDat(iSpecies,iSrc,iz,ix,iy) * factor_in / &
                           & (cnc_air * xSize(i1d) * ySize(i1d))
            end do
          end do
        endif
      elseif(ifPerVolume)then   ! normalised to cell volume, no air density
        if(ifPerDz)then
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              i1d = ix+(iy-1)*nx_dispersion
              datOut(i1d) = datOut(i1d) * factor_out + &
                          & arDat(iSpecies,iSrc,iz,ix,iy) * factor_in / &
                          & (xSize(i1d) * ySize(i1d) * &
                           & (dz_past(i1d)*weight_past + dz_future(i1d)*(1.-weight_past)))
            end do
          end do
        else   ! no dz
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              i1d = ix+(iy-1)*nx_dispersion
              datOut(i1d) = datOut(i1d) * factor_out + &
                          & arDat(iSpecies,iSrc,iz,ix,iy) * factor_in / &
                          & (xSize(i1d) * ySize(i1d))
            end do
          end do
        endif
      elseif(ifPerDensity)then  ! normalized to air density, no volume
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            i1d = ix+(iy-1)*nx_dispersion
            cnc_air = (weight_past*rho_past(i1d) + (1.0-weight_past)*rho_future(i1d)) / &
                    & molecular_weight_air
            datOut(i1d) = datOut(i1d) * factor_out + &
                        & arDat(iSpecies,iSrc,iz,ix,iy) * factor_in / cnc_air
          end do
        end do
      else      ! no normalization
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            i1d = ix+(iy-1)*nx_dispersion
            datOut(i1d) = datOut(i1d) * factor_out + arDat(iSpecies,iSrc,iz,ix,iy) * factor_in
          end do
        end do
      endif  ! type of scaling

    end subroutine merge_data_arrays
    
    !========================================================================

    subroutine merge_cumulative_arrays(arDat, in_extra_factor, &
                                     & xSize, ySize, dz_past, dz_future, rho_past, rho_future, &
                                     & weight_past, iSrc, iSpecies, iz, &
                                     & ifPerVolume, ifPerDensity, ifPerDz, &
                                     & datOut)
      implicit none

      ! Imported parameters
      real, dimension(:,:,:,:,:), pointer :: arDat
      logical, intent(in) :: ifPerVolume, ifPerDensity, ifPerDz
      real, dimension(:), pointer :: datOut
      real, dimension(:), pointer :: xSize, ySize, dz_past, dz_future, rho_past, rho_future
      real, intent(in) :: weight_past
      integer, intent(in) :: iSrc, iSpecies, iz
      real, intent(in) :: in_extra_factor

      type(silja_time) :: tValidOut, tValidTarget, tValidIn, tAccumStartTarget, & 
                        & tAccumStartIn, tAccumStartOut
      real, dimension(:), pointer :: arTmp
      !
      ! we need the internal field to handle 2 cumulatives
      !
      tValidOut = fu_valid_time(idOut)
      tValidTarget = fu_valid_time(targetId)
      tValidIn = fu_valid_time(idIn)
      tAccumStartTarget = fu_accumulation_start_time(targetId)
      !
      ! So far this function does not work for backward direction in time.
      !
      if(tValidTarget < tValidOut)then
        call set_error('Inverse time direction is not supported','merge_cumul_2_cumulative')
        call unset_error('merge_cumul_2_cumulative')
        return
      endif
      !
      ! There will be several cases of mutual disposition of accumulation
      ! periods, and they will be considered one-by-one
      !
      if(tValidIn + fu_validity_length(idIn) <= tAccumStartTarget)then
        !
        ! Before-target-period accumulation - drop it to the internal id
        !
        idOut_internal= idIn
        call set_met_src(idOut_internal, silam_internal_src)
        call merge_data_arrays(arDat, ifPerVolume, ifPerDensity, ifPerDz, &
                             & xSize, ySize, dz_past, dz_future, rho_past, rho_future, &
                             & weight_past, iSrc, iSpecies, iz, &
                             & dataOut_internal, in_extra_factor, 0.)
        !
        ! May be, it is time for output ? Possible, if targetId is void
        ! and has zero accumulation length
        !
        if(tValidIn == tValidTarget)then
          idOut = targetId
          call set_level(idOut, fu_level(idIn))
          call set_validity_length(idOut, fu_validity_length(idIn))
          dataOut(1:fs) = 0.
        endif

      elseif(tValidIn <= tValidTarget)then
        !
        ! idIn covers at least part of the target accumulation period.
        !
        tAccumStartIn = fu_accumulation_start_time(idIn)
        tAccumStartOut = fu_accumulation_start_time(idOut)
        tAccumStartTarget = fu_accumulation_start_time(targetId)

        if(tAccumStartIn == tAccumStartTarget .or. &   !Accum.per. coincide
         & tAccumStartIn == tAccumStartOut .or. &      !Input start same as out
         & tAccumStartIn == fu_valid_time(idOut_internal) .or. &            !Input adds-up to internal
         & tAccumStartIn == fu_accumulation_start_time(idOut_internal))then !Input start same as internal
          !
          ! New field continues the internal one or substitutes the output
          !
          idOut = idIn
          call merge_data_arrays(arDat, ifPerVolume, ifPerDensity, ifPerDz, &
                               & xSize, ySize, dz_past, dz_future, rho_past, rho_future, &
                               & weight_past, iSrc, iSpecies, iz, &
                               & dataOut, in_extra_factor, 0.) !dataOut(1:fs) = dataIn(1:fs)

        elseif(tAccumStartIn == tValidOut)then
          !
          ! New field adds-up to the currently stored accumulation
          !
          call merge_data_arrays(arDat, ifPerVolume, ifPerDensity, ifPerDz, &
                               & xSize, ySize, dz_past, dz_future, rho_past, rho_future, &
                               & weight_past, iSrc, iSpecies, iz, &
                               & dataOut, in_extra_factor, 1.) !dataOut(1:fs) = dataOut(1:fs) + dataIn(1:fs)
          call set_accumulation_length(idOut,(tValidIn - tAccumStartOut))
          call set_valid_time(idOut, tValidIn)
          call set_validity_length(idOut, fu_validity_length(idIn))

        else
          !
          ! New field overlaps or has holes with regard to the current output
          !
          dataOut(1:fs) = 0.
          call set_valid_time(idOut, tValidTarget)
          call report_error('Inconsistent accumulation',.true.)
          return
        endif
        !
        ! May be, it is time for output ?
        ! Then store internal field and, for averaged fields, divide the 
        ! accumulated data by the period
        !
        tValidOut = fu_valid_time(idOut)
        if(tValidOut == tValidTarget)then
          !
          ! For output we may need to subtract the internal field - if its valid
          ! time == start of targetId accumulation and its accumulation start ==
          ! accumulation start of output field
          ! But storing the current output is always useful.
          !
          arTmp => fu_work_array()
          arTmp(1:fs) = dataOut(1:fs)

          if(tAccumStartOut == tAccumStartTarget) then
            !
            ! Do nothing - Out is OK as it is
            !
          elseif(tAccumStartTarget == fu_valid_time(idOut_internal) .and. &
               & tAccumStartOut == fu_accumulation_start_time(idOut_internal))then
            !
            ! Out and Out_internal start accumulation from the same time, prior
            ! to the targetId, which starts at the valid time of Out_internal
            ! Subtraction is justified
            !
            dataOut(1:fs) = dataOut(:) - dataOut_internal(:)
          else
            dataOut(1:fs) = 0.
            call set_valid_time(idOut, fu_valid_time(targetId))
            call free_work_array(arTmp)
            call report_error('Cum/Ave->cum/ave failed at last step',.true.)
            return
          endif
          dataOut_internal(1:fs) = arTmp(1:fs)
          call free_work_array(arTmp)
          idOut_internal = idOut
          call set_met_src(idOut_internal, silam_internal_src)
          idOut = targetId
          call set_level(idOut, fu_level(idIn))
          call set_validity_length(idOut, fu_validity_length(idIn))

        endif  ! time to output: tValidOut == tValidTarget

      elseif(tValidIn > tValidTarget)then
        !
        ! idIn passes through the end of the target accumulation period.
        !
        call report_error('Input accumulation passes through output time', .false.)
        dataOut(1:fs) = 0.
        call set_valid_time(idOut, tValidTarget)
        return
      else
        !
        ! Something strange happend.
        !
        call report_error('Strange times detected', .false.)
        dataOut(1:fs) = 0.
        call set_valid_time(idOut, tValidTarget)
        return
      endif  ! validity of input and target

    end subroutine merge_cumulative_arrays

    !========================================================================
    
    subroutine report_error(chMsg, ifReportInternal)
      implicit none
      character(len=*), intent(in) :: chMsg
      logical, intent(in) :: ifReportInternal
        call msg('')
        call msg('Input field id:')
        call report(idIn)
        call msg('Currently stored field id')
        call report(idOut)
        call msg('Current target field id')
        call report(targetId)
        call report(stackPtr)
        call set_error(chMsg,'merge_cumul_2_instant')
        if(ifReportInternal)then
          call msg('Currently stored intermediate field id')
          call report(idOut_internal)
        endif
    end subroutine report_error
    
  end subroutine merge_mass_map_2_stack


  !*************************************************************************

  subroutine merge_vector_data_to_stack(idIn, dataIn, stackOut, targetId, ifRandomise) !, direction)
    !
    ! Merges the given data array to the stack. The method of merging
    ! is determined by target id and type of field in the stack.
    ! TargetId tells how the field should look like - for example, whether it
    ! should be averaged, accumulated, etc.
    !
    ! Be carefull with target_Id - it has undefined level if the quantity is 3D.
    ! Also, be carefull if the targetId is averaged. Then we have to keep 
    ! 2 fields - one for before-the-integration time period and one for 
    ! during-the-integration time period. They are to be subtracted one 
    ! from another immediately before the output. The before-the-integration 
    ! data are marked with SILAM intermediate met_src.= silam_internal_src
    !
    ! Finally, we allow two sources of input data: MassMap and field (i.e., just 1d vector)
    !
    implicit none

    ! Imported parameters 
    type(silja_field_id), intent(in) :: idIn, targetId
    real, dimension(:), intent(in), target :: dataIn
    type(silja_stack), intent(inout) :: stackOut
    logical, intent(in) :: ifRandomise

    ! Local variables
    !
    type(silja_field_id), pointer :: idTmp, idOut, idOut_internal
    type(silja_field), pointer :: field
    real, dimension(:), pointer :: dataOut, dataOut_internal, dataInOutgrid, datawork
    integer :: fs, i, timeDir, nx_target, ny_target,  nx_in, ny_in
    logical :: found, found_internal, ifTmpGrid
    real :: avreage_seconds, shift_lat, shift_lon
    logical :: selection_done, grid_ok
    type(silja_grid) :: grid_new
    character(len=*), parameter :: sub_name = 'merge_vector_data_to_stack'
    

    !
    ! First of all, find the receiving fields
    !
    call find_receiving_fields_4_merge(idIn, stackOut, targetId, &
                                     & found, found_internal, timeDir, &
                                     & idOut, dataOut, idOut_internal, dataOut_internal)
    if(error .or. .not. found)return

    call grid_dimensions(fu_grid(targetId), nx_target, ny_target)
    call grid_dimensions(fu_grid(idIn), nx_in, ny_in)
    if(error)return

    fs = nx_target * ny_target


    if (.not. (fu_grid(targetId) == fu_grid(idIn))) then
      
#ifdef DEBUG
      call msg('Forcing grid_new into the target grid, forbidding Arakawa')
#endif
      grid_new = fu_grid(targetId)
      
      dataInOutGrid => fu_work_array(fu_number_of_gridpoints(grid_new))
      if(error)return
      call grid_data_horizontal_select(fu_grid(idIn), &
                                     & dataIn, selection_done, grid_new, &
                                     & dataInOutGrid, &
                                     & ifRandomise, &
                                     & iAccuracy=5, &
                                     & method=linear, &
                                     & fMissingValue=real_missing, &
                                     & iOutside=notAllowed)
      if (error) return
      
      grid_ok = fu_grid(targetId) == grid_new
      if (.not. grid_ok .and. &
        & fu_number_of_gridpoints(fu_grid(targetId)) == fu_number_of_gridpoints(grid_new)) then
        call grid_shift_indices(fu_grid(targetId), grid_new, shift_lon, shift_lat)
        if (abs(shift_lat) <= 0.5 .and. abs(shift_lon) <= 0.5) then
          ! some arakawa shift => just ignore in output
          grid_new = fu_grid(targetId)
          grid_ok = .true.
        end if
      end if
      
      if (.not. grid_ok) then
        call msg('Target grid:')
        call report(fu_grid(targetId))
        call msg('Grid_new')
        call report(grid_new)
        call set_error('Can''t get grids to match', sub_name)
        return
      end if
    else
      datainOutgrid => dataIn
    end if




    !-------------------------------------------------------------------
    !
    ! When the run just starts, there is nothing in the stack, therefore the
    ! validity_length is set to missing. Then the first field that comes puts itself
    ! into the data and overwrites the validity_length.
    ! In principle, all times must be missing but then it is nearly impossible to 
    ! get around such field
    !
    if(.not. defined(fu_validity_length(idOut)))then
      dataOut(1:fs) = dataInOutgrid(1:fs)
!      call set_validity_length(idOut, fu_validity_length(idIn))
      call set_validity_length(idOut, zero_interval)
!      if(found_internal) call set_validity_length(idOut_internal, fu_validity_length(idIn))
      
    else

!!!   The stuff is wrong!       
!!!    !
!!!    ! If valid times are already OK - fields are the same.
!!!    ! The only thing, which might eventually be out of sense is validity length
!!!    !
!!!    if(fu_valid_time(idOut) == fu_valid_time(targetId))then
!!!      call set_validity_length(idOut, fu_validity_length(targetId))
!!!      return
!!!    endif

    !-------------------------------------------------------------------
    !
    ! Having both ids known, let's check what should be done with the data,
    ! depending on the averaging type of idIn, idOut and targetId
    !
       select case(fu_field_kind(targetId))

         case(forecast_flag)  
           !
           ! Target ID is instant
           !
           select case(fu_field_kind(idIn))

             case(forecast_flag, averaged_flag) ! , period_valid_flag)
               call merge_inst_and_aver_2_instant(idIn, dataInOutgrid, idOut, dataOut, targetId, fs)

             case(accumulated_flag) !, accumulated_period_valid_flag)
               call merge_cumul_2_instant(idIn, dataInOutgrid, idOut, dataOut, targetId, fs)
               
             case default
               dataOut(1:fs) = 0.
               call set_valid_time(idOut, fu_valid_time(targetId))
               call set_error('Unsupported input field type, instant target',sub_name)
           end select


         case(accumulated_flag, averaged_flag) !, accumulated_period_valid_flag) ! Target ID
           !
           ! TargetID is cumulative or averaged.
           ! Cumulative and averaged fields have to be integrated in time
           ! from beginning of the current accumulation period till its end.
           ! Average field will at the end be divided by length of integration.
           !
           select case(fu_field_kind(idIn))

             case (forecast_flag)  
               !call merge_inst_2_cumulative(idIn, dataInOutgrid, idOut, dataOut, targetId, fs)
               call merge_interval_2_cumulative(idIn, dataInOutgrid, idOut, dataOut, targetId, fs)

             case (accumulated_flag) !, accumulated_period_valid_flag)

               if(found_internal)then 
                 !
                 ! we need the internal field to handle 2 cumulatives
                 !
                 call merge_cumul_2_cumulative(idIn, dataInOutgrid, idOut, dataOut, &
                                            & idOut_internal, dataOut_internal, targetId, fs)

               else
                 !
                 ! May be, the in-field is just what we expect ? If we want it to be 
                 ! just dropped into the stack - the internal field may be not needed
                 !
                 if(idIn == targetId)then
                   idOut = targetId
                   call set_validity_length(idOut, fu_validity_length(idIn))
                   dataOut(1:fs) = dataInOutgrid(1:fs)
                 else
                   call msg("")
                   call msg("")
                   call msg_warning("Can not find internal field for merging",sub_name)
                   call msg('Input ID')
                   call report(idIn)
                   call msg("")
                   call msg('Target ID')
                   call report(targetId)
                   call msg("")
                   call msg('STACK to search for')
                   call report(stackOut)
                   call msg("")
                   call set_error('Can not find internal field for merging',sub_name)
                 endif
               endif     ! .not.found_internal

             case (averaged_flag) !, accumulated_period_valid_flag)
               !
               !
               ! Integrated or averaged - take with care in order not to 
               ! overlap the accumulation periods. Use intermediate field to store
               ! the prior-to-target accumulation.
               !
               if(found_internal)then
                 !
                 ! Convert average field to cumulative one and do not forget the 
                 ! internal field needed to handle 2 cumulatives
                 !
                 dataInOutgrid(1:fs) = dataInOutgrid(1:fs) * fu_sec(fu_accumulation_length(idIn))
                 call merge_cumul_2_cumulative(idIn, dataInOutgrid, idOut, dataOut, &
                                            & idOut_internal, dataOut_internal, targetId, fs)

               else
                 !
                 ! May be, the in-field is just what we expect ? If we want it to be 
                 ! just dropped into the stack - the internal field may be not needed
                 !
                 if(idIn == targetId)then
                   idOut = targetId
                   call set_validity_length(idOut, fu_validity_length(idIn))
                   dataOut(1:fs) = dataInOutgrid(1:fs)
                 else
                   call set_error('Can not find internal field for merging',sub_name)
                   call msg('Input ID')
                   call report(idIn)
                   call msg('STACK to search for')
                   call report(stackOut)
                 endif
               endif     ! .not.found_internal

             case default ! idIn
               !
               ! Types like difference_flag, period_averaging_flag
               !
               call set_error('Non-supported averaging type of input field',sub_name)
           end select

           if(fu_field_kind(targetId) == averaged_flag)then
             !
             ! For averaged field, if this is the last data collection before the
             ! output - divide the integrated stuff by the accumulation length
             !
   !          call msg("Finalizing averaging for "+ fu_quantity_short_string(fu_quantity(targetId)))
             if(fu_valid_time(targetId) == fu_valid_time(idOut) .and. fu_field_kind(IdOut) == accumulated_flag )then

               ! accumulation complete.
               avreage_seconds = fu_sec(fu_accumulation_length(idOut))
               if( avreage_seconds == 0)then
                 dataOut(1:fs) = 0.
               else
                 dataOut(1:fs) = dataOut(1:fs) / abs(avreage_seconds)
               endif
               idOut = targetId
               call set_level(idOut, fu_level(idIn))
               !call set_validity_length(idOut, fu_set_interval_sec(-avreage_seconds))
               !  Hack! 
               call set_validity_length(idOut, zero_interval) 
             endif
           endif

         case default ! target ID
           dataOut(1:fs) = 0.
           call set_valid_time(idOut, fu_valid_time(targetId))
           call set_error('Non-supported averaging type of output field',sub_name)
       end select
    endif ! 
    if (.not. associated(dataInOutGrid, dataIn)) call free_work_array(dataInOutgrid)

!call msg('')
!call msg('The sum for:' + fu_quantity_short_string(fu_quantity(idOut)) + fu_substance_name(idOut),sum(dataOut(1:fs)))
!call report(idOut)

  end subroutine merge_vector_data_to_stack


  !**********************************************************************************

  subroutine merge_inst_and_aver_2_instant(idIn, dataIn, idOut, dataOut, targetId, fs)
    !
    ! Merges the data from all-but-accumulated fields into an instant output field
    ! Instant field is just picked if time is right
    ! Average field will substitute the instant one - at the last time moment
    !
    implicit none

    ! Imported parameters
    type(silja_field_id), pointer :: idOut
    real, dimension(:), pointer :: dataIn, dataOut
    type(silja_field_id), intent(in) ::  idIn,targetId
    integer, intent(in) :: fs

    !
    ! Two options are allowed: (i) valid times are exactly the same,
    ! (ii) target valid time is inside the averaging period of input field
    !
    if(fu_valid_time(idIn) == fu_valid_time(targetId)) then ! Last time ?
      dataOut(1:fs) = dataIn(1:fs)
      idOut = idIn
    elseif(fu_accumulated_id(idIn))then
      if(fu_between_times(fu_valid_time(targetId), &  ! time
                          & fu_valid_time(idIn), &      ! limit 1
                          & fu_valid_time(idIn) + fu_accumulation_length(idIn), &  ! limit 2
                          & .true.)) then                ! if include borders
        dataOut(1:fs) = dataIn(1:fs)
        idOut = idIn
      endif
    endif

  end subroutine merge_inst_and_aver_2_instant


  !**********************************************************************************

  subroutine merge_cumul_2_instant(idIn, dataIn, idOut, dataOut, targetId, fs)
    !
    ! Instant field will be substituted with the mean one between two last
    ! cumulative fields. If before that moment, take the new field
    !
    implicit none

    ! Imported parameters
    type(silja_field_id), pointer :: idOut
    real, dimension(:), pointer :: dataIn, dataOut
    type(silja_field_id), intent(in) ::  idIn,targetId
    integer, intent(in) :: fs

    if(fu_valid_time(idIn) == fu_valid_time(targetId)) then  ! Last time ?

      if(fu_accumulation_start_time(idIn) == fu_accumulation_start_time(idOut))then
        !
        ! Same start of the fields, their difference gives mean over last timestep
        !
        dataOut(1:fs) = (dataIn(1:fs) - dataOut(1:fs)) / &
                      & abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut)))
        idOut = targetId

      elseif(fu_accumulation_start_time(idIn) == fu_valid_time(idOut))then  
        !
        ! Last accumulation covers the needs
        !
        dataOut(1:fs) = dataIn(1:fs) / abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut)))
        idOut = targetId
      else
        call msg('')
        call msg('Input field id:')
        call report(idIn)
        call msg('Currently stored field id')
        call report(idOut)
        call msg('Current target field id')
        call report(targetId)
        dataOut(1:fs) = 0.
        call set_valid_time(idOut, fu_valid_time(targetId))
        call set_error('Can not make instant field from accumulated','merge_cumul_2_instant')
        return
      endif

    else    ! Intermediate collection - just store
      dataOut(1:fs) = dataIn(1:fs)
      idOut = idIn
    endif

  end subroutine merge_cumul_2_instant


!  !**********************************************************************************
!
!  subroutine merge_inst_2_cumulative(idIn, dataIn, idOut, dataOut, targetId, fs)
!    !
!    ! idIn, when targetId==accumulated/averaged
!    !
!    ! Instant field. Integrate from previous valid time till current one,
!    ! if we are inside the accumulation period
!    implicit none
!
!    ! Imported parameters
!    type(silja_field_id), intent(inout) :: idOut
!    real, dimension(:), pointer :: dataIn, dataOut
!    type(silja_field_id), intent(in) ::  idIn,targetId
!    integer, intent(in) :: fs
!
!    real :: fTmp
!    !
!    ! Before integrating the new value let's check the starting time
!    ! of the main field and target ID. For example, it may happen that
!    ! we already should start the new accumulation cycle. That is the case if
!    ! the accumulations start time is later (in the forward run) than the 
!    ! accumulation start time of the output field but the output valid time 
!    ! is the latest
!    !
!    if(fu_between_times(fu_valid_time(targetId), &  ! time
!                      & fu_accumulation_start_time(targetId), &  ! limit 1
!                      & fu_valid_time(idIn), .false.) .or. &  ! limit 2, include borders
!!     & fu_between_times(fu_accumulation_start_time(targetId), &  ! time
!!                      & fu_accumulation_start_time(idOut), &  ! limit 1
!!                      & fu_valid_time(idOut), .false.)) then   ! limit 2, include borders
!     & fu_accumulation_start_time(targetId) /= fu_accumulation_start_time(idOut)) then
!     ! New accumulation
!      call msg('')
!      call msg('')
!      call msg('')
!      call msg('Enter zeroying if: idIn, idOut, target id')
!      if (fu_accumulation_start_time(targetId) /= fu_accumulation_start_time(idOut)) then
!              call msg('Acc start time differ')
!      else
!               call msg('Valid time beyond the target ID')
!      endif
!      call report(idIn)
!      call msg('')
!      call report(idOut)
!      call msg('')
!      call report(targetId)
!      call msg('done')
!      dataOut(1:fs) = 0.
!      call set_accumulation_length(idOut, zero_interval)
!      call set_analysis_time(idOut,fu_accumulation_start_time(targetId))
!      call set_valid_time(idOut,fu_accumulation_start_time(targetId))
!    endif
!    !
!    ! The main action happens if and only if the idIn is in-between the target valid_time
!    ! and target accumulation start time
!    !
!    if(fu_between_times(fu_valid_time(idIn), &  ! time
!                      & fu_accumulation_start_time(targetId), &  ! limit 1
!                      & fu_valid_time(targetId), .true.)) then   ! limit 2, true to include borders
!      !
!      ! Of course, it is not always correct integration - in reality, the new
!      ! field is valid at valid_time, while below it is used from 
!      ! last_integrated time till valid_time. 
!      ! Ideally, there should be symmetrical integration period, but for that
!      ! I would need 3 fields - past, current and future, which is too complicated.
!      ! See also io_server, where this error is somewhat coped in routine write_output.
!      !
!      fTmp = abs(fu_sec(fu_valid_time(idIn) - fu_valid_time(idOut)))
!
!      dataOut(1:fs) = dataOut(1:fs) + dataIn(1:fs) * fTmp
!
!      call set_accumulation_length(idOut, fu_valid_time(idIn) - fu_accumulation_start_time(idOut))
!      call set_valid_time(idOut, fu_valid_time(idIn))
!      call set_validity_length(idOut, fu_validity_length(idIn))
!    endif
!
!  end subroutine merge_inst_2_cumulative
  !**********************************************************************************

  subroutine merge_interval_2_cumulative(idIn, dataIn, idOut, dataOut, targetId, fs)
    !
    ! idIn is period-valid
    ! targetId is accumulated/averaged
    !
    ! No instant fields here. They can not be averaged without interpolation,
    ! which requires two fields.
    ! Cumulative can be backwards, but validity of idIn goes normally.
    ! Some balck magick is done to make the thing work properly with 
    ! averaging -- valid_time of  idOut and targetId is _earliest_ for adjoint
    !              accumulation time is negative, but weights are always  positive.
    !              
    
    implicit none

    ! Imported parameters
    type(silja_field_id), intent(inout) :: idOut
    real, dimension(:), pointer :: dataIn, dataOut
    type(silja_field_id), intent(in) ::  idIn,targetId
    type(silja_time) :: valid_in, valid_in_prev
    integer, intent(in) :: fs

    real :: fTmp
    integer :: iTmp
  
!    iTmp = temperature_2m_flag
!    !    iTmp = large_scale_rain_int_flag
!
!    call msg("Merging to cumulative field for:" + fu_quantity_short_string(fu_quantity(targetId)))
!    if (fu_quantity(targetId) == iTmp ) then
!          ! New accumulation
!            call msg('')
!            call msg('')
!            call msg('Before merging....')
!           call msg('idIn:')
!            call report(idIn)
!            call msg('')
!            call msg('idOut:')
!            call report(idOut)
!            call msg('')
!            call msg('target id:')
!            call report(targetId)
!            call msg('done')
!    endif


    !
    ! First, figure out which of the times of the input fields is the right 
    ! one to use for accumulation  -- latest for forward, and earliest for
    ! adjoint

    if (fu_sec(fu_accumulation_length(targetId)) > 0 ) then !forward
            valid_in_prev = fu_valid_time(idIn)
            valid_in =  valid_in_prev + fu_validity_Length(idIn)
    else
            valid_in = fu_valid_time(idIn)
            valid_in_prev = valid_in + fu_validity_Length(idIn)
    endif


    if (fu_between_times(valid_in_prev, &  ! time
                  & fu_valid_time(IdOut), &  ! limit 1
                  & valid_in, .false.)) then 
           call msg_warning("Discontinuity in averaging!!!", "merge_interval_2_cumulative")
           ! Should we die here????
    endif 

    !Cut validity to current TargetID
    if (fu_between_times(fu_valid_time(targetId), &  ! time
                  & fu_accumulation_start_time(targetId), &  ! limit 1
                  & valid_in, .false.)) then 
      valid_in = fu_valid_time(targetId)
    endif 

    ! Is the averaging complete already???
    ! Should happen with single-time fields
    if (valid_in == fu_valid_time(idOut)) then
      call msg("Can't add anything new for:"+fu_quantity_short_string(fu_quantity(targetId)))
      return
    endif
            
    
    ! Before integrating the new value let's check the starting time
    ! of the main field and target ID. For example, it may happen that
    ! we already should start the new accumulation cycle. That is the case if
    ! the accumulations start time is later (in the forward run) than the 
    ! accumulation start time of the output field but the output valid time 
    ! is the latest
    !
    if(fu_between_times(fu_valid_time(targetId), &  ! time
                      & fu_accumulation_start_time(targetId), &  ! limit 1
                      & valid_in, .false.) .or. &  ! limit 2, not including borders
     & fu_accumulation_start_time(targetId) /= fu_accumulation_start_time(idOut)) then
     ! New accumulation
      call msg('')
      call msg('')
      call msg_warning('Enter zeroying if: idIn, idOut, target id','merge_interval_2_cumulative')
      if (fu_accumulation_start_time(targetId) /= fu_accumulation_start_time(idOut)) then
        call msg('Acc start time differ')
      else
        call msg('Valid time beyond the target ID')
      endif
      call report(idIn)
      call msg('')
      call report(idOut)
      call msg('')
      call report(targetId)
      call msg('done')
      dataOut(1:fs) = 0.
      call set_accumulation_length(idOut, zero_interval)
      call set_analysis_time(idOut,fu_accumulation_start_time(targetId))
      call set_valid_time(idOut,fu_accumulation_start_time(targetId))
    endif


    ! The main action happens if and only if the idIn is in-between the target valid_time
    ! and target accumulation start time
    !
    if(fu_between_times(valid_in, &  ! time
                      & fu_accumulation_start_time(targetId), &  ! limit 1
                      & fu_valid_time(targetId), .true.)) then   ! limit 2, true to include borders
      
      fTmp = abs(fu_sec(valid_in - fu_valid_time(idOut)))

      dataOut(1:fs) = dataOut(1:fs) + dataIn(1:fs) * fTmp

      call set_accumulation_length(idOut, valid_in - fu_accumulation_start_time(idOut))
      call set_field_kind(idOut, accumulated_flag)
      call set_valid_time(idOut, valid_in)
    endif

  end subroutine merge_interval_2_cumulative


  !*********************************************************************************

  subroutine merge_cumul_2_cumulative(idIn, dataIn, idOut, dataOut, &
                                    & idOut_internal, dataOut_internal, targetId, fs)
    !
    ! Merges two cumulative fields.
    !
    implicit none

    ! Imported parameters
    type(silja_field_id), pointer :: idOut, idOut_internal
    real, dimension(:), pointer :: dataIn, dataOut, dataOut_internal
    type(silja_field_id), intent(in) ::  idIn,targetId
    integer, intent(in) :: fs

    ! Internal variables
    real, dimension(:), pointer :: arTmp
    type(silja_time) :: tAccumStartIn, tAccumStartOut, tAccumStartTarget, &
                      & tValidIn, tValidOut, tValidTarget

    tValidOut = fu_valid_time(idOut)
    tValidTarget = fu_valid_time(targetId)

    !
    ! So far this function does not work for abckward direction in time.
    !
    if(tValidTarget < tValidOut)then
      call set_error('Inverse time direction is not supported','merge_cumul_2_cumulative')
      call unset_error('merge_cumul_2_cumulative')
      return
    endif

    tValidIn = fu_valid_time(idIn)
    tAccumStartTarget = fu_accumulation_start_time(targetId)

    !
    ! There will be several cases of mutual disposition of accumulation
    ! periods, and they will be considered one-by-one
    !
    if(tValidIn + fu_validity_length(idIn) <= tAccumStartTarget)then
      !
      ! Before-target-period accumulation - drop it to the internal id
      !
      idOut_internal= idIn
      call set_met_src(idOut_internal, silam_internal_src)
      dataOut_internal(1:fs) = dataIn(1:fs)
      !
      ! May be, it is time for output ? Possible, if targetId is void
      ! and has zero accumulation length
      !
      if(tValidIn == tValidTarget)then
        idOut = targetId
        call set_level(idOut, fu_level(idIn))
        call set_validity_length(idOut, fu_validity_length(idIn))
        dataOut(1:fs) = 0.
      endif

    elseif(tValidIn <= tValidTarget)then
      !
      ! idIn covers at least part of the target accumulation period.
      !
      tAccumStartIn = fu_accumulation_start_time(idIn)
      tAccumStartOut = fu_accumulation_start_time(idOut)
      tAccumStartTarget = fu_accumulation_start_time(targetId)

      if(tAccumStartIn == tAccumStartTarget .or. &   !Accum.per. coincide
       & tAccumStartIn == tAccumStartOut .or. &      !Input start same as out
       & tAccumStartIn == fu_valid_time(idOut_internal) .or. &            !Input adds-up to internal
       & tAccumStartIn == fu_accumulation_start_time(idOut_internal))then !Input start same as internal
        !
        ! New field continues the internal one or substitutes the output
        !
        idOut = idIn
        dataOut(1:fs) = dataIn(1:fs)

      elseif(tAccumStartIn == tValidOut)then
        !
        ! New field adds-up to the currently stored accumulation
        !
        dataOut(1:fs) = dataOut(1:fs) + dataIn(1:fs)
        call set_accumulation_length(idOut,(tValidIn - tAccumStartOut))
        call set_valid_time(idOut, tValidIn)
        call set_validity_length(idOut, fu_validity_length(idIn))

      else
        !
        ! New field overlaps or has holes with regard to the current output
        !
        call msg("")
        call msg("")
        call msg('Inconsistent accumulation')
        call msg("")
        call msg('Input field id:')
        call report(idIn)
        call msg("")
        call msg('Currently stored field id')
        call report(idOut)
        call msg("")
        call msg('Currently stored intermediate field id')
        call report(idOut_internal)
        dataOut(1:fs) = 0.
        call set_valid_time(idOut, tValidTarget)
        call set_error('Inconsistent accumulation','merge_cumul_2_cumulative')
        return
      endif
      !
      ! May be, it is time for output ?
      ! Then store internal field and, for averaged fields, divide the 
      ! accumulated data by the period
      !
      tValidOut = fu_valid_time(idOut)
      if(tValidOut == tValidTarget)then
        !
        ! For output we may need to subtract the internal field - if its valid
        ! time == start of targetId accumulation and its accumulation start ==
        ! accumulation start of output field
        ! But storing the current output is always useful.
        !
        arTmp => fu_work_array()
        arTmp(1:fs) = dataOut(1:fs)

        if(tAccumStartOut == tAccumStartTarget) then
          !
          ! Do nothing - Out is OK as it is
          !
        elseif(tAccumStartTarget == fu_valid_time(idOut_internal) .and. &
             & tAccumStartOut == fu_accumulation_start_time(idOut_internal))then
          !
          ! Out and Out_internal start accumulation from the same time, prior
          ! to the targetId, which starts at the valid time of Out_internal
          ! Subtraction is justified
          !
          dataOut(1:fs) = dataOut(:) - dataOut_internal(:)

        else
          call msg('Input field id:')
          call report(idIn)
          call msg('Target field id')
          call report(targetId)
          call msg('Currently stored field id')
          call report(idOut)
          call msg('Currently stored intermediate field id')
          call report(idOut_internal)
          dataOut(1:fs) = 0.
          call set_valid_time(idOut, fu_valid_time(targetId))
          call set_error('Cum/Ave->cum/ave failed at last step', 'merge_cumul_2_cumulative')
          call free_work_array(arTmp)
          return
        endif

        dataOut_internal(1:fs) = arTmp(1:fs)
        call free_work_array(arTmp)
        idOut_internal = idOut
        call set_met_src(idOut_internal, silam_internal_src)
        idOut = targetId
        call set_level(idOut, fu_level(idIn))
        call set_validity_length(idOut, fu_validity_length(idIn))

      endif

    elseif(tValidIn > tValidTarget)then
      !
      ! idIn passes through the end of the target accumulation period.
      !
      call msg('Input field id:')
      call report(idIn)
      call msg('Target field id')
      call report(targetId)
      dataOut(1:fs) = 0.
      call set_valid_time(idOut, tValidTarget)
      call set_error('Input accumulation passes through output time','merge_cumul_2_cumulative')
      return
    else
      !
      ! Something strange happend.
      !
      call msg('Input field id:')
      call report(idIn)
      call msg('Target field id')
      call report(targetId)
      call msg('Currently stored field id')
      call report(idOut)
      call msg('Currently stored intermediate field id')
      call report(idOut_internal)
      dataOut(1:fs) = 0.
      call set_valid_time(idOut, tValidTarget)
      call set_error('Strange times detected', 'merge_cumul_2_cumulative')
      return
    endif

  end subroutine merge_cumul_2_cumulative


!!  !*****************************************************************************
!!  
!!  subroutine merge_lagrange_mass_2_stack(lpDyn, lpMass, lpStatus, &   ! Coord and masses to merge
!!                                       & nop, &                ! number of particles to look at
!!                                       & idTemplate, &         ! Main parameters of fields to make
!!                                       & verticalTemplate, &   ! ... and their vertical
!!                                       & ifPerVolume, ifPerDensity, ifPerDz, &  ! normalization
!!                                       & iSpecies, iSrc, weight_past, &
!!                                       & data_buf, & ! for xSize, ySize, zSize etc
!!                                       & stackPtr, targetId)
!!    implicit none
!!    !
!!    ! Merges the mass map to the output stack. The sub picks out a requested species and
!!    ! emission source. Normalization allowed is per-volume and per-air-denstiy
!!    ! Does the same thing as merge_vector_data_to_stack below but for mass_map
!!    ! ATTENTION.
!!    ! It is called after Eulerian data collection and the stuff goes into the same fields. Also,
!!    ! Lagrangian particles noramally project many-to-many to the output grid cells. Therefore, 
!!    ! we have to be careful: nullification of the fields may be mandatory, especially if 
!!    ! the run is purely Lagrangian.
!!    !
!!    
!!    ! Imported parameters
!!    real, dimension(:,:), pointer :: lpDyn, lpMass
!!    integer, dimension(:), pointer :: lpStatus
!!    type(silja_field_id), intent(in) :: idTemplate    ! Main parameters of the fields to make
!!    type(silja_field_id), intent(inout) :: targetId   ! Main parameters of the target output fields
!!    type(silam_vertical), intent(in) :: verticalTemplate ! vertical to make
!!    logical, intent(in) :: ifPerVolume, ifPerDensity, ifPerDz ! normalization
!!    real, dimension(:), pointer :: xSize, ySize, zSize
!!!    type(field_3d_data_ptr), pointer :: dz_past_p3d, dz_future_p3d, rho_past_p3d, rho_future_p3d
!!    integer, intent(in) :: iSpecies, iSrc, nop   ! "filters": we need only this source and only up to nop particles
!!    real, intent(in) :: weight_past
!!    type(silja_stack), pointer :: stackPtr
!!    
!!    ! Local variables
!!    type(silja_field_id) :: idTmp
!!    type(silja_field_id), pointer :: idOut_internal
!!    type(silja_field), pointer :: field
!!    type(field_3d_data_ptr) :: field3d
!!    real, dimension(:), pointer :: dataOut_internal
!!    integer :: ix, iy, iz, nz, timeDir, fs
!!    logical :: found, found_internal
!!    real, dimension(:,:,:,:,:), pointer :: arDat
!!    real :: cnc_air
!!
!!    !
!!    ! First of all, find the receiving fields
!!    !
!!    idTmp = idTemplate
!!    nz = fu_NbrOfLevels(verticalTemplate)
!!    do iz = 1, nz
!!      call set_level(idTmp, fu_level(verticalTemplate,iz))
!!      call find_receiving_fields_4_merge(idTmp, stackPtr, targetId, &
!!                                       & found, found_internal, timeDir, &
!!                                       & field3d%p2d(iz)%idPtr, field3d%p2d(iz)%ptr, & ! receiving flds
!!                                       & idOut_internal, dataOut_internal)
!!      if(error .or. .not. found)return
!!    enddo ! iz
!!
!!    fs = fu_number_of_gridpoints(fu_grid(targetId))
!!    !
!!    ! When the run just starts, there is nothing in the stack, therefore the
!!    ! validity_length is set to missing. Then the first field that comes puts itself
!!    ! into the data and overwrites the validity_length.
!!    ! In principle, all times must be missing but then it is nearly impossible to 
!!    ! get around such field.
!!    !
!!    if(.not. defined(fu_validity_length(field3d%p2d(1)%idPtr)))then
!!      call merge_data_arrays(ifPerVolume, ifPerDensity, ifPerDz, field3d, 1.0, 0.0, &  ! replace
!!                           & xSize, ySize, zSize)
!!      do iz = 1, nz
!!        call set_validity_length(field3d%p2d(iz)%idPtr, fu_validity_length(idTemplate))
!!      enddo
!!      return
!!    endif
!!    !
!!    ! If valid times are already OK - fields are the same.
!!    ! The only thing, which might eventually be out of sense is validity length
!!    !
!!    if(fu_valid_time(field3d%p2d(1)%idPtr) == fu_valid_time(targetId))then
!!      do iz = 1, nz
!!        call set_validity_length(field3d%p2d(iz)%idPtr, fu_validity_length(targetId))
!!      end do
!!      return
!!    endif
!!
!!    !-------------------------------------------------------------------
!!    !
!!    ! Having both ids known, let's check what should be done with the data,
!!    ! depending on the averaging type of idIn, idOut and targetId
!!    ! Note that the input stuff is the lagrangian lpSet, i.e. instant fields
!!    !
!!    select case(fu_field_kind(targetId))
!!
!!      case(forecast_flag)  
!!        !
!!        ! Target ID is instant
!!        !
!!        select case(fu_field_kind(idTemplate))
!!
!!          case(forecast_flag) ! instant + aver -> instant
!!            if(fu_valid_time(idTemplate) == fu_valid_time(targetId)) then ! Last time ?
!!              call merge_data_arrays(ifPerVolume, ifPerDensity, ifPerDz, field3d, 1.0, 0.0, &  ! replace
!!                                   & xSize, ySize, zSize)
!!              do iz = 1, nz
!!                field3d%p2d(iz)%idPtr = idTemplate
!!                call set_level(field3d%p2d(iz)%idPtr, fu_level(verticalTemplate,iz))
!!              end do
!!            elseif(fu_accumulated(idTemplate))then
!!              if(fu_between_times(fu_valid_time(targetId), &  ! time
!!                                & fu_valid_time(idTemplate), &      ! limit 1
!!                                & fu_valid_time(idTemplate) + fu_accumulation_length(idTemplate), &  ! limit 2
!!                                & .true.)) then                ! if include borders
!!                call merge_data_arrays(ifPerVolume, ifPerDensity, ifPerDz, field3d, 1.0, 0.0, &  ! replace
!!                                     & xSize, ySize, zSize)
!!                do iz = 1, nz
!!                  field3d%p2d(iz)%idPtr = idTemplate
!!                  call set_level(field3d%p2d(iz)%idPtr, fu_level(verticalTemplate,iz))
!!                end do
!!              endif
!!            endif
!!
!!          case default
!!            do iz = 1, nz
!!              field3d%p2d(iz)%ptr(1:fs) = 0.
!!              call set_valid_time(field3d%p2d(iz)%idPtr, fu_valid_time(targetId))
!!            end do
!!            call report_error('Unsupported input field type, instant target', .false.)
!!            return
!!        end select
!!
!!      case(accumulated_flag, averaged_flag) !, accumulated_period_valid_flag) ! Target ID
!!        !
!!        ! TargetID is cumulative or averaged.
!!        ! Cumulative and averaged fields have to be integrated in time
!!        ! from beginning of the current accumulation period till its end.
!!        ! Average field will at the end be divided by length of integration.
!!        !
!!        select case(fu_field_kind(idTemplate))
!!
!!          case (forecast_flag)  
!!            if(fu_between_times(fu_valid_time(targetId), &  ! time
!!                              & fu_accumulation_start_time(targetId), &  ! limit 1
!!                              & fu_valid_time(idTemplate), .false.) .or. &  ! limit 2, include borders
!!             & fu_accumulation_start_time(targetId) /= fu_accumulation_start_time(field3d%p2d(1)%idPtr)) then
!!              do iz = 1, nz
!!                field3d%p2d(iz)%ptr(1:fs) = 0.
!!                call set_accumulation_length(field3d%p2d(iz)%idPtr, zero_interval)
!!                call set_analysis_time(field3d%p2d(iz)%idPtr, fu_accumulation_start_time(targetId))
!!                call set_valid_time(field3d%p2d(iz)%idPtr,fu_accumulation_start_time(targetId))
!!              enddo
!!            endif
!!            !
!!            ! The main action happens if and only if the idIn is in-between the target valid_time
!!            ! and target accumulation start time
!!            !
!!            if(fu_between_times(fu_valid_time(idTemplate), &  ! time
!!                              & fu_accumulation_start_time(targetId), &  ! limit 1
!!                              & fu_valid_time(targetId), .true.)) then   ! limit 2, true to include borders
!!              !
!!              ! Of course, it is not always correct integration - in reality, the new
!!              ! field is valid at valid_time, while below it is used from 
!!              ! last_integrated time till valid_time. 
!!              ! Ideally, there should be symmetrical integration period, but for that
!!              ! I would need 3 fields - past, current and future, which is too complicated.
!!              ! See also io_server, where this error is somewhat coped in routine write_output.
!!              !
!!              call merge_data_arrays(ifPerVolume, ifPerDensity, ifPerDz, field3d, &
!!                                   & abs(fu_sec(fu_valid_time(idTemplate) - fu_valid_time(field3d%p2d(1)%idPtr))), &
!!                                   & 1., & !      dataOut(1:fs) = dataOut(1:fs) + dataIn(1:fs) * fTmp
!!                                   & xSize, ySize, zSize)
!!              do iz = 1, nz
!!                call set_accumulation_length(field3d%p2d(iz)%idPtr, &
!!                                           & fu_valid_time(idTemplate) - fu_accumulation_start_time(field3d%p2d(1)%idPtr))
!!                call set_valid_time(field3d%p2d(iz)%idPtr, fu_valid_time(idTemplate))
!!                call set_validity_length(field3d%p2d(iz)%idPtr, fu_validity_length(idTemplate))
!!              enddo
!!            endif
!!
!!          case default ! idIn
!!            !
!!            ! Types like difference_flag, period_averaging_flag
!!            !
!!            call set_error('Non-supported averaging type of input field','merge_lagrange_mass_2_stack')
!!            return
!!        end select
!!
!!        if(fu_field_kind(targetId) == averaged_flag)then
!!          !
!!          ! For averaged field, if this is the last data collection before the
!!          ! output - divide the integrated stuff by the accumulation length
!!          !
!!          if(fu_valid_time(targetId) == fu_valid_time(idTemplate))then
!!            do iz = 1, nz
!!              if(fu_accumulation_length(field3d%p2d(1)%idPtr) == zero_interval)then
!!                field3d%p2d(iz)%ptr(1:fs) = 0.
!!              else
!!                field3d%p2d(iz)%ptr(1:fs) = field3d%p2d(iz)%ptr(1:fs) / &
!!                                                       & abs(fu_sec(fu_accumulation_length(field3d%p2d(1)%idPtr)))
!!              endif
!!              field3d%p2d(iz)%idPtr = targetId
!!              call set_level(field3d%p2d(iz)%idPtr, fu_level(verticalTemplate,iz))
!!              call set_validity_length(field3d%p2d(iz)%idPtr, fu_validity_length(idTemplate))
!!            end do
!!          endif
!!        endif
!!
!!      case default ! target ID
!!        do iz = 1, nz
!!          field3d%p2d(iz)%ptr(1:fs) = 0.
!!          call set_valid_time(field3d%p2d(iz)%idPtr, fu_valid_time(targetId))
!!        enddo
!!        call set_error('Non-supported averaging type of output field','merge_lagrange_mass_2_stack')
!!        return
!!    end select
!!    
!!    contains
!!    
!!    !========================================================================
!!    
!!    subroutine merge_data_arrays(ifPerVolume, ifPerDensity, ifPerDz, datOut3d, factor_in, factor_out, &
!!                               & xSize, ySize, zSize)
!!      !
!!      ! Replaces the output data vector with the stuff from mass map, possibly with 
!!      ! scaling with const factors both in and out arrays, with grid cell volume, and 
!!      ! with air density.
!!      ! outNew = outOld * factorOut + in * factorIn (/cellArea) (/cell-thickness) (/air_density)
!!      ! It serves for replacing data, for adding and subtracting them with arbitrary scaling
!!      !
!!      implicit none
!!      ! Imported parameters
!!      logical, intent(in) :: ifPerVolume, ifPerDensity, ifPerDz
!!      type(field_3d_data_ptr), intent(inout) :: datOut3d
!!      real, dimension(:), pointer :: xSize, ySize, zSize
!!      real, intent(in) :: factor_in, factor_out
!!      ! Local variables
!!      integer :: iPart, i1d, ix, iy, iz
!!      real :: cnc_air
!!
!!      if(ifPerVolume .and. ifPerDensity)then    ! normalised to cell volume and air density
!!        call set_error('Lagrange output does not supprto VMR','merge_data_arrays')
!!        return
!!!        if(ifPerDz)then
!!!          !  Full normalization: (dx * dy * dz * cnc_air)
!!!          do iPart = 1, nop
!!!            if(lpStatus(iPart) == int_missing)cycle  ! empty
!!!            if(mod(lpStatus(iPart),100) /= iSrc)cycle  ! wrong source ID
!!!            ix = nint(lpDyn(lp_x,iPart))    ! all in meteo_grid
!!!            iy = nint(lpDyn(lp_y,iPart))
!!!            iz = nint(lpDyn(lp_z,iPart))
!!!            i1d = ix+(iy-1) * nx_meteo
!!!            cnc_air = (weight_past * rho_past_p3d%p2d(iz)%ptr(i1d) + &
!!!                    & (1.0-weight_past) * rho_future_p3d%p2d(iz)%ptr(i1d))  / &
!!!                    & molecular_weight_air
!!!            datOut(i1d) = datOut(i1d) * factor_out + &
!!!                        & lpMass(iSpecies,iPart) * factor_in / &
!!!                        & (cnc_air * xSize(i1d) * ySize(i1d) * &
!!!                         & (dz_past_p3d%p2d(iz)%ptr(i1d) * weight_past + &
!!!                          & dz_future_p3d%p2d(iz)%ptr(i1d) * (1.-weight_past)))
!!!          end do
!!!        else
!!!          ! No dz: (dx * dy * cnc_air)
!!!          do iPart = 1, nop
!!!            if(lpStatus(iPart) == int_missing)cycle  ! empty
!!!            if(mod(lpStatus(iPart),100) /= iSrc)cycle  ! wrong source ID
!!!            ix = nint(lpDyn(lp_x,iPart))    ! all in meteo_grid
!!!            iy = nint(lpDyn(lp_y,iPart))
!!!            iz = nint(lpDyn(lp_z,iPart))
!!!            i1d = ix+(iy-1) * nx_meteo
!!!            cnc_air = (weight_past * rho_past_p3d%p2d(iz)%ptr(i1d) + &
!!!                    & (1.0-weight_past) * rho_future_p3d%p2d(iz)%ptr(i1d))  / &
!!!                    & molecular_weight_air
!!!            datOut(i1d) = datOut(i1d) * factor_out + &
!!!                        & lpMass(iSpecies,iPart) * factor_in / &
!!!                        & (cnc_air * xSize(i1d) * ySize(i1d))
!!!          end do
!!!        endif  ! ifPerDz
!!      elseif(ifPerVolume)then   
!!        ! normalised to cell volume only, no air density
!!        if(ifPerDz)then
!!          ! Full volume: dx*dy*dz. Note that dz is static, 1D array (1:nz)
!!          do iPart = 1, nop
!!            if(lpStatus(iPart) == int_missing)cycle  ! empty
!!            if(mod(lpStatus(iPart),100) /= iSrc)cycle  ! wrong source ID
!!            ix = nint(lpDyn(lp_x,iPart))    ! all in meteo_grid
!!            iy = nint(lpDyn(lp_y,iPart))
!!            iz = nint(lpDyn(lp_z,iPart))
!!            i1d = ix+(iy-1) * nx_meteo
!!            datOut3d%p2d(iz)%ptr(i1d) = datOut3d%p2d(iz)%ptr(i1d) * factor_out + &
!!                        & lpMass(iSpecies,iPart) * factor_in / &
!!                        & (xSize(i1d) * ySize(i1d) * zSize(iz))
!!!                         & (dz_past_p3d%p2d(iz)%ptr(i1d) * weight_past + &
!!!                          & dz_future_p3d%p2d(iz)%ptr(i1d) * (1.-weight_past)))
!!          end do
!!        else
!!          ! only area: dx*dy
!!          do iPart = 1, nop
!!            if(lpStatus(iPart) == int_missing)cycle  ! empty
!!            if(mod(lpStatus(iPart),100) /= iSrc)cycle  ! wrong source ID
!!            ix = nint(lpDyn(lp_x,iPart))    ! all in meteo_grid
!!            iy = nint(lpDyn(lp_y,iPart))
!!            iz = nint(lpDyn(lp_z,iPart))
!!            i1d = ix+(iy-1) * nx_meteo
!!            datOut3d%p2d(iz)%ptr(i1d) = datOut3d%p2d(iz)%ptr(i1d) * factor_out + &
!!                        & lpMass(iSpecies,iPart) * factor_in / &
!!                        & (xSize(i1d) * ySize(i1d))
!!          end do
!!        endif
!!      elseif(ifPerDensity)then  
!!        call set_error('Lagrange output does not support VMR','merge_data_arrays')
!!        return
!!!        ! normalized to air density, no volume
!!!        do iPart = 1, nop
!!!          if(lpStatus(iPart) == int_missing)cycle  ! empty
!!!          if(mod(lpStatus(iPart),100) /= iSrc)cycle  ! wrong source ID
!!!          ix = nint(lpDyn(lp_x,iPart))    ! all in meteo_grid
!!!          iy = nint(lpDyn(lp_y,iPart))
!!!          iz = nint(lpDyn(lp_z,iPart))
!!!          i1d = ix+(iy-1) * nx_meteo
!!!          cnc_air = (weight_past * rho_past_p3d%p2d(iz)%ptr(i1d) + &
!!!                  & (1.0-weight_past) * rho_future_p3d%p2d(iz)%ptr(i1d))  / &
!!!                  & molecular_weight_air
!!!          datOut(i1d) = datOut(i1d) * factor_out + lpMass(iSpecies,iPart) * factor_in / cnc_air
!!!        end do
!!      else      
!!        ! no normalization at all
!!        do iPart = 1, nop
!!          if(lpStatus(iPart) == int_missing)cycle  ! empty
!!          if(mod(lpStatus(iPart),100) /= iSrc)cycle  ! wrong source ID
!!          ix = nint(lpDyn(lp_x,iPart))    ! all in meteo_grid
!!          iy = nint(lpDyn(lp_y,iPart))
!!          iz = nint(lpDyn(lp_z,iPart))
!!          i1d = ix+(iy-1) * nx_meteo
!!          datOut3d%p2d(iz)%ptr(i1d) = datOut3d%p2d(iz)%ptr(i1d) * factor_out + lpMass(iSpecies,iPart) * factor_in
!!        end do
!!      endif  ! type of scaling
!!
!!    end subroutine merge_data_arrays
!!    
!!    !========================================================================
!!    
!!    subroutine report_error(chMsg, ifReportInternal)
!!      implicit none
!!      character(len=*), intent(in) :: chMsg
!!      logical, intent(in) :: ifReportInternal
!!        call msg('')
!!        call msg('Input field id:')
!!        call report(idTemplate)
!!        call msg('Desired vertical:')
!!        call report(verticalTemplate)
!!        call msg('Currently stored field id')
!!        call report(field3d%p2d(1)%idPtr)
!!        call msg('Current target field id')
!!        call report(targetId)
!!        call report(stackPtr)
!!        call set_error(chMsg,'merge_cumul_2_instant')
!!        if(ifReportInternal)then
!!          call msg('Currently stored intermediate field id')
!!          call report(idOut_internal)
!!        endif
!!    end subroutine report_error
!!
!!  end subroutine merge_lagrange_mass_2_stack

end module dispersion_server

