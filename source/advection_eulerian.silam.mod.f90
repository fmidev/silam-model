MODULE advection_eulerian

  ! Description:
  ! Dispatcher for eulerian advections
  !
  ! Authors: Mikhail Sofiev, Julius Vira, Rostislav Kouznetsov, FMI, firstname.lastname@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:

  use depositions
  use ini_boundary_conditions
  use chemistry_manager
  use silam_partitioning
!  use silam_levels
!  use advection_eulerian_v4
  use advection_eulerian_v5

  !$use OMP_LIB

!!!#define DEBUG

  IMPLICIT NONE

  private

  public InitEulerAdvectionFields
  public advect_eulerian_cloud
  public fu_if_bulk_eulerian_advection
  public euler_adv_input_needs

  public  advect_rect
  public  advect_tri 
  public  advect_step

  logical, private, save :: ifXFirst = .true.

 


CONTAINS

  ! ****************************************************************

  subroutine euler_adv_input_needs(adv_method, ifMakeMassfluxes, &
                                 & q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    ! Include the meteorological quantities needed by the eulerian
    ! advection methods.
    implicit none
    INTEGER, INTENT(in) :: adv_method
    logical, intent(in) :: ifMakeMassfluxes
    INTEGER, DIMENSION(:), INTENT(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st
    ! Local variables
    integer :: iTmp
     
    SELECT CASE(adv_method)

      case(no_advection)
      
      CASE(adv_euler_Galperin_3d_bulk)
        iTmp = fu_merge_integer_to_array(Kz_scalar_1m_flag, q_met_dynamic)
        iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dynamic)

!      case(adv_euler_Galperin_v4)
!        call  euler_adv_v4_input_needs( ifIncludeVerticalDependentFlds, &
!                        & ifMakeMassfluxes, &
!                        & q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
      case(adv_euler_Galperin_v5)
        call  euler_adv_v5_input_needs(q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
      CASE DEFAULT
        CALL set_error('Unknown advection method:' + fu_str(adv_method),'euler_adv_input_needs')
        RETURN
    END SELECT

   end subroutine euler_adv_input_needs


  !********************************************************************************

  subroutine InitEulerAdvectionFields(adv_method, adv_variant, smoother_factor, &
      & nsources, nspecies,  npassengers, ifMolecDiffusion, ifSubgridDiffusion)
    !
    ! Calls  InitEulerAdvectionFields_vX.
    !
    implicit none

    ! Imported parameter
    integer, intent(in) :: adv_method, adv_variant
    integer, intent(in) :: nsources, nspecies, npassengers
    real, intent(in) :: smoother_factor
    logical, intent(in) :: ifMolecDiffusion, ifSubgridDiffusion

    ! Local variables
    integer :: iLev,iStatus, nthreads, nxy, nz, ithread, nPoles
    type (silja_level), dimension(max_levels) :: leveltops
    type (silam_vertical) :: vertdispTop
    
    nthreads = 1
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(nthreads)
    !$OMP MASTER
    !$ nthreads = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL
    !
    ! Size of arrays for horizontal advection
    ! If MPI, we have to account for a few extra lines for each MPI border (public from partitioning)
    nxy = max(nx_dispersion + wing_depth_e + wing_depth_w, &
        &     ny_dispersion + wing_depth_n + wing_depth_s) + 1
    !
    ! Stupidity check.
    !
    if(nx_dispersion < 2 .or. ny_dispersion < 2 .or. nx_dispersion > 2000 .or. ny_dispersion > 2000)then
      call msg('nx and ny of dispersion grid:',nx_dispersion, ny_dispersion)
      call msg_warning('Strange dispersion grid parameters','InitEulerAdvectionFields')
    endif
    
    call msg ("InitEulerAdvectionFields method:" + fu_name(adv_method))

    SELECT CASE(adv_method)

      case(no_advection)

      CASE(adv_euler_Galperin_3d_bulk)
        call set_error('bulk advection is broken','InitEulerAdvectionFields')
      case(adv_euler_Galperin_v5)
         call InitEulerAdvectionFields_v5(nsources, nspecies,  npassengers, nxy, nz_dispersion, & 
           & nthreads, adv_variant, smoother_factor, ifMolecDiffusion, ifSubgridDiffusion)
      CASE DEFAULT
        CALL set_error('Unknown advection method:' + fu_str(adv_method),'InitEulerAdvectionFields')
    END SELECT
  end subroutine InitEulerAdvectionFields

  logical function fu_if_bulk_eulerian_advection(which_type) result(if_bulk) 
    ! Check if we have specieswise or bulk advection. This information
    ! is needed for the emission and boundary conditions.
    implicit none
    integer, intent(in) :: which_type
    
    if_bulk = (which_type == adv_euler_Galperin_3d_bulk) 
  end function fu_if_bulk_eulerian_advection


  !===============================================================================
  !===============================================================================
  !===============================================================================
  !
  !    Eulerian advection.
  !
  !===============================================================================
  !===============================================================================
  !===============================================================================

  !*********************************************************************************

  subroutine advect_eulerian_cloud(advection_method, &
                                 & mapConc, mapPx_conc, mapPy_conc, mapPz_conc, mapAerosol, &
                                 & interpCoefMeteo2DispHoriz, interpCoefMeteo2DispVert, &
                                 & ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp, &
                                 & IniBoundaryRules, pBoundaryBuffer, &
                                 & met_buf, disp_buf, &
                                 & now, timestep_sec, weight_past, &
                                 & garbage_mass, mapDryDep, mapCnc2m, &
                                 & x0_mass, xM_mass, y0_mass, yM_mass, &
                                 & bottom_mass, top_mass, &
                                 & chem_rules, have_negatives, wdr)
    !
    ! Depending on the type of advection requested, selects the appropriate subroutine(s)
    !
    implicit none
    !
    ! Imported parameters
    !
    integer, intent(in) :: advection_method
    type(TMass_map), intent(inout) :: mapConc, mapPx_conc, mapPy_conc, mapPz_conc, mapAerosol, &
                              & mapDryDep, mapCnc2m
    type(THorizInterpStruct), pointer :: interpCoefMeteo2DispHoriz
    type(TVertInterpStruct), pointer :: interpCoefMeteo2DispVert
    logical, intent(in) :: ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp
    type(Tini_boundary_rules), intent(in) :: IniBoundaryRules
    type(TboundaryBuffer), pointer :: pBoundaryBuffer
    real, dimension(:,:), pointer :: garbage_mass                          ! (nSrc,nSpecies)
    real, dimension(:,:,:,:), pointer :: x0_mass, xM_mass, y0_mass, yM_mass  ! (nSrc,nSpecies,2,nz)
    real, dimension(:,:,:), pointer :: bottom_mass, top_mass !  (nSrc,nSpecies,2)
    TYPE(silja_time), INTENT(in) :: now
    real, intent(in) :: timestep_sec, weight_past
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    type(Tchem_rules), intent(in) :: chem_rules
    logical, intent(in) :: have_negatives
    TYPE(silja_wdr), INTENT(in) :: wdr
    logical :: ifMoments
    logical, save :: ifFirst = .true.

    ! Local declarations:
    character(len=20) :: counter_name
    integer :: iy

    !
    ! The main selector of the advection type
    !
    ! Vertical interpolation coefficients are already refined in
    ! dispersion supplementary along with computation of the dispersion vertical wind.
    !

    if (debug_level > 0) call msg('Sum of mapConc before transport', sum(mapConc%arm))


    ifMoments = .false. ! Centers of mass in mass map


    select case(advection_method)

      case(no_advection)
        return
      
      case(adv_euler_Galperin_3D_bulk)
        !
        ! Bulk 3D advection: all 3 axes are handled in this sub
        !
         call set_error('3D bulk advection is broken. Revive or use v4','advect_eulerian_cloud')



      case(adv_euler_Galperin_v5)
        !
        ! Semi-analytic advection, version 5
        !

        ifMoments = .false. ! Centers of mass in mass map
        if (.not. ifXfirst) then
           counter_name = 'Vertical transport'
           call start_count(chCounterNm=counter_name)
           call adv_diffusion_vertical_v5(mapConc, &
                             & mapPx_conc, mapPy_conc, mapPz_conc, &
                             & mapAerosol, &
                             & interpCoefMeteo2DispHoriz, interpCoefMeteo2DispVert, &
                             & ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp, &
                             & met_buf, disp_buf, &
                             & timestep_sec, weight_past, &
                             & garbage_mass, &
                             & bottom_mass, top_mass, &
                             & mapDryDep, mapCnc2m, &
                             & pBoundaryBuffer, chem_rules, &
                             & have_negatives, wdr, ifMoments, ifXfirst, &
                             & .false.)   !now >= fu_set_time_utc(1980,2,11,8,0,0.))
            if (error) return
            call stop_count(chCounterNm=counter_name)
        endif

        !Some estimate of load imbalance
        if (smpi_global_tasks > 1) then
          counter_name = 'Before Hadv barrier'
          call start_count(chCounterNm = counter_name)
          call smpi_advection_barrier()
          call stop_count(chCounterNm = counter_name)
        endif

        counter_name = 'Horizontal transport'
        call start_count(chCounterNm=counter_name)
        call adv_euler_Galp_xy_v5(mapConc, mapPx_conc, &
                                & mapPy_conc, mapPz_conc, &
                                & disp_buf, &
                                & timestep_sec, weight_past, &
                                & garbage_mass, &
                                & x0_mass, xM_mass, y0_mass, yM_mass, &
                                & chem_rules, ifMoments, have_negatives, &
                                & pBoundaryBuffer, &
                                & ifXfirst)
        call stop_count(chCounterNm=counter_name)
        if (error) return
        if (ifXfirst) then 
           counter_name = 'Vertical transport'
           call start_count(chCounterNm=counter_name)
           call adv_diffusion_vertical_v5(mapConc, &
                             & mapPx_conc, mapPy_conc, mapPz_conc, &
                             & mapAerosol, &
                             & interpCoefMeteo2DispHoriz, interpCoefMeteo2DispVert, &
                             & ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp, &
                             & met_buf, disp_buf, &
                             & timestep_sec, weight_past, &
                             & garbage_mass, &
                             & bottom_mass, top_mass, &
                             & mapDryDep, mapCnc2m, &
                             & pBoundaryBuffer, chem_rules, &
                             & have_negatives, wdr, ifMoments, ifXfirst, &
                             & .false.)    !now >= fu_set_time_utc(1980,2,11,7,0,0.))
            if (error) return
            call stop_count(chCounterNm=counter_name)
        endif
#ifdef DEBUG                
       if(any(pBoundaryBuffer%iBoundaryType(1:5) /= zero_boundary_type))then
         call report_inout_mass_stuff_v5(incoming, mapConc, pBoundaryBuffer,.false.) !Not by levels
         call report_inout_mass_stuff_v5(outgoing, mapConc, pBoundaryBuffer,.false.) 
       endif
#endif

      case default
        call set_error('Unknown advection method:' + &
                     & fu_name(advection_method), 'advect_eulerian_cloud')
    end select  !  advection


    if (debug_level > 0) call msg('Sum of mapConc after transport', sum(mapConc%arm))

    if (ifMoments) then
       call flip_cm(mapConc, mapPx_conc, mapPy_conc, mapPz_conc, have_negatives)
       ifMoments = .false.
    endif

   ifXfirst = .not. ifXfirst
        
   ifFirst = .false.


  contains
    
    subroutine flip_cm(map_c, map_m1,  map_m2,  map_m3, have_negatives)
       !Moments to CM
      implicit none
      type(Tmass_map), intent(inout) :: map_c, map_m1,  map_m2,  map_m3
      logical, intent(in) :: have_negatives
      integer :: ix,iy,iz,iSrc,isp
      real :: fTmp

      !$OMP PARALLEL default(shared) private(ix,iy,iz,iSrc,isp, ftmp)
      !$OMP DO collapse(4)
      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          do iz = 1,nz_dispersion
           do isrc= 1, map_c%nSrc
            do isp = 1, map_c%nSpecies
              if (abs(map_c%arm(isp,iSrc,iz,ix,iy)) > 0.) then 
                 fTmp = map_m1%arm(isp,iSrc,iz,ix,iy) / map_c%arm(isp,iSrc,iz,ix,iy)
                 if (have_negatives) fTmp = max(0.5,min(-0.5, fTmp))
                 map_m1%arm(isp,iSrc,iz,ix,iy) = fTmp

                 fTmp = map_m2%arm(isp,iSrc,iz,ix,iy) / map_c%arm(isp,iSrc,iz,ix,iy)
                 if (have_negatives) fTmp = max(0.5,min(-0.5, fTmp))
                 map_m2%arm(isp,iSrc,iz,ix,iy) = fTmp 

                 fTmp = map_m3%arm(isp,iSrc,iz,ix,iy) / map_c%arm(isp,iSrc,iz,ix,iy)
                 if (have_negatives) fTmp = max(0.5,min(-0.5, fTmp))
                 map_m3%arm(isp,iSrc,iz,ix,iy) = fTmp 
               else
                  map_m1%arm(isp,iSrc,iz,ix,iy) = 0.
                  map_m2%arm(isp,iSrc,iz,ix,iy) = 0.
                  map_m3%arm(isp,iSrc,iz,ix,iy) = 0.
               endif
            end do
           end do 
          end do
        end do
      end do
      !$OMP ENDDO
      !$OMP END PARALLEL 

    end subroutine flip_cm


  end subroutine advect_eulerian_cloud

  

END MODULE advection_eulerian

