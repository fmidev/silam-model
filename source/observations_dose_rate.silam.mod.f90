module observations_dose_rate
  ! This module is built over observations_in_situ by Sebastian Heinonen
  ! The following comments are for the observations_in_situ and should in some time get replaced:
  ! 
  ! This module defines subroutines and datatypes for time-series type observations of
  ! either in-situ or column-integral type. The observations act on the model fields
  ! either in forward (observe) or adjoint (inject) mode. In the forward mode, the
  ! modelled observations are stored, in the adjoint mode, the discrepancy y - Hx between
  ! modelled and real observations is injected as a forcing.
  !
  ! A numerical trick: the injected mass is multiplied with a number >>1 and divided when
  ! the gradient is collected. This is because the adjoint runs have no source terms, and
  ! hence the much lower concentrations would end up negligible compared to the low-mass
  ! threshold. This is needed even though the observations and variances given in units
  ! /m3 or /m2 are converted into cell or cell-area integrated values.
  !
  ! Sebastian Heinonen 2017
  
  use dispersion_server
  use optical_density
  use observations_in_situ

  implicit none

  public inject_dose_rate
  public observe_dose_rate
  public fu_init_dose_rate_obs
  public fu_init_dose_rate_addition
  public test_dose_rate
  public restart_dose_rate

  ! PRIVATE subroutines start here
  private fu_get_dose_rate_at_station
  private fu_get_relative_air_dose
  private fu_get_ground_air_ratio
  private fu_get_relative_ground_dose
  private cart_to_sph
  private bruteforce
  private cloudshine_laguerre
  private cloudshine_legendre
  private legendre30_params
  private legendre_params
  private laguerre_params
  private cerc_params

  
  !***********************************************************************************************

  type t_dose_rate_obs_addition
     !Type to give additional info needed to calculate dose rate.
     !Points to observation.
     type(inSituObservation), pointer :: obs => null()
     !Obs_map contains the weights for in function of space and species.
     real, dimension(:,:,:,:), allocatable :: obs_map          !x,y,z,spc
     !Relative dose has fraction of dose in space and species (total is 1 each ind_obs).
     real, dimension(:,:,:,:,:), allocatable :: relative_dose  !x,y,z,spc,obs_ind
     !Groundshine_multiplier contains the value that has to be multiplied by 
     !the concentration in area to get the dose from the radiation from ground.
     real, dimension(:), allocatable :: groundshine_multiplier !spc
     !Cloud_ground_ratio gives the percentage of total dose from groundshine (rest from cloudshine).
     !Relative wet and dry deposition doses
     real, dimension(:,:), allocatable :: cloud_ground_ratio, &  !wet,dry;ind_obs => (2,n_obs)
                                        & relative_wet,relative_dry !spc,ind_obs
     !The minimum and maximum indices that define the size in x, y and z directions of the obs_map etc.
     integer :: x_min,x_max,y_min,y_max,z_min,z_max
     !Include ground shine
     logical :: ifGroundToo
     !Defined when initialized.
     logical :: defined = .false.
  end type t_dose_rate_obs_addition

contains
  
  
  !***********************************************************************************************
  !
  !     Initializations
  !
  !***********************************************************************************************
  
  
  ! Initialization. Returns a new observation object with given data,
  ! locations, and time windows. The returned object is ready for use
  ! at any time.  Array pointers endTimes and durations should be
  ! allocated and dipersionGrid should be defined in advance.
  ! For dose rate only species that have gamma radiation will be used.

  function fu_init_dose_rate_addition(ifGroundToo, obs, vertical, transport_species) result(a)
    !called when observation and station are already initialized
    implicit none
    !input
    logical, intent(in) :: ifGroundToo
    type(inSituObservation), intent(in), target :: obs
    type(silam_vertical), intent(in) :: vertical
    type(silam_species), dimension(:), intent(in) :: transport_species
    
    !return initialized dose rate observation addition
    type(t_dose_rate_obs_addition) :: a
    
    !local
    real :: x,y
    real :: Rmax, deltax, deltay, deltah, dx,dy,dh, xpos, ypos, fract
    real, dimension(6,14) :: cerc
    real, dimension(14) :: exp_par
    real, dimension(:), allocatable :: heights
    real, dimension(:,:,:), allocatable :: distances
    real, dimension(:,:,:,:), allocatable :: mapping
    real, dimension(:,:,:,:), allocatable :: obs_map_tabE
    real, parameter, dimension(13) :: CERC_energies = (/12.5, 17.5, 25.0, 40.0, 57.5, 82.5, 150.0, &
                                                      & 350.0, 750.0, 1250.0, 1750.0, 3000.0, 7000.0/)
    integer :: nx, ny, offx, offy, gnx, gny, i, j, k, nEnergy, n, x0, y0
    integer :: ix_min, ix_max, jy_min, jy_max, kz_min, kz_max
    integer :: isp_obs, isp_transp, ilev
    integer :: ispecies_transp, ispecies_obs, iE, itabE, iCercEne
    !local pointers
    real(r4k), dimension(:,:), pointer :: pEnergies
    type(silam_material), pointer :: material
    type(silam_nuclide), pointer :: nuclide
    
    
    a%obs => obs                   !appoint observation
    a%ifGroundToo = ifGroundToo    !define if groundshine is included
    
    xpos = a%obs%station%xLoc_dispGrid      ! local variables defined for xpos,
    ypos = a%obs%station%yLoc_dispGrid      ! ypos and
    ilev = a%obs%station%ind_lev   ! ilev to make the code more readable
    
    !smpi decomposition subroutine brings the necessary offsets plus global and local maximum indices
    call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
    !project point to grid brings necessary wholeMPIdispersiongrid and indices in it
    call project_point_to_grid(a%obs%station%lon, a%obs%station%lat, &
                              & wholeMPIdispersion_grid, x, y)
    y0 = nint(y) !center indices into nearest integers
    x0 = nint(x)
    
    cerc = cerc_params()      !downloads the params from cerc paper
    exp_par = exp_int_par()   ! and for exponential integral values
    nEnergy = size(cerc(1,:)) !number of energies
    Rmax = -log(0.001)/(minval(cerc(2,:))) !area chosen so that the smallest attennuation is diminished in that area into 0.1%
    
    
    
    
    
    !******************************************************************************************
    !extent inside local grid in indices searched so that all grid points are inside Rmax in all directions
    !X
    !The loops assume that indices are allowed without wholeMpiDispersionGrid reductions (offx and offy)
    ! to vary from 1 to nx/ny and with reductions from 1 to gnx/gny
    if (xpos < 0.5) then !X AXIS
       !check the negative x axis reach
       deltax = -xpos & 
              & * fu_dx_cell_m(wholeMPIdispersion_grid, x0, y0)
       i = 0
       !Find the border (similar for all directions..): either cell center is further than Rmax,
       !index out of global indices OR out of local MPI indices
       do while ((deltax > -Rmax) .and. (x0-i >= 1) .and. (x0-i-offx >= 1))
          i = i + 1
          deltax = deltax - 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0-i, y0) &
                 & - 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0-i-1, y0)
       end do
       ix_min = x0-i
       
       !check the positive x axis reach
       deltax = (0.5-xpos) & 
              & * fu_dx_cell_m(wholeMPIdispersion_grid, x0, y0) &
              & + 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0+1, y0)
       i = 1
       do while ((deltax < Rmax) .and. (x0+i <= gnx) .and. (x0+i-offx <= nx))
          i = i + 1
          deltax = deltax + 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0+i, y0) &
                 & + 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0+i-1, y0)
       end do
       ix_max = x0+i
    else
       !check the negative x axis reach
       deltax = -(xpos-0.5) & 
              & * fu_dx_cell_m(wholeMPIdispersion_grid, x0, y0) &
              & - 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0-1, y0)
       i = 0
       do while ((deltax > -Rmax) .and. (x0-i >= 1) .and. (x0-i-offx >= 1))
          i = i + 1
          deltax = deltax - 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0-i, y0) &
                 & - 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0-i-1, y0)
       end do
       ix_min = x0-i+1
       
       !check the positive x axis reach
       deltax = (1.0-xpos) & 
              & * fu_dx_cell_m(wholeMPIdispersion_grid, x0, y0)
       i = 1
       do while ((deltax < Rmax) .and. (x0+i <= gnx) .and. (x0+i-offx <= nx))
          i = i + 1
          deltax = deltax + 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0+i, y0) &
                 & + 0.5*fu_dx_cell_m(wholeMPIdispersion_grid, x0+i-1, y0)
       end do
       ix_max = x0+i
    end if
    !Y
    if (ypos < 0.5) then !Y AXIS
       !check the negative y axis reach
       deltay = -ypos & 
              & *fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0)
       j = 0
       do while ((deltay > -Rmax) .and. (y0-j >= 1) .and. (y0-j-offy >= 1))
          j = j + 1
          deltay = deltay - 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0-j) &
                 & - 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0-j-1)
       end do
       jy_min = y0-j
       
       !check the positive y axis reach
       deltay = (0.5-ypos) & 
              & * fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0) &
              & + 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0+1)
       j = 1
       do while ((deltay < Rmax) .and. (y0+j <= gny) .and. (y0+j-offy <= ny))
          j = j + 1
          deltay = deltay + 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0+j) &
                 & + 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0+j-1)
       end do
       jy_max = y0+j
    else
       !check the negative y axis reach
       deltay = -(ypos-0.5) & 
              & * fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0) &
              & - 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0-1)
       j = 0
       do while ((deltay > -Rmax) .and. (y0-j >= 1) .and. (y0-j-offy >= 1))
          j = j + 1
          deltay = deltay - 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0-j) &
                 & - 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0-j-1)
       end do
       jy_min = y0-j+1
       
       !check the positive y axis reach
       deltay = (1.0-ypos) & 
              & * fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0)
       j = 1
       do while ((deltay < Rmax) .and. (y0+j <= gny) .and. (y0+j-offy <= ny))
          j = j + 1
          deltay = deltay + 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0+j) &
                 & + 0.5*fu_dy_cell_m(wholeMPIdispersion_grid, x0, y0+j-1)
       end do
       jy_max = y0+j
    end if
    !Z
    !check negative z axis reach
    if (ilev > 1) then
       k = 1
       deltah = -fu_layer_thickness_m(fu_level(vertical, ilev - k))
       do while ((deltah > -Rmax) .and. (ilev - (k+1) > 0))
          k = k + 1
          deltah = deltah - fu_layer_thickness_m(fu_level(vertical, ilev - k))
       end do
       kz_min = ilev - k
    else
       kz_min = ilev
    end if
    
    !check the positive z axis reach
    deltah = fu_layer_thickness_m(fu_level(vertical, ilev))
    k = 0
    do while (deltah < Rmax) !error if Rmax extends over the top of silam. assumed to be impossible
       k = k + 1
       deltah = deltah + fu_layer_thickness_m(fu_level(vertical, ilev + k))
    end do
    kz_max = ilev + k
    !end find extent
    !***************************************************************************************************
    
    
    
    
    !store the indices with reduced grid offsets
    a%x_min = ix_min - offx
    a%x_max = ix_max - offx
    a%y_min = jy_min - offy
    a%y_max = jy_max - offy
    a%z_min = kz_min
    a%z_max = kz_max
    
    
    !***************************************************************************************************
    !ALLOCATIONS
    !Allocating empty layer around the physical reach
    allocate(distances(a%x_min-1:a%x_max+1,a%y_min-1:a%y_max+1,2))
    allocate(heights(a%z_min-1:a%z_max+1))
    allocate(mapping(a%x_min-1:a%x_max+1,a%y_min-1:a%y_max+1,a%z_min-1:a%z_max+1,5)) !weight, distx, disty, distz, vol
    allocate(obs_map_tabE(a%x_min-1:a%x_max+1, &      !x
                        & a%y_min-1:a%y_max+1, &      !y
                        & a%z_min-1:a%z_max+1, &      !z
                        & nEnergy))                   !energy
    allocate(a%obs_map(a%x_min:a%x_max, &             !x
                     & a%y_min:a%y_max, &             !y
                     & a%z_min:a%z_max, &             !z
                     & a%obs%num_obs_species))        !spc
    allocate(a%groundshine_multiplier(a%obs%num_obs_species))
    allocate(a%relative_dose(a%x_min:a%x_max, &       !x
                           & a%y_min:a%y_max, &       !y
                           & a%z_min:a%z_max, &       !z
                           & a%obs%num_obs_species, & !spc
                           & size(a%obs%obsData)))       !ind_obs
    allocate(a%cloud_ground_ratio(2,size(a%obs%obsData)))!ind_obs
    allocate(a%relative_wet(a%obs%num_obs_species, &  !spc
                          & size(a%obs%obsData)))        !ind_obs
    allocate(a%relative_dry(a%obs%num_obs_species, &  !spc
                          & size(a%obs%obsData)))        !ind_obs
    
    !END ALLOCATIONS
    !***************************************************************************************
    
    
    
    
    !complimentary loop for the next one
    do i = a%x_min-1,a%x_max+1
       do j = a%y_min-1,a%y_max+1
          distances(i,j,1) = fu_dx_cell_m(wholeMPIdispersion_grid, i, j)
          distances(i,j,2) = fu_dy_cell_m(wholeMPIdispersion_grid, i, j)
       end do
    end do
    do k = a%z_min-1,a%z_max+1
       if (k == a%z_min-1) then
          heights(k) = fu_layer_thickness_m(fu_level(vertical, k+1))
       else
          heights(k) = fu_layer_thickness_m(fu_level(vertical, k))
       end if
    end do
    !uses distances defined above, thickness and cell dx and dy in m to define observation mapping
    !that has distances of needed grid points from the station center
    do i = a%x_min-1,a%x_max+1
       do j = a%y_min-1,a%y_max+1
          do k = a%z_min-1,a%z_max+1
             if (i == a%x_min-1 .or. i == a%x_max+1 .or. &
               & j == a%y_min-1 .or. j == a%y_max+1 .or. &
               & k == a%z_min-1 .or. k == a%z_max+1) then
                !WEIGHT OF EACH POINT INITIALIZED AS 1.0 except in the outer layer it is 0.0
                mapping(i,j,k,1) = 0.0
             else
                mapping(i,j,k,1) = 1.0
             end if
             !POSITION IN X AXIS
             if (xpos < 0.5) then
                if (i <= x0) then
                   if (i == x0) then
                      mapping(i,j,k,2) = -xpos &
                                       & * distances(i,j,1)
                   else
                      mapping(i,j,k,2) = -xpos &
                                       & * distances(x0,j,1) &
                                       & - 0.5*sum(distances(i+1:y0,j,1)) &
                                       & - 0.5*sum(distances(i:(x0-1),j,1))
                   end if
                else !i > x0
                   if (i == x0+1) then
                      mapping(i,j,k,2) = 0.5 * distances(i,j,1) &
                                       & + (0.5-xpos) &
                                       & * distances(i-1,j,1)
                   else
                      mapping(i,j,k,2) = 0.5 * distances((x0+1),j,1) &
                                       & + (0.5-xpos) &
                                       & * distances(x0,j,1) &
                                       & + 0.5*sum(distances((x0+1):i-1,j,1)) &
                                       & + 0.5*sum(distances((x0+2):i,j,1))
                   end if
                end if
             else ! xpos > 0.5
                if (i <= x0) then
                   if (i == x0) then
                      mapping(i,j,k,2) = -0.5 * distances(i,j,1) &
                                       & - (xpos - 0.5) & 
                                       & * distances(i+1,j,1)
                   else
                      mapping(i,j,k,2) = -0.5 * distances(x0,j,1) &
                                       & - (xpos - 0.5) &
                                       & * distances(x0+1,j,1) &
                                       & - 0.5*sum(distances((i+1):x0,j,1)) &
                                       & - 0.5*sum(distances(i:(x0-1),j,1))
                   end if
                else !i > x0
                   if (i == x0+1) then
                      mapping(i,j,k,2) = (1.0 - xpos) &
                                       & * distances(i,j,1)
                   else
                      mapping(i,j,k,2) = (1.0 - xpos) &
                                       & * distances(x0+1,j,1) &
                                       & + 0.5*sum(distances((x0+1):i-1,j,1)) &
                                       & + 0.5*sum(distances((x0+2):i,j,1))
                   end if
                end if
             end if
             !POSITION IN Y AXIS
             if (ypos < 0.5) then
                if (j <= y0) then
                   if (j == y0) then
                      mapping(i,j,k,3) = -ypos & 
                                       & * distances(i,j,2)
                   else
                      mapping(i,j,k,3) = -ypos & 
                                       & * distances(i,y0,2) &
                                       & - 0.5*sum(distances(i,j+1:y0,2)) &
                                       & - 0.5*sum(distances(i,j:(y0-1),2))
                   end if
                else !j > y0
                   if (j == y0+1) then
                      mapping(i,j,k,3) = 0.5 * distances(i,j,2) &
                                       & + (0.5-ypos) & 
                                       & * distances(i,j-1,2)
                   else
                      mapping(i,j,k,3) = 0.5 * distances(i,(y0+1),2) &
                                       & + (0.5-ypos) & 
                                       & * distances(i,y0,2) &
                                       & + 0.5*sum(distances(i,(y0+1):j-1,2)) &
                                       & + 0.5*sum(distances(i,(y0+2):j,2))
                   end if
                end if
             else ! ypos > 0.5
                if (j <= y0) then
                   if (j == y0) then
                      mapping(i,j,k,3) = -0.5 * distances(i,j,2) &
                                       & - (ypos - 0.5) & 
                                       & * distances(i,j+1,2)
                   else
                      mapping(i,j,k,3) = -0.5 * distances(i,y0,2) &
                                       & - (ypos - 0.5) & 
                                       & * distances(i,y0+1,2) &
                                       & - 0.5*sum(distances(i,j+1:y0,2)) &
                                       & - 0.5*sum(distances(i,j:(y0-1),2))
                   end if
                else !j > y0
                   if (j == y0+1) then
                      mapping(i,j,k,3) = (1.0 - ypos) &
                                       & * distances(i,j,2)
                   else
                      mapping(i,j,k,3) = (1.0 - ypos) & 
                                       & * distances(i,y0+1,2) &
                                       & + 0.5*sum(distances(i,(y0+1):j-1,2)) &
                                       & + 0.5*sum(distances(i,(y0+2):j,2))
                   end if
                end if
             end if
             !NEGATIVE IN Z AXIS
             if (k == ilev) then
                mapping(i,j,k,4) = -0.5*heights(ilev)
             else if (k < ilev) then
                mapping(i,j,k,4) = -0.5*sum(heights(k:ilev)) &
                                 & - 0.5*sum(heights(k:ilev-1))
             else if (k == ilev+1) then !POSITIVE Z
                mapping(i,j,k,4) = 0.5*heights(ilev)
             else !positive larger that 1
                mapping(i,j,k,4) = 0.5*sum(heights(ilev:k)) &
                                 & + 0.5*sum(heights(ilev+1:k))
             end if
             !VOLUME OF EACH POINT
             mapping(i,j,k,5) = distances(i,j,1) * distances(i,j,2) * heights(k)
          end do
       end do
    end do
    !from cartesian coordinates just got, use function to get the corresponding spherical coordinates
    mapping(:,:,:,2:4) = cart_to_sph(mapping(:,:,:,2:4))
    if (minval(mapping(:,:,:,5)) == 0.0) call set_error('Volume of cell zero','fu_init_dose_rate_addition')
    if (error) return !Cell volume cant be 0
    do n = 1,nEnergy
       obs_map_tabE(:,:,:,n) &
              & = cloudshine_laguerre(mapping(:,:,:,:4),n) & !get the weights from laguerre approximate polynomial integral for all energies
                & * cerc(5,n) * cerc(6,n) &                  !multiplied by energy, absorption coeff and Sv/Gy =dose/amount of rad coeff to get dose rate
                & / mapping(:,:,:,5)                         !divided by the volume of the cell
    end do
    
    do ispecies_obs = 1,a%obs%num_obs_species
       ispecies_transp = a%obs%ind_obs2transp(ispecies_obs)
       fract = a%obs%scale_transp2obs(ispecies_obs)
       material => fu_material(transport_species(ispecies_transp)) !function fu_material defined in chemical setup
       if (error) return
       nuclide => fu_nuclide(material)
       if (error) return !only nuclides
       pEnergies => fu_decay_enegries(nuclide) !spell_error in function name
       do iE = 1,size(pEnergies(1,:))
          itabE = 14
          do iCercEne = 1, size(CERC_energies)
             if(pEnergies(1,iE) .lt. CERC_energies(iCercEne))then
                itabE = iCercEne
                exit
             endif
          enddo
          !to get map in species
          a%obs_map(:,:,:,ispecies_obs) = obs_map_tabE(a%x_min:a%x_max, & !x
                                                     & a%y_min:a%y_max, & !y
                                                     & a%z_min:a%z_max, & !z
                                                     & itabE) &           !tabulated energy
                                      & * fract &                         !unit correction 
                                      & * pEnergies(2,iE)                 !emitting probability
          !assumed that deposition is in units of Bq/m2 and measurement device is 2m above ground
          !http://cdn.intechweb.org/pdfs/21651.pdf eq. 16
          a%groundshine_multiplier(ispecies_obs) = exp_par(itabE) &
                                               & * pEnergies(2,iE) &
                                               & * cerc(5,itabE) &
                                               & * cerc(6,itabE) &
                                               & / 2.0
       end do
    end do
    !Deallocate temporary memories
    deallocate(distances,heights,mapping,obs_map_tabE)
    
    if (any(isnan(a%groundshine_multiplier))) call set_error('Multiplier cant be NaN','fu_init_dose_rate_addition')
    if (any(isnan(a%obs_map))) call set_error('Obs map cant be NaN','fu_init_dose_rate_addition')
    if (any(a%obs_map < 0.0)) call set_error('Obs map cant be negative','fu_init_dose_rate_addition')
    
    !initialize the undefined to zero
    a%relative_dose = 0.0
    a%cloud_ground_ratio = 0.0
    a%relative_wet = 0.0
    a%relative_dry = 0.0
    
    !finally the addition is defined
    a%defined = .true.
  end function fu_init_dose_rate_addition

  
  !******************************************************************************************

  function fu_init_dose_rate_obs(endTimes, durations, obsData, variance, &
                                  & n_values, variance_flag, &
                                  & transport_species, &
                                  & obs_unit, station, &
                                  & level, vertical, tag) &
                                  & result(newObservation)
    implicit none
    ! input
    type(silja_time), dimension(:), intent(in)  :: endTimes
    type(silja_interval), dimension(:), intent(in) :: durations
    real, dimension(:), intent(in) :: obsData, variance
    integer, intent(in) :: variance_flag, n_values
    ! the durations are currently assumed to be given in hours. True,
    ! silja_interval would make sense in this context.
    type(observationStation), intent(in), target:: station
    type(silam_species), dimension(:), intent(in) :: transport_species
    character(len=*), intent(in) :: obs_unit
    type(silam_vertical), intent(in) :: vertical
    type(silja_level), intent(in) :: level
    character(len=*), intent(in) :: tag
    
    type(inSituObservation) :: newObservation

    ! local
    integer :: status, nx, ny, i, isp_obs, no_gamma_spc
    real :: x,y, weight, conversion ! x and y are in grid coordinates
    type(silam_material), pointer :: material
    type(silam_nuclide), pointer :: nuclide
    real(r4k), dimension(:,:), pointer :: pEnergies
    type(chemical_adaptor) :: adaptor
    logical :: isDuplicate
    integer, dimension(:), allocatable :: obs2trn
    real, dimension(:), allocatable :: scales

    allocate(newobservation%endTimes(n_values), newObservation%durations(n_values), &
            & newObservation%obsData(n_values), newObservation%variance(n_values), stat=status)
    if (status /= 0) then
      call set_error('Allocate failed', 'fu_init_dose_rate_obs')
      return
    end if

    newObservation%tag = tag
    newObservation%endTimes = endTimes(1:n_values)
    newObservation%durations = durations(1:n_values)
    newObservation%obsData = obsData(1:n_values)
    newObservation%dataLength = n_values
    newObservation%station => station
    allocate(newObservation%modelData(n_values), stat=status)
    if (status /= 0) then
      call set_error('Cannot initialize observation: allocate failed.','fu_init_dose_rate_obs')
      return
    end if
    newObservation%modelData = 0.0

    select case(variance_flag)
    case (variable_variance)
      newObservation%variance = variance(1:n_values)
    case (constant_variance)
      newObservation%variance = variance(1)
    case default
      call set_error('Strange variance_flag', 'fu_init_dose_rate_obs')
      return
    end select

    ! Chemicals. Have a list of contributing species + unit (ie. kg).
    ! - in this context we are looking at all transport species and picking nuclides that emit gamma from them
    if (error) return
    
    !First allocate all species as observed
    allocate(obs2trn(size(transport_species)), &
           & scales(size(transport_species)), stat=status)
    if (fu_fails(status == 0, 'Allocate failed', 'fu_init_dose_rate_obs')) return
    
    no_gamma_spc = 0
    do isp_obs = 1, size(transport_species)
       isDuplicate = .false.
       if (isp_obs > 1) then !Not for the first spc.
          do i = 1,isp_obs-1 !Loop over species that have already gone..
             if (transport_species(i) == transport_species(isp_obs)) then !..then if species is the same..
                if (fu_material(transport_species(i)) == fu_material(transport_species(isp_obs))) then !..and if the material is the same..
                   isDuplicate = .true. !..then this species is a duplicate..
                end if
             end if
          end do
       end if
       if (isDuplicate) return !..and thus skipped.
       material => fu_material(transport_species(isp_obs))
       if (error) return
       nuclide => fu_nuclide(material)
       if (error) return !only nuclides
       pEnergies => fu_decay_enegries(nuclide)
       if (error) return
       if (pEnergies(1,1) > 0.0) then !SEGFAULT???
          ! From basic unit to the unit of observation - preparation to the obs operator.
          conversion = fu_conversion_factor(fu_basic_mass_unit(material), obs_unit, material)
          if (fu_fails(conversion > 0.0, 'Bad conversion factor to obs unit', 'fu_init_dose_rate_obs')) return
          no_gamma_spc = no_gamma_spc + 1
          obs2trn(no_gamma_spc) = isp_obs
          scales(no_gamma_spc) = conversion
          if (error) return
          if (debug_level > 0) then
             call msg('Observed species', isp_obs)
             call msg('Factor', conversion)
          end if
       end if
    end do
    newObservation%num_obs_species = no_gamma_spc
    allocate(newObservation%ind_obs2transp(no_gamma_spc), &
           & newObservation%scale_transp2obs(no_gamma_spc), stat=status)
    !Cut the extra
    newObservation%ind_obs2transp = obs2trn(:no_gamma_spc)
    newObservation%scale_transp2obs = scales(:no_gamma_spc)
    ! We assume the measurement is on the first model layer. These go
    ! with the defaults.
    !
    if (defined(level)) then
       call set_error('No support for defined levels yet', 'fu_init_dose_rate_obs')
       return
    end if
    
    newObservation%cell_volume &  
         & = newObservation%station%cell_area & 
         &   * fu_layer_thickness_m(fu_level(vertical, 1)) !function fu_layer_thickness defined in silja levels
    
    newObservation%obsData = obsData(1:n_values)
  end function fu_init_dose_rate_obs
 
  
  !**********************************************************************************
  !
  ! Inject and observe
  !
  !**********************************************************************************

  subroutine inject_dose_rate(addition, obs, map_c, map_px, map_py, map_pz, wetdep, drydep, now_, timestep)
    implicit none
    
    type(t_dose_rate_obs_addition), intent(inout) :: addition
    type(inSituObservation), intent(inout) :: obs
    type(Tmass_map), intent(inout) :: map_c, map_px, map_py, map_pz, wetdep, drydep
    
    type(silja_time), intent(in) :: now_
    type(silja_interval), intent(in) :: timestep
    
    !internal
    integer :: ind_obs, ind_src, isp_obs, i, j, k, isp_transp, x, y
    real :: weight, inj_air_dose, inj_wet_ground_dose, inj_dry_ground_dose !injects in dose
    real :: incr_air_mass, incr_wet_ground, incr_dry_ground !increments in mass/concetration
    real :: moment_multiplier, moment_divident, time_fract, infty
    type(silja_interval) :: overlap
    type(silja_time) :: time_start, time_end, step_start, now
    
    if ((.not. defined(map_c)) .or. (.not. defined(wetdep)) .or. (.not. defined(drydep))) then
       call set_error('Concentration or deposition maps not defined','inject_dose_rate')
    else if ((.not. defined(map_px)) .or. (.not. defined(map_py)) .or. (.not. defined(map_pz))) then
       call set_error('Moment maps not defined','inject_dose_rate')
    end if
    if (error) return
    
    !Change time direction if needed
    if (obs%time_direction == forwards) then
      obs%observationTimestep = obs%datalength
      obs%time_direction = backwards
    end if
    
    !Define time variables from obs and subroutine input
    now = now_
    step_start = now + timestep
    ind_obs = obs%observationTimestep
    time_end = obs%endTimes(ind_obs)
    time_start = time_end - obs%durations(ind_obs)
    !Go backwards in time
    do while (time_start > now .and. ind_obs > 1)
      call advance()
    end do
    if (time_start >= now) return ! used all observations
    
    !Find station nearest gridpoint
    x = addition%obs%station%ix_dispersion
    y = addition%obs%station%iy_dispersion
    do while (time_end > step_start)
       overlap = fu_time_overlap(time_start, time_end, step_start, now)  
       if (fu_fails(.not.(overlap == zero_interval), 'Impossible error', 'inject_dose_rate')) return
       time_fract = fu_sec(overlap) / fu_sec(obs%durations(ind_obs))
       !Data assimilation rules are different for negative and positive values
       if (obs%obsData(ind_obs) .eps. 0.0) then
          ind_src = DA_ZERO
          weight = 1.0
       else if (obs%obsData(ind_obs) > obs%modelData(ind_obs)) then
          ind_src = DA_NEGATIVE
          weight = DA_NEGT_COEF_OBS
       else
          ind_src = DA_POSITIVE
          weight = 1.0
       end if
       if (ind_src == int_missing) then
          call set_error('ind_src is set int_missing, probably no 4Dvar DA so inject shouldnt be used', &
                       & 'inject_dose_rate')
       end if
       if (error) return
       if (addition%ifGroundToo .eqv. .false.) then !inject only to air
          inj_air_dose = time_fract &
                       & * (obs%modelData(ind_obs) &!difference between model and
                          & - obs%obsData(ind_obs)) &  !observation data
                       & * observation_scaling &    !numerical trick see the beginning of this module
                       & / obs%variance(ind_obs)    !divide by variance
          !Same logic in the below formulas for inject doses where 
          !ground is getting injected too plus the cloud_ground_ratio 
          !proportianating the dose between air, wet and dry.
          !Not injected in ground
          inj_wet_ground_dose = 0.0
          inj_dry_ground_dose = 0.0
       else !Inject proportionally to air and ground (DEFAULT)
          !air
          inj_air_dose = time_fract &
                       & * (obs%modelData(ind_obs) & 
                          & - obs%obsData(ind_obs)) &
                       & * (1.0 - addition%cloud_ground_ratio(1,ind_obs) & 
                          & - addition%cloud_ground_ratio(2,ind_obs)) & 
                       & * observation_scaling &
                       & / obs%variance(ind_obs)
          !ground wet
          inj_wet_ground_dose = time_fract &
                              & * (obs%modelData(ind_obs) &
                                 & - obs%obsData(ind_obs)) &
                              & * addition%cloud_ground_ratio(1,ind_obs) &
                              & * observation_scaling &
                              & / obs%variance(ind_obs)
          !ground dry
          inj_dry_ground_dose = time_fract &
                              & * (obs%modelData(ind_obs) &
                                 & - obs%obsData(ind_obs)) &
                              & * addition%cloud_ground_ratio(2,ind_obs) & 
                              & * observation_scaling &
                              & / obs%variance(ind_obs)
       end if
       if (isnan(inj_air_dose)) call set_error('Injecting NaNs','inject_dose_rate')
       if (isnan(inj_wet_ground_dose)) call set_error('Injecting NaNs','inject_dose_rate')
       if (isnan(inj_dry_ground_dose)) call set_error('Injecting NaNs','inject_dose_rate')
       if (error) return
       do isp_obs = 1, obs%num_obs_species
          isp_transp = obs%ind_obs2transp(isp_obs)
          do i = addition%x_min,addition%x_max
             do j = addition%y_min,addition%y_max
                do k = addition%z_min,addition%z_max
                   !Makes the concentration valid in these indices
                   map_c%ifColumnValid(ind_src, i, j) = .true.
                   map_c%ifGridValid(1, ind_src) = .true.
                   
                   if (addition%obs_map(i,j,k,isp_obs) == 0.0) then
                      incr_air_mass = 0.0
                   else
                      incr_air_mass &
                           & = weight * inj_air_dose &                         !DA (rules) weight * dose
                           & * addition%relative_dose(i,j,k,isp_obs,ind_obs) & !total dose * relative dose(i,j,k,ispc) = dose(i,j,k,ispc)
                           & / addition%obs_map(i,j,k,isp_obs)                 !dose(i,j,k,ispc)/weight(i,j,k,ispc) = C(i,j,k,ispc)
                   end if
                   
                   if (isnan(incr_air_mass)) call set_error('NaN increment mass','inject_dose_rate')
                   if (abs(incr_air_mass) > huge(infty)) call set_error('Infinity increment','inject_dose_rate')
                   if (error) return
                   
                   map_c%arm(isp_transp,ind_src,k,i,j) = map_c%arm(isp_transp,ind_src,k,i,j) + incr_air_mass
                   
                   !Note: moments not injected currently
                   !Moments are the weighted average of previous moments and incr with moment 0
                   !moment_divident = map_c%arm(isp_transp, ind_src, k, i, j) + incr_air_mass
                   !if (moment_divident .ne. 0.0) then
                   !   moment_multiplier = map_c%arm(isp_transp, ind_src, k, i, j) / moment_divident
                   !   map_px%arm(isp_transp,ind_src,k,i,j) = map_px%arm(isp_transp,ind_src,k,i,j)*moment_multiplier
                   !   map_py%arm(isp_transp,ind_src,k,i,j) = map_py%arm(isp_transp,ind_src,k,i,j)*moment_multiplier
                   !   map_pz%arm(isp_transp,ind_src,k,i,j) = map_pz%arm(isp_transp,ind_src,k,i,j)*moment_multiplier
                   !end if
                end do
             end do
          end do
          if (addition%ifGroundToo .eqv. .true.) then
             !to save time: the ground increments were already set to 0 if injection only to air
             if (addition%groundshine_multiplier(isp_obs) .ne. 0.0) then
                !Logic here same as in above injection into air
                incr_wet_ground = weight * inj_wet_ground_dose & 
                                & * addition%relative_wet(isp_obs,ind_obs) &
                                & / addition%groundshine_multiplier(isp_obs)
                incr_dry_ground = weight * inj_dry_ground_dose &
                                & * addition%relative_dry(isp_obs,ind_obs) &
                                & / addition%groundshine_multiplier(isp_obs)
                wetdep%arm(isp_transp,ind_src,1,x,y) = wetdep%arm(isp_transp,ind_src,1,x,y) + incr_wet_ground
                drydep%arm(isp_transp,ind_src,1,x,y) = drydep%arm(isp_transp,ind_src,1,x,y) + incr_dry_ground
                if (isnan(drydep%arm(isp_transp,ind_src,1,x,y))) call set_error('NaNs in injected deposition data','inject_dose_rate')
                if (isnan(wetdep%arm(isp_transp,ind_src,1,x,y))) call set_error('NaNs in injected deposition data','inject_dose_rate')
                if (error) return
             end if
          end if
       end do
       if (time_start < step_start) then
          exit
       else if (ind_obs == 1) then
          exit
       else
          call advance()
       end if
       
    end do
    if (any(isnan(map_c%arm))) call set_error('NaNs in concentration','inject_dose_rate')
    if (any(map_c%arm > abs(huge(infty)))) call set_error('Infinity in concentration','inject_dose_rate')
    if (error) return
    
    obs%observationTimestep = ind_obs
  contains
    
    subroutine advance()
      implicit none
      ind_obs = ind_obs - 1
      time_end = obs%endTimes(ind_obs)
      time_start = time_end - obs%durations(ind_obs)
    end subroutine advance
      
  end subroutine inject_dose_rate
  
  !************************************************************************

  subroutine observe_dose_rate(addition, obs, map_c, map_wet, map_dry, now, timestep)
    !Most of observe_in_situ preserved, only type addition added to
    !make observation and injection of dose rate possible.
    implicit none
    type(inSituObservation), intent(inout) :: obs
    type(t_dose_rate_obs_addition), intent(inout) :: addition
    type(Tmass_map), intent(in) :: map_c, map_wet, map_dry
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    
    integer :: ind_obs
    real :: time_fract, v
    type(silja_time) :: time_start, time_end, step_end
    type(silja_interval) :: overlap
    
    if ((.not. defined(map_c)) .or. (.not. defined(map_wet)) .or. (.not. defined(map_dry))) then
       call set_error('Concentration or deposition maps not defined','observe_dose_rate')
    end if
    if (error) return
    
    if (obs%time_direction == backwards) then
      obs%observationTimestep = 1
      obs%time_direction = forwards
    end if
    
    ind_obs = obs%observationTimestep
    
    time_end = obs%endTimes(ind_obs)
    time_start = time_end - obs%durations(ind_obs)
    
    do while (time_end <= now .and. ind_obs < obs%datalength)
      call advance()
    end do
    
    if (time_end <= now) return
    
    step_end = now + timestep
    do while (time_start < step_end)
      overlap = fu_time_overlap(time_start, time_end, now, step_end)
      if (fu_fails(.not.(overlap == zero_interval), 'Impossible error', 'observe_dose_rate')) return
      time_fract = fu_sec(overlap) / fu_sec(obs%durations(ind_obs))
      !Obs model data is the actual dose rate at the station and addition model data
      !also contains how it is distributed and weighted in physical space and nuclides.
      obs%modeldata(ind_obs) = obs%modeldata(ind_obs) & 
                             & + time_fract &
                               & * fu_get_dose_rate_at_station(addition, map_c, map_wet, map_dry)
      !Above is the actual dose rate observation and
      !these below are used as help in injection to achieve proportionality
      addition%relative_dose(:,:,:,:,ind_obs) = addition%relative_dose(:,:,:,:,ind_obs) &
                                              & + time_fract &
                                                & * fu_get_relative_air_dose(addition, map_c)
      
      addition%cloud_ground_ratio(:,ind_obs) = addition%cloud_ground_ratio(:,ind_obs) &
                                             & + time_fract & 
                                               & * fu_get_ground_air_ratio(addition, map_c, map_wet, map_dry)
                                             
      addition%relative_wet(:,ind_obs) = time_fract * fu_get_relative_ground_dose(addition, map_wet)
      addition%relative_dry(:,ind_obs) = time_fract * fu_get_relative_ground_dose(addition, map_dry)
      
      !Error triggers originally for testing
      if ((sum(addition%relative_wet(:,ind_obs)) < 0.99999) &
        & .and. (any(addition%relative_wet(:,ind_obs) /= 0.0))) call set_error('total ratio cant be under 1','observe_dose_rate')
      if ((sum(addition%relative_dry(:,ind_obs)) < 0.99999) &
        & .and. (any(addition%relative_dry(:,ind_obs) /= 0.0))) call set_error('total ratio cant be under 1','observe_dose_rate')
      if (sum(addition%relative_wet(:,ind_obs)) > 1.00001) call set_error('total ratio cant be over 1','observe_dose_rate')
      if (sum(addition%relative_dry(:,ind_obs)) > 1.00001) call set_error('total ratio cant be over 1','observe_dose_rate')
      if (sum(addition%cloud_ground_ratio(:,ind_obs)) > 1.00001) call set_error('total ratio cant be over 1','observe_dose_rate')
      if (any(addition%cloud_ground_ratio(:,ind_obs) > 1.00001)) call set_error('ratio cant be over 1','observe_dose_rate')
      if (any(addition%cloud_ground_ratio(:,ind_obs) < 0.0)) call set_error('ratio cant be under 0','observe_dose_rate')
      if (error) return
      
      if (time_end > step_end) then
        !timestep ends before observation -> leave the remainder for the next time step.
        exit
      else if (ind_obs < obs%datalength) then
        call advance()
      else
        ! all obs times processed
        exit
      end if
    end do
    obs%observationTimestep = ind_obs
    
  contains
    
    subroutine advance()
      implicit none
      ind_obs = ind_obs + 1
      time_end = obs%endTimes(ind_obs)
      time_start = time_end - obs%durations(ind_obs)
      !Values are reset for each new observation index
      obs%modeldata(ind_obs) = 0.0
      addition%relative_dose(:,:,:,:,ind_obs) = 0.0
      addition%cloud_ground_ratio(:,ind_obs) = 0.0
      addition%relative_wet(:,ind_obs) = 0.0
      addition%relative_dry(:,ind_obs) = 0.0
    end subroutine advance
    
  end subroutine observe_dose_rate


  !************************************************************************************
    
  
  subroutine test_dose_rate(vert, timestep, obs_len, &
                              & obs, observed_data, injected_data, &
                              & spc, map_c, wetdep, drydep, map_px, map_py, map_pz)
    implicit none
    type(silam_vertical), intent(in) :: vert
    type(silja_interval), intent(in) :: timestep
    integer, intent(in) :: obs_len
    type(inSituObservation), intent(out) :: obs
    
    type(tmass_map), pointer :: map_c, wetdep, drydep, map_px, map_py, map_pz
    type(silam_species), dimension(:), intent(in) :: spc
    integer :: ii
    type(silja_interval) :: shift
    real, dimension(obs_len) :: obs_data, variance
    type(silja_time), dimension(obs_len) :: obs_times_start
    type(silja_interval), dimension(obs_len) :: obs_durations
    type(observationStation) :: station
    type(t_dose_rate_obs_addition) :: addition
    type(silja_time) :: now = time_missing, time_end=time_missing, time_start=time_missing
    real :: x, y, injected_data, observed_data
    integer :: ispecies_obs, ispecies_transp
    !To test dose rate as I did without having actual observations but by calling
    ! the observe and inject + necessary initializations do the following:
    ! -> Set call the injectAll at dispersion_supplementary module inside and in the 
    !    beginning of the loop over time:
    !       call injectAll(obs_ptrs, cloud, meteo_ptr, &
    !                    & simRules%chemicalRules, simRules%dynamicsRules, &
    !                    & simrules%periodToCompute, now)
    !    ps. the injectAll wont be reached otherwise without correct 4dvar ini-file
    ! -> Include calling this subroutine test_dose_rate inside the injectAll. 
    !    First set some values and point to maps and then call the routine:
    !       type(inSituObservation) :: obs
    !       real :: observed_data, injected_data
    !       type(Tmass_map), pointer :: map_px, map_py, map_pz
    !    Then:
    !       map_px => fu_advection_moment_X_MM_ptr(cloud)
    !       map_py => fu_advection_moment_Y_MM_ptr(cloud)
    !       map_pz => fu_advection_moment_Z_MM_ptr(cloud)
    !       call test_dose_rate(dispersion_vertical, &
    !                         & timestep, 1, &
    !                         & obs, &
    !                         & observed_data, injected_data, &
    !                         & fu_species_transport(cloud), mapConc, mapWetdep, mapDrydep, map_px, map_py, map_pz)
    !       print*,'the obs data:',observed_data
    !       print*,'the inj data:',injected_data
    !    Printing just to get the observed and injected dose rates 
    !    to compare what the injection did to the value.
    ! -> Also set ind_src to one in inject as it isnt set if the run isnt adjoint:
    !       ind_src = 1.0
    
    time_start = fu_set_time_utc(2015, 7, 22, 06, 00, 0.0)
    shift = fu_set_interval_min(20)
    
    time_end = time_start + shift + one_hour*obs_len + shift
    do ii = 1,obs_len
       obs_times_start(ii) = time_start + shift + (one_hour*ii)
    end do
    obs_durations = one_hour
    obs_data = 1.0   !are actually 
    variance = 100.0 !observed
    call project_point_to_grid(0.0, 0.0, wholeMPIdispersion_grid,x,y)
    station = fu_initObservationStation('S1', 'S1', 0.0, 0.0, 1000.0, wholeMPIdispersion_grid)
    obs = fu_init_dose_rate_obs(obs_times_start, obs_durations, obs_data, variance, &
                              & obs_len, constant_variance, spc, &
                              & 'mole', station, level_missing, vert, 'test')
    addition = fu_init_dose_rate_addition(.true., obs, vert, spc)
    
    now = time_start
    do while (now <= time_end) 
      call observe_dose_rate(addition, obs, map_c, wetdep, drydep, now, timestep)
      if (error) return
      now = now + timestep
    end do
    obs%obsData = obs%modelData(1)*0.9
    observed_data = obs%modelData(1)

    do while (now >= time_start)
      call inject_dose_rate(addition, obs, map_c, map_px, map_py, map_pz, wetdep, drydep, now, fu_opposite(timestep))
      if (error) return
      now = now - timestep
    end do
    do while (now <= time_end) 
      call observe_dose_rate(addition, obs, map_c, wetdep, drydep, now, timestep)
      if (error) return
      now = now + timestep
    end do

    injected_data = obs%modelData(1)
  end subroutine test_dose_rate



  !**********************************************************************************
  !
  ! Auxiliary routines
  ! laguerre, legendre and bruteforce ways to calc integration weights, only one used
  !
  !**********************************************************************************

  function fu_get_dose_rate_at_station(addition, mapconc, wetdep, drydep) result(v) !needs deposition concentration in ground
    !gets a value corresponding to the  observation from the massmap

    implicit none
    real :: v
    type(Tmass_map), intent(in)  :: mapconc,wetdep,drydep
    type(t_dose_rate_obs_addition), intent(in) :: addition

    ! local
    integer :: ispecies_transp,ispecies_obs,x,y,i,j,k
    
    !get the dose rate at station from cloudshine
    v = 0.0
    do ispecies_obs = 1,addition%obs%num_obs_species
       ispecies_transp = addition%obs%ind_obs2transp(ispecies_obs)
       do i = addition%x_min,addition%x_max
          do j = addition%y_min,addition%y_max
             do k = addition%z_min,addition%z_max
                v = v + mapconc%arm(ispecies_transp,1,k,i,j) &
                      & * addition%obs_map(i,j,k,ispecies_obs)
             end do
          end do
       end do
    end do
    !get the dose rate at station from groundshine (closest grid point represents the ground source)
    x = addition%obs%station%ix_dispersion
    y = addition%obs%station%iy_dispersion
    do ispecies_obs = 1,addition%obs%num_obs_species
       ispecies_transp = addition%obs%ind_obs2transp(ispecies_obs)
       v = v + (wetdep%arm(ispecies_transp,1,1,x,y) + drydep%arm(ispecies_transp,1,1,x,y)) &
             & * addition%groundshine_multiplier(ispecies_obs)
    end do
    
  end function fu_get_dose_rate_at_station
  
  
  function fu_get_relative_air_dose(addition, mapconc) result(weights)
    !gets dose rate impact mapping for each nuclide
    ! so that inject knows which nuclides have what relative impact on the dose rate measurement
    
    type(Tmass_map), intent(in)  :: mapconc
    type(t_dose_rate_obs_addition), intent(in) :: addition
    real, dimension(:,:,:,:), allocatable :: weights

    ! local
    integer :: ispecies_transp, ispecies_obs, i, j, k

    allocate(weights(addition%x_min:addition%x_max, &
                   & addition%y_min:addition%y_max, &
                   & addition%z_min:addition%z_max, &
                   & addition%obs%num_obs_species))
    
    do i = addition%x_min,addition%x_max
       do j = addition%y_min,addition%y_max
          do k = addition%z_min,addition%z_max
             do ispecies_obs = 1,addition%obs%num_obs_species
                ispecies_transp = addition%obs%ind_obs2transp(ispecies_obs)
                weights(i,j,k,ispecies_obs) = mapconc%arm(ispecies_transp,1,k,i,j) &
                                            & * addition%obs_map(i,j,k,ispecies_obs)
             end do
          end do
       end do
    end do
    if (sum(weights) .ne. 0.0) weights = weights/sum(weights) !to relative value
  end function fu_get_relative_air_dose
  
  
  
  
  function fu_get_ground_air_ratio(addition, mapconc, wetdep, drydep) result(ratio)
    !get the ratio between groundshine and cloudshine doses
    !if cloudshine dose is 50 and groundshine 30, the ratio will be 0.375 
    ! as 37.5% of total dose is from the ground
    type(Tmass_map), intent(in) :: mapconc, wetdep, drydep
    type(t_dose_rate_obs_addition), intent(in) :: addition
    real, dimension(2) :: ratio
    real :: cloud, ground_dry, ground_wet, total
    integer :: ispecies_transp, ispecies_obs, x, y,i,j,k
    !get the dose rate at station from cloudshine
    cloud = 0.0
    do ispecies_obs = 1,addition%obs%num_obs_species
       do i = addition%x_min,addition%x_max
          do j = addition%y_min,addition%y_max
             do k = addition%z_min,addition%z_max
                ispecies_transp = addition%obs%ind_obs2transp(ispecies_obs)
                cloud = cloud + mapconc%arm(ispecies_transp,1,k,i,j) &
                              & * addition%obs_map(i,j,k,ispecies_obs)
             end do
          end do
       end do
    end do
    !get the dose rate at station from groundshine
    ground_wet = 0.0
    ground_dry = 0.0
    x = addition%obs%station%ix_dispersion
    y = addition%obs%station%iy_dispersion
    do ispecies_obs = 1,addition%obs%num_obs_species
       ispecies_transp = addition%obs%ind_obs2transp(ispecies_obs)
       ground_wet = ground_wet + wetdep%arm(ispecies_transp,1,1,x,y) &
                               & * addition%groundshine_multiplier(ispecies_obs)
       ground_dry = ground_dry + drydep%arm(ispecies_transp,1,1,x,y) &
                               & * addition%groundshine_multiplier(ispecies_obs)
    end do
    total = cloud + ground_wet + ground_dry
    if (total == 0.0) then
       ratio = 0.0
    else !division may resulty to small negative/unphysical ratios
       ratio(1) = max(0.0,ground_wet/total)
       ratio(2) = max(0.0,ground_dry/total)
    end if
    
  end function fu_get_ground_air_ratio
  
  function fu_get_relative_ground_dose(addition, depmap) result(weights)
    type(Tmass_map), intent(in) :: depmap
    type(t_dose_rate_obs_addition), intent(in) :: addition
    real, dimension(:), allocatable :: weights
    integer :: ispecies_transp, ispecies_obs, x, y
    allocate(weights(addition%obs%num_obs_species))
    x = addition%obs%station%ix_dispersion
    y = addition%obs%station%iy_dispersion
    do ispecies_obs = 1,addition%obs%num_obs_species
       ispecies_transp = addition%obs%ind_obs2transp(ispecies_obs)
       weights(ispecies_obs) = depmap%arm(ispecies_transp,1,1,x,y) &
                             & * addition%groundshine_multiplier(ispecies_obs)
    end do
    if (sum(weights) .ne. 0.0) weights = weights/sum(weights) 
  end function fu_get_relative_ground_dose
  
  
  
  !Resets dose rate obs and addition to the beginning, nullifies on request
  subroutine restart_dose_rate(obs, addition, ifNullifyModel)
    implicit none
    type(inSituObservation), intent(inout) :: obs
    type(t_dose_rate_obs_addition), intent(inout) :: addition
    logical, intent(in) :: ifNullifyModel
    
    obs%observationTimestep = 1
    obs%time_direction = forwards
    if(ifNullifyModel) then
       obs%modeldata = 0.0
       addition%relative_dose = 0.0
       addition%cloud_ground_ratio = 0.0
       addition%relative_wet = 0.0
       addition%relative_dry = 0.0
    end if
  end subroutine restart_dose_rate
  


  function cart_to_sph(C)
    !takes cartesian concentration and makes the field spherical
    !input in x,y,z
    !output in theta,phi,r
    implicit none
    real, dimension(:,:,:,:) :: C
    real, dimension(:,:,:,:), allocatable :: cart_to_sph, C_cart
    integer :: i, j, k, nlat, nlon, nz, npar
    real :: theta, phi, r
    nlat = int(size(C(:,1,1,1)))
    nlon = int(size(C(1,:,1,1)))
    nz = int(size(C(1,1,:,1)))
    npar = int(size(C(1,1,1,:)))
    allocate(cart_to_sph(nlat,nlon,nz,npar))
    C_cart = C
    do k = 1,nz
       do j = 1,nlon
          do i = 1,nlat
             r = sqrt(C_cart(i,j,k,1)**2.0 &
                    & + C_cart(i,j,k,2)**2.0 &
                    & + C_cart(i,j,k,3)**2.0)
             theta = atan(C_cart(i,j,k,2)/C_cart(i,j,k,1)) + 0.5*pi
             if (C_cart(i,j,k,1) == 0.0) then
                theta = pi
             end if
             phi = acos(C_cart(i,j,k,3)/r)
             if (r < 1.0) then
                phi = 0.0
             end if
             cart_to_sph(i,j,k,:) = (/ theta, phi, r /)
          end do
       end do
    end do
  end function cart_to_sph

  function bruteforce(C,pss)
    !divides the vicinity of the origin into point source size (pss) cubes and that each radiate to the origin
    implicit none
    real, dimension(:), allocatable :: xyz
    real, dimension(:,:,:), allocatable :: x,y,z,bruteforce
    real, dimension(:,:,:,:) :: C
    real, dimension(6,14) :: cerc
    real :: pss,Rmax,d,dist,dist_to_cell,shortest,volume,B,minx,maxx,miny,maxy,minz,maxz,R
    integer :: n,i,j,k,ii,jj,kk,q,w,e,i1,j1,k1,nlat,nlon,nz,npar
    nlat = size(C(:,1,1,1))
    nlon = size(C(1,:,1,1))
    nz = size(C(1,1,:,1))
    npar = size(C(1,1,1,:))
    cerc = cerc_params()
    allocate(x(nlat,nlon,nz))
    allocate(y(nlat,nlon,nz))
    allocate(z(nlat,nlon,nz))
    allocate(bruteforce(nlat,nlon,nz))
    x = C(:,:,:,2) !coordinates into temporary arrays to simplify equations
    y = C(:,:,:,3)
    z = C(:,:,:,4)
    Rmax = log(0.001)/(-cerc(2,12)) !maximum distance where 0.1% of incident radiation emitted is captured due to attenuation
    n = ceiling(2.0*Rmax/pss) !number of cells for the extent with the predetermined point source size
    d = 2.0*Rmax/n
    allocate(xyz(n+1))
    do i = 1,n
       xyz(i) = i*d-Rmax
    end do
    volume = d**3.0 !volume of cells calculated
    minx = minval(x) !boundaries found before loop to get them only once
    maxx = maxval(x)
    miny = minval(y)
    maxy = maxval(y)
    minz = minval(z)
    maxz = maxval(z)
    bruteforce = 0.0
    do kk = 1,n
       do jj = 1,n
          do ii = 1,n
             dist = max(1.0,sqrt(xyz(ii)**2.0 + xyz(jj)**2.0 & !if distance becomes less that unity then dose will approach infinity
                                & + xyz(kk)**2.0))             !=unphysical
             if (dist < Rmax .and.  & !if distance is within limits of real field and wanted integration limit
                & xyz(ii) > minx .and. xyz(ii) < maxx .and. &
                & xyz(jj) > miny .and. xyz(jj) < maxy .and. &
                & xyz(kk) > minz .and. xyz(kk) < maxz) then
                shortest = earth_radius
                do k = 1,nz
                   do j = 1,nlon
                      do i = 1,nlat
                         R = sqrt((x(i,j,k))**2.0 + (y(i,j,k))**2.0 + (z(i,j,k))**2.0)
                         dist_to_cell = sqrt((x(i,j,k)-xyz(ii))**2.0 &
                                            & + (y(i,j,k)-xyz(jj))**2.0 &
                                            & + (z(i,j,k)-xyz(kk))**2.0)
                         if (dist_to_cell < shortest) then
                            shortest = dist_to_cell
                            i1 = i
                            j1 = j
                            k1 = k
                         end if
                      end do
                   end do
                end do
                B = 1.0 + cerc(3,12) * cerc(2,12) * dist & 
                          & * exp(cerc(4,12) * cerc(2,12) * dist)
                bruteforce(i1,j1,k1) = bruteforce(i1,j1,k1) &
                                     & + B*exp(-cerc(2,12)*dist)/(dist**2.0)
             end if
          end do
       end do
    end do
    bruteforce = bruteforce*volume/(4.0*pi)
  end function bruteforce
  
  
  function cloudshine_legendre(C_leg,nEnergy)
    implicit none
    !Quite straightforwardly explained in: 'http://www.cerc.co.uk/environmental-software/assets/data/doc_techspec/CERC_ADMS5_P20_01.pdf'
    !given concentrations in theta,phi,r field in measuring devices spherical coordinates and radius where legendre integration is extended gives the total cloud shine weights to all real world grid cells that can then be multiplied cell by cell and taken a sum to obtain the total
    real, dimension(:,:,:,:) :: C_leg
    real, dimension(:,:,:), allocatable :: cloudshine_legendre
    real, dimension(:,:), allocatable :: leg
    real, dimension(:), allocatable :: thetaj, phik, ri
    real, dimension(2,15) :: leg15
    real, dimension(2,30) :: leg30
    real, dimension(6,14) :: cerc
    real :: B,thetamin,thetamax,phimin,phimax,shortest,dist,Rmax,Rmin
    integer :: i,j,k,l,m,n,ll,mm,nn,nlat,nlon,nz,npar,nEnergy
    nlat = size(C_leg(:,1,1,1))
    nlon = size(C_leg(1,:,1,1))
    nz = size(C_leg(1,1,:,1))
    npar = size(C_leg(1,1,1,:))
    leg30 = legendre30_params()
    leg15 = legendre_params()
    allocate(leg(size(leg15(:,1)),size(leg15(1,:)))) !change these two lines from leg15 to leg30 to implement
    leg = leg15                                      !twice the integration points for additional accuracy
    cerc = cerc_params()
    allocate(cloudshine_legendre(nlat,nlon,nz))
    !theta, phi and r changed because legendre integral is between -1 and 1 (abscissae)
    allocate(thetaj(size(leg(1,:))))
    allocate(phik(size(leg(1,:))))
    allocate(ri(size(leg(1,:))))
    thetamax = pi
    thetamin = 0.0
    phimax = 2.0*pi
    phimin = 0.0
    thetaj = 0.5*(thetamax-thetamin)*leg(2,:) + 0.5*(thetamax+thetamin)
    phik = 0.5*(phimax-phimin)*leg(2,:) + 0.5*(phimax+phimin)
    !integral extended to radius where incident radiation is reduced to 0.1% from original
    Rmax = -log(0.001)/(cerc(2,nEnergy))
    !Rmin can be set to f.e. 30m
    Rmin = 0.0
    ri = 0.5*(Rmax-Rmin)*leg(2,:) + 0.5*(Rmax+Rmin)
    cloudshine_legendre = 0.0
    !looping over legendre integral points
    do k = 1,size(phik)
       do j = 1,size(thetaj)
          do i = 1,size(ri)
             !find the closest points to predetermined legendre points
             shortest = earth_radius
             !looping over the grid coordinates and finding the nearest of grid coords to legendre integral coords
             do n = 1,nz
                do m = 1,nlon
                   do l = 1,nlat
                      dist = sqrt(C_leg(l,m,n,4)**2.0 &
                                 & + ri(i)**2.0 &
                                 & - 2.0*C_leg(l,m,n,4)*ri(i) &
                                   & * (sin(C_leg(l,m,n,2))*sin(thetaj(j)) &
                                   & * cos(C_leg(l,m,n,3)-phik(k)) &
                                   & + cos(C_leg(l,m,n,2))*cos(thetaj(j))))
                      if (dist < shortest) then
                         shortest = dist
                         ll = l
                         mm = m
                         nn = n
                      end if
                   end do
                end do
             end do
             B = 1.0 + cerc(3,nEnergy) * cerc(2,nEnergy) * ri(i) &
                       & * exp(cerc(4,nEnergy) * cerc(2,nEnergy) * ri(i))
             cloudshine_legendre(ll,mm,nn) = cloudshine_legendre(ll,mm,nn) & 
                                           & + leg(1,i)*leg(1,j)*leg(1,k)*B &
                                             & * exp(-cerc(2,nEnergy)*ri(i))*sin(thetaj(j))
          end do
       end do
    end do
    cloudshine_legendre = cloudshine_legendre*pi*(Rmax-Rmin)/16.0
  end function cloudshine_legendre


  function cloudshine_laguerre(C_lag,nEnergy)
    implicit none
    !cloudshine_legendre more commented, otherwise similar but r integral in laguerre and not legendre style
    !r in laguerre assumes integral is from 0 to inf
    !laguerre abscissae and weights differ and integral is 2 part with some changing of variables in the integrals
    !current form achieved from this cerc paper: 'http://www.cerc.co.uk/environmental-software/assets/data/doc_techspec/CERC_ADMS5_P20_01.pdf'
    !given concentrations in theta,phi,r field in measuring devices spherical coordinates and radius where legendre for theta and phi..
    !and laguerre integration for r is extended gives the total cloud shine to that point
    real, dimension(:,:,:,:) :: C_lag
    real, dimension(:,:,:), allocatable :: cloudshine_laguerre
    real, dimension(:,:), allocatable :: leg
    real, dimension(:), allocatable :: thetaj, phik
    real, dimension(2,15) :: leg15,lag
    real, dimension(2,30) :: leg30
    real, dimension(15) :: ri1, ri2
    real, dimension(6,14) :: cerc
    real :: maxr,term,thetamax,thetamin,phimax,phimin,shortest1,shortest2,dist1,dist2,dc
    integer :: i,j,k,l,m,n,ll1,ll2,mm1,mm2,nn1,nn2,nlat,nlon,nz,npar,nEnergy
    nlat = size(C_lag(:,1,1,1))
    nlon = size(C_lag(1,:,1,1))
    nz = size(C_lag(1,1,:,1))
    npar = size(C_lag(1,1,1,:))
    leg30 = legendre30_params()
    leg15 = legendre_params()
    allocate(leg(size(leg15(:,1)),size(leg15(1,:))))
    leg = leg15
    lag = laguerre_params()
    cerc = cerc_params()
    allocate(cloudshine_laguerre(nlat,nlon,nz))
    !theta and phi changed because legendre integral is between -1 and 1
    !laguerre integral from 0 to inf assuming exp(-x) dependence
    allocate(thetaj(size(leg(1,:))))
    allocate(phik(size(leg(1,:))))
    thetamax = pi
    thetamin = 0.0
    phimax = 2.0*pi
    phimin = 0.0
    thetaj = 0.5*(thetamax-thetamin)*leg(2,:) + 0.5*(thetamax+thetamin)
    phik = 0.5*(phimax-phimin)*leg(2,:) + 0.5*(phimax+phimin)
    !integral over r in 2 parts
    ri1 = lag(2,:)/cerc(2,nEnergy)
    ri2 = lag(2,:)/(cerc(2,nEnergy)*(1.0-cerc(4,nEnergy)))
    maxr = max(maxval(ri1),maxval(ri2))
    
    
    cloudshine_laguerre = 0.0
    do k = 1,size(phik)
       do j = 1,size(thetaj)
          do i = 1,size(ri1)
             !find the closest points to predetermined legendre points
             shortest1 = earth_radius !starting to search for closest distance with this large distance
             shortest2 = earth_radius
             do n = 1,nz
                do m = 1,nlon
                   do l = 1,nlat
                      !term calculated separately for simplicity
                      term = sin(C_lag(l,m,n,2))*sin(thetaj(j))*cos(C_lag(l,m,n,3)-phik(k)) &
                           & + cos(C_lag(l,m,n,2))*cos(thetaj(j))
                      dist1 = sqrt(C_lag(l,m,n,4)**2.0 &
                                  & + ri1(i)**2.0 &
                                  & - 2.0*C_lag(l,m,n,4)*ri1(i)*term)
                      dist2 = sqrt(C_lag(l,m,n,4)**2.0 &
                                  & + ri2(i)**2.0 &
                                  & - 2.0*C_lag(l,m,n,4)*ri2(i)*term)
                      if (dist1 < shortest1) then
                         shortest1 = dist1
                         ll1 = l
                         mm1 = m
                         nn1 = n
                      end if
                      if (dist2 < shortest2) then
                         shortest2 = dist2
                         ll2 = l
                         mm2 = m
                         nn2 = n
                      end if
                   end do
                end do
             end do
             cloudshine_laguerre(ll1,mm1,nn1) = cloudshine_laguerre(ll1,mm1,nn1) &
                                              & + lag(1,i)*leg(1,j)*leg(1,k) &
                                                & * (1.0/cerc(2,nEnergy))*sin(thetaj(j))
             cloudshine_laguerre(ll2,mm2,nn2) = cloudshine_laguerre(ll2,mm2,nn2) &
                                              & + lag(1,i)*leg(1,j)*leg(1,k) &
                                                & * (ri2(i)*cerc(3,nEnergy)/(1.0-cerc(4,nEnergy)))*sin(thetaj(j))
          end do
       end do
    end do
    
    cloudshine_laguerre = cloudshine_laguerre*pi/8.0
  end function cloudshine_laguerre
  

  !*******************************************************************************
  !
  ! Functions that download parameters
  !
  !*******************************************************************************


  function legendre_params()
    !returns 2*15 Legendre Gaussian quadrature params with (i,j) -> i(1) = weights i(2) = abscissae, j = n
    !copied from 'https://pomax.github.io/bezierinfo/legendre-gauss.html#n10'
    real :: legendre1d(30)
    real :: legendre_params(2,15)
    legendre1d(1:30) = (/ 0.2025782419255613, 0.0000000000000000, 0.1984314853271116, -0.2011940939974345, 0.1984314853271116, &
         0.2011940939974345, 0.1861610000155622, -0.3941513470775634, 0.1861610000155622, 0.3941513470775634, &
         0.1662692058169939, -0.5709721726085388, 0.1662692058169939, 0.5709721726085388, 0.1395706779261543, &
         -0.7244177313601701, 0.1395706779261543, 0.7244177313601701, 0.1071592204671719, -0.8482065834104272, &
         0.1071592204671719, 0.8482065834104272, 0.0703660474881081, -0.9372733924007060, 0.0703660474881081, &
         0.9372733924007060, 0.0307532419961173, -0.9879925180204854, 0.0307532419961173, 0.9879925180204854 /)
    legendre_params(1,:) = legendre1d(1:29:2)!w
    legendre_params(2,:) = legendre1d(2:30:2)!x
  end function legendre_params

  function legendre30_params()
    !2*30 30n legendre params from 'https://pomax.github.io/bezierinfo/legendre-gauss.html#n10'
    real :: legendre1d(60)
    real :: legendre30_params(2,30)
    legendre1d(1:60) = (/ 0.1028526528935588, -0.0514718425553177, &
         0.1028526528935588, 0.0514718425553177, &
         0.1017623897484055, -0.1538699136085835, &
         0.1017623897484055, 0.1538699136085835, &
         0.0995934205867953, -0.2546369261678899, &
         0.0995934205867953, 0.2546369261678899, &
         0.0963687371746443, -0.3527047255308781, &
         0.0963687371746443, 0.3527047255308781, &
         0.0921225222377861, -0.4470337695380892, &
         0.0921225222377861, 0.4470337695380892, &
         0.0868997872010830, -0.5366241481420199, &
         0.0868997872010830, 0.5366241481420199, &
         0.0807558952294202, -0.6205261829892429, &
         0.0807558952294202, 0.6205261829892429, &
         0.0737559747377052, -0.6978504947933158, &
         0.0737559747377052, 0.6978504947933158, &
         0.0659742298821805, -0.7677774321048262, &
         0.0659742298821805, 0.7677774321048262, &
         0.0574931562176191, -0.8295657623827684, &
         0.0574931562176191, 0.8295657623827684, &
         0.0484026728305941, -0.8825605357920527, &
         0.0484026728305941, 0.8825605357920527, &
         0.0387991925696271, -0.9262000474292743, &
         0.0387991925696271, 0.9262000474292743, &
         0.0287847078833234, -0.9600218649683075, &
         0.0287847078833234, 0.9600218649683075, &
         0.0184664683110910, 0.9836681232797472, &
         0.0184664683110910, 0.9836681232797472, &
         0.0079681924961666, -0.9968934840746495, &
         0.0079681924961666, 0.9968934840746495  /)
    legendre30_params(1,:) = legendre1d(1:59:2)!w
    legendre30_params(2,:) = legendre1d(2:60:2)!x
  end function legendre30_params

  function laguerre_params()
    !returns 2*15 Laguerre Gaussian quadrature params with (i,j) -> i(1) = weights i(2) = abscissae, j = n
    !copied from 'http://keisan.casio.com/exec/system/1281279441'
    real :: laguerre1d(30)
    real :: laguerre_params(2,15)
    laguerre1d(1:30) = (/ 0.0933078120172818047629, 0.2182348859400868899, 0.49269174030188390896, 0.342210177922883329639, &
         1.21559541207094946373, 0.26302757794168009741, 2.269949526203743202474, 0.12642581810593053584, &
         3.667622721751437277249, 0.040206864921000914842, 5.425336627413553165344, 0.008563877803611838364, &
         7.565916226613067860497, 0.00121243614721425207622, 10.12022856801911273479, 1.1167439234425194199E-4, &
         13.1302824821757235641, 6.4599267620229009247E-6, 16.65440770832995782252, 2.22631690709627263033E-7, &
         20.77647889944876677292, 4.2274303849793650074E-9, 25.62389422672878014459, 3.92189726704108929038E-11, &
         31.40751916975393851524, 1.456515264073126406333E-13, 38.53068330648600941625, &
         1.48302705111330133546E-16, 48.02608557268579434657, 1.60059490621113323105E-20 /)
    laguerre_params(2,:) = laguerre1d(1:29:2)
    laguerre_params(1,:) = laguerre1d(2:30:2)
  end function laguerre_params


  function cerc_params()
    !cerc params so that (i,j) i(1)=E(MeV), i(2)=mu, i(3)=a, i(4)=b, i(5)=myE(Gym+2), i(6)=Cb(Sv/Gy)
    !copied from 'http://www.cerc.co.uk/environmental-software/assets/data/doc_techspec/CERC_ADMS5_P20_01.pdf'
    real :: cerc1d(84)
    real :: cerc_params(6,14)
    integer :: i
    cerc1d(1:84) = (/ 0.01,  0.623, 0.025, -0.0464, 7.43E-16, 0.00296, &
                    & 0.015, 0.187, 0.0947, -0.0484, 3.12E-16, 0.0183, &
                    & 0.02, 0.0893, 0.2652, -0.0463, 1.68E-16, 0.0543, &
                    & 0.03, 0.0411, 1.055, -0.0192, 0.721E-16, 0.191, &
                    & 0.05, 0.0253, 3.498, 0.0729, 0.323E-16, 0.557, &
                    & 0.065, 0.0226, 4.209, 0.1169, 0.278E-16, 0.63, &
                    & 0.1, 0.0195, 4.033, 0.1653, 0.371E-16, 0.765, &
                    & 0.2, 0.0159, 2.678, 0.1678, 0.856E-16, 0.703, &
                    & 0.5, 0.0112, 1.748, 0.1014, 2.38E-16, 0.689, &
                    & 1.0, 0.00821, 1.269, 0.0559, 4.47E-16, 0.732, &
                    & 1.5, 0.00668, 1.040, 0.0338, 6.12E-16, 0.765, &
                    & 2.0, 0.00574, 0.891, 0.0215, 7.50E-16, 0.791, &
                    & 4.0, 0.00398, 0.5879, 0.0022, 12.0E-16, 0.850, &
                    & 10.0, 2.65E-3, 0.3113, -0.0194, 23.1E-16, 0.935 /)
    do i = 1,6
       cerc_params(i,:) = cerc1d(i:(84-6+i):6)
    end do
  end function cerc_params
  
  function exp_int_par()
    !exponential integral parameters corresponding the energy bins in cerc params
    !for attenuation coefficient multiplied by 2m height (assumed height from ground for the
    ! measurement device)
    !values fetched from the online calculator by feeding it with the corresponding mu * 2m
    !http://keisan.casio.com/exec/system/1180573423
    real :: exp_int_par(14)
    exp_int_par(:) = (/ 0.14733, 0.748, 1.316, 2.002, 2.457, &
                      & 2.564, 2.7056, 2.90262, 3.24375, 3.54839, &
                      & 3.75159, 3.90138, 4.26405, 4.66813 /)
  end function exp_int_par


end module Observations_dose_rate
