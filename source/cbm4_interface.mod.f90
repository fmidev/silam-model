
module cbm4_interface
  use cbm4_Parameters
  use cbm4_Integrator, only : rosenbrock, init_solver, a
  use cbm4_Global
  use cbm4_Precision
  use cbm4_adj_Integrator, only : rosenbrockAdj, init_solver_adj
  use cocktail_basic
  !$use omp_lib
  implicit none

  private
  
  ! public subroutines - but exported only one level up!
  public rosenbrock
  public set_rates_cbm4
  public set_const_rates_cbm4
  public set_fixed_cbm4
  public register_species_int
  public inventory_int
  public init_solver

  public init_solver_adj
  public rosenbrockAdj

  public a

  ! only this will be exported globally
  public init_chemicals_cbm4


  integer, parameter, private :: num_species = nvar ! from _Parameters module
  integer, parameter, public :: num_species_cbm4 = num_species ! just for laziness
  integer, parameter, public :: num_react_cbm4 = nreact
  integer, parameter, public :: num_fixed_cbm4 = nfix
  integer, parameter, public :: precision_cbm4 = sp ! from Precision module
  character(len=substNmLen), dimension(num_species), parameter, private :: &
       & subst_names = (/ &
                      & 'O1D ','H2O2','PAN ','CRO ','TOL ','N2O5', &
& 'XYL ','XO2N','HONO','PNA ','TO2 ','HNO3', &
& 'ROR ','CRES','MGLY','CO  ','ETH ','XO2 ', &
& 'OPEN','PAR ','HCHO','C5H8','OLE ','ALD2', &
& 'O3  ','NO2 ','OH  ','HO2 ','O   ','NO3 ', &
& 'NO  ','C2O3' &
                      & /)
  type(silam_species), dimension(num_species), public, save, target :: species_cbm4
  
  
contains

 
  subroutine set_rates_cbm4(temp, press, cos_theta, cloud_att)
    implicit none
    real, intent(in) :: temp, press, cos_theta, cloud_att
    
        RCONST(1) = (photo(1.108E-02,0.397,0.183,cos_theta)*cloud_att)
    RCONST(2) = (ARR2(1.4E+3,1175.0))
    RCONST(3) = (ARR2(1.8E-12,-1370.0))
    ! RCONST(4) = constant rate coefficient
    RCONST(5) = (ARR2(1.6E-13,687.0))
    RCONST(6) = (ARR2(2.2E-13,602.0))
    RCONST(7) = (ARR2(1.2E-13,-2450.0))
    RCONST(8) = (photo(5.219e-4,0.322,0.079,cos_theta)*cloud_att)
    RCONST(9) = (photo(8.978e-5,1.436,0.936,cos_theta)*cloud_att)
    RCONST(10) = (ARR2(1.9E+8,390.0))
    ! RCONST(11) = constant rate coefficient
    RCONST(12) = (ARR2(1.6E-12,-940.0))
    RCONST(13) = (ARR2(1.4E-14,-580.0))
    RCONST(14) = (photo(1.853e-1,0.189,0.112,cos_theta)*cloud_att)
    RCONST(15) = (ARR2(1.3E-11,250.0))
    RCONST(16) = (ARR2(2.5E-14,-1230.0))
    RCONST(17) = (ARR2(5.3E-13,256.0))
    ! RCONST(18) = constant rate coefficient
    RCONST(19) = (ARR2(3.5E+14,-10897.0))
    RCONST(20) = (ARR2(1.8E-20,530.0))
    ! RCONST(21) = constant rate coefficient
    RCONST(22) = (ARR2(4.5E-13,806.0))
    RCONST(23) = (1.511e-03*cos_theta*cloud_att)
    ! RCONST(24) = constant rate coefficient
    ! RCONST(25) = constant rate coefficient
    RCONST(26) = (ARR2(1.0E-12,713.0))
    RCONST(27) = (ARR2(5.1E-15,1000.0))
    RCONST(28) = (ARR2(3.7E-12,240.0))
    RCONST(29) = (ARR2(1.2E-13,749.0))
    RCONST(30) = (ARR2(4.8E+13,-10121.0))
    RCONST(31) = (ARR2(1.3E-12,380.0))
    RCONST(32) = (ARR2(5.9E-14,1150.0))
    RCONST(33) = (ARR2(2.2E-38,5800.0))
    RCONST(34) = (photo(1.057e-5,0.8,0.243,cos_theta)*cloud_att)
    RCONST(35) = (ARR2(3.1E-12,-187.0))
    ! RCONST(36) = constant rate coefficient
    ! RCONST(37) = constant rate coefficient
    RCONST(38) = (photo(4.866e-5,0.781,0.349,cos_theta)*cloud_att)
    RCONST(39) = (photo(6.79e-5,0.565,0.275,cos_theta)*cloud_att)
    RCONST(40) = (ARR2(3.0E-11,-1550.0))
    ! RCONST(41) = constant rate coefficient
    RCONST(42) = (ARR2(1.2E-11,-986.0))
    RCONST(43) = (ARR2(7.0E-12,250.0))
    ! RCONST(44) = constant rate coefficient
    RCONST(45) = (4.00E-06*cos_theta*cloud_att)
    RCONST(46) = (ARR2(3.51e-11,-180.0))
    RCONST(47) = (ARR2(2.62e-12,380.0))
    RCONST(48) = (ARR2(2.0e16,-13500.0))
    ! RCONST(49) = constant rate coefficient
    ! RCONST(50) = constant rate coefficient
    RCONST(51) = (ARR2(1.1E+2,-1710.0))
    ! RCONST(52) = constant rate coefficient
    RCONST(53) = (ARR2(1.0E+15,-8000.0))
    ! RCONST(54) = constant rate coefficient
    ! RCONST(55) = constant rate coefficient
    RCONST(56) = (ARR2(1.2E-11,-324.0))
    RCONST(57) = (ARR2(5.2E-12,504.0))
    RCONST(58) = (ARR2(1.4E-14,-2105.0))
    ! RCONST(59) = constant rate coefficient
    RCONST(60) = (ARR2(1.0E-11,-792.0))
    RCONST(61) = (ARR2(2.0E-12,411.0))
    RCONST(62) = (ARR2(1.3E-14,-2633.0))
    RCONST(63) = (ARR2(2.1E-12,322.0))
    ! RCONST(64) = constant rate coefficient
    ! RCONST(65) = constant rate coefficient
    ! RCONST(66) = constant rate coefficient
    ! RCONST(67) = constant rate coefficient
    ! RCONST(68) = constant rate coefficient
    RCONST(69) = (ARR2(1.7E-11,116.0))
    ! RCONST(70) = constant rate coefficient
    RCONST(71) = ((5.334E-05*cos_theta*cloud_att))
    RCONST(72) = (ARR2(5.4E-17,-500.0))
    ! RCONST(73) = constant rate coefficient
    RCONST(74) = ((1.654E-04*cos_theta*cloud_att))
    ! RCONST(75) = constant rate coefficient
    ! RCONST(76) = constant rate coefficient
    ! RCONST(77) = constant rate coefficient
    ! RCONST(78) = constant rate coefficient
    ! RCONST(79) = constant rate coefficient
    RCONST(80) = (ARR2(1.7E-14,1300.0))
    ! RCONST(81) = constant rate coefficient

  contains

    REAL(kind=sp) FUNCTION ARR2( A0,B0 )
      REAL A0,B0           
      ARR2 =  DBLE(A0) * EXP( DBLE(B0)/TEMP )              
    END FUNCTION ARR2

    real function photo(l, m, n, cos_theta)
      implicit none
      real, intent(in) :: l, n, m, cos_theta
      if (cos_theta < 1e-5) then
        photo = 0.0
        return
      end if
      photo = l * cos_theta**m * exp(-n / cos_theta)
    end function photo

  end subroutine set_rates_cbm4

  subroutine set_const_rates_cbm4()
    implicit none
        RCONST(4) = 9.3e-12
    RCONST(11) = 2.2e-10
    RCONST(18) = 1.3e-21
    RCONST(21) = 4.39999e-40
    RCONST(24) = 6.6e-12
    RCONST(25) = 1e-20
    RCONST(36) = 2.2e-13
    RCONST(37) = 1e-11
    RCONST(41) = 6.3e-16
    RCONST(44) = 2.5e-15
    RCONST(49) = 2e-12
    RCONST(50) = 6.5e-12
    RCONST(52) = 8.1e-13
    RCONST(54) = 1600
    RCONST(55) = 1.5e-11
    RCONST(59) = 7.7e-15
    RCONST(64) = 8.1e-12
    RCONST(65) = 4.2
    RCONST(66) = 4.1e-11
    RCONST(67) = 2.2e-11
    RCONST(68) = 1.4e-11
    RCONST(70) = 3e-11
    RCONST(73) = 1.7e-11
    RCONST(75) = 1.8e-11
    RCONST(76) = 9.6e-11
    RCONST(77) = 1.2e-17
    RCONST(78) = 3.2e-13
    RCONST(79) = 8.1e-12
    RCONST(81) = 6.8e-13
  end subroutine set_const_rates_cbm4

  subroutine set_fixed_cbm4(cnc_H2O)
    implicit none
    real, intent(in) :: cnc_H2O
    
    fix(1) = cnc_H2O

  end subroutine set_fixed_cbm4

  subroutine init_chemicals_cbm4()
    implicit none

    integer :: i

    do i = 1, num_species_cbm4
      call set_species(species_cbm4(i), fu_get_material_ptr(subst_names(i)), in_gas_phase)
      if (error) return
    end do

  end subroutine init_chemicals_cbm4
  
  subroutine inventory_int(flag_transf, &
                         & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                         & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                         & iClaimedSpecies)
    implicit none
    integer, intent(in) :: flag_transf ! the transformation id defined in an upper module
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: iEmis

    call addSpecies(speciesTransp, nSpeciesTransp, species_cbm4, num_species, .true.)
    if(error)return

    !
    ! Presence of emission species in the list of transported cbm species means ownership
    !
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_cbm4, num_species) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('cbm4 owns:' + fu_substance_name(speciesEmis(iEmis)))
          iClaimedSpecies(iEmis) = flag_transf
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(iEmis))
          call set_error('Cannot claim ownership because someone claimed it already',&
                       & 'inventory_int')
          return
        endif
      endif  ! emission species belongs to my transport
    end do  ! emission species

  end subroutine inventory_int
  
    !************************************************************************************

  subroutine register_species_int(speciesTransp, speciesShortlived, speciesAerosol,&
                              & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    !
    ! The integrator requires fixed indices for the species. At this
    ! point, we just check that the transport species follow this ordering.
    implicit none
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    
    ! Local variable
    integer :: i

    do i = 1, num_species
      if (.not. fu_index(species_cbm4(i), speciesTransp) == i) then
        call msg('Looking for:')
        call report(species_cbm4(i))
        call msg('Index in cbm4, transport:', i, fu_index(species_cbm4(i), speciesTransp))
        call set_error('Incompatible ordering of transport species', &
                     & 'register_species_int')
        return
      end if
    end do  ! species
    
  end subroutine register_species_int

end module cbm4_interface
