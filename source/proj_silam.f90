module proj_silam
  !
  !  Wrapper for proj4
  !

  !$ use OMP_LIB
  use, intrinsic :: ISO_C_BINDING, only: c_int, c_char, c_ptr, c_null_ptr, c_associated, c_size_t, c_loc, c_funptr, c_null_funptr
  use toolbox

  implicit none
  private

  public lonalt2proj
  public lonalt2proj_pt
  public test_proj

  private projFromStr
!
#ifdef USE_PROJ4
  interface
    function pj_init_plus(projstr) bind(C,name='pj_init_plus')
      USE ISO_C_BINDING
      IMPLICIT NONE
      type(C_PTR) :: pj_init_plus
      character(kind=C_CHAR) :: projstr(*)
    end function pj_init_plus

    function pj_transform(src, dst, point_count, point_offset, &
              & x, y, z) bind(C,name='pj_transform')
      USE ISO_C_BINDING
      IMPLICIT NONE
      integer(C_INT) :: pj_transform
      type(C_PTR), value :: src, dst
      integer(C_LONG), value :: point_count
      integer(C_INT), value :: point_offset
      type(C_PTR), value :: x, y, z
    end function pj_transform

  end interface
#endif

  character (len=*), public, parameter :: lonlat_proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" !EPSG:4326




  !!!! Registry of projections
  integer, parameter :: max_nbr_of_projections = 10
  integer :: nProjectionsKnown = int_missing
  type(c_ptr), dimension(max_nbr_of_projections) ::  known_projection_ptr = c_null_ptr
  character (len=fnlen), dimension(max_nbr_of_projections) :: known_projection_strings 
  
  

!  +o_lat_p=<latitude>
!Latitude of the North pole of the unrotated source CRS, expressed in the rotated geographic CRS.
!+o_lon_p=<longitude>
!Longitude of the North pole of the unrotated source CRS, expressed in the rotated geographic CRS.
!+lon_0=<value>
!Longitude of projection center. (Defaults to 0.0)
 

  contains

  integer function projFromStr(projstr)
    !returns number of projection

    character(len=*), intent(in) :: projstr
    character (len=*), parameter :: sub_name='projFromStr'


    integer :: iTmp

#ifdef USE_PROJ4
    projFromStr = int_missing
    if (nProjectionsKnown < 1) then
      !$OMP CRITICAL (proj_lock)
      if (nProjectionsKnown < 1) then ! Ues, same check inside CRITICAL
        known_projection_ptr(1) = pj_init_plus(lonlat_proj4//achar(0))
        known_projection_strings(1) = lonlat_proj4
        nProjectionsKnown = 1
      endif
      !$OMP END CRITICAL (proj_lock)
      if (.not. c_associated(known_projection_ptr(1))) then
        call set_error("pj_init_plus failed", sub_name)
        return
      endif
    endif

    do iTmp =  1, nProjectionsKnown
      if (known_projection_strings(iTmp) == projstr) exit
    enddo


    if (iTmp <=  nProjectionsKnown) then
      projFromStr = iTmp
      return
    endif

    !$OMP CRITICAL (proj_lock)
      do iTmp =  iTmp, nProjectionsKnown !!Could appear by now
        if (known_projection_strings(iTmp) == projstr) exit
      enddo

      if (iTmp <=  nProjectionsKnown) then !! Appeared while we were waiting for proj_lock
        projFromStr = iTmp
      elseif (iTmp > max_nbr_of_projections) then

        call msg_warning ("max_nbr_of_projections exceeded"//fu_str(max_nbr_of_projections), &
                                                                & sub_name)
        call msg("Projections registered so far:")
        do iTmp =  1, nProjectionsKnown
          call msg('"'//trim(known_projection_strings(iTmp))//'"')
        enddo
        call set_error ("max_nbr_of_projections exceeded", sub_name)

      else
        known_projection_ptr(iTmp) = pj_init_plus(trim(projstr)//achar(0))
        if (c_associated(known_projection_ptr(iTmp))) then
          known_projection_strings(iTmp) = trim(projstr)
          nProjectionsKnown = iTmp
          projFromStr = iTmp
        else
          call set_error("pj_init_plus failed", sub_name)
        endif
      endif
    !$OMP END CRITICAL (proj_lock)
#else
  call set_error("compiled without USE_PROJ4", sub_name)
#endif    
  end function projFromStr


  !******************************************************************

  subroutine lonalt2proj(projstr, x, y, point_count, direction)

    !
    ! Projects in-place lon and lat to projection
    !

    character(len=*), intent(in) :: projstr 
    integer, intent(in) :: point_count !! length of x and y
    integer, intent(in) :: direction
    real(8), dimension (:), target, intent(inout) :: x, y !! projection x,y on input, lon,lat on output
    character (len=*), parameter :: sub_name='lonalt2proj'

    integer :: projno, iStat
    integer(8) :: lcount

    lcount = point_count


#ifdef USE_PROJ4 
    projno = projFromStr(projstr)
    if (error) then
      call msg("projFromStr failed for:")
      call msg('"'//trim(projstr)//'"')
      call set_error("projFromStr failed", sub_name)
      return
    endif

    if (direction == forwards) then
      iStat = pj_transform(known_projection_ptr(1), known_projection_ptr(projno), &
                  &   lcount, 1, C_LOC(x(1)), C_LOC(y(1)), C_NULL_PTR)
    elseif (direction == backwards) then
      iStat = pj_transform(known_projection_ptr(projno), known_projection_ptr(1), &
                  &   lcount, 1, C_LOC(x(1)), C_LOC(y(1)), C_NULL_PTR)
    else
      call set_error("Unknown direction: "//fu_str(direction),sub_name)
      return
    endif
    if (iStat /= 0) then
      call set_error("pj_transform failed",sub_name)
    endif
#else
  call set_error("compiled without USE_PROJ4", sub_name)
  !
  ! FIXME Handling of rotated lon-lat should be here
  !
#endif

  end subroutine lonalt2proj

  subroutine lonalt2proj_pt(projstr, x, y, direction)
    !! HUOM!!! pj_transform treats all angles in radians!

    !
    ! Projects in-place lon and lat to projection
    !

    character(len=*), intent(in) :: projstr 
    integer, intent(in) :: direction
    real(8), intent(inout), target :: x, y !! projection x,y on input, lon,lat on output
    character (len=*), parameter :: sub_name='lonalt2proj'
    integer(8), parameter :: one = 1

    integer :: projno, iStat


#ifdef USE_PROJ4 
    projno = projFromStr(projstr)
    if (error) then
      call msg("projFromStr failed for:")
      call msg('"'//trim(projstr)//'"')
      call set_error("projFromStr failed", sub_name)
      return
    endif

    if (direction == forwards) then
      iStat = pj_transform(known_projection_ptr(1), known_projection_ptr(projno), &
                  &   one, 1, C_LOC(x), C_LOC(y), C_NULL_PTR)
    elseif (direction == backwards) then
      iStat = pj_transform(known_projection_ptr(projno), known_projection_ptr(1), &
                  &   one, 1, C_LOC(x), C_LOC(y), C_NULL_PTR)
    else
      call set_error("Unknown direction: "//fu_str(direction),sub_name)
      return
    endif
    if (iStat /= 0) then
      call set_error("pj_transform failed",sub_name)
    endif
#else
  call set_error("compiled without USE_PROJ4", sub_name)
  !
  ! FIXME Handling of rotated lon-lat bypassing proj should be here
  !
#endif

  end subroutine lonalt2proj_pt



  subroutine test_proj()
    implicit none

    integer, parameter :: np = 4
    real(8), dimension(np) :: lats, lons
    integer :: i

    character (len=*), parameter :: hirlam_rll_proj4 = &
           & "+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0"
    character (len=*), parameter :: lcc_example="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"

    !! corners of EU meteo in projected lon-lat
    lons(1:np) = (/-26., -26., 40., 40./)*ddeg_to_rad
    lats(1:np) = (/-35., 22.5, 22.5, -35./)*ddeg_to_rad

    call msg("rotated-pole coordinates lon, lat")
    call msg(hirlam_rll_proj4)
    do i=1,np
      call msg(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
    enddo

    call lonalt2proj(hirlam_rll_proj4, lons, lats, np, backwards)

    call msg("Same in WGS84 coordinates lon, lat")
    call msg(lonlat_proj4)
    do i=1,np
      call msg(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
    enddo

    call lonalt2proj(hirlam_rll_proj4, lons, lats, np, forwards)

    call msg("And back! lon, lat")
    call msg(hirlam_rll_proj4)
    do i=1,np
      call msg(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
    enddo
    call msg("Done")


      


    

    


  end subroutine test_proj



end module proj_silam
