module proj_silam
  !
  !  Wrapper for proj4
  !

  !$ use OMP_LIB
  use, intrinsic :: ISO_C_BINDING, only: c_int, c_char,c_double, c_ptr, c_null_ptr, c_associated, c_size_t, c_loc, c_funptr, c_null_funptr,c_sizeof
  use toolbox

  implicit none
  private

  !public proj_trans!(generalized)
  public ll2proj
  public ll2proj_pt
  public test_proj

  !private projFromStr
 
#ifdef USE_PROJ6
   interface
     function proj_context_create() bind(C,name='proj_context_create')
       use iso_c_binding
       implicit none
       type(c_ptr) :: proj_context_create
     end function
     subroutine proj_context_destroy(ctx) bind(C,name='proj_context_destroy')
       use iso_c_binding
       implicit none
       type(c_ptr),value :: ctx
     end subroutine
     function proj_create_crs_to_crs(ctx,s_srs,t_srs,area) bind(C,name='proj_create_crs_to_crs')
       use iso_c_binding
       implicit none
       type(C_PTR)            :: proj_create_crs_to_crs
       character(kind=C_CHAR) :: s_srs(*)
       character(kind=C_CHAR) :: t_srs(*)
       type(C_PTR), value     :: ctx, area
     end function
     subroutine proj_destroy(PJ) bind(C,name='proj_destroy')
       use iso_c_binding
       implicit none
       type(c_ptr), value :: PJ
     end subroutine
     function proj_trans_generic(pj,dir,x,sx,nx,y,sy,ny,z,sz,nz,t,st,nt) bind(C,name='proj_trans_generic')
       use iso_c_binding
       implicit none
       integer(c_int)       :: proj_trans_generic
       type(c_ptr), value   :: pj
       integer(c_int),value :: dir
       type(c_ptr),value    :: x,y,z,t
       integer(c_int),value :: sx,sy,sz,st
       integer(c_int),value :: nx,ny,nz,nt
     end function
   end interface
#else
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
       integer(C_INT), value :: point_count
      ! integer(C_LONG), value :: point_count
       integer(C_INT), value :: point_offset
       type(C_PTR), value :: x, y, z
     end function pj_transform
  
   end interface
#endif
#endif

character (len=*), public, parameter :: lonlat_proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" !EPSG:4326
  
contains

subroutine proj_trans(s_srs,t_srs, x, y, np, direction)
    !
    !"Any 2 any" coordinate system transofmration:
    ! Projects IN-PLACE arrays of coordinates from one projection to another.
    !
    implicit none
    character(len=*), intent(in) :: s_srs,t_srs           !source & target spatial reference system (srs) (proj4str, epsgs, wkt, etc.)
    integer, intent(in) :: np                             ! number of points (length of x or y)
    integer, intent(in) :: direction                      ! forward / backwards
    real(8), dimension (:), target, intent(inout) :: x, y ! coordinates to transform
    character (len=*), parameter :: sub_name='proj_trans'
    
#ifdef USE_PROJ6
  type(c_ptr) :: CTX = c_null_ptr
  type(c_ptr) ::  PJ = c_null_ptr
  type(c_ptr) ::  area= c_null_ptr
  integer     :: dir, iStat
  integer     :: sp              !size (in bytes of coordinates datatype)

  print*,"Usando PROJ6.."  

  CTX=proj_context_create()
  PJ=proj_create_crs_to_crs(CTX, trim(s_srs)//achar(0),           &
                           &     trim(t_srs)//achar(0), C_NULL_PTR)

  if      ( direction == forwards ) then
     dir=1;
  else if ( direction == backwards) then
     dir=-1;
  else
     call set_error("Unknown direction: "//fu_str(direction),sub_name)
     return
  endif

  sp=c_sizeof(x(1))         !Step length (bytes) between consecutive elements of the array
  iStat=proj_trans_generic( PJ, dir,             &
                          & c_loc(x)   , sp, np, &
                          & c_loc(y)   , sp, np, &
                          & C_NULL_PTR ,  0, 0 , &
                          & C_NULL_PTR ,  0, 0  )
  if (iStat /= np) then                             
    print*,"Number of transformations failed = ",np-iStat
    call set_error("proj_trans_generic failed",sub_name)
  endif

  !Clean:
  call proj_destroy(PJ)
  call proj_context_destroy(CTX)

#else
#ifdef USE_PROJ4 
   integer     :: iStat
   type(c_ptr) :: PJ1 = c_null_ptr
   type(c_ptr) :: PJ2 = c_null_ptr

   print*,"Usando PROJ4.."  

   x=x*ddeg_to_rad  !convert degrees to radians
   y=y*ddeg_to_rad  !convert degrees to radians

   PJ1= pj_init_plus(trim(s_srs)//achar(0))
   PJ2= pj_init_plus(trim(t_srs)//achar(0))

   if (direction == forwards) then
      iStat = pj_transform(PJ1, PJ2, np, 1, C_LOC(x), C_LOC(y), C_NULL_PTR)
   elseif (direction == backwards) then
      iStat = pj_transform(PJ2, PJ1, np, 1, C_LOC(x), C_LOC(y), C_NULL_PTR)
   else
     call set_error("Unknown direction: "//fu_str(direction),sub_name)
     return
   endif
   if (iStat /= 0) then
     call set_error("pj_transform failed",sub_name)
   endif
   x=x*drad_to_deg     !return to degrees!
   y=y*drad_to_deg     !return to degrees!

#else
     call set_error("compiled without USE_PROJ4 nor USE_PROJ6", sub_name)
     ! FIXME Handling of rotated lon-lat should be here
#endif
#endif
end subroutine



!!WRAPERS!!
! Subroutines to handle latlon (ll) transformations and single value transformation

!ll2proj:
  subroutine ll2proj(t_srs, x, y, np, direction)
    !
    ! Projects in-place lon and lat to any other projection
    implicit none
    character(len=*), intent(in) :: t_srs !source & target spatial reference system (srs) (proj4str, epsgs, wkt, etc.)
    integer, intent(in) :: np             ! number of points (length of x or y)
    integer, intent(in) :: direction      ! forward / backwards
    real(8), dimension (:), target, intent(inout) :: x, y ! coordinates to transform
    character (len=*), parameter :: sub_name='ll2proj'
    call proj_trans(lonlat_proj4,t_srs,x,y,np,direction)
  end subroutine

!ll2proj_pt:
  subroutine ll2proj_pt(t_srs, x, y, direction)
    !      
    ! same that ll2proj but for a single value
    implicit none
    character(len=*), intent(in) :: t_srs  !source & target spatial reference system (srs) (proj4str, epsgs, wkt, etc.)
    integer, intent(in) :: direction       !forward / backwards
    real(8), target, intent(inout) :: x, y !coordinates to transform
    real(8), dimension(1) :: xa, ya        !coordinates to transform (rank-1 array form)
    character (len=*), parameter :: sub_name='ll2proj_pt'
    xa=x;ya=y
    call proj_trans(lonlat_proj4,t_srs,xa,ya,1,direction)
    x=xa(1);y=ya(1)
  end subroutine
!ll2proj_pt:
  subroutine proj_trans_pt(s_srs, t_srs, x, y, direction)
    !      
    ! same that ll2proj but for a single value
    implicit none
    character(len=*), intent(in) :: s_srs, t_srs  !source & target spatial reference system (srs) (proj4str, epsgs, wkt, etc.)
    integer, intent(in) :: direction       !forward / backwards
    real(8), target, intent(inout) :: x, y !coordinates to transform
    real(8), dimension(1) :: xa, ya        !coordinates to transform (rank-1 array form)
    character (len=*), parameter :: sub_name='ll2proj_pt'
    xa=x;ya=y
    call proj_trans(s_srs,t_srs,xa,ya,1,direction)
    x=xa(1);y=ya(1)
  end subroutine

!TEST:
  subroutine test_proj()
     implicit none

     integer, parameter :: np = 4
     real(8), dimension(np) :: lats, lons
     integer :: i
     real(8) :: x, y
     character (len=*), parameter :: hirlam_rll_proj4 ="+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0"
     character (len=*), parameter ::       lcc_example="+proj=lcc  +lon_0=-90 +lat_1=33 +lat_2=45"
     character (len=*), parameter ::      merc_example="+proj=merc +lon_0=-90 +lat_1=33 +lat_2=45"

     !! corners of EU meteo in projected lon-lat
     lons(1:np) = (/-26., -26., 40., 40./)  
     lats(1:np) = (/-35., 22.5, 22.5, -35./)

     call msg("rotated-pole coordinates lon, lat")
     call msg(hirlam_rll_proj4)
     do i=1,np
       call msg(fu_str(i), lons(i),lats(i))
     enddo

     call ll2proj(hirlam_rll_proj4, lons, lats, np, backwards)

     call msg("Same in WGS84 coordinates lon, lat")
     call msg(lonlat_proj4)
     do i=1,np
       call msg(fu_str(i), lons(i),lats(i))
     enddo

     call ll2proj(hirlam_rll_proj4, lons, lats, np, forwards)

     call msg("And back! lon, lat")
     call msg(hirlam_rll_proj4)
     do i=1,np
       call msg(fu_str(i), lons(i),lats(i))
     enddo
     call msg("Done")

     !Aditional test: single value
     print*,"(t) Aditional test 1: single value transformation.."
     x=lons(1); y=lats(1)
     print*,"original lons(1),y=lats(1) ",lons(1),lats(1)
     call ll2proj_pt(hirlam_rll_proj4,x,y,backwards)
     print*,"transformed to x,y         ",x,y
     call ll2proj_pt(hirlam_rll_proj4,x,y,forwards)
     print*,"back to lons(1),y=lats(1): ",x,y
     call msg("Done")

     !Aditional test: any 2 any
     print*,"(t) Aditional test 2: any to any transformation.."
     
     print*,"Original lons(1),y=lats(1) ",lons(1),lats(1)
     x=lons(1); y=lats(1)
     call proj_trans_pt(hirlam_rll_proj4, lcc_example, x, y, forwards )
     print*,"transformed to lcc         ",x,y
     call proj_trans_pt(hirlam_rll_proj4, lcc_example, x, y, backwards)
     print*,"back to lons(1),y=lats(1): ",x,y
     call msg("Done")

  end subroutine test_proj



end module proj_silam
