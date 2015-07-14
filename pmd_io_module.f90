!*******************************************************************************
!  _____  __  __ _____         _____ ____  
! |  __ \|  \/  |  __ \       |_   _/ __ \
! | |__) | \  / | |  | |        | || |  | |
! |  ___/| |\/| | |  | |        | || |  | |
! | |    | |  | | |__| |       _| || |__| |
! |_|    |_|  |_|_____/ _____ |_____\____/ _module
!                      |_____|
!
! Input/Output module for the Porta Model Data (PMD)
!
!*******************************************************************************
!   Fortran 90
!   Flavio Calvo:               Geneva, Locarno
!   2014-09-24
!*******************************************************************************
!
!------*************------------------------------------------------------------
module pmd_io_module
!-------------------------------------------------------------------------------
! NAME:
!   pmd_io_module ('porta_model_data_input_output_module')
!
! PURPOSE:
!   Provide Fortran I/O routines for the PMD type specifications
!
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   use pmd_io_module
!
! VARIABLES (public):
!   None
!
! TYPES:
!   pmd_header_type:            Description             Default value
!                               -----------             -------------
!     magic_str                 Magic string            'portapmd'
!     endianness                0: little endian        0
!                               1: big endian
!     int_sz                    Size of integers        4
!                               (in bytes)
!     db_sz                     Size of doubles         8
!                               (in bytes)
!     pmd_version               PMD version             2
!     localtime                 Local time (ymdhms)     (/0,0,0,0,0,0/)
!                               year since 1900
!                               month (0-11)
!                               day (1-31)
!                               hour (0-23)
!                               minute (0-59)
!                               second (0-59)
!     periodic                  0: non-periodic         (/0,0/)
!                               1: periodic
!                               (in X and Y)
!     domain_sz                 Domain size             (/0.,0.,0./)
!     domain_origin             Domain origin           (/0.,0.,0./)
!     dimensions                Dimensions              (/0,0,0/)
!     xc1                       X node coordinates      (array of) 0.
!                               (increasing order)
!     xc2                       Y node coordinates      (array of) 0.
!                               (increasing order)
!     xc3                       Z node coordinates      (array of) 0.
!                               (increasing order)
!     n_radtheta                Number of inclination   0
!                               angles per octant
!     n_radphi                  Number of azimuthal     0
!                               angles per octant
!     atomic_module             Atomic module name      ''
!     pmd_comment               Model comment           ''
!     module_hd_sz              Size of module header   0
!     unused_var                Undefined               0
!
! SUBROUTINES (high-level only, public, '<' / '>' indicates input/output
!              parameters):
!   pmd_openrd(<unit, <file, >pmd_header, >stat=stat)
!   pmd_openwr(<unit, <file, <pmd_header, >stat=stat)
!   pmd_openap(<unit, <file, >pmd_header, >stat=stat)
!   pmd_rd_box(<unit, >box)
!   pmd_append_box(<unit, <box)
!   pmd_close(<unit)
!
! SUBROUTINES (low-level and high-level, public, '<' / '>' as above,
!              '*' mandatory parameter, with documentation):
!   pmd_openrd(<unit*, <file*, >pmd_header, >stat)
!     -> Open file read-only, read pmd_header, output number of boxes in the
!        node data region (stat variable). If stat<0, then -stat is the output
!        iostat value from the Fortran open function. Position file at the
!        beginning of the first box.
!        Additional stat/error codes:
!          * -500: File format is not PMD
!          * -501: Porta header is not properly formatted
!          * -502: Module header is not properly formatted
!          * -503: Node data is not properly formatted
!        /!\ File is not closed after one of these exceptions is raised.
!
!   pmd_openwr(<unit*, <file*, <pmd_header, <convert, >stat)
!     -> Open file write-only, write header if pmd_header specified, stat as
!        above. If specified, convert is passed to Fortran open statement.
!        If not, convert will be choosen according to the specified endianness
!        in pmd_header. If neither pmd_header nor convert are specified,
!        convert will take 'native' as default value.
!
!   pmd_openap(<unit*, <file*, >pmd_header, >stat)
!     -> Open file read-write and position at the end of the file. The rest is
!        as in pmd_openrd. File must exist.
!        Same additional stat/error codes as for pmd_openrd.
!
!   pmd_types_sz(<unit*, int_sz*, db_sz*, stat)
!     -> Return the size of integers and doubles in opened file. File position
!        remains unchanged, int_sz and db_sz are of type integer*1.
!        stat is either the value returned by fseek if there is an error or
!        by read in iostat.
!        This works even if the file/header is (partially) corrupted.
!
!   pmd_rd_header(<unit*, >pmd_header*)
!     -> Read header from opened file. Position is unchanged.
!
!   pmd_wr_header(<unit*, <pmd_header*)
!     -> Write header to opened file. Position is unchanged.
!
!   pmd_rd_box(<unit*, >box*, >iostat)
!     -> Read box to flexible type variable (array of reals/double precision).
!        File should be opened and at appropriate position.
!        iostat is the corresponding value returned by Fortran read statement.
!        This allows to trigger EOF.
!
!   pmd_append_box(<unit*, <box*)
!     -> Write box from flexible type variable (array of reals/double
!        precision).
!        Box is appened at the end and position remains unchanged.
!
!   pmd_set_sz(<pmd_int_sz*, <pmd_db_sz*)
!     -> Set size of integers and doubles. pmd_rd and pmd_wr will read/write
!        data accordingly. Calls to pmd_rd_header and pmd_wr_header will
!        automatically change the value of those variables according to the
!        content of pmd_header.
!
!   pmd_rd(<unit*, >value*, >iostat)
!     -> Read variable of flexible type from opened file. 3D arrays of doubles
!        also accept iostat parameter, which corresponds to the output of the
!        Fortran read function (in order to detect EOF). One needs to make
!        sure pmd_int_sz and pmd_db_sz were properly set.
!
!   pmd_wr(<unit*, <value*)
!     -> Write variable of flexible type to opened file. One need to make sure
!        pmd_int_sz and pmd_db_sz were properly set.
!
!   pmd_seek(<unit*, <pmd_header*, <pos*)
!     -> Go to some position specified by pos in opened file:
!          pos = 0: module header
!          pos > 0: box #pos in node data section
!          pos < 0: box #(N-|pos|+1) in node data section, with N the total
!                   number of boxes
!
!   pmd_set_localtime(<pmd_header*)
!     -> Set localtime in pmd_header when provided value is not valid.
!
!   (integer) pmd_header_sz(<pmd_header*)
!     -> Returns size of Porta header in bytes.
!
!   (integer) pmd_header_sz(<unit*)
!     -> Returns size of Porta header in bytes, even if header is corrupted
!        and/or does not actually have this size.
!        This function can be used to check that header indeed has the
!        appropriate size.
!
!   (integer) pmd_box_sz(<pmd_header*)
!     -> Returns size of a box (in node data section) in bytes.
!
! MODIFICATION HISTORY:
!   2014-09-24 (F. C. Locarno) Written
!-------------------------------------------------------------------------------
!
implicit none
private
!
integer*1                               :: integer_size=4, double_size=8
character(len=8)                        :: pmd_magic_str='portapmd'
!
type, public                            :: pmd_header_type
  character(len=8)                      :: magic_str='portapmd'
  integer*1                             :: endianness=0
  integer*1                             :: int_sz=4, db_sz=8
  integer                               :: pmd_version=2
  integer, dimension(6)                 :: localtime=(/0,0,0,0,0,0/)
  integer*1, dimension(2)               :: periodic=(/0,0/)
  double precision, dimension(3)        :: domain_sz=(/0.,0.,0./), &
                                           domain_origin=(/0.,0.,0./)
  integer, dimension(3)                 :: dimensions=(/0,0,0/)
  double precision, dimension(8192)     :: xc1=0., xc2=0., xc3=0.
  integer                               :: n_radtheta=0, n_radphi=0
  character(len=1023)                   :: atomic_module=''
  character(len=4096)                   :: pmd_comment=''
  integer                               :: module_hd_sz=0
  integer                               :: unused_var=0
end type pmd_header_type
!
interface pmd_wr
  module procedure pmd_wr_int8, pmd_wr_char, pmd_wr_int, &
                   pmd_wr_real, pmd_wr_db, &
                   pmd_wr_int8_1d, pmd_wr_char_1d, pmd_wr_int_1d, &
                   pmd_wr_real_1d, pmd_wr_db_1d, &
                   pmd_wr_int8_2d, pmd_wr_char_2d, pmd_wr_int_2d, &
                   pmd_wr_real_2d, pmd_wr_db_2d, &
                   pmd_wr_int8_3d, pmd_wr_char_3d, pmd_wr_int_3d, &
                   pmd_wr_real_3d, pmd_wr_db_3d
end interface pmd_wr
!
interface pmd_rd
  module procedure pmd_rd_int8, pmd_rd_char, pmd_rd_int, &
                   pmd_rd_real, pmd_rd_db, &
                   pmd_rd_int8_1d, pmd_rd_char_1d, pmd_rd_int_1d, &
                   pmd_rd_real_1d, pmd_rd_db_1d, &
                   pmd_rd_int8_2d, pmd_rd_char_2d, pmd_rd_int_2d, &
                   pmd_rd_real_2d, pmd_rd_db_2d, &
                   pmd_rd_int8_3d, pmd_rd_char_3d, pmd_rd_int_3d, &
                   pmd_rd_real_3d, pmd_rd_db_3d
end interface pmd_rd
!
interface pmd_append_box
  module procedure pmd_append_box_real_3d, pmd_append_box_db_3d
end interface pmd_append_box
!
interface pmd_rd_box
  module procedure pmd_rd_box_real_3d, pmd_rd_box_db_3d
end interface pmd_rd_box
!
interface pmd_header_sz
  module procedure pmd_header_sz_hd, pmd_header_sz_unit
end interface
!
public pmd_openwr, pmd_openap, pmd_openrd, pmd_close, &
       pmd_set_sz, pmd_wr, pmd_rd, pmd_wr_header, pmd_rd_header, &
       pmd_append_box, pmd_rd_box, pmd_types_sz, pmd_set_localtime, &
       pmd_header_sz, pmd_box_sz, pmd_seek, pmd_fseek, pmd_ftell, &
       pmd_seek_node, pmd_append_node, pmd_rd_node
!
contains
!
subroutine pmd_fseek(unit, pos, whence)
  !
  implicit none
  !
  integer, intent(in)                   :: unit, pos, whence
  character(len=80)                     :: act
  integer                               :: pos0
  !
  select case (whence)
    case (0)
      pos0 = pos+1
    case (1)
      inquire(unit, pos=pos0)
      pos0 = pos+pos0
    case (2)
      inquire(unit, size=pos0)
      pos0 = pos0+1-pos
    case default
      pos0 = pos+1
  end select
  !
  inquire(unit, action=act)
  if(trim(act) == 'READ' .or. trim(act) == 'read') then
    if(pos0 == 0) then
      rewind(unit)
    else
      read(unit, pos=pos0)
    endif
  else
    if(pos0 == 0) then
      rewind(unit)
    else
      write(unit, pos=pos0)
    endif
  endif
end subroutine pmd_fseek
!
subroutine pmd_ftell(unit, pos)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(out)                  :: pos
  !
  inquire(unit, pos=pos)
  pos=pos-1
end subroutine pmd_ftell
!
subroutine pmd_set_sz(pmd_int_sz, pmd_db_sz)
  !
  implicit none
  !
  integer, intent(in)                   :: pmd_int_sz, pmd_db_sz
  !
  integer_size = pmd_int_sz
  double_size = pmd_db_sz
end subroutine pmd_set_sz
!
subroutine pmd_wr_int8(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, intent(in)                 :: i
  !
  write(unit), i
end subroutine pmd_wr_int8
!
subroutine pmd_wr_char(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), intent(in)              :: c
  integer*1                             :: i
  integer                               :: pos
  !
  do pos=1,len(c)
    i=int(ichar(c(pos:pos)),1)
    write(unit), i
  end do
end subroutine pmd_wr_char
!
subroutine pmd_wr_int(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(in)                   :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  !
  select case (integer_size)
    case (1)
      i1=int(u,1)
      write(unit), i1
    case (2)
      i2=int(u,2)
      write(unit), i2
    case (4)
      i4=int(u,4)
      write(unit), i4
    case (8)
      i8=int(u,8)
      write(unit), i8
    case default
      i4=int(u,4)
      write(unit), i4
  end select
end subroutine pmd_wr_int
!
subroutine pmd_wr_real(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, intent(in)                      :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  !
  select case (double_size)
    case (4)
      f4=real(r,4)
      write(unit), f4
    case (8)
      f8=real(r,8)
      write(unit), f8
    case (16)
      f16=real(r,16)
      write(unit), f16
    case default
      f8=real(r,8)
      write(unit), f8
  end select
end subroutine pmd_wr_real
!
subroutine pmd_wr_db(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, intent(in)          :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  !
  select case (double_size)
    case (4)
      f4=real(d,4)
      write(unit), f4
    case (8)
      f8=real(d,8)
      write(unit), f8
    case (16)
      f16=real(d,16)
      write(unit), f16
    case default
      f8=real(d,8)
      write(unit), f8
  end select
end subroutine pmd_wr_db
!
subroutine pmd_wr_int8_1d(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, dimension(:), intent(in)   :: i
  !
  write(unit), i
end subroutine pmd_wr_int8_1d
!
subroutine pmd_wr_char_1d(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), dimension(:), &
                             intent(in) :: c
  integer*1                             :: i
  integer                               :: item, pos
  !
  do item=1,size(c,1)
    do pos=1,len(c(item))
      i=int(ichar(c(item)(pos:pos)),1)
      write(unit), i
    end do
  end do
end subroutine pmd_wr_char_1d
!
subroutine pmd_wr_int_1d(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, dimension(:), intent(in)     :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  integer                               :: item
  !
  select case (integer_size)
    case (1)
      do item=1,size(u,1)
        i1=int(u(item),1)
        write(unit), i1
      end do
    case (2)
      do item=1,size(u,1)
        i2=int(u(item),2)
        write(unit), i2
      end do
    case (4)
      do item=1,size(u,1)
        i4=int(u(item),4)
        write(unit), i4
      end do
    case (8)
      do item=1,size(u,1)
        i8=int(u(item),8)
        write(unit), i8
      end do
    case default
      do item=1,size(u,1)
        i4=int(u(item),4)
        write(unit), i4
      end do
  end select
end subroutine pmd_wr_int_1d
!
subroutine pmd_wr_real_1d(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:), intent(in)        :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: item
  !
  select case (double_size)
    case (4)
      do item=1,size(r,1)
        f4=real(r(item),4)
        write(unit), f4
      end do
    case (8)
      do item=1,size(r,1)
        f8=real(r(item),8)
        write(unit), f8
      end do
    case (16)
      do item=1,size(r,1)
        f16=real(r(item),16)
        write(unit), f16
      end do
    case default
      do item=1,size(r,1)
        f8=real(r(item),8)
        write(unit), f8
      end do
  end select
end subroutine pmd_wr_real_1d
!
subroutine pmd_wr_db_1d(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:), &
                             intent(in) :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: item
  !
  select case (double_size)
    case (4)
      do item=1,size(d,1)
        f4=real(d(item),4)
        write(unit), f4
      end do
    case (8)
      do item=1,size(d,1)
        f8=real(d(item),8)
        write(unit), f8
      end do
    case (16)
      do item=1,size(d,1)
        f16=real(d(item),16)
        write(unit), f16
      end do
    case default
      do item=1,size(d,1)
        f8=real(d(item),8)
        write(unit), f8
      end do
  end select
end subroutine pmd_wr_db_1d
!
subroutine pmd_wr_int8_2d(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, dimension(:,:), intent(in) :: i
  !
  write(unit), i
end subroutine pmd_wr_int8_2d
!
subroutine pmd_wr_char_2d(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), dimension(:,:), &
                             intent(in) :: c
  integer*1                             :: i
  integer                               :: j1, j2, pos
  !
  do j2=1,size(c,2)
    do j1=1,size(c,1)
      do pos=1,len(c(j1,j2))
        i=int(ichar(c(j1,j2)(pos:pos)),1)
        write(unit), i
      end do
    end do
  end do
end subroutine pmd_wr_char_2d
!
subroutine pmd_wr_int_2d(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, dimension(:,:), intent(in)   :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  integer                               :: j1, j2
  !
  select case (integer_size)
    case (1)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          i1=int(u(j1,j2),1)
          write(unit), i1
        end do
      end do
    case (2)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          i2=int(u(j1,j2),2)
          write(unit), i2
        end do
      end do
    case (4)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          i4=int(u(j1,j2),4)
          write(unit), i4
        end do
      end do
    case (8)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          i8=int(u(j1,j2),8)
          write(unit), i8
        end do
      end do
    case default
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          i4=int(u(j1,j2),4)
          write(unit), i4
        end do
      end do
  end select
end subroutine pmd_wr_int_2d
!
subroutine pmd_wr_real_2d(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:,:), intent(in)      :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2
  !
  select case (double_size)
    case (4)
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          f4=real(r(j1,j2),4)
          write(unit), f4
        end do
      end do
    case (8)
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          f8=real(r(j1,j2),8)
          write(unit), f8
        end do
      end do
    case (16)
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          f16=real(r(j1,j2),16)
          write(unit), f16
        end do
      end do
    case default
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          f8=real(r(j1,j2),8)
          write(unit), f8
        end do
      end do
  end select
end subroutine pmd_wr_real_2d
!
subroutine pmd_wr_db_2d(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:), &
                             intent(in) :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2
  !
  select case (double_size)
    case (4)
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          f4=real(d(j1,j2),4)
          write(unit), f4
        end do
      end do
    case (8)
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          f8=real(d(j1,j2),8)
          write(unit), f8
        end do
      end do
    case (16)
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          f16=real(d(j1,j2),16)
          write(unit), f16
        end do
      end do
    case default
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          f8=real(d(j1,j2),8)
          write(unit), f8
        end do
      end do
  end select
end subroutine pmd_wr_db_2d
!
subroutine pmd_wr_int8_3d(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, dimension(:,:,:), &
                             intent(in) :: i
  !
  write(unit), i
end subroutine pmd_wr_int8_3d
!
subroutine pmd_wr_char_3d(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), dimension(:,:,:), &
                             intent(in) :: c
  integer*1                             :: i
  integer                               :: j1, j2, j3, pos
  !
  do j3=1,size(c,3)
    do j2=1,size(c,2)
      do j1=1,size(c,1)
        do pos=1,len(c(j1,j2,j3))
          i=int(ichar(c(j1,j2,j3)(pos:pos)),1)
          write(unit), i
        end do
      end do
    end do
  end do
end subroutine pmd_wr_char_3d
!
subroutine pmd_wr_int_3d(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, dimension(:,:,:), intent(in) :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  integer                               :: j1, j2, j3
  !
  select case (integer_size)
    case (1)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            i1=int(u(j1,j2,j3),1)
            write(unit), i1
          end do
        end do
      end do
    case (2)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            i2=int(u(j1,j2,j3),2)
            write(unit), i2
          end do
        end do
      end do
    case (4)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            i4=int(u(j1,j2,j3),4)
            write(unit), i4
          end do
        end do
      end do
    case (8)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            i8=int(u(j1,j2,j3),8)
            write(unit), i8
          end do
        end do
      end do
    case default
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            i4=int(u(j1,j2,j3),4)
            write(unit), i4
          end do
        end do
      end do
  end select
end subroutine pmd_wr_int_3d
!
subroutine pmd_wr_real_3d(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:,:,:), intent(in)    :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2, j3
  !
  select case (double_size)
    case (4)
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            f4=real(r(j1,j2,j3),4)
            write(unit), f4
          end do
        end do
      end do
    case (8)
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            f8=real(r(j1,j2,j3),8)
            write(unit), f8
          end do
        end do
      end do
    case (16)
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            f16=real(r(j1,j2,j3),16)
            write(unit), f16
          end do
        end do
      end do
    case default
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            f8=real(r(j1,j2,j3),8)
            write(unit), f8
          end do
        end do
      end do
  end select
end subroutine pmd_wr_real_3d
!
subroutine pmd_wr_db_3d(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:,:), &
                             intent(in) :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2, j3
  !
  select case (double_size)
    case (4)
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            f4=real(d(j1,j2,j3),4)
            write(unit), f4
          end do
        end do
      end do
    case (8)
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            f8=real(d(j1,j2,j3),8)
            write(unit), f8
          end do
        end do
      end do
    case (16)
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            f16=real(d(j1,j2,j3),16)
            write(unit), f16
          end do
        end do
      end do
    case default
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            f8=real(d(j1,j2,j3),8)
            write(unit), f8
          end do
        end do
      end do
  end select
end subroutine pmd_wr_db_3d
!
subroutine pmd_rd_int8(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, intent(out)                :: i
  !
  read(unit), i
end subroutine pmd_rd_int8
!
subroutine pmd_rd_char(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), intent(inout)           :: c
  integer*1                             :: i
  integer                               :: pos
  !
  do pos=1,len(c)
    read(unit), i
    c(pos:pos)=achar(int(i))
  end do
end subroutine pmd_rd_char
!
subroutine pmd_rd_int(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(out)                  :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  !
  select case (integer_size)
    case (1)
      read(unit), i1
      u=int(i1)
    case (2)
      read(unit), i2
      u=int(i2)
    case (4)
      read(unit), i4
      u=int(i4)
    case (8)
      read(unit), i8
      u=int(i8)
    case default
      read(unit), i4
      u=int(i4)
  end select
end subroutine pmd_rd_int
!
subroutine pmd_rd_real(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, intent(out)                     :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  !
  select case (double_size)
    case (4)
      read(unit), f4
      r=real(f4)
    case (8)
      read(unit), f8
      r=real(f8)
    case (16)
      read(unit), f16
      r=real(f16)
    case default
      read(unit), f8
      r=real(f8)
  end select
end subroutine pmd_rd_real
!
subroutine pmd_rd_db(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, intent(out)         :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  !
  select case (double_size)
    case (4)
      read(unit), f4
      d=dble(f4)
    case (8)
      read(unit), f8
      d=dble(f8)
    case (16)
      read(unit), f16
      d=dble(f16)
    case default
      read(unit), f8
      d=dble(f8)
  end select
end subroutine pmd_rd_db
!
subroutine pmd_rd_int8_1d(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, dimension(:), intent(out)  :: i
  !
  read(unit), i
end subroutine pmd_rd_int8_1d
!
subroutine pmd_rd_char_1d(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), dimension(:), &
                          intent(inout) :: c
  integer*1                             :: i
  integer                               :: item, pos
  !
  do item=1,size(c,1)
    do pos=1,len(c(item))
      read(unit), i
      c(item)(pos:pos)=achar(i)
    end do
  end do
end subroutine pmd_rd_char_1d
!
subroutine pmd_rd_int_1d(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, dimension(:), intent(inout)  :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  integer                               :: item
  !
  select case (integer_size)
    case (1)
      do item=1,size(u,1)
        read(unit), i1
        u(item)=int(i1)
      end do
    case (2)
      do item=1,size(u,1)
        read(unit), i2
        u(item)=int(i2)
      end do
    case (4)
      do item=1,size(u,1)
        read(unit), i4
        u(item)=int(i4)
      end do
    case (8)
      do item=1,size(u,1)
        read(unit), i8
        u(item)=int(i8)
      end do
    case default
      do item=1,size(u,1)
        read(unit), i4
        u(item)=int(i4)
      end do
  end select
end subroutine pmd_rd_int_1d
!
subroutine pmd_rd_real_1d(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:), intent(inout)     :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: item
  !
  select case (double_size)
    case (4)
      do item=1,size(r,1)
        read(unit), f4
        r(item)=real(f4)
      end do
    case (8)
      do item=1,size(r,1)
        read(unit), f8
        r(item)=real(f8)
      end do
    case (16)
      do item=1,size(r,1)
        read(unit), f16
        r(item)=real(f16)
      end do
    case default
      do item=1,size(r,1)
        read(unit), f8
        r(item)=real(f8)
      end do
  end select
end subroutine pmd_rd_real_1d
!
subroutine pmd_rd_db_1d(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:), &
                          intent(inout) :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: item
  !
  select case (double_size)
    case (4)
      do item=1,size(d,1)
        read(unit), f4
        d(item)=dble(f4)
      end do
    case (8)
      do item=1,size(d,1)
        read(unit), f8
        d(item)=dble(f8)
      end do
    case (16)
      do item=1,size(d,1)
        read(unit), f16
        d(item)=dble(f16)
      end do
    case default
      do item=1,size(d,1)
        read(unit), f8
        d(item)=dble(f8)
      end do
  end select
end subroutine pmd_rd_db_1d
!
subroutine pmd_rd_int8_2d(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, dimension(:,:), &
                            intent(out) :: i
  !
  read(unit), i
end subroutine pmd_rd_int8_2d
!
subroutine pmd_rd_char_2d(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), dimension(:,:), &
                          intent(inout) :: c
  integer*1                             :: i
  integer                               :: j1, j2, pos
  !
  do j2=1,size(c,2)
    do j1=1,size(c,1)
      do pos=1,len(c(j1,j2))
        read(unit), i
        c(j1,j2)(pos:pos)=achar(i)
      end do
    end do
  end do
end subroutine pmd_rd_char_2d
!
subroutine pmd_rd_int_2d(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, dimension(:,:), &
                          intent(inout) :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  integer                               :: j1, j2
  !
  select case (integer_size)
    case (1)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          read(unit), i1
          u(j1,j2)=int(i1)
        end do
      end do
    case (2)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          read(unit), i2
          u(j1,j2)=int(i2)
        end do
      end do
    case (4)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          read(unit), i4
          u(j1,j2)=int(i4)
        end do
      end do
    case (8)
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          read(unit), i8
          u(j1,j2)=int(i8)
        end do
      end do
    case default
      do j2=1,size(u,2)
        do j1=1,size(u,1)
          read(unit), i4
          u(j1,j2)=int(i4)
        end do
      end do
  end select
end subroutine pmd_rd_int_2d
!
subroutine pmd_rd_real_2d(unit, r)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:,:), intent(inout)   :: r
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2
  !
  select case (double_size)
    case (4)
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          read(unit), f4
          r(j1,j2)=real(f4)
        end do
      end do
    case (8)
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          read(unit), f8
          r(j1,j2)=real(f8)
        end do
      end do
    case (16)
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          read(unit), f16
          r(j1,j2)=real(f16)
        end do
      end do
    case default
      do j2=1,size(r,2)
        do j1=1,size(r,1)
          read(unit), f8
          r(j1,j2)=real(f8)
        end do
      end do
  end select
end subroutine pmd_rd_real_2d
!
subroutine pmd_rd_db_2d(unit, d)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:), &
                          intent(inout) :: d
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2
  !
  select case (double_size)
    case (4)
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          read(unit), f4
          d(j1,j2)=dble(f4)
        end do
      end do
    case (8)
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          read(unit), f8
          d(j1,j2)=dble(f8)
        end do
      end do
    case (16)
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          read(unit), f16
          d(j1,j2)=dble(f16)
        end do
      end do
    case default
      do j2=1,size(d,2)
        do j1=1,size(d,1)
          read(unit), f8
          d(j1,j2)=dble(f8)
        end do
      end do
  end select
end subroutine pmd_rd_db_2d
!
subroutine pmd_rd_int8_3d(unit, i)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, dimension(:,:,:), &
                            intent(out) :: i
  !
  read(unit), i
end subroutine pmd_rd_int8_3d
!
subroutine pmd_rd_char_3d(unit, c)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), dimension(:,:,:), &
                          intent(inout) :: c
  integer*1                             :: i
  integer                               :: j1, j2, j3, pos
  !
  do j3=1,size(c,3)
    do j2=1,size(c,2)
      do j1=1,size(c,1)
        do pos=1,len(c(j1,j2,j3))
          read(unit), i
          c(j1,j2,j3)(pos:pos)=achar(i)
        end do
      end do
    end do
  end do
end subroutine pmd_rd_char_3d
!
subroutine pmd_rd_int_3d(unit, u)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, dimension(:,:,:), &
                          intent(inout) :: u
  integer*1                             :: i1
  integer*2                             :: i2
  integer*4                             :: i4
  integer*8                             :: i8
  integer                               :: j1, j2, j3
  !
  select case (integer_size)
    case (1)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            read(unit), i1
            u(j1,j2,j3)=int(i1)
          end do
        end do
      end do
    case (2)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            read(unit), i2
            u(j1,j2,j3)=int(i2)
          end do
        end do
      end do
    case (4)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            read(unit), i4
            u(j1,j2,j3)=int(i4)
          end do
        end do
      end do
    case (8)
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            read(unit), i8
            u(j1,j2,j3)=int(i8)
          end do
        end do
      end do
    case default
      do j3=1,size(u,3)
        do j2=1,size(u,2)
          do j1=1,size(u,1)
            read(unit), i4
            u(j1,j2,j3)=int(i4)
          end do
        end do
      end do
  end select
end subroutine pmd_rd_int_3d
!
subroutine pmd_rd_real_3d(unit, r, iostat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:,:,:), intent(inout) :: r
  integer, optional, intent(out)        :: iostat
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2, j3, iostat0
  !
  select case (double_size)
    case (4)
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            read(unit, iostat=iostat0), f4
            r(j1,j2,j3)=real(f4)
          end do
        end do
      end do
    case (8)
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            read(unit, iostat=iostat0), f8
            r(j1,j2,j3)=real(f8)
          end do
        end do
      end do
    case (16)
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            read(unit, iostat=iostat0), f16
            r(j1,j2,j3)=real(f16)
          end do
        end do
      end do
    case default
      do j3=1,size(r,3)
        do j2=1,size(r,2)
          do j1=1,size(r,1)
            read(unit, iostat=iostat0), f8
            r(j1,j2,j3)=real(f8)
          end do
        end do
      end do
  end select
  !
  if (present(iostat)) iostat=iostat0
end subroutine pmd_rd_real_3d
!
subroutine pmd_rd_db_3d(unit, d, iostat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:,:), &
                          intent(inout) :: d
  integer, optional, intent(out)        :: iostat
  real*4                                :: f4
  real*8                                :: f8
  real*16                               :: f16
  integer                               :: j1, j2, j3, iostat0
  !
  select case (double_size)
    case (4)
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            read(unit, iostat=iostat0), f4
            d(j1,j2,j3)=dble(f4)
          end do
        end do
      end do
    case (8)
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            read(unit, iostat=iostat0), f8
            d(j1,j2,j3)=dble(f8)
          end do
        end do
      end do
    case (16)
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            read(unit, iostat=iostat0), f16
            d(j1,j2,j3)=dble(f16)
          end do
        end do
      end do
    case default
      do j3=1,size(d,3)
        do j2=1,size(d,2)
          do j1=1,size(d,1)
            read(unit, iostat=iostat0), f8
            d(j1,j2,j3)=dble(f8)
          end do
        end do
      end do
  end select
  !
  if (present(iostat)) iostat=iostat0
end subroutine pmd_rd_db_3d
!
subroutine pmd_wr_header(unit, pmd_header)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  type(pmd_header_type), intent(in)     :: pmd_header
  integer                               :: pos
  character(len=1023)                   :: atomic_module
  character(len=4096)                   :: pmd_comment
  !
  call pmd_ftell(unit, pos)
  call pmd_fseek(unit, 0, 0)
  !
  integer_size  = pmd_header%int_sz
  double_size   = pmd_header%db_sz
  atomic_module = repeat(achar(0), 1023)
  pmd_comment   = repeat(achar(0), 4096)
  atomic_module(1:len(trim(pmd_header%atomic_module))) &
        = trim(pmd_header%atomic_module)
  pmd_comment(1:len(trim(pmd_header%pmd_comment))) &
        = trim(pmd_header%pmd_comment)
  !
  call pmd_wr(unit, pmd_header%magic_str)
  call pmd_wr(unit, pmd_header%endianness)
  call pmd_wr(unit, pmd_header%int_sz)
  call pmd_wr(unit, pmd_header%db_sz)
  call pmd_wr(unit, pmd_header%pmd_version)
  call pmd_wr(unit, pmd_header%localtime)
  call pmd_wr(unit, pmd_header%periodic)
  call pmd_wr(unit, pmd_header%domain_sz)
  call pmd_wr(unit, pmd_header%domain_origin)
  call pmd_wr(unit, pmd_header%dimensions)
  call pmd_wr(unit, pmd_header%xc1)
  call pmd_wr(unit, pmd_header%xc2)
  call pmd_wr(unit, pmd_header%xc3)
  call pmd_wr(unit, pmd_header%n_radtheta)
  call pmd_wr(unit, pmd_header%n_radphi)
  call pmd_wr(unit, atomic_module)
  call pmd_wr(unit, pmd_comment)
  call pmd_wr(unit, pmd_header%module_hd_sz)
  call pmd_wr(unit, pmd_header%unused_var)
  !
  call pmd_fseek(unit, pos, 0)
end subroutine pmd_wr_header
!
subroutine pmd_rd_header(unit, pmd_header)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  type(pmd_header_type), intent(out)    :: pmd_header
  integer                               :: pos
  !
  call pmd_ftell(unit, pos)
  call pmd_fseek(unit, 0, 0)
  !
  call pmd_rd(unit, pmd_header%magic_str)
  call pmd_rd(unit, pmd_header%endianness)
  call pmd_rd(unit, pmd_header%int_sz)
  call pmd_rd(unit, pmd_header%db_sz)
  !
  integer_size = pmd_header%int_sz
  double_size  = pmd_header%db_sz
  !
  call pmd_rd(unit, pmd_header%pmd_version)
  call pmd_rd(unit, pmd_header%localtime)
  call pmd_rd(unit, pmd_header%periodic)
  call pmd_rd(unit, pmd_header%domain_sz)
  call pmd_rd(unit, pmd_header%domain_origin)
  call pmd_rd(unit, pmd_header%dimensions)
  call pmd_rd(unit, pmd_header%xc1)
  call pmd_rd(unit, pmd_header%xc2)
  call pmd_rd(unit, pmd_header%xc3)
  call pmd_rd(unit, pmd_header%n_radtheta)
  call pmd_rd(unit, pmd_header%n_radphi)
  call pmd_rd(unit, pmd_header%atomic_module)
  call pmd_rd(unit, pmd_header%pmd_comment)
  call pmd_rd(unit, pmd_header%module_hd_sz)
  call pmd_rd(unit, pmd_header%unused_var)
  !
  call pmd_fseek(unit, pos, 0)
end subroutine pmd_rd_header
!
subroutine pmd_types_sz(unit, int_sz, db_sz, stat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer*1, intent(out)                :: int_sz, db_sz
  integer, optional, intent(out)        :: stat
  integer                               :: pos, stat0=0
  !
  call pmd_ftell(unit, pos)
  call pmd_fseek(unit, 9, 0) !, status=stat0)
  !
  if (stat0 == 0) read(unit, iostat=stat0), int_sz
  if (stat0 == 0) read(unit, iostat=stat0), db_sz
  !
  call pmd_fseek(unit, pos, 0)
  !
  if (present(stat)) stat=stat0
end subroutine pmd_types_sz
!
function pmd_header_sz_unit(unit)
  !
  implicit none
  !
  integer                               :: pmd_header_sz_unit
  integer, intent(in)                   :: unit
  integer*1                             :: int_sz, db_sz
  integer                               :: stat0
  !
  call pmd_types_sz(unit, int_sz, db_sz, stat=stat0)
  if (stat0 == 0) then
    pmd_header_sz_unit = 5132+14*int(int_sz)+24582*int(db_sz)
  else
    pmd_header_sz_unit = 0
  endif
end function pmd_header_sz_unit
!
function pmd_header_sz_hd(pmd_header)
  !
  implicit none
  !
  integer                               :: pmd_header_sz_hd
  type(pmd_header_type), intent(in)     :: pmd_header
  !
  pmd_header_sz_hd = 5132+14*pmd_header%int_sz+24582*pmd_header%db_sz
end function pmd_header_sz_hd
!
function pmd_box_sz(pmd_header)
  !
  implicit none
  !
  integer                               :: pmd_box_sz
  type(pmd_header_type), intent(in)     :: pmd_header
  !
  pmd_box_sz = pmd_header%dimensions(1) * pmd_header%dimensions(2) &
               * pmd_header%dimensions(3) * pmd_header%db_sz
end function pmd_box_sz
!
function pmd_node_sz(unit, pmd_header)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer                               :: pmd_node_sz
  type(pmd_header_type), intent(in)     :: pmd_header
  integer                               :: sz, header_sz, mod_hd_sz, &
                                           ncells
  !
  inquire(unit, size=sz)
  header_sz = pmd_header_sz(pmd_header)
  mod_hd_sz = pmd_header%module_hd_sz
  ncells = pmd_header%dimensions(1) * pmd_header%dimensions(2) &
           * pmd_header%dimensions(3)
  pmd_node_sz = int((sz-header_sz-mod_hd_sz)/ncells)
end function pmd_node_sz
!
subroutine pmd_seek(unit, pmd_header, pos)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  type(pmd_header_type), intent(in)     :: pmd_header
  integer, intent(in)                   :: pos
  integer                               :: hd_sz, tobox0, box_sz
  !
  hd_sz = pmd_header_sz(pmd_header)
  tobox0 = hd_sz+pmd_header%module_hd_sz
  box_sz = pmd_box_sz(pmd_header)
  if (pos == 0) call pmd_fseek(unit, hd_sz, 0)
  if (pos > 0) call pmd_fseek(unit, tobox0+(pos-1)*box_sz, 0)
  if (pos < 0) call pmd_fseek(unit, -pos*box_sz, 2)
end subroutine pmd_seek
!
subroutine pmd_seek_node(unit, pmd_header, nodeX, nodeY, nodeZ)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  type(pmd_header_type), intent(in)     :: pmd_header
  integer, intent(in)                   :: nodeX, nodeY, nodeZ
  integer                               :: hd_sz, tonode0, node_sz, &
                                           dimX, dimY, dimZ, node
  !
  hd_sz   = pmd_header_sz(pmd_header)
  tonode0 = hd_sz+pmd_header%module_hd_sz
  node_sz = pmd_node_sz(unit, pmd_header)
  dimX    = pmd_header%dimensions(0)
  dimY    = pmd_header%dimensions(1)
  dimZ    = pmd_header%dimensions(2)
  node    = nodeX + nodeY * dimX + nodeZ * dimX*dimY
  if (nodeX>=0.and.nodeX<dimX.and.nodeY>=0.and.nodeY<dimY &
      .and.nodeZ>=0.and.nodeZ<dimZ) then
    call pmd_fseek(unit, tonode0+node*node_sz, 0)
  endif
end subroutine pmd_seek_node
!
subroutine pmd_append_box_real_3d(unit, box)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:,:,:), intent(in)    :: box
  integer                               :: pos
  !
  call pmd_types_sz(unit, integer_size, double_size)
  call pmd_ftell(unit, pos)
  call pmd_fseek(unit, 0, 2)
  call pmd_wr(unit, box)
  call pmd_fseek(unit, pos, 0)
end subroutine pmd_append_box_real_3d
!
subroutine pmd_append_box_db_3d(unit, box)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:,:), &
                             intent(in) :: box
  integer                               :: pos
  !
  call pmd_types_sz(unit, integer_size, double_size)
  call pmd_ftell(unit, pos)
  call pmd_fseek(unit, 0, 2)
  call pmd_wr(unit, box)
  call pmd_fseek(unit, pos, 0)
end subroutine pmd_append_box_db_3d
!
subroutine pmd_rd_box_real_3d(unit, box, iostat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  real, dimension(:,:,:), intent(inout) :: box
  integer, optional, intent(out)        :: iostat
  integer                               :: iostat0
  !
  call pmd_types_sz(unit, integer_size, double_size)
  call pmd_rd(unit, box, iostat=iostat0)
  !
  if (present(iostat)) iostat=iostat0
  !
end subroutine pmd_rd_box_real_3d
!
subroutine pmd_rd_box_db_3d(unit, box, iostat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:,:), &
                          intent(inout) :: box
  integer, optional, intent(out)        :: iostat
  integer                               :: iostat0
  !
  call pmd_types_sz(unit, integer_size, double_size)
  call pmd_rd(unit, box, iostat=iostat0)
  !
  if (present(iostat)) iostat=iostat0
  !
end subroutine pmd_rd_box_db_3d
!
subroutine pmd_append_node(unit, node)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character, dimension(:), intent(in)   :: node
  integer                               :: pos
  !
  call pmd_types_sz(unit, integer_size, double_size)
  call pmd_ftell(unit, pos)
  call pmd_fseek(unit, 0, 2)
  call pmd_wr(unit, node)
  call pmd_fseek(unit, pos, 0)
end subroutine pmd_append_node
!
subroutine pmd_rd_node(unit, node)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character, dimension(:), intent(out)  :: node
  !
  call pmd_types_sz(unit, integer_size, double_size)
  call pmd_rd(unit, node)
  !
end subroutine pmd_rd_node
!
subroutine pmd_openwr(unit, file, pmd_header, convert, stat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), intent(in)              :: file
  type(pmd_header_type), optional, &
                            intent(in)  :: pmd_header
  character(*), optional, intent(in)    :: convert
  integer, optional, intent(out)        :: stat
  character(80)                         :: convert0
  integer                               :: iostat0, i
  !
  if (present(convert)) then
    convert0 = convert
  else
    if (present(pmd_header)) then
      if (pmd_header%endianness == 0) then
        convert0 = 'little_endian'
      elseif (pmd_header%endianness == 1) then
        convert0 = 'big_endian'
      else
        convert0 = 'native'
      endif
    else
      convert0 = 'native'
    endif
  endif
  !
  open(unit, file=file, action='readwrite', access='stream', &
       status='replace', convert=trim(convert0), iostat=iostat0)
  !
  if (iostat0 == 0) then
  if (present(pmd_header)) then
    call pmd_wr_header(unit, pmd_header)
    call pmd_fseek(unit, 0, 2)
    do i=1,pmd_header%module_hd_sz
      call pmd_wr(unit, char(0))
    end do
  endif
  endif ! Open
  !
  if (present(stat)) stat = iostat0
  !
end subroutine pmd_openwr
!
subroutine pmd_openap(unit, file, pmd_header, stat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), intent(in)              :: file
  type(pmd_header_type), optional, &
                            intent(out) :: pmd_header
  integer, optional, intent(out)        :: stat
  type(pmd_header_type)                 :: pmd_header0
  integer                               :: iostat0, stat0, sz, hd_sz, data_sz
  integer*1                             :: endianness
  character(len=8)                      :: magic_str
  !
  open(unit, file=file, action='read', access='stream', &
       status='old', convert='native', iostat=iostat0)
  !
  if (iostat0 > 0) then
    stat0 = -iostat0
  else
    stat0 = iostat0
  endif
  !
  if (stat0 == 0) then
  inquire(unit, size=sz)
  if (sz < 8) then
    if (present(stat)) stat = -500
    return
  else if (sz < 9) then
    if (present(stat)) stat = -501
    return
  endif
  call pmd_fseek(unit, 8, 0)
  call pmd_rd(unit, endianness)
  close(unit)
  if (endianness == 0) then
    open(unit, file=file, action='readwrite', access='stream', &
         status='old', convert='little_endian', iostat=iostat0)
  else
    open(unit, file=file, action='readwrite', access='stream', &
         status='old', convert='big_endian', iostat=iostat0)
  endif
  !
  if (iostat0 > 0) then
    stat0 = -iostat0
  else
    stat0 = iostat0
  endif
  !
  if (stat0 == 0) then
  inquire(unit, size=sz)
  if(sz >= 8) then
    call pmd_fseek(unit, 0, 0)
    read(unit), magic_str
    if (magic_str /= pmd_magic_str) then
      if (present(stat)) stat = -500
      return
    endif
  else
    if (present(stat)) stat = -500
    return
  endif
  hd_sz=pmd_header_sz(unit)
  if (sz<hd_sz) then
    if (present(stat)) stat = -501
    return
  endif
  call pmd_rd_header(unit, pmd_header0)
  if (sz<hd_sz+pmd_header0%module_hd_sz) then
    stat0 = -502
  endif
  data_sz=sz-hd_sz-pmd_header0%module_hd_sz
  if (pmd_box_sz(pmd_header0) > 0 .and. stat0 == 0) then
    if (mod(data_sz,pmd_box_sz(pmd_header0)) /= 0) then
      stat0 = -503
    endif
  endif
  if (pmd_box_sz(pmd_header0) > 0 .and. stat0 == 0) then
    !sz = int(data_sz / pmd_box_sz(pmd_header0))
    sz = pmd_node_sz(unit, pmd_header)
  else
    sz = 0
  endif
  call pmd_fseek(unit, 0, 2)
  !
  endif ! Open 2nd time
  endif ! Open 1st time
  if (stat0 == 0) stat0 = sz
  !
  if (present(stat)) stat = stat0
  if (present(pmd_header)) pmd_header = pmd_header0
  !
end subroutine pmd_openap
!
subroutine pmd_openrd(unit, file, pmd_header, stat)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character(*), intent(in)              :: file
  type(pmd_header_type), optional, &
                            intent(out) :: pmd_header
  integer, optional, intent(out)        :: stat
  type(pmd_header_type)                 :: pmd_header0
  integer                               :: iostat0, stat0, sz, hd_sz, data_sz
  integer*1                             :: endianness
  character(len=8)                      :: magic_str
  !
  open(unit, file=file, action='read', access='stream', &
       status='old', convert='native', iostat=iostat0)
  !
  if (iostat0 > 0) then
    stat0 = -iostat0
  else
    stat0 = iostat0
  endif
  !
  if (stat0 == 0) then
  inquire(unit, size=sz)
  if (sz < 8) then
    if (present(stat)) stat = -500
    return
  else if (sz < 9) then
    if (present(stat)) stat = -501
    return
  endif
  call pmd_fseek(unit, 8, 0)
  call pmd_rd(unit, endianness)
  close(unit)
  if (endianness == 0) then
    open(unit, file=file, action='read', access='stream', &
         status='old', convert='little_endian', iostat=iostat0)
  else
    open(unit, file=file, action='read', access='stream', &
         status='old', convert='big_endian', iostat=iostat0)
  endif
  !
  if (iostat0 > 0) then
    stat0 = -iostat0
  else
    stat0 = iostat0
  endif
  !
  if (stat0 == 0) then
  inquire(unit, size=sz)
  if (sz >= 8) then
    call pmd_fseek(unit, 0, 0)
    read(unit), magic_str
    if (magic_str /= pmd_magic_str) then
      if (present(stat)) stat = -500
      return
    endif
  else
    if (present(stat)) stat = -500
    return
  endif
  hd_sz=pmd_header_sz(unit)
  if (sz<hd_sz) then
    if (present(stat)) stat = -501
    return
  endif
  call pmd_rd_header(unit, pmd_header0)
  if (sz<hd_sz+pmd_header0%module_hd_sz) then
    stat0 = -502
  endif
  data_sz=sz-hd_sz-pmd_header0%module_hd_sz
  if (pmd_box_sz(pmd_header0) > 0 .and. stat0 == 0) then
    if (mod(data_sz,pmd_box_sz(pmd_header0)) /= 0) then
      stat0 = -503
    endif
  endif
  if (pmd_box_sz(pmd_header0) > 0 .and. stat0 == 0) then
    !sz = int(data_sz / pmd_box_sz(pmd_header0))
    sz = pmd_node_sz(unit, pmd_header0)
  else
    sz = 0
  endif
  !
  if (present(pmd_header)) then
    if (sz>0) then
      call pmd_seek(unit, pmd_header0, 1)
    else
      call pmd_fseek(unit, 0, 0)
    endif
    pmd_header = pmd_header0
  else
    call pmd_fseek(unit, 0, 0)
  endif
  endif ! Open 2nd time
  endif ! Open 1st time
  if (stat0 == 0) stat0 = sz
  !
  if (present(stat)) stat = stat0
  !
end subroutine pmd_openrd
!
subroutine pmd_close(unit)
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  !
  close(unit)
end subroutine pmd_close
!
subroutine pmd_set_localtime(pmd_header)
  !
  implicit none
  !
  type(pmd_header_type), intent(inout)  :: pmd_header
  integer, dimension(9)                 :: time_values
  !
  call date_and_time(values=time_values)
  pmd_header%localtime(1)       = time_values(1)-1900
  pmd_header%localtime(2)       = time_values(2)-1
  pmd_header%localtime(3)       = time_values(3)
  pmd_header%localtime(4)       = time_values(5)
  pmd_header%localtime(5)       = time_values(6)
  pmd_header%localtime(6)       = time_values(7)
end subroutine pmd_set_localtime
!
end module pmd_io_module
!
!*******************************************************************************
!
module pypmd
!
implicit none
!
contains
!
! From http://fortranwiki.org/fortran/show/newunit
!
! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.
integer function newunit()
! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=1000
  logical :: opened
  integer :: lun
! begin
  newunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      newunit=lun
      exit
    end if
  end do
end function newunit
!
subroutine openrd(unit, file, magic_str, endianness, int_sz, db_sz, &
                  pmd_version, localtime, periodic, domain_sz, domain_origin, &
                  dimensions, xc1, xc2, xc3, n_radtheta, n_radphi, &
                  atomic_module, pmd_comment, module_hd_sz, unused_var, stat)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(out)                  :: unit
  character(*), intent(in)              :: file
  character(len=8), intent(out)         :: magic_str
  integer*1, intent(out)                :: endianness
  integer*1, intent(out)                :: int_sz, db_sz
  integer, intent(out)                  :: pmd_version
  integer, dimension(6), intent(out)    :: localtime
  integer*1, dimension(2), intent(out)  :: periodic
  double precision, dimension(3), &
                      intent(out)       :: domain_sz, domain_origin
  integer, dimension(3), intent(out)    :: dimensions
  double precision, dimension(8192), &
                         intent(out)    :: xc1, xc2, xc3
  integer, intent(out)                  :: n_radtheta, n_radphi
  character(len=1023), intent(out)      :: atomic_module
  character(len=4096), intent(out)      :: pmd_comment
  integer, intent(out)                  :: module_hd_sz
  integer, intent(out)                  :: unused_var
  integer, intent(out)                  :: stat
  type(pmd_header_type)                 :: pmd_header
  !
  unit = newunit()
  call pmd_openrd(unit, file, pmd_header, stat=stat)
  !
  magic_str     = pmd_header%magic_str
  endianness    = pmd_header%endianness
  int_sz        = pmd_header%int_sz
  db_sz         = pmd_header%db_sz
  pmd_version   = pmd_header%pmd_version
  localtime     = pmd_header%localtime
  periodic      = pmd_header%periodic
  domain_sz     = pmd_header%domain_sz
  domain_origin = pmd_header%domain_origin
  dimensions    = pmd_header%dimensions
  xc1           = pmd_header%xc1
  xc2           = pmd_header%xc2
  xc3           = pmd_header%xc3
  n_radtheta    = pmd_header%n_radtheta
  n_radphi      = pmd_header%n_radphi
  atomic_module = pmd_header%atomic_module
  pmd_comment   = pmd_header%pmd_comment
  module_hd_sz  = pmd_header%module_hd_sz
  unused_var    = pmd_header%unused_var
end subroutine openrd
!
subroutine openwr(unit, ltime, file, magic_str, endianness, int_sz, db_sz, &
                  pmd_version, localtime, periodic, domain_sz, domain_origin, &
                  dimensions, xc1, xc2, xc3, n_radtheta, n_radphi, &
                  atomic_module, pmd_comment, module_hd_sz, unused_var, stat)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(out)                  :: unit
  integer, dimension(6), intent(out)    :: ltime
  character(*), intent(in)              :: file
  character(len=8), intent(in)          :: magic_str
  integer*1, intent(in)                 :: endianness
  integer*1, intent(in)                 :: int_sz, db_sz
  integer, intent(in)                   :: pmd_version
  integer, dimension(6), intent(in)     :: localtime
  integer*1, dimension(2), intent(in)   :: periodic
  double precision, dimension(3), &
                      intent(in)        :: domain_sz, domain_origin
  integer, dimension(3), intent(in)     :: dimensions
  double precision, dimension(8192), &
                         intent(in)     :: xc1, xc2, xc3
  integer, intent(in)                   :: n_radtheta, n_radphi
  character(len=1023), intent(in)       :: atomic_module
  character(len=4096), intent(in)       :: pmd_comment
  integer, intent(in)                   :: module_hd_sz
  integer, intent(in)                   :: unused_var
  integer, intent(out)                  :: stat
  type(pmd_header_type)                 :: pmd_header
  !
  pmd_header%magic_str          = magic_str
  pmd_header%endianness         = endianness
  pmd_header%int_sz             = int_sz
  pmd_header%db_sz              = db_sz
  pmd_header%pmd_version        = pmd_version
  pmd_header%localtime          = localtime
  pmd_header%periodic           = periodic
  pmd_header%domain_sz          = domain_sz
  pmd_header%domain_origin      = domain_origin
  pmd_header%dimensions         = dimensions
  pmd_header%xc1                = xc1
  pmd_header%xc2                = xc2
  pmd_header%xc3                = xc3
  pmd_header%n_radtheta         = n_radtheta
  pmd_header%n_radphi           = n_radphi
  pmd_header%atomic_module      = atomic_module
  pmd_header%pmd_comment        = pmd_comment
  pmd_header%module_hd_sz       = module_hd_sz
  pmd_header%unused_var         = unused_var
  !
  if(pmd_header%localtime(1)<0 &
     .or.pmd_header%localtime(2)<1.or.pmd_header%localtime(2)>31 &
     .or.pmd_header%localtime(3)<0.or.pmd_header%localtime(3)>23 &
     .or.pmd_header%localtime(4)<0.or.pmd_header%localtime(4)>59 &
     .or.pmd_header%localtime(5)<0.or.pmd_header%localtime(5)>61) &
     call pmd_set_localtime(pmd_header)
  ltime = pmd_header%localtime
  !
  unit = newunit()
  call pmd_openwr(unit, file, pmd_header, stat=stat)
end subroutine openwr
!
subroutine openap(unit, file, magic_str, endianness, int_sz, db_sz, &
                  pmd_version, localtime, periodic, domain_sz, domain_origin, &
                  dimensions, xc1, xc2, xc3, n_radtheta, n_radphi, &
                  atomic_module, pmd_comment, module_hd_sz, unused_var, stat)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(out)                  :: unit
  character(*), intent(in)              :: file
  character(len=8), intent(out)         :: magic_str
  integer*1, intent(out)                :: endianness
  integer*1, intent(out)                :: int_sz, db_sz
  integer, intent(out)                  :: pmd_version
  integer, dimension(6), intent(out)    :: localtime
  integer*1, dimension(2), intent(out)  :: periodic
  double precision, dimension(3), &
                            intent(out) :: domain_sz, domain_origin
  integer, dimension(3), intent(out)    :: dimensions
  double precision, dimension(8192), &
                         intent(out)    :: xc1, xc2, xc3
  integer, intent(out)                  :: n_radtheta, n_radphi
  character(len=1023), intent(out)      :: atomic_module
  character(len=4096), intent(out)      :: pmd_comment
  integer, intent(out)                  :: module_hd_sz
  integer, intent(out)                  :: unused_var
  integer, intent(out)                  :: stat
  type(pmd_header_type)                 :: pmd_header
  !
  unit = newunit()
  call pmd_openap(unit, file, pmd_header, stat=stat)
  !
  magic_str     = pmd_header%magic_str
  endianness    = pmd_header%endianness
  int_sz        = pmd_header%int_sz
  db_sz         = pmd_header%db_sz
  pmd_version   = pmd_header%pmd_version
  localtime     = pmd_header%localtime
  periodic      = pmd_header%periodic
  domain_sz     = pmd_header%domain_sz
  domain_origin = pmd_header%domain_origin
  dimensions    = pmd_header%dimensions
  xc1           = pmd_header%xc1
  xc2           = pmd_header%xc2
  xc3           = pmd_header%xc3
  n_radtheta    = pmd_header%n_radtheta
  n_radphi      = pmd_header%n_radphi
  atomic_module = pmd_header%atomic_module
  pmd_comment   = pmd_header%pmd_comment
  module_hd_sz  = pmd_header%module_hd_sz
  unused_var    = pmd_header%unused_var
end subroutine openap
!
subroutine rd_box(unit, box, iostat)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:,:), &
                          intent(inout) :: box
  integer, intent(out)                  :: iostat
  !
  call pmd_rd_box(unit, box, iostat=iostat)
end subroutine rd_box
!
subroutine rd_nbox(unit, box, iostat, nx, ny, nz)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(in)                   :: nx, ny, nz
  double precision, dimension(nx,ny,nz), &
                          intent(inout) :: box
  integer, intent(out)                  :: iostat
  !
  call pmd_rd_box(unit, box, iostat=iostat)
end subroutine rd_nbox
!
subroutine append_box(unit, box)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  double precision, dimension(:,:,:), &
                             intent(in) :: box
  !
  call pmd_append_box(unit, box)
end subroutine append_box
!
subroutine append_nbox(unit, box, nx, ny, nz)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(in)                   :: nx, ny, nz
  double precision, dimension(nx,ny,nz), &
                             intent(in) :: box
  !
  call pmd_append_box(unit, box)
end subroutine append_nbox
!
subroutine rd_node(unit, node, sz)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(in)                   :: sz
  character, dimension(sz), intent(out) :: node
  !
  call pmd_rd_node(unit, node)
end subroutine rd_node
!
subroutine append_node(unit, node)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character, dimension(:), intent(in)   :: node
  !
  call pmd_append_node(unit, node)
end subroutine append_node
!
subroutine rd_module_hd(unit, module_header, sz)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  integer, intent(in)                   :: sz
  character, dimension(sz), intent(out) :: module_header
  integer                               :: hd_sz, pos
  !
  if (sz <= 0) return
  call pmd_ftell(unit, pos)
  hd_sz = pmd_header_sz(unit)
  call pmd_fseek(unit, hd_sz, 0)
  call pmd_rd(unit, module_header)
  call pmd_fseek(unit, pos, 0) 
end subroutine rd_module_hd
!
subroutine wr_module_hd(unit, module_header)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  character, dimension(:), intent(in)   :: module_header
  integer                               :: hd_sz, pos
  !
  call pmd_ftell(unit, pos)
  hd_sz = pmd_header_sz(unit)
  call pmd_fseek(unit, hd_sz, 0)
  call pmd_wr(unit, module_header)
  call pmd_fseek(unit, pos, 0) 
end subroutine wr_module_hd
!
subroutine close(unit)
  !
  use pmd_io_module
  !
  implicit none
  !
  integer, intent(in)                   :: unit
  !
  call pmd_close(unit)
end subroutine close
!
end module pypmd
!
!*******************************************************************************
!
!program test
!  !
!  use pmd_io_module
!  !
!  implicit none
!  !
!  type(pmd_header_type)                 :: header, header_check
!  integer                               :: stat
!  integer, dimension(2,3)               :: arr
!  character, dimension(12)              :: letters
!  character*2, dimension(1,2,3)         :: words
!  integer                               :: sz, i
!  double precision                      :: db
!  integer                               :: pos
!  !
!  header%periodic       = (/1,1/)
!  header%domain_sz      = (/1.,1.,1./)
!  header%dimensions     = (/2,3,4/)
!  header%n_radtheta     = 2
!  header%n_radphi       = 2
!  header%atomic_module  = 'The atomic module'
!  header%pmd_comment    = 'No comment'
!  header%module_hd_sz   = 4
!  !
!  call pmd_set_localtime(header)
!  call pmd_openwr(1, 'test.bin', header, stat=stat)
!  !call pmd_openwr(1, 'test.bin', stat=stat)
!  print*, 'Openwr test.bin', stat
!  !write(1), 'portapmd'
!  !inquire(1, size=sz)
!  !print*, sz
!  print*, 'box_sz', pmd_box_sz(header)
!  call pmd_seek(1, header, 0)
!  call pmd_wr(1, '1234')
!  do i=1,24
!    db = i
!    call pmd_wr(1, db)
!  end do
!  inquire(1, size=sz)
!  print*, sz
!  close(1)
!  call pmd_openrd(1, 'test.bin', header_check, stat=stat)
!  call pmd_ftell(1, pos)
!  !inquire(1,pos=pos)
!  print*, 't', pos
!  call pmd_fseek(1,1,1)
!  !read(1,pos=pos+2)
!  call pmd_ftell(1,pos)
!  !inquire(1,pos=pos)
!  print*, 't', pos
!  call pmd_fseek(1,pos,0)
!  !read(1,pos=pos+2)
!  call pmd_ftell(1,pos)
!  !inquire(1,pos=pos)
!  print*, 't', pos
!  print*, 'Openrd test.bin', stat
!  inquire(1, size=sz)
!  print*, sz
!  close(1)
!  call pmd_openwr(1, 'test_check.bin', header_check, stat=stat)
!  print*, 'Openwr test_check.bin', stat
!  close(1)
!  call pmd_openap(1, 'test.bin', stat=stat)
!  print*, 'Openap test.bin', stat
!  close(1)
!  arr=reshape((/1,2,3,4,5,6/),shape(arr))
!  letters=(/'a','b','c','d','e','f','g','h','i','j','k','l'/)
!  words=reshape((/'ab','cd','ef','gh','ij','kz'/),shape(words))
!  call pmd_wr(1, arr)
!  call pmd_wr(1, letters)
!  call pmd_wr(1, words)
!  !write(1), letters, words, arr
!  print*, header%int_sz, header%db_sz, pmd_header_sz(header), pmd_box_sz(header)
!end program test
