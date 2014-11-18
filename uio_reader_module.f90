!*******************************************************************************
!
!  _    _ _____ ____                            _           
! | |  | |_   _/ __ \                          | |          
! | |  | | | || |  | |       _ __ ___  __ _  __| | ___ _ __ 
! | |  | | | || |  | |      | '__/ _ \/ _` |/ _` |/ _ \ '__|
! | |__| |_| || |__| |      | | |  __/ (_| | (_| |  __/ |   
!  \____/|_____\____/ _____ |_|  \___|\__,_|\__,_|\___|_|   _module
!                    |_____|
!
! Reader module for arbitrary UIO files
!
!*******************************************************************************
!   Fortran 90
!   Flavio Calvo:               Geneva, Locarno
!   2014-08-26
!*******************************************************************************
!
!------*************------------------------------------------------------------
module integer_ll_module
include "linkedlistofint.f90"
end module integer_ll_module
!-------------------------------------------------------------------------------
module uio_reader_module
!-------------------------------------------------------------------------------
! NAME:
!   uio_reader_module ('universal_input_output_reader_module')
!
! PURPOSE:
!   Read UIO file
!
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   use uio_reader_module
!
! VARIABLES (public):
!   record:     pointer to callback function
!   inte*, real*, comp*, char*, string1D, string2D, string3D, string4D
!   (* stands for '0', '1D', '2D', '3D' and '4D' respectively)
!   tab:        uio_table
!   name0:      name of variable
!   unit0:      unit of variable
!   label:      corresponds to 'identifier' in other uio_modules
!   qtype:      type of variable (interger, real, character,...)
!   use_dq:     logical, whether or not derived quantities should be computed
!   r_modelfile:model filename
!   r_parfile:  parfile filename
!   r_nc:       unit in which modelfile is opened
!   r_nmodel:   number of snapshots in modelfile
!   r_imodel:   current snapshot
!
! TYPES:
!   None
!
! SUBROUTINES (public):
!   uio_dq_wrapper:    identifies known quantities and calls rhd_dq_init
!                      and uio_rd accordingly (to use only when use_dq=.true.)
!   uio_data_rd:       identifies data type and reads it, unless uio_dq_wrapper
!                      recognizes the quantity and reads it first
!   uio_label_count:   counts number of some kind of labels in a box
!   uio_struct_rd:     reads UIO files calling uio_data_rd when data is found
!
! HOW IT WORKS:
!   1. uio_struct_rd starts reading some general UIO file. uio_struct_read is
!      the Fortran equivalent of the corresponding uio_struct_read IDL function.
!   2. When data is found, uio_data_rd is called. When use_dq=.true.,
!      uio_data_rd calls uio_dq_wrapper
!   3. uio_dq_wrapper is aware of the quantities that are expected in an
!      STA/END/FULL file. The rhd_dq_module is called to allocate memory for
!      those quantities and uio_rd is called to actually read them. If the
!      quantity is not recognized, "processed" variable is set to "False"
!   4. If use_dq was set to False or uio_dq_wrapper did not recognize some
!      quantity, uio_data_rd will read the quantity on his own, using uio_rd
!
! HOW IS DATA PASSED TO A FORTRAN PROGRAM OR TO PYTHON?
!   There are two possible situations:
!     * use_dq is set to True and data is recognized by uio_dq_wrapper
!       In this case the data will be potentially needed by the rhd_dq_module
!       to compute derived quantities. Memory has been allocated by
!       rhd_dq_module and rhd_dq_module owns this data. Then the callback
!       functions is called with the argument 'link', and it is the
!       responsability of the Fortran/Python program to access the data
!       directly in the rhd_dq_module
!     * use_dq is set to False or data is not recognized by uio_dq_wrapper
!       In this case data is stored in a temporary variable inte*, real*,
!       comp*, char*, string1D, string2D, string3D or string4D depending on
!       its type. The callback function is called with 'copy' argument and
!       it is the responsability of the Fortran/Python program to get the data
!       and make a copy of it (the original data will be immediately destroyed
!       after the callback function is called).
!       For tables it is a bit more complicated (the callback function is
!       called with 'table' argument a first time and it is repeatedly called
!       with 'table:...' argument for each column of the table).
!
!  MODIFICATION HISTORY:
!    2014-08-26 (F. C. Locarno) Written
!
!-------------------------------------------------------------------------------
!
use uio_table_module
use rhd_dq_module
use integer_ll_module
!
implicit none
!
procedure(), pointer                    :: record
!
integer                                 :: inte0
integer, allocatable, dimension(:)      :: inte1D
integer, allocatable, dimension(:,:)    :: inte2D
integer, allocatable, dimension(:,:,:)  :: inte3D
integer, allocatable, dimension(:,:,:,:)        :: inte4D
real                                    :: real0
real, allocatable, dimension(:)         :: real1D
real, allocatable, dimension(:,:)       :: real2D
real, allocatable, dimension(:,:,:)     :: real3D
real, allocatable, dimension(:,:,:,:)   :: real4D
complex                                 :: comp0
complex, allocatable, dimension(:)      :: comp1D
complex, allocatable, dimension(:,:)    :: comp2D
complex, allocatable, dimension(:,:,:)  :: comp3D
complex, allocatable, dimension(:,:,:,:)        :: comp4D
character*(let)                                 :: char0
character, allocatable, dimension(:,:)          :: char1D
character, allocatable, dimension(:,:,:)        :: char2D
character, allocatable, dimension(:,:,:,:)      :: char3D
character, allocatable, dimension(:,:,:,:,:)    :: char4D
character*(let), allocatable, dimension(:)      :: string1D
character*(let), allocatable, dimension(:,:)    :: string2D
character*(let), allocatable, dimension(:,:,:)  :: string3D
character*(let), allocatable, dimension(:,:,:,:)        :: string4D
type(uio_table_type)                    :: tab
character*(let)                         :: name0, unit0
character*(let)                         :: label, qtype
integer                                 :: ndim
integer                                 :: m1, n1, m2, n2, m3, n3, m4, n4
logical                                 :: use_dq
character*(let)                         :: r_modelfile, r_parfile
integer                                 :: r_nc, r_nmodel, r_imodel
type(integers_ll)                       :: uio_index
!
contains
!
subroutine uio_dq_wrapper(nc, termt, ntt, ident, name, unit, processed, &
                          norec, outstr, ierr)
use uio_bulk_module
!use uio_base_module
!use uio_siz_module
use uio_nam_module
use rhd_dq_module
use rhd_box_module
implicit none
!
! --- I/O parameters ---
integer,                  intent(in)    :: nc
logical,       optional,  intent(in)    :: norec
character*(*),         intent(inout)    :: termt(2,nttmx)
integer,               intent(inout)    :: ntt
character*(*),           intent(out)    :: ident
character*(*), optional, intent(out)    :: name, unit
logical,                  intent(out)   :: processed
character*(*), optional, intent(out)    :: outstr
integer,       optional, intent(out)    :: ierr
!
! --- Local variables
!
logical                                 :: norec0
character*(let)                         :: outstr0
integer                                 :: ierr0
!
outstr0='done'
ierr0=0
!
! --- Additional keywords ---
if (present(norec)) then
        norec0=norec
else
        norec0=.False.
endif
!
processed=.True.
if (qtype == realna .and. ndim == 3) then
        !print*, 'Inside wrapper'
        select case (trim(label))
          case ('rho')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'rho')
                  call uio_rd(nc, termt, ntt, dq_rho, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('ei')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'ei')
                  call uio_rd(nc, termt, ntt, dq_ei, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('v1')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'v1')
                  call uio_rd(nc, termt, ntt, dq_v1, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('v2')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'v2')
                  call uio_rd(nc, termt, ntt, dq_v2, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('v3')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'v3')
                  call uio_rd(nc, termt, ntt, dq_v3, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('bb1')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'Bb1')
                  call uio_rd(nc, termt, ntt, dq_Bb1, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('bb2')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'Bb2')
                  call uio_rd(nc, termt, ntt, dq_Bb2, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('bb3')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'Bb3')
                  call uio_rd(nc, termt, ntt, dq_Bb3, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('xb1')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'xb1')
                  call uio_rd(nc, termt, ntt, dq_xb1, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('xb2')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'xb2')
                  call uio_rd(nc, termt, ntt, dq_xb2, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('xb3')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'xb3')
                  call uio_rd(nc, termt, ntt, dq_xb3, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('xc1')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'xc1')
                  call uio_rd(nc, termt, ntt, dq_xc1, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('xc2')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'xc2')
                  call uio_rd(nc, termt, ntt, dq_xc2, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case ('xc3')
                  call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, 'xc3')
                  call uio_rd(nc, termt, ntt, dq_xc3, &
                              ident, mode=3, name=name, unit=unit, &
                              outstr=outstr0, ierr=ierr0)
                  if (.not. norec0) call record('link')
          case default
                  processed=.False.
        end select
else if (qtype == intena .and. ndim == 2 .and. trim(label) == 'dimension') then
        allocate(inte2D(m1:n1,m2:n2))
        call uio_rd(nc, termt, ntt, inte2D, &
                    ident, mode=3, name=name, unit=unit, &
                    outstr=outstr0, ierr=ierr0)
        dq_m1 = inte2D(m1,m2)
        dq_n1 = inte2D(m1+1,m2)
        dq_m2 = inte2D(m1,m2+1)
        dq_n2 = inte2D(m1+1,m2+1)
        dq_m3 = inte2D(m1,m2+2)
        dq_n3 = inte2D(m1+1,m2+2)
        if (.not. norec0) call record('link')
        deallocate(inte2D)
else if (qtype == realna .and. ndim == 1 .and. n1<=m1 &
         .and. trim(label) == 'time') then
        call uio_rd(nc, termt, ntt, dq_time, &
                    ident, name=name, unit=unit, &
                    outstr=outstr0, ierr=ierr0)
        if (.not. norec0) call record('link')
else if (qtype == intena .and. ndim == 1 .and. n1<=m1 &
         .and. trim(label) == 'itime') then
        call uio_rd(nc, termt, ntt, dq_itime, &
                    ident, name=name, unit=unit, &
                    outstr=outstr0, ierr=ierr0)
        if (.not. norec0) call record('link')
else
        processed=.False.
end if
!
if (present(outstr)) outstr=outstr0
if (present(ierr)) ierr=ierr0
!print*, 'End wrapper'
end subroutine uio_dq_wrapper
subroutine uio_data_rd(nc, termt, ntt, ident, name, unit, norec, outstr, ierr)
use uio_bulk_module
use uio_base_module
use uio_siz_module
use uio_nam_module
!
implicit none
!
! --- I/O parameters ---
integer,                  intent(in)    :: nc
logical,       optional,  intent(in)    :: norec
character*(*),         intent(inout)    :: termt(2,nttmx)
integer,               intent(inout)    :: ntt
character*(*),           intent(out)    :: ident
character*(*), optional, intent(out)    :: name, unit
character*(*), optional, intent(out)    :: outstr
integer,       optional, intent(out)    :: ierr
!
! --- Local variables
!
character*(let)                         :: outstr0
integer                                 :: ierr0
logical                                 :: norec0, processed
character*(let)                         :: dimstr
integer, dimension(4)                   :: ilow, iup
integer                                 :: ncol, nlin
integer                                 :: i, j
!------------------------------------------------------------------------------
outstr0='done'
ierr0=0
qtype=termt(1,1)
ident=termt(1,2)
label=ident
!print*, 'Ident: ', ident
!print*, termt
!
! --- Additional keywords ---
call uio_exkeyw(termt,ntt, namena,name0)
call uio_exkeyw(termt,ntt, unitna,unit0)
if (present(name)) name=trim(name0)
if (present(unit)) unit=trim(unit0)
if (present(norec)) then
        norec0=norec
else
        norec0=.False.
endif
!
call uio_exkeyw(termt, ntt, dimna, dimstr)
!print*, 'dimstr: ', dimstr
call uio_st2dim(dimstr, ilow, iup, ndim=ndim)
m1=ilow(1)
m2=ilow(2)
m3=ilow(3)
m4=ilow(4)
n1=iup(1)
n2=iup(2)
n3=iup(3)
n4=iup(4)
!print*, 'Infos: ', qtype, ilow, iup, ndim
!print*, m1, m2, m3, n1, n2, n3
if (use_dq) then
        call uio_dq_wrapper(nc, termt, ntt, ident, name, unit, processed, &
                            norec=norec0, outstr=outstr0, ierr=ierr0)
        if (processed) return
endif
if (qtype == tabna) then
  call uio_rdtab(nc, termt, ntt, tab, ident, &
                 outstr=outstr0, ierr=ierr0)
  call uio_exkeyw(termt,ntt,namena,name0,mode=10)
  call record('table')
  call uio_tabaskhead(tab, icol=1, ncol=ncol, outstr=outstr0, ierr=ierr0)
  do j=2,ncol
        call uio_tabaskhead(tab, icol=j, type=qtype, ident=label, &
                            ncol=ncol, nlin=nlin, unit=unit0, &
                            termt=termt, ntt=ntt, outstr=outstr0, ierr=ierr0)
        call uio_exkeyw(termt,ntt,namena,name0,mode=10)
        select case (qtype)
          case (intena)
                allocate(inte1D(nlin))
                call uio_tabget(tab, label, inte1D, termt=termt, ntt=ntt, &
                                outstr=outstr0, ierr=ierr0)
                if (.not. norec0) call record('table:'//ident)
                deallocate(inte1D)
          case (realna)
                allocate(real1D(nlin))
                call uio_tabget(tab, label, real1D, termt=termt, ntt=ntt, &
                                outstr=outstr0, ierr=ierr0)
                if (.not. norec0) call record('table:'//ident)
                deallocate(real1D)
          !case (compna)
          !      allocate(comp1D(nlin))
          !      call uio_tabget(tab, label, comp1D, termt=termt, ntt=ntt, &
          !                      outstr=outstr0, ierr=ierr0)
          !      if (.not. norec0) call record('table:'//ident)
          !      deallocate(comp1D)
          case (charna)
                allocate(char1D(nlin,let))
                allocate(string1D(nlin))
                call uio_tabget(tab, label, string1D, termt=termt, ntt=ntt, &
                                outstr=outstr0, ierr=ierr0)
                do i=1,let
                    char1D(:,i)=string1D(:)(i:i)
                end do
                if (.not. norec0) call record('table:'//ident)
                deallocate(char1D, string1D)
        end select
  end do
  !call uio_tabaskhead(tab, icol=2, type=qtype, ident=label, &
  !                    ncol=ncol, nlin=nlin, unit=unit0, &
  !                    termt=termt, ntt=ntt, outstr=outstr0, ierr=ierr0)
  !call uio_exkeyw(termt,ntt,namena,name0,mode=10)
  !allocate(inte1D(nlin))
  !call uio_tabget(tab, ident, inte1D, termt=termt, ntt=ntt, &
  !                outstr=outstr0, ierr=ierr0)
  !call record('table:'//ident)
  !deallocate(inte1D)
else
  select case (ndim)
    case (1)
            if (n1 <= m1) then
                    select case (qtype)
                      case (intena)
                              call uio_rd(nc, termt, ntt, inte0, &
                                          ident, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, inte0
                      case (realna)
                              call uio_rd(nc, termt, ntt, real0, &
                                          ident, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, real0
                      case (compna)
                              call uio_rd(nc, termt, ntt, comp0, &
                                          ident, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, comp0
                      case (charna)
                              call uio_rd(nc, termt, ntt, char0, &
                                          ident, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, char0
                              if (trim(label) == 'box_id') then
                                      label=char0
                                      call record('rename')
                                      label=ident
                              endif
                      case default
                              call uio_skipda(nc, termt, ntt, &
                                              outstr=outstr0, ierr=ierr0)
                    end select
                    if (.not. norec0) call record('copy')
            else
                    select case (qtype)
                      case (intena)
                              allocate(inte1D(m1:n1))
                              call uio_rd(nc, termt, ntt, inte1D, &
                                          ident, mode=3, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, inte1D
                              if (.not. norec0) call record('copy')
                              deallocate(inte1D)
                      case (realna)
                              allocate(real1D(m1:n1))
                              call uio_rd(nc, termt, ntt, real1D, &
                                          ident, mode=3, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, real1D
                              if (.not. norec0) call record('copy')
                              deallocate(real1D)
                      case (compna)
                              allocate(comp1D(m1:n1))
                              call uio_rd(nc, termt, ntt, comp1D, &
                                          ident, mode=3, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, comp1D
                              if (.not. norec0) call record('copy')
                              deallocate(comp1D)
                      case (charna)
                              allocate(char1D(m1:n1,let))
                              allocate(string1D(m1:n1))
                              call uio_rd(nc, termt, ntt, string1D, &
                                          ident, mode=3, name=name, unit=unit, &
                                          outstr=outstr0, ierr=ierr0)
                              !print*, string1D
                              do i=1,let
                                  char1D(:,i)=string1D(:)(i:i)
                              end do
                              if (.not. norec0) call record('copy')
                              deallocate(char1D, string1D)
                      case default
                              call uio_skipda(nc, termt, ntt, &
                                              outstr=outstr0, ierr=ierr0)
                    end select
            endif
    case (2)
            select case (qtype)
              case (intena)
                      allocate(inte2D(m1:n1,m2:n2))
                      call uio_rd(nc, termt, ntt, inte2D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, inte2D
                      if (.not. norec0) call record('copy')
                      deallocate(inte2D)
              case (realna)
                      allocate(real2D(m1:n1,m2:n2))
                      call uio_rd(nc, termt, ntt, real2D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, real2D
                      if (.not. norec0) call record('copy')
                      deallocate(real2D)
              case (compna)
                      allocate(comp2D(m1:n1,m2:n2))
                      call uio_rd(nc, termt, ntt, comp2D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, comp2D
                      if (.not. norec0) call record('copy')
                      deallocate(comp2D)
              case (charna)
                      allocate(char2D(m1:n1,m2:n2,let))
                      allocate(string2D(m1:n1,m2:n2))
                      call uio_rd(nc, termt, ntt, string2D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, string2D
                      do i=1,let
                          char2D(:,:,i)=string2D(:,:)(i:i)
                      end do
                      if (.not. norec0) call record('copy')
                      deallocate(char2D, string2D)
              case default
                      call uio_skipda(nc, termt, ntt, &
                                      outstr=outstr0, ierr=ierr0)
            end select
    case (3)
            select case (qtype)
              case (intena)
                      allocate(inte3D(m1:n1,m2:n2,m3:n3))
                      call uio_rd(nc, termt, ntt, inte3D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, inte3D
                      if (.not. norec0) call record('copy')
                      deallocate(inte3D)
              case (realna)
                      allocate(real3D(m1:n1,m2:n2,m3:n3))
                      call uio_rd(nc, termt, ntt, real3D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, real3D
                      if (.not. norec0) call record('copy')
                      deallocate(real3D)
              case (compna)
                      allocate(comp3D(m1:n1,m2:n2,m3:n3))
                      call uio_rd(nc, termt, ntt, comp3D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, comp3D
                      if (.not. norec0) call record('copy')
                      deallocate(comp3D)
              case (charna)
                      allocate(char3D(m1:n1,m2:n2,m3:n3,let))
                      allocate(string3D(m1:n1,m2:n2,m3:n3))
                      call uio_rd(nc, termt, ntt, string3D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, string3D
                      do i=1,let
                          char3D(:,:,:,i)=string3D(:,:,:)(i:i)
                      end do
                      if (.not. norec0) call record('copy')
                      deallocate(char3D, string3D)
              case default
                      call uio_skipda(nc, termt, ntt, &
                                      outstr=outstr0, ierr=ierr0)
            end select
    case (4)
            select case (qtype)
              case (intena)
                      allocate(inte4D(m1:n1,m2:n2,m3:n3,m4:n4))
                      call uio_rd(nc, termt, ntt, inte4D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, inte4D
                      if (.not. norec0) call record('copy')
                      deallocate(inte4D)
              case (realna)
                      allocate(real4D(m1:n1,m2:n2,m3:n3,m4:n4))
                      call uio_rd(nc, termt, ntt, real4D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, real4D
                      if (.not. norec0) call record('copy')
                      deallocate(real4D)
              case (compna)
                      allocate(comp4D(m1:n1,m2:n2,m3:n3,m4:n4))
                      call uio_rd(nc, termt, ntt, comp4D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, comp4D
                      if (.not. norec0) call record('copy')
                      deallocate(comp4D)
              case (charna)
                      allocate(char4D(m1:n1,m2:n2,m3:n3,m4:n4,let))
                      allocate(string4D(m1:n1,m2:n2,m3:n3,m4:n4))
                      call uio_rd(nc, termt, ntt, string4D, &
                                  ident, mode=3, name=name, unit=unit, &
                                  outstr=outstr0, ierr=ierr0)
                      !print*, string4D
                      do i=1,let
                          char4D(:,:,:,:,i)=string4D(:,:,:,:)(i:i)
                      end do
                      if (.not. norec0) call record('copy')
                      deallocate(char4D, string4D)
              case default
                      call uio_skipda(nc, termt, ntt, &
                                      outstr=outstr0, ierr=ierr0)
            end select
    case default
            call uio_skipda(nc, termt, ntt, &
                            outstr=outstr0, ierr=ierr0)
  end select
endif
!
if (present(outstr)) outstr=outstr0
if (present(ierr)) ierr=ierr0
end subroutine uio_data_rd
!
subroutine uio_lbl_count(filename, channel, startlabel, endlabel, nlabel, &
                         maxlabel, outstr, ierr)
!use uio_base_module
use uio_bulk_module
!use uio_siz_module
use uio_nam_module
!use uio_table_module
!
implicit none
!
! --- I/O parameters ---
character*(*),           intent(in)     :: filename
integer, optional,       intent(in)     :: channel
character*(*),           intent(in)     :: startlabel, endlabel
integer, optional,       intent(in)     :: maxlabel
integer,                 intent(out)    :: nlabel
character*(*), optional, intent(out)    :: outstr
integer, optional,       intent(out)    :: ierr
!
! --- Local variables ---
character*(let)                         :: outstr0
integer                                 :: ierr0
integer                                 :: nc
integer                                 :: maxlabel0
logical                                 :: eof_flag
integer                                 :: ntt
character*(let)                         :: termt(2,nttmx)
character*(let)                         :: ident
!------------------------------------------------------------------------------
outstr0='done'
ierr0=0
!
! --- Open file ---
if (present(channel)) then
        nc=channel
else
        print*, 'Reopen file.......!'
        call uio_openrd(nc, filename, outstr=outstr0, ierr=ierr0)
endif
!
if (present(maxlabel)) then
        maxlabel0=maxlabel
else
        maxlabel0=-1
endif
!
! --- Count labels unitl EOF reached or maxlabel reached ---
nlabel=0
eof_flag=.false.
do while (.not. eof_flag .and. (maxlabel0 < 0 .or. nlabel < maxlabel0))
        !if (maxlabel0 >= 0 .and. nlabel >= maxlabel0) exit
        call uio_rdhd(nc, termt, ntt, outstr=outstr0, ierr=ierr0)
        if (outstr0 == 'eof') then
                eof_flag=.true.
        else
                qtype=termt(1,1)
                if (qtype == lablna) then
                        call uio_rdlabl(nc, termt, ntt, ident, &
                                        outstr=outstr0, ierr=ierr0)
                        if (ident == startlabel) then
                                nlabel = nlabel+1
                                if (maxlabel0 >= 0 .and. &
                                    nlabel >= maxlabel0) exit
                        endif
                else
                        call uio_skipda(nc, termt, ntt, &
                                        outstr=outstr0, ierr=ierr0)
                endif
        endif
end do
!
! --- Close file again ---
if (.not. present(channel)) call uio_closrd(nc, outstr=outstr0, ierr=ierr0)
!
if (present(outstr)) outstr=outstr0
if (present(ierr)) ierr=ierr0
end subroutine uio_lbl_count
recursive subroutine uio_struct_rd(filename, channel, headread, &
                                   startlabel, endlabel, nlabel, &
                                   substartlabel, subendlabel, outstr, ierr, &
                                   name_flag, unit_flag)
!use uio_base_module
use uio_bulk_module
use uio_siz_module
use uio_nam_module
use uio_table_module
!
implicit none
!
! --- I/O parameters ---
character*(*),           intent(in)     :: filename
integer, optional,       intent(in)     :: channel
logical, optional,       intent(in)     :: headread
character*(*), optional, intent(in)     :: startlabel, endlabel
integer, optional,       intent(in)     :: nlabel
character*(*), optional, intent(in)     :: substartlabel, subendlabel
character*(*), optional, intent(out)    :: outstr
integer, optional,       intent(out)    :: ierr
logical, optional,       intent(in)     :: name_flag, unit_flag
!
! --- Local variables ---
character*(let)                         :: outstr0
integer                                 :: ierr0
logical                                 :: name_flag0, unit_flag0
character*(let)                         :: startlabel0, endlabel0, &
                                           substartlabel0, subendlabel0
integer                                 :: nc
logical                                 :: eof_flag, found_flag, id_flag, &
                                           his_flag, des_flag, ver_flag
integer                                 :: position
integer                                 :: ntt
character*(let)                         :: termt(2,nttmx)
character*(let)                         :: ident
integer                                 :: nlabel0, number_of_labels
!------------------------------------------------------------------------------
outstr0='done'
ierr0=0
!
if (present(name_flag)) then
        name_flag0=name_flag
else
        name_flag0=.false.
endif
if (present(unit_flag)) then
        unit_flag0=unit_flag
else
        unit_flag0=.false.
endif
!
if (present(startlabel)) then
        startlabel0=startlabel
else
        startlabel0=''
endif
if (present(endlabel)) then
        endlabel0=endlabel
else
        endlabel0=''
endif
if (present(substartlabel)) then
        substartlabel0=substartlabel
else
        substartlabel0=''
endif
if (present(subendlabel)) then
        subendlabel0=subendlabel
else
        subendlabel0=''
endif
!
! --- Open file ---
if (present(channel)) then
        nc=channel
else
        call uio_openrd(nc, filename, outstr=outstr0, ierr=ierr0)
endif
!
! --- Read header ---
!print*, 'Read header'
if (present(headread)) then
        ! --- Flags set to zero ---
        eof_flag=.false.      ! --- Becomes true if EOF is reached
        found_flag=.false.    ! --- Becomes true if label is reached
        id_flag=.false.
        his_flag=.false.
        des_flag=.false.
        ver_flag=.false.
        !!! Create head structure
        label='head'
        call record('down')
        !print*, 'Down to head'
        ! --- Search until label found or EOF reached ---
        do while ((.not. eof_flag) .and. (.not. found_flag))
                call uio_chpos(nc, 'get', posout=position, &
                               outstr=outstr0, ierr=ierr0)
                if (ierr0>0) stop
                call uio_rdhd(nc, termt, ntt, outstr=outstr0, ierr=ierr0)
                if (outstr0 == 'eof') then
                        eof_flag=.true.
                else
                        qtype=termt(1,1)
                        if (qtype == lablna) then
                                found_flag=.true.
                                call uio_chpos(nc, 'goto', position, &
                                               outstr=outstr0, ierr=ierr0)
                        else if (qtype == fifona) then
                                call uio_skipda(nc, termt, ntt, &
                                                outstr=outstr0, ierr=ierr0)
                        else
                                !print*, 'Call to uio_data_rd in head'
                                call uio_data_rd(nc, termt, ntt, &
                                                ident, norec=.True., &
                                                outstr=outstr0, ierr=ierr0)
                                !print*, 'uio_data_rd okay in head'
                                !if (qtype == tabna) then
                                !        call uio_rdtab(nc, termt, ntt, tab, &
                                !                       ident, outstr=outstr0, &
                                !                       ierr=ierr0)
                                !else
                                !        call uio_data_rd(nc, termt, ntt, &
                                !                        ident, outstr=outstr0, &
                                !                        ierr=ierr0)
                                        !call uio_skipda(nc, termt, ntt, &
                                        !                outstr=outstr0, &
                                        !                ierr=ierr0)
                                        !print*, 'Unknown qtype :', qtype, ' ', &
                                        !termt(1,2)
                                        !!! Read data
                                        !call uio_rd(nc, termt, ntt, data, &
                                        !            ident, &
                                        !            name=name, &
                                        !            unit=unit, &
                                        !            outstr=outstr0, ierr=ierr0)
                                !endif
                                ! print*, qtype, ' ', ident
                                !
                                if (.not. (((ident == 'file_id')     &
                                            .and. id_flag ) .or. &
                                         ((ident == 'description') &
                                            .and. des_flag) .or. &
                                         ((ident == 'history')     &
                                            .and. his_flag) .or. &
                                         ((ident == 'version')     &
                                            .and. ver_flag))) then
                                        !!! Create structure
                                        !head_structure=create_struct( &
                                        !    head_structure, ident, data)
                                        call uio_chpos(nc, 'goto', position, &
                                               outstr=outstr0, ierr=ierr0)
                                        call uio_rdhd(nc, termt, ntt, &
                                                outstr=outstr0, ierr=ierr0)
                                        call uio_data_rd(nc, termt, ntt, &
                                                ident, outstr=outstr0, &
                                                ierr=ierr0)
                                endif
                                !
                                id_flag =id_flag  .or. (ident == 'file_id')
                                des_flag=des_flag .or. (ident == 'description')
                                his_flag=his_flag .or. (ident == 'history')
                                ver_flag=ver_flag .or. (ident == 'version')
                        endif
                endif
        end do
!print*, 'Requesting up'
call record('up')
!print*, 'Up okay'
if (present(outstr)) outstr=outstr0
if (present(ierr)) ierr=ierr0
endif
!
! --- Search until start label with specified number is reached ---
if (present(startlabel)) then
        !print*, 'Look for startlabel'
        ! --- Default value for label number ---
        if (present(nlabel)) then
                !print*, 'Present nlabel!!!'
                nlabel0=nlabel
                if (nlabel == 0) nlabel0 = 1
                if (nlabel < 0) then
                        call uio_lbl_count(filename, nc, startlabel, &
                                           endlabel, nlabel0, &
                                           outstr=outstr0, ierr=ierr0)
                        nlabel0 = nlabel0+nlabel+1
                endif
        else
                nlabel0=1
        endif
        ! --- End flags and counter set to zero ---
        eof_flag=.false.      ! --- Becomes true if EOF is reached ---
        found_flag=.false.    ! --- Becomes true if start label is found ---
        number_of_labels=0
        ! --- Search unitl label found or EOF reached ---
        do while ((.not. eof_flag) .and. (.not. found_flag))
                call uio_rdhd(nc, termt, ntt, outstr=outstr0, ierr=ierr0)
                if (outstr0 == 'eof') then
                        eof_flag=.true.
                else
                        qtype=termt(1,1)
                        if (qtype == lablna) then
                                call uio_rdlabl(nc, termt, ntt, ident, &
                                                outstr=outstr0, ierr=ierr0)
                                if (ident == startlabel0) then
                                        number_of_labels = number_of_labels+1
                                        if (number_of_labels >= nlabel0) then
                                                found_flag=.true.
                                        endif
                                endif
                        else
                                call uio_skipda(nc, termt, ntt, &
                                                outstr=outstr0, ierr=ierr0)
                        endif
                endif
        end do
        if (eof_flag) then
                ! --- Error ---
                if (.not. present(channel)) then
                        call uio_closrd(nc, outstr=outstr0, ierr=ierr0)
                endif
                outstr0='start label not found'
                ierr0=2
                if (.not. present(channel)) then
                        call uio_err(outstr=outstr0, ierr=ierr0, nc=nc, &
                                     routine='uio_struct_rd', errlevel=5)

                else
                        call uio_err(outstr=outstr0, ierr=ierr0, nc=nc, &
                                     routine='uio_struct_rd', errlevel=0)
                endif
                r_imodel = -1
        endif
!print*, 'Finished with startlabel'
endif
!
! --- Read all file data entries except labels and fileform entries ---
if(present(endlabel)) then
!print*, 'Read all file data entries'
eof_flag=.false.      ! --- Becomes true if EOF is reached ---
found_flag=.false.    ! --- Becomes true if start label is found ---
!!! Create structure
!if (n_elements(head_structure)) then begin
!  data_structure=create_struct(data_structure, 'head', head_structure)
!endif
do while ((.not. eof_flag) .and. (.not. found_flag))
        call uio_chpos(nc, 'get', posout=position, outstr=outstr0, ierr=ierr0)
        !print*, 'position: ', position
        call uio_rdhd(nc, termt, ntt, outstr=outstr0, ierr=ierr0)
        !print*, 'Termt12: ', termt(1,2)
        !
        if (outstr0 == 'eof') then
                eof_flag=.true.
        else
                qtype=termt(1,1)
                !print*, qtype
                if (qtype == fifona) then
                        call uio_skipda(nc, termt, ntt, &
                                        outstr=outstr0, ierr=ierr0)
                else if (qtype == lablna) then
                        !print*, 'Ident before: ', ident
                        call uio_rdlabl(nc, termt, ntt, ident, &
                                        outstr=outstr0, ierr=ierr0)
                        !print*, 'Ident after: ', ident
                        found_flag=(ident == endlabel)
                        if (.not. found_flag) then
                                ! --- If no end label is found, test if entry
                                !     is new start label ---
                                found_flag=(ident == startlabel)
                                if (found_flag) then
                                        call uio_chpos(nc, 'goto', position, &
                                                       outstr=outstr0, &
                                                       ierr=ierr0)
                                endif
                        endif
                        if ((.not. found_flag) .and. &
                              present(substartlabel)) then
                            if (ident == substartlabel) then
                                    ! --- Load box ---
                                    !print*, 'chpos'
                                    !print*, nc, position
                                    call uio_chpos(nc, 'goto', &
                                      position, outstr=outstr0, &
                                      ierr=ierr0)
                                    !print*, 'Request down in '//ident
                                    label=ident
                                    call record('down')
                                    !print*, 'Down ok, read'
                                    call uio_struct_rd( &
                                      filename, channel=nc, &
                                      startlabel=substartlabel, &
                                      endlabel=subendlabel, &
                                      name_flag=name_flag0, &
                                      unit_flag=unit_flag0)
                                    !print*, 'Read ok, request up'
                                    call record('up')
                                    !print*, 'Up ok'
                                    !!! Feed structure...
          !if ( (where(tag_names(data) eq 'box_id'))(0) ne 0) then begin
          !  ident=data.box_id
          !endif else begin
          !  ident='box'
          !endelse
          !;
          !data_structure=create_struct(data_structure, ident, data)
          !if name_flag ne 0 then $ 
          !  data_structure=create_struct(data_structure, ident+'__name', name) 
          !if unit_flag ne 0 then $ 
          !  data_structure=create_struct(data_structure, ident+'__unit', unit  )
                            endif
                    endif
            else
                    call uio_data_rd(nc, termt, ntt, &
                                    ident, outstr=outstr0, &
                                    ierr=ierr0)
                    !if (qtype == tabna) then
                    !        call uio_rdtab(nc, termt, ntt, tab, &
                    !                       ident, outstr=outstr0, &
                    !                       ierr=ierr0)
                    !else
                    !        call uio_data_rd(nc, termt, ntt, &
                    !                        ident, outstr=outstr0, &
                    !                        ierr=ierr0)
                            !call uio_skipda(nc, termt, ntt, &
                            !                outstr=outstr0, &
                            !                ierr=ierr0)
                            !!! Read data
                            !print*, 'Unknown qtype: ', qtype, ' ', termt(1,2)
                            !call uio_rd(nc, termt, ntt, data, &
                            !            ident, name=name, &
                            !            unit=unit, &
                            !            outstr=outstr0, ierr=ierr0)
                    !endif
                    !print*, qtype, ' ', ident
                    !if (qtype == 'character') print*, '  ', data
                    if (scan(ident, '-') /= 0) then
                            !print*, 'Misshaped identifier: ', qtype, &
                            !' ', data
                    else
                            !!! Populate structure...
          !data_structure=create_struct(data_structure, ident, data)
          !if name_flag ne 0 then $ 
          !  data_structure=create_struct(data_structure, ident+'__name', name) 
          !if unit_flag ne 0 then $ 
          !  data_structure=create_struct(data_structure, ident+'__unit', unit  )
                    endif
            endif
    endif
end do
endif
!
if (ierr0 > 0) then
        call uio_err(outstr=outstr0, ierr=ierr0, nc=nc, &
                     routine='uio_struc_rd', errlevel=5)
endif
!
! --- Close file again ---
if (.not. present(channel)) call uio_closrd(nc, outstr=outstr0, ierr=ierr0)
!
if (present(outstr)) outstr=outstr0
if (present(ierr)) ierr=ierr0
!
end subroutine uio_struct_rd
end module uio_reader_module
