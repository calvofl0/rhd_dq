module visit_vars
!
implicit none
!
type visit_data
  character(len=80)                     :: name="", unit="", mesh=""
  real, pointer, dimension(:)           :: p
  type(visit_data), pointer             :: head, tail, next
  logical                               :: is_head=.false., is_tail=.true.
end type visit_data
type visit_mesh
  character(len=80)                     :: name="", unitX="", unitY="", unitZ=""
  integer                               :: dims=3
  real, pointer, dimension(:)           :: xb1, xb2, xb3
  type(visit_mesh), pointer             :: head, tail, next
  logical                               :: is_head=.false., is_tail=.true.
end type visit_mesh
type visit_sim
  character(len=80)                     :: name="", comment="", path=""
end type visit_sim
!
type(visit_sim)                         :: simulation
type(visit_mesh), target                :: meshes
type(visit_data), target                :: scalars, vectors
end module visit_vars
!
subroutine pv_reset()
use visit_vars
!
!
implicit none
!
type(visit_mesh), pointer               :: mesh, mesh0
type(visit_data), pointer               :: scalar, scalar0, vector, vector0
logical                                 :: resume
!
simulation%name = ""
simulation%comment = ""
simulation%path = ""
mesh=>meshes%next
resume = .false.
if (meshes%is_tail .eqv. .false.) then; do
  mesh0=>mesh%next
  if (mesh%is_tail .eqv. .true.) resume=.true.
  deallocate(mesh)
  if (resume) exit
  mesh=>mesh0
  nullify(mesh0)
end do; end if
meshes%name = ""
meshes%unitX = ""
meshes%unitY = ""
meshes%unitZ = ""
meshes%dims = 3
nullify(meshes%xb1)
nullify(meshes%xb2)
nullify(meshes%xb3)
nullify(meshes%head)
nullify(meshes%tail)
nullify(meshes%next)
meshes%is_head = .false.
meshes%is_tail = .true.
scalar=>scalars%next
resume = .false.
if (scalars%is_tail .eqv. .false.) then; do
  scalar0=>scalar%next
  if (scalar%is_tail .eqv. .true.) resume=.true.
  deallocate(scalar)
  if (resume) exit
  scalar=>scalar0
  nullify(scalar0)
end do; end if
scalars%name = ""
scalars%unit = ""
scalars%mesh = ""
nullify(scalars%p)
nullify(scalars%head)
nullify(scalars%tail)
nullify(scalars%next)
scalars%is_head = .false.
scalars%is_tail = .true.
resume = .false.
if (vectors%is_tail .eqv. .false.) then; do
  vector0=>vector%next
  if (vector%is_tail .eqv. .true.) resume=.true.
  deallocate(vector)
  if (resume) exit
  vector=>vector0
  nullify(vector0)
end do; end if
vectors%name = ""
vectors%unit = ""
vectors%mesh = ""
nullify(vectors%p)
nullify(vectors%head)
nullify(vectors%tail)
nullify(vectors%next)
vectors%is_head = .false.
vectors%is_tail = .true.
end subroutine pv_reset
!
subroutine pv_regSim(name, comment, path)
use visit_vars
!
implicit none
!
character(len=*)                        :: name, comment, path
!
simulation%name = name
simulation%comment = comment
simulation%path = path
end subroutine pv_regSim
!
subroutine pv_regMesh3D(name, xb1, xb2, xb3, unitX, unitY, unitZ)
use visit_vars
!
implicit none
!
character(len=*), intent(in)            :: name, unitX, unitY, unitZ
real, target, dimension(:), intent(in)  :: xb1, xb2, xb3
type(visit_mesh), pointer               :: mesh, new_mesh
!
if (meshes%is_head .eqv. .true.) then
  allocate(new_mesh)
  new_mesh%name = name
  new_mesh%unitX = unitX
  new_mesh%unitY = unitY
  new_mesh%unitZ = unitZ
  new_mesh%dims = 3
  new_mesh%xb1 => xb1
  new_mesh%xb2 => xb2
  new_mesh%xb3 => xb3
  new_mesh%head => meshes
  new_mesh%next => new_mesh
  new_mesh%tail => new_mesh
  mesh => meshes
  do
    mesh%tail => new_mesh
    mesh => mesh%next
    if (mesh%is_tail .eqv. .true.) exit
  end do
  mesh%next => new_mesh
  mesh%is_tail = .false.
else
  meshes%is_head = .true.
  meshes%name = name
  meshes%unitX = unitX
  meshes%unitY = unitY
  meshes%unitZ = unitZ
  meshes%xb1 => xb1
  meshes%xb2 => xb2
  meshes%xb3 => xb3
  meshes%head => meshes
  meshes%next => meshes
  meshes%tail => meshes
end if
end subroutine pv_regMesh3D
!
subroutine pv_regMesh2D(name, xb1, xb2, unitX, unitY)
use visit_vars
!
implicit none
!
character(len=*), intent(in)            :: name, unitX, unitY
real, target, dimension(:), intent(in)  :: xb1, xb2
type(visit_mesh), pointer               :: mesh, new_mesh
!
!print*, "Name: ", name, ", unitX: ", unitX, ", unitY: ", unitY
if (meshes%is_head .eqv. .true.) then
  allocate(new_mesh)
  new_mesh%name = name
  new_mesh%unitX = unitX
  new_mesh%unitY = unitY
  new_mesh%dims = 2
  new_mesh%xb1 => xb1
  new_mesh%xb2 => xb2
  new_mesh%head => meshes
  new_mesh%next => new_mesh
  new_mesh%tail => new_mesh
  mesh => meshes
  do
    mesh%tail => new_mesh
    mesh => mesh%next
    if (mesh%is_tail .eqv. .true.) exit
  end do
  mesh%next => new_mesh
  mesh%is_tail = .false.
else
  meshes%is_head = .true.
  meshes%name = name
  meshes%unitX = unitX
  meshes%unitY = unitY
  meshes%dims = 2
  meshes%xb1 => xb1
  meshes%xb2 => xb2
  meshes%head => meshes
  meshes%next => meshes
  meshes%tail => meshes
end if
mesh => meshes
do
  mesh => mesh%next
  if (mesh%is_tail .eqv. .true.) exit
end do
end subroutine pv_regMesh2D
!
subroutine pv_regScalar(name, mesh, p, unit)
use visit_vars
!
implicit none
!
character(len=*), intent(in)            :: name, unit, mesh
real, target, dimension(:), intent(in)  :: p
type(visit_data), pointer               :: scalar, new_scalar
!
if (scalars%is_head .eqv. .true.) then
  allocate(new_scalar)
  new_scalar%name = name
  new_scalar%unit = unit
  new_scalar%mesh = mesh
  new_scalar%p => p
  new_scalar%head => scalars
  new_scalar%next => new_scalar
  new_scalar%tail => new_scalar
  scalar => scalars
  do
    scalar%tail => new_scalar
    scalar => scalar%next
    if (scalar%is_tail .eqv. .true.) exit
  end do
  scalar%next => new_scalar
  scalar%is_tail = .false.
else
  scalars%is_head = .true.
  scalars%name = name
  scalars%unit = unit
  scalars%mesh = mesh
  scalars%p => p
  scalars%head => scalars
  scalars%next => scalars
  scalars%tail => scalars
end if
!
end subroutine pv_regScalar
!
subroutine pv_regVector(name, mesh, p, unit)
use visit_vars
!
implicit none
!
character(len=*), intent(in)            :: name, unit, mesh
real, target, dimension(:), intent(in)  :: p
type(visit_data), pointer               :: vector, new_vector
!
if (vectors%is_head .eqv. .true.) then
  allocate(new_vector)
  new_vector%name = name
  new_vector%unit = unit
  new_vector%mesh = mesh
  new_vector%p => p
  new_vector%head => vectors
  new_vector%next => new_vector
  new_vector%tail => new_vector
  vector => vectors
  do
    vector%tail => new_vector
    vector => vector%next
    if (vector%is_tail .eqv. .true.) exit
  end do
  vector%next => new_vector
  vector%is_tail = .false.
else
  vectors%is_head = .true.
  vectors%name = name
  vectors%unit = unit
  vectors%mesh = mesh
  vectors%p => p
  vectors%head => vectors
  vectors%next => vectors
  vectors%tail => vectors
end if
!
end subroutine pv_regVector
!
function elapsed_time(t0)
!
implicit none
!
real                                    :: elapsed_time
integer, intent(in)                     :: t0
integer                                 :: t, count_rate, count_max
!
call system_clock(t, count_rate, count_max)
if (t>=t0) elapsed_time=real(t-t0)/real(count_rate)
if (t<t0) elapsed_time=(real(count_max)-real(t0-t))/real(count_rate)
end function elapsed_time
!
subroutine pv_runsim()
use visit_vars
!
implicit none
!
include "visitfortransimV2interface.inc"
!
interface
function elapsed_time(t0)
real                                    :: elapsed_time
integer, intent(in)                     :: t0
end function elapsed_time
end interface
!
! VisIt variables
!
integer visitstate, result, runflag, blocking, err
!
integer                                 :: t0
!
call system_clock(t0)
!
err = visitsetoptions("-cli -nowin -nosplash -nowindowmetrics", 38)
err = visitsetupenv()
err = visitinitializesim(trim(simulation%name), len(trim(simulation%name)), &
        trim(simulation%comment), len(trim(simulation%comment)), &
        trim(simulation%path), len(trim(simulation%path)), &
        VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN, &
        VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN, &
        VISIT_F77NULLSTRING, VISIT_F77NULLSTRINGLEN)
!
runflag = 1
!
! Initialization
!
MAIN_LOOP: do
!
if(runflag.eq.1) then
  blocking = 0
  if (elapsed_time(t0)>60) exit
else
  blocking = 1
endif
!
visitstate = visitdetectinput(blocking, -1)     
!
if (visitstate.lt.0) then
  exit MAIN_LOOP
elseif (visitstate.eq.0) then
  !call simulate_one_timestep()
  call sleep(1)
elseif (visitstate.eq.1) then
  runflag = 0
  result = visitattemptconnection()
  !if (result.eq.1) then
  !  write (6,*) 'VisIt connected!'
  !else
  !  write (6,*) 'VisIt did not connect!'
  !endif
elseif (visitstate.eq.2) then
  runflag = 0
  if (visitprocessenginecommand().eq.0) then
    result = visitdisconnect()
    runflag = 1
    exit
  endif
endif
end do MAIN_LOOP
!
end subroutine pv_runsim

subroutine simulate_one_timestep()
!write (6,*) 'Simulating time step'
call sleep(1)
end subroutine simulate_one_timestep

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! These functions must be defined to satisfy the visitfortransimV2interface lib.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!---------------------------------------------------------------------------
! visitcommandcallback
!---------------------------------------------------------------------------
subroutine visitcommandcallback (cmd, lcmd, args, largs)
implicit none
character*8 cmd, args
integer     lcmd, largs
end subroutine

!---------------------------------------------------------------------------
! visitbroadcastintfunction
!---------------------------------------------------------------------------
integer function visitbroadcastintfunction(value, sender)
implicit none
integer value, sender
!     REPLACE WITH MPI COMMUNICATION IF SIMULATION IS PARALLEL
visitbroadcastintfunction = 0
end function

!---------------------------------------------------------------------------
! visitbroadcaststringfunction
!---------------------------------------------------------------------------
integer function visitbroadcaststringfunction(str, lstr, sender)
implicit none
character*8 str
integer     lstr, sender
!     REPLACE WITH MPI COMMUNICATION IF SIMULATION IS PARALLEL
visitbroadcaststringfunction = 0
end function

!---------------------------------------------------------------------------
! visitslaveprocesscallback
!---------------------------------------------------------------------------
subroutine visitslaveprocesscallback ()
implicit none
!     REPLACE WITH MPI COMMUNICATION IF SIMULATION IS PARALLEL
end subroutine

!---------------------------------------------------------------------------
! visitactivatetimestep
!---------------------------------------------------------------------------
integer function visitactivatetimestep()
implicit none
include "visitfortransimV2interface.inc"
visitactivatetimestep = VISIT_OKAY
end function

!---------------------------------------------------------------------------
! visitgetmetadata
!---------------------------------------------------------------------------
integer function visitgetmetadata()
use visit_vars
implicit none
include "visitfortransimV2interface.inc"
! SIMSTATE common block
integer runflag, simcycle
real simtime
type(visit_mesh), pointer               :: mesh
type(visit_data), pointer               :: scalar, vector
integer                                 :: i
! Local variables
integer md, m3d, vmd, err
print *, 'Inside visitgetmetadata'
mesh=>meshes
scalar=>scalars
vector=>vectors
!
if(visitmdsimalloc(md).eq.VISIT_OKAY) then
  err = visitmdsimsetcycletime(md, simcycle, simtime)
  if(runflag.eq.1) then
    err = visitmdsimsetmode(md, VISIT_SIMMODE_RUNNING)
  else
    err = visitmdsimsetmode(md, VISIT_SIMMODE_STOPPED)
  endif
endif
!
! Set real mesh's properties
do
  if(visitmdmeshalloc(m3d).eq.VISIT_OKAY) then
    err = visitmdmeshsetname(m3d, trim(mesh%name), len(trim(mesh%name)))
    err = visitmdmeshsetmeshtype(m3d, VISIT_MESHTYPE_RECTILINEAR)
    err = visitmdmeshsettopologicaldim(m3d, mesh%dims)
    err = visitmdmeshsetspatialdim(m3d, mesh%dims)
    err = visitmdmeshsetxunits(m3d, trim(mesh%unitX), len(trim(mesh%unitX)))
    err = visitmdmeshsetxlabel(m3d, "Length", 6)
    if(mesh%dims>1) then
      err = visitmdmeshsetyunits(m3d, trim(mesh%unitY), len(trim(mesh%unitY)))
      err = visitmdmeshsetylabel(m3d, "Width", 5)
    endif
    if(mesh%dims>2) then
      err = visitmdmeshsetzunits(m3d, trim(mesh%unitZ), len(trim(mesh%unitZ)))
      err = visitmdmeshsetzlabel(m3d, "Depth", 5)
    endif
    err = visitmdsimaddmesh(md, m3d)
  endif
  if (mesh%is_tail .eqv. .true.) exit
  mesh => mesh%next
end do
!
! Add scalar variables on mesh3d
if (scalar%is_head .eqv. .true.) then; do
  if(visitmdvaralloc(vmd).eq.VISIT_OKAY) then
    err = visitmdvarsetname(vmd, trim(scalar%name), len(trim(scalar%name)))
    err = visitmdvarsetmeshname(vmd, trim(scalar%mesh), len(trim(scalar%mesh)))
    err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)
    err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_ZONE)

    err = visitmdsimaddvariable(md, vmd)
  endif
  if (scalar%is_tail .eqv. .true.) exit
  scalar => scalar%next
end do; end if
!
! Add vector variables on mesh3d
if (vector%is_head .eqv. .true.) then; do
  if(visitmdvaralloc(vmd).eq.VISIT_OKAY) then
    err = visitmdvarsetname(vmd, trim(vector%name), len(trim(vector%name)))
    err = visitmdvarsetmeshname(vmd, trim(vector%mesh), len(trim(vector%mesh)))
    err = visitmdvarsettype(vmd, VISIT_VARTYPE_VECTOR)
    err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_ZONE)

    err = visitmdsimaddvariable(md, vmd)
  endif
  if (vector%is_tail .eqv. .true.) exit
  vector => vector%next
end do; end if
!
visitgetmetadata = md
end function

!---------------------------------------------------------------------------
! visitgetmesh
!---------------------------------------------------------------------------
integer function visitgetmesh(domain, name, lname)
use visit_vars
implicit none
character*8 name
integer     domain, lname
include "visitfortransimV2interface.inc" 
type(visit_mesh), pointer               :: mesh
integer h, x, y, z, err, i
!
print *, 'Inside visitgetmesh'
mesh=>meshes
h = VISIT_INVALID_HANDLE
!
do
  if(visitstrcmp(name, lname, trim(mesh%name), len(trim(mesh%name))).eq.0) then
   if(visitrectmeshalloc(h).eq.VISIT_OKAY) then
     print *, 'Starting allocation...'
     err = visitvardataalloc(x)
     print *, 'Allocation finished.'
     err = visitvardatasetf(x,VISIT_OWNER_SIM,1,size(mesh%xb1),mesh%xb1)
     if(mesh%dims>1) then
       print *, 'Starting allocation...'
       err = visitvardataalloc(y)
       print *, 'Allocation finished.'
       err = visitvardatasetf(y,VISIT_OWNER_SIM,1,size(mesh%xb2),mesh%xb2)
     endif
     if(mesh%dims>2) then
       print *, 'Starting allocation...'
       err = visitvardataalloc(z)
       print *, 'Allocation finished.'
       err = visitvardatasetf(z,VISIT_OWNER_SIM,1,size(mesh%xb3),mesh%xb3)
     endif
     print *, 'Data sent.'

     ! 1D meshes are not supported by VisIt
     !if(mesh%dims==1) err = visitrectmeshsetcoordsx(h, x)
     if(mesh%dims==2) err = visitrectmeshsetcoordsxy(h, x, y)
     if(mesh%dims==3) err = visitrectmeshsetcoordsxyz(h, x, y, z)
     print *, 'Mesh sent.'
   endif
  endif
  if (mesh%is_tail .eqv. .true.) exit
  mesh => mesh%next
end do
visitgetmesh = h
end function

!---------------------------------------------------------------------------
! visitgetvariable
!---------------------------------------------------------------------------
integer function visitgetvariable(domain, name, lname)
use visit_vars
implicit none
character*8 name
integer     domain, lname
include "visitfortransimV2interface.inc"
! Local vars
type(visit_data), pointer               :: scalar, vector
integer h, nvals, err, i

print *, 'Inside visitgetvariable'
scalar=>scalars
vector=>vectors
h = VISIT_INVALID_HANDLE
!
do
  if(visitstrcmp(name, lname, trim(scalar%name), len(trim(scalar%name))).eq.0) then
    if(visitvardataalloc(h).eq.VISIT_OKAY) then
      nvals = size(scalar%p)
      print*, loc(scalar%p)
      err = visitvardatasetf(h, VISIT_OWNER_VISIT,1,nvals,scalar%p)
      print *, 'All seems ok...'
    endif
  endif
  if (scalar%is_tail .eqv. .true.) exit
  scalar => scalar%next
end do
!
do
  if(visitstrcmp(name, lname, trim(vector%name), len(trim(vector%name))).eq.0) then
    if(visitvardataalloc(h).eq.VISIT_OKAY) then
      nvals = int(size(vector%p)/3)
      print*, loc(vector%p)
      err = visitvardatasetf(h, VISIT_OWNER_VISIT,3,nvals,vector%p)
      print *, 'All seems ok...'
    endif
  endif
  if (vector%is_tail .eqv. .true.) exit
  vector => vector%next
end do
!
visitgetvariable = h
end function

!---------------------------------------------------------------------------
! visitgetcurve
!---------------------------------------------------------------------------
integer function visitgetcurve(name, lname)
implicit none
character*8 name
integer     lname
include "visitfortransimV2interface.inc"
visitgetcurve = VISIT_INVALID_HANDLE
end function

!---------------------------------------------------------------------------
! visitgetdomainlist
!---------------------------------------------------------------------------
integer function visitgetdomainlist(name, lname)
implicit none
character*8 name
integer     lname
include "visitfortransimV2interface.inc"
visitgetdomainlist = VISIT_INVALID_HANDLE
end function

!---------------------------------------------------------------------------
! visitgetdomainbounds
!---------------------------------------------------------------------------
integer function visitgetdomainbounds(name, lname)
implicit none
character*8 name
integer     lname
include "visitfortransimV2interface.inc"
visitgetdomainbounds = VISIT_INVALID_HANDLE
end function

!---------------------------------------------------------------------------
! visitgetdomainnesting
!---------------------------------------------------------------------------
integer function visitgetdomainnesting(name, lname)
implicit none
character*8 name
integer     lname
include "visitfortransimV2interface.inc"
visitgetdomainnesting = VISIT_INVALID_HANDLE
end function

!---------------------------------------------------------------------------
! visitgetmaterial
!---------------------------------------------------------------------------
integer function visitgetmaterial(domain, name, lname)
implicit none
character*8 name
integer     domain, lname
include "visitfortransimV2interface.inc"
visitgetmaterial = VISIT_INVALID_HANDLE
end function
