!**************************************************************************************************
!
!   #####   #    #  #####          #####    ####
!   #    #  #    #  #    #         #    #  #    #
!   #    #  ######  #    #         #    #  #    #
!   #####   #    #  #    #         #    #  #    #
!   #   #   #    #  #    #         #    #  #    #
!   #    #  #    #  #####  ######  #####    ###### _module
!
!   Modules for Radiation-Hydrodynamics Code
!
!**************************************************************************************************
!   Fortran90
!   Flavio Calvo:          Geneva, Locarno, Freiburg
!   2014-02-01
!
!   Conventions, external modules, and rhd_mhd_Bb2Bc subroutine are from CO5BOLD MHD module.
!**************************************************************************************************
! Modules:
!   rhd_dq_module:  Derived quantities external routines for RHD
!**************************************************************************************************
!
!------*************-------------------------------------------------------------------------------
module rhd_dq_module
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_module ('rhd_derived_quantities_module')
!
! PURPOSE:
!   Derived quantities external routines for RHD
!
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   use rhd_dq_module
!
! VARIABLES:
!   None
!
! MODULES:
!   rhd_dq_box:        Box variables (when PYBOLD not defined)
!   box:               Box variables (when PYBOLD defined)
!   gasinter_module:
!   tabinter_module:
!   uio_bulk_module:
!   opta_module:
!   rhd_gl_module:
!   rhd_sub_module:
!   rhd_action_module:
!   str_module:
!
! SUBROUTINES: (contained)
!   rhd_dq_init:       Extract relevant hydrodynamic variables from FULL file
!   rhd_dq_vars:       Set HD vars and compute derived quantities (T, P, rho, ...)
!   rhd_eosopa_RdData: Read EOS and OPA using given parameter file
!   rhd_keep_only:     Deallocate unused variables (free memory)
!   rhd_int_tau:       Integrate opacities to get optical depth
!   diff3D:            Differentiate cell-centered quantity
!   resample1D:        Resample a 1D quantity
!   resample3D:        Resample a 3D quantity
!
! MODIFICATION HISTORY:
!   2014-02-01 (F.C. Freiburg) Written
!--------------------------------------------------------------------------------------------------
!
#ifndef PYBOLD
use rhd_dq_box
#else
use box
#endif
!
implicit none
!
public :: rhd_dq_init
public :: rhd_dq_vars
!
contains

!----------**************--------------------------------------------------------------------------
subroutine rhd_eosopa_RdData(parfile, nc_p)
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_eosopa_RdData ('rhd_eosopa_Read_Data')
!
! PURPOSE:
!   Read eos and opa from location specified in parameter file.
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_eosopa_RdData(parfile, nc_p)
!
! INPUT:
!   parfile:        ('parameter_file') character, name of parameter file
!   nc_p:           ('number_channel_print') integer, logical number for text output
!
! OUTPUT:
!
! ROUTINES:
!   tabinter_reset:        Set all quantities in EOS table to zero
!   uio_init:              Initialization procedure for input/output routines
!   rhd_par_Read:          Read parameters from parameter file and put them into global variable
!   rhd_par_Toaction:      Prepare action variable from input parameters
!   str_pathandfile:       Merge path name and file name
!   tabinter_rdcoeff:      Read coefficients for bicubic interpolation of EOS
!   dfopta:                Read opacity table
!   
!
! MODULES:
!   tabinter_module:       Routines for handling of interpolation coefficients
!   uio_bulk_module:       Main set of uio-routines
!   opta_module:           Opacity interpolation routines
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   rhd_sub_module:        Collection of additional subroutines
!   rhd_action_module:     Routines to handle the control structure 'action'
!   str_module:            String handling routines
!
! SIDE EFFECTS:
!   Open file, read data from it, close file.
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Initialize I/O, open par file (in uio form), extract eos and opa location, 
!   open eos/opa file, read eos/opa, close eos/opa file
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-02-01 (F.C. Freiburg) Written
!--------------------------------------------------------------------------------------------------
use tabinter_module
use uio_bulk_module
use opta_module
use rhd_gl_module
use rhd_sub_module
use rhd_action_module
use str_module
!
implicit none
!
! --- I/O parameters ---
character(len=*),        intent(in)     :: parfile
integer,                 intent(in)     :: nc_p
!
! --- Local variables ---
character(len=256)                      :: filename
character(len=80)                       :: outstr
integer                                 :: ncopa, ierr, nband
type(action_type)                       :: action
!--------------------------------------------------------------------------------------------------
!
! === Initialization ==============================================================================
!
call tabinter_reset(1, nc_p)
#ifndef PYBOLD
call uio_init(progrm='rhd_dq_module')
#endif
!
! === Read PAR file ===============================================================================
!
call rhd_par_Read(parfile, nc_p)
! --- Prepare action variable ---
call rhd_par_Toaction(action)
!
! === Read EOS file ===============================================================================
!
call str_pathandfile(par%eospath, par%eosfile, filename)
call tabinter_rdcoeff(filename, nc_p, outstr, ierr)
! === Read OPA file ===============================================================================
!
call str_pathandfile(par%opapath, par%opafile, filename)
call uio_chunit('get', nc1out=ncopa)
open(ncopa, file=filename, form='formatted', status='old', iostat=ierr)
call dfopta(2, nc_p, ncopa, nband)
close(ncopa)
call uio_chunit('return', nc1in=ncopa)
end subroutine rhd_eosopa_RdData
!
!----------**************--------------------------------------------------------------------------
#ifdef PYBOLD
subroutine rhd_dq_init(parfile, file, nmodel, imodel, nc_p)
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_init ('rhd_derived_quantities_initialization')
!
! PURPOSE:
!   Prepare to use FULL file to extract all relevant hydrodynamic variables
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_dq_init(parfile, file, nc_p, imodel=imodel)
!
! INPUT:
!   parfile:        ('parameter_file') character, name of parameter file
!   file:           ('file') character, name of FULL file
!   nmodel:         ('number_model') number of model datasets in file
!   nc_p:           ('number_channel_print') integer, logical number for text output
!   imodel:         ('time_index') integer, time index of requested model
!
! OUTPUT:
!
! ROUTINES:
!   rhd_eosopa_RdData:     Read EOS and OPA files
!   rhd_box_Read:          Read box from file in UIO form
!   rhd_keep_only:         Remove unused variables
!
! MODULES:
!   const_module:          Global physical and mathematical constants, units
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   !rhd_box_module:        Box-handling routines
!   !rhd_sub_module:        Collection of additional subroutines
!   rhd_action_module:     Routines to handle the control structure 'action'
!
! SIDE EFFECTS:
!   !Open file, read data from it, close file.
!
! RESTRICTIONS:
!
! PROCEDURE:
!   If not present, determine imodel. !Read full file, allocate memory, copy needed
!   variables, deallocate memory.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-08-17 (F.C. Locarno) Written
!--------------------------------------------------------------------------------------------------
use const_module
use rhd_gl_module
!use rhd_box_module
!use rhd_sub_module
use rhd_action_module
!
implicit none
!
! --- I/O parameters ---
character(len=*),        intent(in)     :: parfile
character(len=*),        intent(in)     :: file
integer,                 intent(in)     :: nmodel
integer,                 intent(in)     :: nc_p
integer, optional,       intent(in)     :: imodel
!
! --- Local variables ---
integer                                 :: ncin, imodel0
!type(box_type), pointer                 :: model
real, allocatable, dimension(:,:,:)     :: dTde
logical                                 :: Bb_flag
real, dimension(1)                      :: tmp_array1
!--------------------------------------------------------------------------------------------------
! === Extract model ===============================================================================
!
! Reload parfile and data file if those have changed
!
if(dq_status.ne.0) then
  if((trim(parfile) /= trim(dq_parfile)).or.(trim(file) /= trim(dq_file))) then
    call rhd_keep_only(base=.false.)
  endif
endif
if(dq_status.eq.0) then
  call rhd_eosopa_RdData(trim(parfile), nc_p)
  dq_C_radHtautop = par%C_radHtautop
  dq_grav = par%grav
  dq_Bb_flag = ((par%A%hdscheme == hdscheme_RoeMHD) .or. &
                              (par%A%hdscheme == hdscheme_HLLMHD) .or. &
                              (par%A%hdscheme == hdscheme_CenMHD))
endif
!
! We allow for negative values for imodel, which are interpretd as counting
! backwards from last available model. Ex.: imodel=-1 means last available
! model.
!
if (present(imodel)) then
  imodel0 = imodel
else
  if(dq_status.eq.0) then
    imodel0 = par%istep_in_start
  else
    imodel0 = dq_imodel
  endif
endif
!!! Set nmodel
if (imodel0 < 0) then
  imodel0 = nmodel+imodel0+1
  if(dq_status.eq.0) then
    dq_imodel = imodel0
  endif
endif
if ((imodel0.ne.dq_imodel).and.(dq_imodel.ge.0)) then
  call rhd_keep_only(base=.false.)
endif
dq_imodel=imodel0
if (dq_status.eq.0) then
  dq_parfile=parfile
  dq_file=file
  dq_nc_p=nc_p
  dq_nmodel = nmodel
endif
end subroutine rhd_dq_init
!
!----------**************--------------------------------------------------------------------------
subroutine rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, var)
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_vars_init ('rhd_derived_quantities_variables_initialization')
!
! PURPOSE:
!   Initialize box variables
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_dq_vars_init(m1, n1, m2, n2, m3, n3, var)
!
! INPUT:
!   m1:              ('m1')
!   n1:              ('n1')
!   m2:              ('m2')
!   n2:              ('n2')
!   m3:              ('m3')
!   n3:              ('n3')
!   !parfile:        ('parameter_file') character, name of parameter file
!   !file:           ('file') character, name of FULL file
!   !nc_p:           ('number_channel_print') integer, logical number for text output
!   !imodel:         ('time_index') integer, time index of requested model
!
! ROUTINES:
!   !rhd_eosopa_RdData:     Read EOS and OPA files
!   !rhd_box_Read:          Read box from file in UIO form
!   !rhd_keep_only:         Remove unused variables
!
! MODULES:
!   !const_module:          Global physical and mathematical constants, units
!   !rhd_gl_module:         Global type definitions and parameters for RHD
!   !rhd_box_module:        Box-handling routines
!   !rhd_sub_module:        Collection of additional subroutines
!   !rhd_action_module:     Routines to handle the control structure 'action'
!
! SIDE EFFECTS:
!   !Open file, read data from it, close file.
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Compute dq_deltax, dq_deltay, dq_deltaz and change magnetic field units
!   !If not present, determine imodel. !Read full file, allocate memory, copy needed
!   !variables, deallocate memory.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-08-17 (F.C. Locarno) Written
!--------------------------------------------------------------------------------------------------
implicit none
!
integer                                 :: m1, m2, m3, n1, n2, n3
! --- I/O parameters ---
character*(*)                           :: var
!--------------------------------------------------------------------------------------------------
if (dq_status.eq.0) then
  select case (var)
    case ('ei')
          allocate(dq_ei(m1:n1, m2:n2, m3:n3))
    case ('rho')
          allocate(dq_rho(m1:n1, m2:n2, m3:n3))
    case ('v1')
          allocate(dq_v1(m1:n1, m2:n2, m3:n3))
    case ('v2')
          allocate(dq_v2(m1:n1, m2:n2, m3:n3))
    case ('v3')
          allocate(dq_v3(m1:n1, m2:n2, m3:n3))
    case ('xb1')
          allocate(dq_xb1(m1:n1))
    case ('xb2')
          allocate(dq_xb2(m2:n2))
    case ('xb3')
          allocate(dq_xb3(m3:n3))
    case ('xc1')
          allocate(dq_xc1(m1:n1))
    case ('xc2')
          allocate(dq_xc2(m2:n2))
    case ('xc3')
          allocate(dq_xc3(m3:n3))
    case ('Bb1')
          allocate(dq_Bb1(m1:n1, m2:n2, m3:n3))
    case ('Bb2')
          allocate(dq_Bb2(m1:n1, m2:n2, m3:n3))
    case ('Bb3')
          allocate(dq_Bb3(m1:n1, m2:n2, m3:n3))
  end select
endif
!if(present(dimension)) then
!  dq_m1=dimension(1,1)
!  dq_n1=dimension(2,1)
!  dq_m2=dimension(1,2)
!  dq_n2=dimension(2,2)
!  dq_m3=dimension(1,3)
!  dq_n3=dimension(2,3)
!endif
end subroutine rhd_dq_vars_init
!
!----------**************--------------------------------------------------------------------------
subroutine rhd_dq_end_init()
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_end_init ('rhd_derived_quantities_end_initialization')
!
! PURPOSE:
!   End initialization; change magnetic field units, define dq_deltax,...
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_dq_end_init(parfile, file, nc_p, imodel=imodel)
!
! INPUT:
!   !parfile:        ('parameter_file') character, name of parameter file
!   !file:           ('file') character, name of FULL file
!   !nc_p:           ('number_channel_print') integer, logical number for text output
!   !imodel:         ('time_index') integer, time index of requested model
!
! OUTPUT:
!
! ROUTINES:
!   !rhd_eosopa_RdData:     Read EOS and OPA files
!   !rhd_box_Read:          Read box from file in UIO form
!   !rhd_keep_only:         Remove unused variables
!
! MODULES:
!   !const_module:          Global physical and mathematical constants, units
!   !rhd_gl_module:         Global type definitions and parameters for RHD
!   !rhd_box_module:        Box-handling routines
!   !rhd_sub_module:        Collection of additional subroutines
!   !rhd_action_module:     Routines to handle the control structure 'action'
!
! SIDE EFFECTS:
!   !Open file, read data from it, close file.
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Compute dq_deltax, dq_deltay, dq_deltaz and change magnetic field units
!   !If not present, determine imodel. !Read full file, allocate memory, copy needed
!   !variables, deallocate memory.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-08-17 (F.C. Locarno) Written
!--------------------------------------------------------------------------------------------------
use const_module
!use rhd_gl_module
!use rhd_box_module
!use rhd_sub_module
!use rhd_action_module
!
implicit none
!
! --- I/O parameters ---
!character(len=*),        intent(in)     :: parfile
!character(len=*),        intent(in)     :: file
!integer,                 intent(in)     :: nc_p
!integer, optional,       intent(in)     :: imodel
!
! --- Local variables ---
!integer                                 :: ncin, nmodel, imodel0
!type(box_type), pointer                 :: model
!real, allocatable, dimension(:,:,:)     :: dTde
!logical                                 :: Bb_flag
real, dimension(1)                      :: tmp_array1
!--------------------------------------------------------------------------------------------------
! === End initialization ==========================================================================
if (dq_status.eq.0) then
  dq_status=1
  !!! Set dq_mi and dq_ni
  dq_nx = dq_n1-dq_m1+1
  dq_ny = dq_n2-dq_m2+1
  dq_nz = dq_n3-dq_m3+1
  !!! Set dq_time
  !tmp_array1 = (dq_xb1(dq_n1+1)-dq_xb1(dq_m1))/dq_nx
  !dq_deltax = tmp_array1(1)
  !tmp_array1 = (dq_xb2(dq_n2+1)-dq_xb2(dq_m2))/dq_ny
  !dq_deltay = tmp_array1(1)
  !tmp_array1 = (dq_xb3(dq_n3+1)-dq_xb3(dq_m3))/dq_nz
  !dq_deltaz = tmp_array1(1)
  if(dq_Bb_flag) then
    write(dq_nc_p,*) 'Bb_flag present.'
    dq_Bb1 = dq_Bb1*sqrt(fourpi_drk)
    dq_Bb2 = dq_Bb2*sqrt(fourpi_drk)
    dq_Bb3 = dq_Bb3*sqrt(fourpi_drk)
  else
    allocate(dq_Bb1(dq_m1:dq_n1+1, dq_m2:dq_n2, dq_m3:dq_n3))
    allocate(dq_Bb2(dq_m1:dq_n1, dq_m2:dq_n2+1, dq_m3:dq_n3))
    allocate(dq_Bb3(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3+1))
    dq_Bb1 = 0.0
    dq_Bb2 = 0.0
    dq_Bb3 = 0.0
  endif
endif
end subroutine rhd_dq_end_init
!--------------------------------------------------------------------------------------------------
#else
subroutine rhd_dq_init(parfile, file, nc_p, imodel)
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_init ('rhd_derived_quantities_initialization')
!
! PURPOSE:
!   Use FULL file to extract all relevant hydrodynamic variables
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_dq_init(parfile, file, nc_p, imodel=imodel)
!
! INPUT:
!   parfile:        ('parameter_file') character, name of parameter file
!   file:           ('file') character, name of FULL file
!   nc_p:           ('number_channel_print') integer, logical number for text output
!   imodel:         ('time_index') integer, time index of requested model
!
! OUTPUT:
!
! ROUTINES:
!   rhd_eosopa_RdData:     Read EOS and OPA files
!   rhd_box_Read:          Read box from file in UIO form
!   rhd_keep_only:         Remove unused variables
!
! MODULES:
!   const_module:          Global physical and mathematical constants, units
!   rhd_gl_module:         Global type definitions and parameters for RHD
!   rhd_box_module:        Box-handling routines
!   rhd_sub_module:        Collection of additional subroutines
!   rhd_action_module:     Routines to handle the control structure 'action'
!
! SIDE EFFECTS:
!   Open file, read data from it, close file.
!
! RESTRICTIONS:
!
! PROCEDURE:
!   If not present, determine imodel. Read full file, allocate memory, copy needed
!   variables, deallocate memory.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-02-01 (F.C. Freiburg) Written
!--------------------------------------------------------------------------------------------------
use const_module
use rhd_gl_module
use rhd_box_module
use rhd_sub_module
use rhd_action_module
!
implicit none
!
! --- I/O parameters ---
character(len=*),        intent(in)     :: parfile
character(len=*),        intent(in)     :: file
integer,                 intent(in)     :: nc_p
integer, optional,       intent(in)     :: imodel
!
! --- Local variables ---
integer                                 :: ncin, nmodel, imodel0
type(box_type), pointer                 :: model
real, allocatable, dimension(:,:,:)     :: dTde
logical                                 :: Bb_flag
real, dimension(1)                      :: tmp_array1
!--------------------------------------------------------------------------------------------------
! === Extract model ===============================================================================
!
! Reload parfile and data file if those have changed
!
if(dq_status.ne.0) then
  if((trim(parfile) /= trim(dq_parfile)).or.(trim(file) /= trim(dq_file))) then
    call rhd_keep_only(base=.false.)
  endif
endif
if(dq_status.eq.0) then
  call rhd_eosopa_RdData(trim(parfile), nc_p)
  dq_C_radHtautop = par%C_radHtautop
  dq_grav = par%grav
  dq_Bb_flag = ((par%A%hdscheme == hdscheme_RoeMHD) .or. &
                              (par%A%hdscheme == hdscheme_HLLMHD) .or. &
                              (par%A%hdscheme == hdscheme_CenMHD))
endif
!
! We allow for negative values for imodel, which are interpretd as counting
! backwards from last available model. Ex.: imodel=-1 means last available
! model.
!
if (present(imodel)) then
  imodel0 = imodel
else
  if(dq_status.eq.0) then
    imodel0 = par%istep_in_start
  else
    imodel0 = dq_imodel
  endif
endif
call rhd_box_Read('open', ncin, nmodel=nmodel, file=trim(file))
if (imodel0 < 0) then
  imodel0 = nmodel+imodel0+1
  if(dq_status.eq.0) then
    dq_imodel = imodel0
  endif
endif
if ((imodel0.ne.dq_imodel).and.(dq_imodel.ge.0)) then
  call rhd_keep_only(base=.false.)
endif
dq_imodel=imodel0
if (dq_status.eq.0) then
  dq_status=1
  dq_parfile=parfile
  dq_file=file
  dq_nc_p=nc_p
  allocate(model)
  call rhd_box_Init(model)
  call rhd_box_Read('model,close', ncin, &
                    box=model, Bb_flag=dq_Bb_flag, file=trim(file), imodel=imodel0)
  dq_m1 = model%m1
  dq_n1 = model%n1
  dq_nx = dq_n1-dq_m1+1
  dq_m2 = model%m2
  dq_n2 = model%n2
  dq_ny = dq_n2-dq_m2+1
  dq_m3 = model%m3
  dq_n3 = model%n3
  dq_nz = dq_n3-dq_m3+1
  dq_time = model%time
  dq_itime = model%itime
  dq_nmodel = nmodel
  allocate(dq_ei(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_ei = model%ei
  allocate(dq_rho(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_rho = model%rho
  allocate(dq_v1(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_v1 = model%v1
  allocate(dq_v2(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_v2 = model%v2
  allocate(dq_v3(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_v3 = model%v3
  allocate(dq_xb1(dq_m1:dq_m1+1))
  dq_xb1 = model%xb1
  tmp_array1 = (dq_xb1(dq_n1+1)-dq_xb1(dq_m1))/dq_nx
  dq_deltax = tmp_array1(1)
  allocate(dq_xb2(dq_m2:dq_m2+1))
  dq_xb2 = model%xb2
  tmp_array1 = (dq_xb2(dq_n2+1)-dq_xb2(dq_m2))/dq_ny
  dq_deltay = tmp_array1(1)
  allocate(dq_xb3(dq_m3:dq_m3+1))
  dq_xb3 = model%xb3
  tmp_array1 = (dq_xb3(dq_n3+1)-dq_xb3(dq_m3))/dq_nz
  dq_deltaz = tmp_array1(1)
  allocate(dq_xc1(dq_m1:dq_n1))
  dq_xc1 = model%xc1
  allocate(dq_xc2(dq_m2:dq_n2))
  dq_xc2 = model%xc2
  allocate(dq_xc3(dq_m3:dq_n3))
  dq_xc3 = model%xc3
  allocate(dq_Bb1(dq_m1:dq_n1+1, dq_m2:dq_n2, dq_m3:dq_n3))
  allocate(dq_Bb2(dq_m1:dq_n1, dq_m2:dq_n2+1, dq_m3:dq_n3))
  allocate(dq_Bb3(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3+1))
  if(dq_Bb_flag) then
    write(dq_nc_p,*) 'Bb_flag present.'
    dq_Bb1 = model%Bb1*sqrt(fourpi_drk)
    dq_Bb2 = model%Bb2*sqrt(fourpi_drk)
    dq_Bb3 = model%Bb3*sqrt(fourpi_drk)
  else
    dq_Bb1 = 0.0
    dq_Bb2 = 0.0
    dq_Bb3 = 0.0
  endif
  call rhd_box_Delete(model)
  deallocate(model)
else
  call rhd_box_Read('close', ncin, file=trim(file))
endif
end subroutine rhd_dq_init
#endif
!
!----------**************--------------------------------------------------------------------------
subroutine rhd_dq_vars(imodel, resample, &
                       rho, ei, v1, v2, v3, xb1, xb2, xb3, xc1, xc2, xc3, &
                       T, P, s, j1, j2, j3, jabs, kappa, Bb1, Bb2, Bb3, &
                       Bc1, Bc2, Bc3, divB, vabs, vh, ekin, plin, vmflux, &
                       g1, g3, cs, mach, mmu, tau, Babs, Bh, B2, emag, cA, &
                       beta, csca, nx, ny, nz, deltax, deltay, deltaz)
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_vars ('rhd_derived_quantities_variables')
!
! PURPOSE:
!   Use FULL file to extract all relevant hydrodynamic variables
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_dq_vars(imodel, resample, &
!                      rho=rho, ei=ei, v1=v1, v2=v2, v3=v3, &
!                      xb1=xb1, xb2=xb2, xb3=xb3, xc1=xc1, xc2=xc2, xc3=xc3, &
!                      T=T, P=P, s=s, j1=j1, j2=j2, j3=j3, jabs=jabs, &
!                      kappa=kappa, Bc1=Bc2, Bc2=Bc2, Bc3=Bc3, divB=divB, &
!                      vabs=vabs, vh=vh, ekin=ekin, plin=plin, vmflux=vmflux, &
!                      g1=g1, g3=g3, cs=cs, mach=mach, mmu=mmu, &
!                      tau=tau, Babs=Babs, Bh=Bh, B2=B2, emag=emag, cA=cA, &
!                      beta=beta, csca=csca, nx=nx, ny=ny, nz=nz, &
!                      deltax=deltax, deltay=deltay, deltaz=deltaz)
!
! INPUT:
!   imodel:         ('time_index') integer, time index of requested model
!   resample:       ('resample') logical, logical number for resampling
!
! OUTPUT:
!   rho:      ('rho') real, dimension(:,:,:), density
!   ei:       ('ei') real, dimension(:,:,:), energy
!   v1:       ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!   v2:       ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!   v3:       ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!   xb1:      ('x_boundary_1') real, 1. coordinate of cell boundaries
!   xb2:      ('x_boundary_2') real, 2. coordinate of cell boundaries
!   xb3:      ('x_boundary_3') real, 3. coordinate of cell boundaries
!   xc1:      ('x_center_1') real, 1. coordinate of cell centers
!   xc2:      ('x_center_2') real, 2. coordinate of cell centers
!   xc3:      ('x_center_3') real, 3. coordinate of cell centers
!   T:        ('T') real, dimension(:,:,:), temperature
!   P:        ('P') real, dimension(:,:,:), pressure
!   s:        ('P') real, dimension(:,:,:), entropy
!   j1:       ('j_center_1') real, dimension(:,:,:), cell center current density, 1. component
!   j2:       ('j_center_2') real, dimension(:,:,:), cell center current density, 2. component
!   j3:       ('j_center_3') real, dimension(:,:,:), cell center current density, 3. component
!   jabs:     ('j_abs') real, dimension(:,:,:), cell center absolute current density
!   kappa:    ('kappa') real, dimension(:,:,:), cell center opacity
!   Bb1:      ('B_boundary_1') real, dimension(:,:,:), cell boundary magnetic field, 1. component
!   Bb2:      ('B_boundary_2') real, dimension(:,:,:), cell boundary magnetic field, 2. component
!   Bb3:      ('B_boundary_3') real, dimension(:,:,:), cell boundary magnetic field, 3. component
!   Bc1:      ('B_center_1') real, dimension(:,:,:), cell center magnetic field, 1. component
!   Bc2:      ('B_center_2') real, dimension(:,:,:), cell center magnetic field, 2. component
!   Bc3:      ('B_center_3') real, dimension(:,:,:), cell center magnetic field, 3. component
!   divB:     ('div_B') real, dimension(:,:,:), cell center magnetic field divergence
!   vabs:     ('v_absolute') real, dimension(:,:,:), cell center speed
!   vh:       ('v_horizontal') real, dimension(:,:,:), cell center horizontal velocity
!   ekin:     ('kinetic_energy') real, dimension(:,:,:), cell center kinetic energy
!   plin:     ('linear_momentum') real, dimension(:,:,:), cell center linear momentum
!   vmflux    ('vertical_mass_flux') real, dimension(:,:,:), cell center vertical mass flux
!   g1:       ('gamma_1') real, dimension(:,:,:), cell center first adiabatic coefficient
!   g3:       ('gamma_3') real, dimension(:,:,:), cell center third adiabatic coefficient
!   cs:       ('sound_speed') real, dimension(:,:,:), cell center sound speed
!   mach:     ('mach_number') real, dimension(:,:,:), cell center mach number
!   mmu:      ('mean_molecular_weight') real, dimension(:,:,:), cell center mean molecular weight
!   tau:      ('optical_depth') real, dimension(:,:,:), cell center optical depth
!   Babs:     ('B_absolute') real, dimension(:,:,:), cell center absolute magnetic field
!   Bh:       ('B_horizontal') real, dimension(:,:,:), cell center horizontal magnetic field
!   B2:       ('B_signed') real, dimension(:,:,:), cell center signed absolute magnetic field
!   emag:     ('magnetic_energy') real, dimension(:,:,:), cell center magnetic energy
!   cA:       ('alfven_speed') real, dimension(:,:,:), cell center alfven speed
!   beta:     ('plasma_beta') real, dimension(:,:,:), cell center plasma beta
!   csca:     ('cs_ca') real, dimension(:,:,:), cell center sound speed to Alfven speed ratio
!   nx        ('number_x') integer, 1. dimension, size
!   ny        ('number_y') integer, 2. dimension, size
!   nz        ('number_z') integer, 3. dimension, size
!   deltax    ('delta_x') real, 1. dimension, pixel mean physical length
!   deltay    ('delta_y') real, 2. dimension, pixel mean physical length
!   deltaz    ('delta_z') real, 3. dimension, pixel mean physical length
!
! ROUTINES:
!   rhd_dq_init:           Extract relevant hydrodynamic variables from FULL file
!   eosinter3:             Evaluate equation of state (via interpolation in table)
!   rhd_int_tau:           Integrate opacities to get optical depth
!   resample1D:            Resample a 1D quantity
!   resample3D:            Resample a 3D quantity
!   diff3D:                Differentiate cell-centered quantity
!
! MODULES:
!   gasinter_module:       Routines for interpolating GAS or EOS quantities
!   const_module:          Global physical and mathematical constants, units
!   rhd_box_module:        Box-handling routines
!   opta_module:           Opacity interpolation routines
!
! SIDE EFFECTS:
!   Allocate memory
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Sort dependencies of requested quantities, compute them and allocate
!   memory to store them.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-02-01 (F.C. Freiburg) Written
!--------------------------------------------------------------------------------------------------
use gasinter_module
use const_module
use rhd_box_module
use opta_module
#ifdef _OPENMP
use omp_lib
#else
#define omp_get_thread_num() 0
#endif
!
implicit none
!
! --- I/O parameters ---
integer, optional,       intent(in)     :: imodel
logical, optional,       intent(in)     :: resample
integer, optional,       intent(out)    :: nx, ny, nz
real, optional,          intent(out)    :: deltax, deltay, deltaz
!
! Fundamental quantities
real, pointer, dimension(:,:,:), &
               optional, intent(out)    :: rho, ei, v1, v2, v3, Bb1, Bb2, Bb3
real, pointer, dimension(:), &
               optional, intent(out)    :: xb1, xb2, xb3, xc1, xc2, xc3
! Derived quantities
real, pointer, dimension(:,:,:), &
               optional, intent(out)    :: T, P, s, j1, j2, j3, kappa, &
                                           Bc1, Bc2, Bc3, divB
! Additional quantities
real, pointer, dimension(:,:,:), &
               optional, intent(out)    :: vabs, vh, ekin, plin, vmflux, g1, g3, &
                                           cs, mach, mmu, tau, Bh, Babs, B2, &
                                           emag, jabs, cA, beta, csca
!
! --- Local variables ---
integer                                 :: ncin, nmodel, imodel0, &
                                           i1, i2, i3, &
                                           needP=0, needT=0, needdP=0, &
                                           needkappa=0, needcs=0, needca=0, &
                                           needj=0, needB=0, numthreads
real, allocatable, dimension(:,:,:)     :: dTde
logical                                 :: Bb_flag, resample0=.false.
real, dimension(1)                      :: tmp_array1
!--------------------------------------------------------------------------------------------------
#ifndef PYBOLD
if(present(imodel)) then
  call rhd_dq_init(dq_parfile, dq_file, dq_nc_p, imodel)
else
  call rhd_dq_init(dq_parfile, dq_file, dq_nc_p)
endif
#endif
if (present(resample)) then
  resample0=resample
endif
! === Sort dependences ============================================================================
if(present(mmu).or.present(beta).or. &
   present(G1).or.present(kappa).or. &
   present(tau).or.present(cs).or. &
   present(csca).or.present(mach)) then
  needP=1
endif
if(present(mmu).or.present(kappa).or.present(tau)) then
  needT=1
endif
if(present(G1).or.present(G3).or. &
   present(cs).or.present(csca).or. &
   present(mach)) then
  needdP=1
endif
if(present(jabs)) then
  needj=1
endif
if(present(Bc1).or.present(Bc2).or. &
   present(Bc3).or.present(Bh).or. &
   present(Babs).or.present(B2).or. &
   present(ca).or.present(csca).or. &
   present(j1).or.present(j2).or.present(j3).or. &
   present(jabs).or.present(beta).or.present(divB)) then
  needB=1
endif
if(present(tau)) then
  needkappa=1
endif
if(present(mach)) then
  needcs=1
endif
! === Compute requested quantities ================================================================
#ifdef _OPENMP
numthreads = omp_get_num_threads()
#endif
#ifdef PYBOLD
#ifdef _OPENMP
call omp_set_num_threads(1)
#endif
#endif
!$OMP PARALLEL
!$OMP BARRIER
!$OMP SECTIONS
!$OMP SECTION
if((present(P).or.needP.eq.1).and..not.allocated(dq_P)) then
  allocate(dq_P(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  allocate(dq_dPdei(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  allocate(dq_dPdrho(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  call eosinter3(dq_rho, dq_ei, dq_P, dq_dPdrho, dq_dPdei)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'P         [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(T).or.needT.eq.1).and..not.allocated(dq_T)) then
  allocate(dq_T(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  allocate(dTde(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  call eosinter3(dq_rho, dq_ei, T=dq_T, dTde=dTde)
  deallocate(dTde)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'T         [ok]'
!OMP END CRITICAL(DQWRITE)
endif
if((present(s)).and..not.allocated(dq_s)) then
  allocate(dq_s(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  call eosinter3(dq_rho, dq_ei, s=dq_s)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 's         [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if(needB.eq.1.and..not.allocated(dq_Bc1)) then
  allocate(dq_Bc1(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  allocate(dq_Bc2(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  allocate(dq_Bc3(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_Bc1(:,:,:) = 0.5*(dq_Bb1(dq_m1:dq_n1,:,:)+dq_Bb1(dq_m1+1:dq_n1+1,:,:))
  dq_Bc2(:,:,:) = 0.5*(dq_Bb2(:,dq_m2:dq_n2,:)+dq_Bb2(:,dq_m2+1:dq_n2+1,:))
  dq_Bc3(:,:,:) = 0.5*(dq_Bb3(:,:,dq_m3:dq_n3)+dq_Bb3(:,:,dq_m3+1:dq_n3+1))
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'B         [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(divB)).and..not.allocated(dq_divB)) then
  allocate(dq_divB(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  do i3=dq_m3,dq_n3
    do i2=dq_m2,dq_n2
      do i1=dq_m1,dq_n1
        dq_divB(i1,i2,i3) &
          = ( (-dq_Bb1(i1,i2,i3)+dq_Bb1(i1+1,i2,i3))*(dq_xb2(i2+1)-dq_xb2(i2))*(dq_xb3(i3+1)-dq_xb3(i3)) &
           +(-dq_Bb2(i1,i2,i3)+dq_Bb2(i1,i2+1,i3))*(dq_xb3(i3+1)-dq_xb3(i3))*(dq_xb1(i1+1)-dq_xb1(i1)) &
           +(-dq_Bb3(i1,i2,i3)+dq_Bb3(i1,i2,i3+1))*(dq_xb1(i1+1)-dq_xb1(i1))*(dq_xb2(i2+1)-dq_xb2(i2)) ) &
           / ((dq_xb2(i2+1)-dq_xb2(i2))*(dq_xb3(i3+1)-dq_xb3(i3))*(dq_xb1(i1+1)-dq_xb1(i1)))
      end do !i1
    end do !i2
  end do !i3
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'divB      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP END SECTIONS
!$OMP BARRIER
!$OMP END PARALLEL
#ifdef PYBOLD
#ifdef _OPENMP
call omp_set_num_threads(numthreads)
#endif
#endif
if(((present(kappa)).or.needkappa.eq.1).and..not.allocated(dq_kappa)) then
  allocate(dq_kappa(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
#ifdef XKAROS_INTERFACE
  dq_kappa=xkaros(dq_P, dq_T, 0)
#else
  call sub_xkaros_logPT2kappa(log(dq_P), dq_T, 0, dq_kappa)
#endif
  write(dq_nc_p,*) 'kappa     [ok]'
endif
#ifdef PYBOLD
#ifdef _OPENMP
call omp_set_num_threads(1)
#endif
#endif
!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
if(((present(j1)).or.needj.eq.1).and..not.allocated(dq_j1)) then
  allocate(dq_j1(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_j1=c_light_drk/fourpi_drk*(diff3D(dq_Bc3, dq_xc2, dq_xb2,2)-diff3D(dq_Bc2, dq_xc3, dq_xb3,3))
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'j1        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if(((present(j2)).or.needj.eq.1).and..not.allocated(dq_j2)) then
  allocate(dq_j2(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_j2=c_light_drk/fourpi_drk*(diff3D(dq_Bc1, dq_xc3, dq_xb3,3)-diff3D(dq_Bc3, dq_xc1, dq_xb1,1))
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'j2        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if(((present(j3)).or.needj.eq.1).and..not.allocated(dq_j3)) then
  allocate(dq_j3(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_j3=c_light_drk/fourpi_drk*(diff3D(dq_Bc2, dq_xc1, dq_xb1,1)-diff3D(dq_Bc1, dq_xc2, dq_xb2,2))
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'j3        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(cs)).and..not.allocated(dq_cs)) then
  allocate(dq_cs(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_cs=sqrt(dq_P/dq_rho*dq_dPdei/dq_rho+dq_dPdrho)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'cs        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP END SECTIONS
!$OMP SECTIONS
!$OMP SECTION
if((present(vabs)).and..not.allocated(dq_vabs)) then
  allocate(dq_vabs(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_vabs=sqrt(dq_v1**2+dq_v2**2+dq_v3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'vabs      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(vh)).and..not.allocated(dq_vh)) then
  allocate(dq_vh(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_vh=sqrt(dq_v1**2+dq_v2**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'vh        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(ekin)).and..not.allocated(dq_ekin)) then
  allocate(dq_ekin(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_ekin=.5*dq_rho*(dq_v1**2+dq_v2**2+dq_v3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'ekin      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(plin)).and..not.allocated(dq_plin)) then
  allocate(dq_plin(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_plin=dq_rho*sqrt(dq_v1**2+dq_v2**2+dq_v3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'plin      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(vmflux)).and..not.allocated(dq_vmflux)) then
  allocate(dq_vmflux(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_vmflux=dq_rho*dq_v3
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'vmflux    [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(g1)).and..not.allocated(dq_g1)) then
  allocate(dq_g1(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_g1=dq_dPdrho*dq_rho/dq_P+dq_dPdei/dq_rho
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'g1        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(g3)).and..not.allocated(dq_g3)) then
  allocate(dq_g3(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_g3=dq_dPdei/dq_rho+1.0
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'g3        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(mach)).and..not.allocated(dq_mach)) then
  allocate(dq_mach(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_mach=sqrt(dq_v1**2+dq_v2**2+dq_v3**2) / &
          sqrt(dq_P/dq_rho*dq_dPdei/dq_rho+dq_dPdrho)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'mach      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(mmu)).and..not.allocated(dq_mmu)) then
  allocate(dq_mmu(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_mmu=8.314e7*dq_rho*dq_T/dq_P   ! R = 8.314e7 [erg/g/K]
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'mmu       [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(tau)).and..not.allocated(dq_tau)) then
  allocate(dq_tau(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  call rhd_int_tau()
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'tau       [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(Bh)).and..not.allocated(dq_Bh)) then
  allocate(dq_Bh(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_Bh=sqrt(dq_Bc1**2+dq_Bc2**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'Bh        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(Babs)).and..not.allocated(dq_Babs)) then
  allocate(dq_Babs(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_Babs=sqrt(dq_Bc1**2+dq_Bc2**2+dq_Bc3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'Babs      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(B2)).and..not.allocated(dq_B2)) then
  allocate(dq_B2(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_B2=dq_Bc1*abs(dq_Bc1)+dq_Bc2*abs(dq_Bc2)+dq_Bc3*abs(dq_Bc3)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'B2        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(emag)).and..not.allocated(dq_emag)) then
  allocate(dq_emag(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_emag=(dq_Bc1**2+dq_Bc2**2+dq_Bc3**2)/(2.0*fourpi_drk)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'emag      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(jabs)).and..not.allocated(dq_jabs)) then
  allocate(dq_jabs(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_jabs=sqrt(dq_j1**2+dq_j2**2+dq_j3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'jabs      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(ca)).and..not.allocated(dq_ca)) then
  allocate(dq_ca(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_ca=sqrt(dq_Bc1**2+dq_Bc2**2+dq_Bc3**2)/sqrt(fourpi_drk*dq_rho)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'ca        [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(beta)).and..not.allocated(dq_beta)) then
  allocate(dq_beta(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_beta=2.0*fourpi_drk*dq_P/(dq_Bc1**2+dq_Bc2**2+dq_Bc3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'beta      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP SECTION
if((present(csca)).and..not.allocated(dq_csca)) then
  allocate(dq_csca(dq_m1:dq_n1, dq_m2:dq_n2, dq_m3:dq_n3))
  dq_csca=sqrt(dq_P/dq_rho*dq_dPdei/dq_rho+dq_dPdrho)*sqrt(fourpi_drk*dq_rho)/sqrt(dq_Bc1**2+dq_Bc2**2+dq_Bc3**2)
!OMP CRITICAL(DQWRITE)
  write(dq_nc_p,*) 'csca      [ok]'
!OMP END CRITICAL(DQWRITE)
endif
!$OMP END SECTIONS
!$OMP BARRIER
!$OMP END PARALLEL
#ifdef PYBOLD
#ifdef _OPENMP
call omp_set_num_threads(numthreads)
#endif
#endif
! === Set pointers and variables to requested quantities ==========================================
if (present(rho)) then
  if (resample0) then
    call resample3D(dq_rho, rho)
  else
    rho => dq_rho
  endif
endif
if (present(ei)) then
  if (resample0) then
    call resample3D(dq_ei, ei)
  else
    ei => dq_ei
  endif
endif
if (present(v1)) then
  if (resample0) then
    call resample3D(dq_v1, v1)
  else
    v1 => dq_v1
  endif
endif
if (present(v2)) then
  if (resample0) then
    call resample3D(dq_v2, v2)
  else
    v2 => dq_v2
  endif
endif
if (present(v3)) then
  if (resample0) then
    call resample3D(dq_v3, v3)
  else
    v3 => dq_v3
  endif
endif
if (present(Bc1)) then
  if (resample0) then
    call resample3D(dq_Bc1, Bc1)
  else
    Bc1 => dq_Bc1
  endif
endif
if (present(Bc2)) then
  if (resample0) then
    call resample3D(dq_Bc2, Bc2)
  else
    Bc2 => dq_Bc2
  endif
endif
if (present(Bc3)) then
  if (resample0) then
    call resample3D(dq_Bc3, Bc3)
  else
    Bc3 => dq_Bc3
  endif
endif
if (present(T)) then
  if (resample0) then
    call resample3D(dq_T, T)
  else
    T => dq_T
  endif
endif
if (present(P)) then
  if (resample0) then
    call resample3D(dq_P, P)
  else
    P => dq_P
  endif
endif
if (present(s)) then
  if (resample0) then
    call resample3D(dq_s, s)
  else
    s => dq_s
  endif
endif
if (present(j1)) then
  if (resample0) then
    call resample3D(dq_j1, j1)
  else
    j1 => dq_j1
  endif
endif
if (present(j2)) then
  if (resample0) then
    call resample3D(dq_j2, j2)
  else
    j2 => dq_j2
  endif
endif
if (present(j3)) then
  if (resample0) then
    call resample3D(dq_j3, j3)
  else
    j3 => dq_j3
  endif
endif
if (present(jabs)) then
  if (resample0) then
    call resample3D(dq_jabs, jabs)
  else
    jabs => dq_jabs
  endif
endif
if (present(kappa)) then
  if (resample0) then
    call resample3D(dq_kappa, kappa)
  else
    kappa => dq_kappa
  endif
endif
if (present(vabs)) then
  if (resample0) then
    call resample3D(dq_vabs, vabs)
  else
    vabs => dq_vabs
  endif
endif
if (present(vh)) then
  if (resample0) then
    call resample3D(dq_vh, vh)
  else
    vh => dq_vh
  endif
endif
if (present(ekin)) then
  if (resample0) then
    call resample3D(dq_ekin, ekin)
  else
    ekin => dq_ekin
  endif
endif
if (present(plin)) then
  if (resample0) then
    call resample3D(dq_plin, plin)
  else
    plin => dq_plin
  endif
endif
if (present(vmflux)) then
  if (resample0) then
    call resample3D(dq_vmflux, vmflux)
  else
    vmflux => dq_vmflux
  endif
endif
if (present(g1)) then
  if (resample0) then
    call resample3D(dq_g1, g1)
  else
    g1 => dq_g1
  endif
endif
if (present(g3)) then
  if (resample0) then
    call resample3D(dq_g3, g3)
  else
    g3 => dq_g3
  endif
endif
if (present(cs)) then
  if (resample0) then
    call resample3D(dq_cs, cs)
  else
    cs => dq_cs
  endif
endif
if (present(mach)) then
  if (resample0) then
    call resample3D(dq_mach, mach)
  else
    mach => dq_mach
  endif
endif
if (present(mmu)) then
  if (resample0) then
    call resample3D(dq_mmu, mmu)
  else
    mmu => dq_mmu
  endif
endif
if (present(tau)) then
  if (resample0) then
    call resample3D(dq_tau, tau)
  else
    tau => dq_tau
  endif
endif
if (present(Babs)) then
  if (resample0) then
    call resample3D(dq_Babs, Babs)
  else
    Babs => dq_Babs
  endif
endif
if (present(Bh)) then
  if (resample0) then
    call resample3D(dq_Bh, Bh)
  else
    Bh => dq_Bh
  endif
endif
if (present(B2)) then
  if (resample0) then
    call resample3D(dq_B2, B2)
  else
    B2 => dq_B2
  endif
endif
if (present(emag)) then
  if (resample0) then
    call resample3D(dq_emag, emag)
  else
    emag => dq_emag
  endif
endif
if (present(cA)) then
  if (resample0) then
    call resample3D(dq_cA, cA)
  else
    cA => dq_cA
  endif
endif
if (present(beta)) then
  if (resample0) then
    call resample3D(dq_beta, beta)
  else
    beta => dq_beta
  endif
endif
if (present(csca)) then
  if (resample0) then
    call resample3D(dq_csca, csca)
  else
    csca => dq_csca
  endif
endif
if (present(xb1)) then
  if (resample0) then
    call resample1D(dq_xb1, xb1)
  else
    xb1 => dq_xb1
  endif
endif
if (present(xb2)) then
  if (resample0) then
    call resample1D(dq_xb2, xb2)
  else
    xb2 => dq_xb2
  endif
endif
if (present(xb3)) then
  if (resample0) then
    call resample1D(dq_xb3, xb3)
  else
    xb3 => dq_xb3
  endif
endif
if (present(xc1)) then
  if (resample0) then
    call resample1D(dq_xc1, xc1)
  else
    xc1 => dq_xc1
  endif
endif
if (present(xc2)) then
  if (resample0) then
    call resample1D(dq_xc2, xc2)
  else
    xc2 => dq_xc2
  endif
endif
if (present(xc3)) then
  if (resample0) then
    call resample1D(dq_xc3, xc3)
  else
    xc3 => dq_xc3
  endif
endif
if (present(nx)) then
  nx = dq_nx
endif
if (present(ny)) then
  ny = dq_ny
endif
if (present(nz)) then
  nz = dq_nz
endif
if (present(deltax)) then
  deltax = dq_deltax
endif
if (present(deltay)) then
  deltay = dq_deltay
endif
if (present(deltaz)) then
  deltaz = dq_deltaz
endif
end subroutine rhd_dq_vars
!
!----------**************--------------------------------------------------------------------------
subroutine rhd_keep_only(base, rho, ei, v1, v2, v3, xb1, xb2, xb3, xc1, xc2, xc3, &
                         T, P, s, j1, j2, j3, jabs, kappa, Bc1, Bc2, Bc3, divB, &
                         vabs, vh, ekin, plin, vmflux, g1, g3, cs, mach, mmu, &
                         tau, Babs, Bh, b2, emag, cA, beta, csca, Bb1, Bb2, Bb3)
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_keep_only ('rhd_keep_only')
!
! PURPOSE:
!   Deallocate variables
!   
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call rhd_keep_only(base=base, rho=rho, ei=ei, v1=v1, v2=v2, v3=v3, &
!                      xb1=xb1, xb2=xb2, xb3=xb3, xc1=xc1, xc2=xc2, xc3=xc3, &
!                      T=T, P=P, s=s, j1=j1, j2=j2, j3=j3, jabs=jabs, &
!                      kappa=kappa, Bc1=Bc2, Bc2=Bc2, Bc3=Bc3, vabs=vabs, &
!                      vh=vh, ekin=ekin, plin=plin, vmflux=vmflux, &
!                      g1=g1, g3=g3, cs=cs, mach=mach, mmu=mmu, &
!                      tau=tau, Babs=Babs, Bh=Bh, B2=B2, emag=emag, cA=cA, &
!                      beta=beta, csca=csca)
!
! INPUT:
!   base:     ('base') logical, keep all base variables (ei, rho, vi, Bbi, xbi, xci)
!   rho:      ('rho') real, dimension(:,:,:), density
!   ei:       ('ei') real, dimension(:,:,:), energy
!   v1:       ('v_1') real, dimension(:,:,:), velocity, 1. coordinate
!   v2:       ('v_2') real, dimension(:,:,:), velocity, 2. coordinate
!   v3:       ('v_3') real, dimension(:,:,:), velocity, 3. coordinate
!   xb1:      ('x_boundary_1') real, 1. coordinate of cell boundaries
!   xb2:      ('x_boundary_2') real, 2. coordinate of cell boundaries
!   xb3:      ('x_boundary_3') real, 3. coordinate of cell boundaries
!   xc1:      ('x_center_1') real, 1. coordinate of cell centers
!   xc2:      ('x_center_2') real, 2. coordinate of cell centers
!   xc3:      ('x_center_3') real, 3. coordinate of cell centers
!   T:        ('T') real, dimension(:,:,:), temperature
!   P:        ('P') real, dimension(:,:,:), pressure
!   s:        ('P') real, dimension(:,:,:), entropy
!   j1:       ('j_center_1') real, dimension(:,:,:), cell center current density, 1. component
!   j2:       ('j_center_2') real, dimension(:,:,:), cell center current density, 2. component
!   j3:       ('j_center_3') real, dimension(:,:,:), cell center current density, 3. component
!   jabs:     ('j_abs') real, dimension(:,:,:), cell center absolute current density
!   kappa:    ('kappa') real, dimension(:,:,:), cell center opacity
!   Bc1:      ('B_center_1') real, dimension(:,:,:), cell center magnetic field, 1. component
!   Bc2:      ('B_center_2') real, dimension(:,:,:), cell center magnetic field, 2. component
!   Bc3:      ('B_center_3') real, dimension(:,:,:), cell center magnetic field, 3. component
!   divB:     ('div_B') real, dimension(:,:,:), cell center magnetic field divergence
!   vabs:     ('v_absolute') real, dimension(:,:,:), cell center speed
!   vh:       ('v_horizontal') real, dimension(:,:,:), cell center horizontal velocity
!   ekin:     ('kinetic_energy') real, dimension(:,:,:), cell center kinetic energy
!   plin:     ('linear_momentum') real, dimension(:,:,:), cell center linear momentum
!   vmflux    ('vertical_mass_flux') real, dimension(:,:,:), cell center vertical mass flux
!   g1:       ('gamma_1') real, dimension(:,:,:), cell center first adiabatic coefficient
!   g3:       ('gamma_3') real, dimension(:,:,:), cell center third adiabatic coefficient
!   cs:       ('sound_speed') real, dimension(:,:,:), cell center sound speed
!   mach:     ('mach_number') real, dimension(:,:,:), cell center mach number
!   mmu:      ('mean_molecular_weight') real, dimension(:,:,:), cell center mean molecular weight
!   tau:      ('optical_depth') real, dimension(:,:,:), cell center optical depth
!   Babs:     ('B_absolute') real, dimension(:,:,:), cell center absolute magnetic field
!   Bh:       ('B_horizontal') real, dimension(:,:,:), cell center horizontal magnetic field
!   B2:       ('B_signed') real, dimension(:,:,:), cell center signed absolute magnetic field
!   emag:     ('magnetic_energy') real, dimension(:,:,:), cell center magnetic energy
!   cA:       ('alfven_speed') real, dimension(:,:,:), cell center alfven speed
!   beta:     ('plasma_beta') real, dimension(:,:,:), cell center plasma beta
!   csca:     ('cs_ca') real, dimension(:,:,:), cell center sound speed to Alfven speed ratio
!
! OUTPUT:
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!   Deallocate memory
!
! RESTRICTIONS:
!
! PROCEDURE:
!   For each uneded variable check if it is allocated. If it is, deallocate.
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-02-20 (F.C. Geneva) Written
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O parameters ---
!
! Fundamental quantities
real, pointer, dimension(:,:,:), &
               optional, intent(in)     :: rho, ei, v1, v2, v3
real, pointer, dimension(:), &
               optional, intent(in)     :: xb1, xb2, xb3, xc1, xc2, xc3
! Derived quantities
real, pointer, dimension(:,:,:), &
               optional, intent(in)     :: T, P, s, j1, j2, j3, kappa, &
                                           Bb1, Bb2, Bb3, Bc1, Bc2, Bc3, divB
! Additional quantities
real, pointer, dimension(:,:,:), &
               optional, intent(in)     :: vabs, vh, ekin, plin, vmflux, g1, g3, &
                                           cs, mach, mmu, tau, Bh, Babs, B2, &
                                           emag, jabs, cA, beta, csca
!
logical, optional, intent(in)           :: base
!
! --- Local variables ---
logical                                 :: base0=.true.
!--------------------------------------------------------------------------------------------------
! === Deallocate variables and pointers ===========================================================
!
if(present(base)) then
  base0=base
endif
if (.not.present(ei).and.allocated(dq_ei).and..not.base0) then
  deallocate(dq_ei)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(rho).and.allocated(dq_rho).and..not.base0) then
  deallocate(dq_rho)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(v1).and.allocated(dq_v1).and..not.base0) then
  deallocate(dq_v1)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(v2).and.allocated(dq_v2).and..not.base0) then
  deallocate(dq_v2)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(v3).and.allocated(dq_v3).and..not.base0) then
  deallocate(dq_v3)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(Bb1).and.allocated(dq_Bb1).and..not.base0) then
  deallocate(dq_Bb1)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(Bb2).and.allocated(dq_Bb2).and..not.base0) then
  deallocate(dq_Bb2)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(Bb3).and.allocated(dq_Bb3).and..not.base0) then
  deallocate(dq_Bb3)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(xb1).and.allocated(dq_xb1).and..not.base0) then
  deallocate(dq_xb1)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(xb2).and.allocated(dq_xb2).and..not.base0) then
  deallocate(dq_xb2)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(xb3).and.allocated(dq_xb3).and..not.base0) then
  deallocate(dq_xb3)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(xc1).and.allocated(dq_xc1).and..not.base0) then
  deallocate(dq_xc1)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(xc2).and.allocated(dq_xc2).and..not.base0) then
  deallocate(dq_xc2)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(xc3).and.allocated(dq_xc3).and..not.base0) then
  deallocate(dq_xc3)
  dq_imodel=-1
  dq_status=0
endif
if (.not.present(T).and.allocated(dq_T)) then
  deallocate(dq_T)
endif
if (.not.present(P).and.allocated(dq_P)) then
  deallocate(dq_P)
endif
if (.not.present(P).and.allocated(dq_dPdei)) then
  deallocate(dq_dPdei)
endif
if (.not.present(P).and.allocated(dq_dPdrho)) then
  deallocate(dq_dPdrho)
endif
if (.not.present(j1).and.allocated(dq_j1)) then
  deallocate(dq_j1)
endif
if (.not.present(j2).and.allocated(dq_j2)) then
  deallocate(dq_j2)
endif
if (.not.present(j3).and.allocated(dq_j3)) then
  deallocate(dq_j3)
endif
if (.not.present(jabs).and.allocated(dq_jabs)) then
  deallocate(dq_jabs)
endif
if (.not.present(kappa).and.allocated(dq_kappa)) then
  deallocate(dq_kappa)
endif
if (.not.present(Bc1).and.allocated(dq_Bc1)) then
  deallocate(dq_Bc1)
endif
if (.not.present(Bc2).and.allocated(dq_Bc2)) then
  deallocate(dq_Bc2)
endif
if (.not.present(Bc3).and.allocated(dq_Bc3)) then
  deallocate(dq_Bc3)
endif
if (.not.present(divB).and.allocated(dq_divB)) then
  deallocate(dq_divB)
endif
if (.not.present(vabs).and.allocated(dq_vabs)) then
  deallocate(dq_vabs)
endif
if (.not.present(vh).and.allocated(dq_vh)) then
  deallocate(dq_vh)
endif
if (.not.present(ekin).and.allocated(dq_ekin)) then
  deallocate(dq_ekin)
endif
if (.not.present(plin).and.allocated(dq_plin)) then
  deallocate(dq_plin)
endif
if (.not.present(vmflux).and.allocated(dq_vmflux)) then
  deallocate(dq_vmflux)
endif
if (.not.present(g1).and.allocated(dq_g1)) then
  deallocate(dq_g1)
endif
if (.not.present(g3).and.allocated(dq_g3)) then
  deallocate(dq_g3)
endif
if (.not.present(cs).and.allocated(dq_cs)) then
  deallocate(dq_cs)
endif
if (.not.present(mach).and.allocated(dq_mach)) then
  deallocate(dq_mach)
endif
if (.not.present(mmu).and.allocated(dq_mmu)) then
  deallocate(dq_mmu)
endif
if (.not.present(tau).and.allocated(dq_tau)) then
  deallocate(dq_tau)
endif
if (.not.present(Babs).and.allocated(dq_Babs)) then
  deallocate(dq_Babs)
endif
if (.not.present(Bh).and.allocated(dq_Bh)) then
  deallocate(dq_Bh)
endif
if (.not.present(B2).and.allocated(dq_B2)) then
  deallocate(dq_B2)
endif
if (.not.present(emag).and.allocated(dq_emag)) then
  deallocate(dq_emag)
endif
if (.not.present(cA).and.allocated(dq_cA)) then
  deallocate(dq_cA)
endif
if (.not.present(beta).and.allocated(dq_beta)) then
  deallocate(dq_beta)
endif
if (.not.present(csca).and.allocated(dq_csca)) then
  deallocate(dq_csca)
endif
end subroutine rhd_keep_only
!
!----------*************---------------------------------------------------------------------------
subroutine rhd_int_tau
!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_int_tau ('rhd_integrate_tau')
!
! PURPOSE:
!   Integrate tau
!   
! CATEGORY:
!   General purpose
!
! CALLING SEQUENCE:
!   call rhd_int_tau()
!
! INPUT:
!
! OUTPUT:
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!   Integrates dq_kappa and stores the result in dq_tau.
!
! RESTRICTIONS:
!
! PROCEDURE:
!   Integrate with the same algorithm used in CAT
!
! EXAMPLE:
!
! MODIFICATION HISTORY:
!   2014-03-20 (F.C. Geneva) Written
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O parameters ---
!
! --- Local variables ---
integer                                 :: i1, i2, i3, k0, k1, k2, &
                                           m1, m2, m3, n1, n2, n3
real, allocatable, dimension(:,:)       :: s3, s4, s5
real, allocatable, dimension(:,:,:)     :: dkds, dxds
real, allocatable, dimension(:)         :: dz
!
!--------------------------------------------------------------------------------------------------
!
m1 = dq_m1
m2 = dq_m2
m3 = dq_m3
n1 = dq_n1
n2 = dq_n2
n3 = dq_n3
k0 = dq_n3
k1 = dq_n3-1
k2 = dq_n3-2
allocate(dkds(m1:n1,m2:n2,m3:n3-1))
allocate(dxds(m1:n1,m2:n2,m3:n3))
allocate(s3(m1:n1,m2:n2))
allocate(s4(m1:n1,m2:n2))
allocate(s5(m1:n1,m2:n2))
allocate(dz(m3:n3))
dz = dq_xb3(m3+1:n3+1)-dq_xb3(m3:n3)
if (dq_C_radHtautop>=0) then
  dq_tau(:,:,n3) = dq_kappa(:,:,n3)*dq_rho(:,:,n3)*dq_C_radHtautop
else
  dq_tau(:,:,n3) = -dq_C_radHtautop*dq_kappa(:,:,n3)*dq_P(:,:,n3)/dq_grav
endif
do i3=m3,k1
  dkds(:,:,i3) = -( dq_kappa(:,:,i3+1)*dq_rho(:,:,i3+1) - &
                    dq_kappa(:,:,i3+1)*dq_rho(:,:,i3+1) ) / dz(i3)
end do
!
! TOP
!
s3 = ( dkds(:,:,k1)*dz(k2) + dkds(:,:,k2)*dz(k1) ) / &
     ( 2.0*(dz(k1)+dz(k2)) )
s4 = min(s3, dkds(:,:,k1), dkds(:,:,k2))
s5 = max(s3, dkds(:,:,k1), dkds(:,:,k2))
dxds(:,:,k0) = 1.5*dkds(:,:,k1) - (max(s4,0.)+min(s5,0.))
!
! INTERIOR
!
do i3=n3-1,m3+1,-1
  k0 = i3
  k1 = i3
  k2 = i3-1
  s3 = ( dkds(:,:,k1)*dz(k2) + dkds(:,:,k2)*dz(k1) ) / &
       ( 2.0*(dz(k1)+dz(k2)) )
  s4 = min(s3, dkds(:,:,k1), dkds(:,:,k2))
  s5 = max(s3, dkds(:,:,k1), dkds(:,:,k2))
  dxds(:,:,k1) = 1.5*dkds(:,:,k1) - (max(s4,0.)+min(s5,0.))
end do
!
! BOTTOM
!
  dxds(:,:,m3) = 1.5*dkds(:,:,m3) - 0.5*dxds(:,:,m3+1)
!
do i3=n3-1,m3,-1
  dq_tau(:,:,i3) = dq_tau(:,:,i3+1) + dz(i3)*( 0.5*( &
    ( dq_kappa(:,:,i3+1)*dq_rho(:,:,i3+1) + dq_kappa(:,:,i3+1)*dq_rho(:,:,i3+1) ) &
    ) ) + dz(i3)*( dxds(:,:,i3+1)-dxds(:,:,i3) )/12.0
end do
!
deallocate(dkds, dxds, s3, s4, s5, dz)
end subroutine rhd_int_tau
!
!----------*************---------------------------------------------------------------------------
function diff3D(qc, vc, vb, iv)
!--------------------------------------------------------------------------------------------------
! NAME:
!   diff3D ('differentiate_3D')
!
! PURPOSE:
!   Compute cell centered derivatives from a boundary centerd quantity.
!
! CATEGORY:
!   General purpose
!
! CALLING SEQUENCE:
!   diff3D = diff3D(qc, vc, vb, iv)
!
! INPUT:
!
!   qc:          ('q_c'), real, dimension(:,:,:), quantity to differentiate
!   vc:          ('v_c'), real, dimension(:), differentiate with respect to
!                         this grid axis, values at cell centre
!   vb:          ('v_b'), real, dimension(:), differentiate with respect to
!                         this grid axis, values at cell boundaries
!   iv:          ('i_v'), integer, index of axis
!
! INPUT/OUTPUT:
!
! OUTPUT:
!   diff3D:      ('rhd_diff_3D'), real, dimension(:,:,:), derivative of 'q_b'
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE: Linearly interpolate boundary centered quantities from cell centered
!            quantities (at box boundary linearly extrapolate from the two nearer
!            values) and compute the derivatives as finite differences from boundary
!            centerd values.
!
! EXAMPLE:
!
! TO DO:
!
! MODIFICATION HISTORY:
!   2014-03-21 (F.C. Geneva) Written
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O parameters ---
real, dimension(dq_m1:dq_n1,dq_m2:dq_n2,dq_m3:dq_n3) &
                                        :: diff3D
real, allocatable, dimension(:,:,:)     :: qc
real, allocatable, dimension(:)         :: vc, vb
integer                                 :: iv
!
! --- Local variables ---
!
integer                                 :: m1, m2, m3, n1, n2, n3
real, allocatable, dimension(:,:,:)     :: qb
!--------------------------------------------------------------------------------------------------
!
m1 = dq_m1
m2 = dq_m2
m3 = dq_m3
n1 = dq_n1
n2 = dq_n2
n3 = dq_n3
!
select case (iv)
  case (1)
    allocate(qb(m1:n1+1,m2:n2,m3:n3))
    qb(m1,:,:)=qc(m1,:,:)+(qc(m1+1,:,:)-qc(m1,:,:))/(vc(m1+1)-vc(m1))*(vb(m1)-vc(m1))
    qb(n1+1,:,:)=qc(n1,:,:)+(qc(n1,:,:)-qc(n1-1,:,:))/(vc(n1)-vc(n1-1))*(vb(n1+1)-vc(n1))
    qb(m1+1:n1,:,:)=0.5*(qc(m1:n1-1,:,:)+qc(m1+1:n1,:,:))
    diff3D = (qb(m1+1:n1+1,:,:)-qb(m1:n1,:,:)) / (vb(2)-vb(1))
  case (2)
    allocate(qb(m1:n1,m2:n2+1,m3:n3))
    qb(:,m2,:)=qc(:,m2,:)+(qc(:,m2+1,:)-qc(:,m2,:))/(vc(m2+1)-vc(m2))*(vb(m2)-vc(m2))
    qb(:,n2+1,:)=qc(:,n2,:)+(qc(:,n2,:)-qc(:,n2-1,:))/(vc(n2)-vc(n2-1))*(vb(n2+1)-vc(n2))
    qb(:,m2+1:n2,:)=0.5*(qc(:,m2:n2-1,:)+qc(:,m2+1:n2,:))
    diff3D = (qb(:,m2+1:n2+1,:)-qb(:,m2:n2,:)) / (vb(2)-vb(1))
  case (3)
    allocate(qb(m1:n1,m2:n2,m3:n3+1))
    qb(:,:,m3)=qc(:,:,m3)+(qc(:,:,m3+1)-qc(:,:,m3))/(vc(m3+1)-vc(m3))*(vb(m3)-vc(m3))
    qb(:,:,n3+1)=qc(:,:,n3)+(qc(:,:,n3)-qc(:,:,n3-1))/(vc(n3)-vc(n3-1))*(vb(n3+1)-vc(n3))
    qb(:,:,m3+1:n3)=0.5*(qc(:,:,m3:n3-1)+qc(:,:,m3+1:n3))
    diff3D = (qb(:,:,m3+1:n3+1)-qb(:,:,m3:n3)) / (vb(2)-vb(1))
end select
deallocate(qb)
end function diff3D
!
!----------*************---------------------------------------------------------------------------
subroutine resample1D(qf, qr)
!--------------------------------------------------------------------------------------------------
! NAME:
!   resample3D ('resample_3D')
!
! PURPOSE:
!   Resample a given quantity to fit a given array.
!
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call resample1D(qf, qr)
!
! INPUT:
!   qf:          ('q_f'), real, dimension(:,:,:), quantity to resample
!
! INPUT/OUTPUT:
!
! OUTPUT:
!   qr:          ('q_r'), real, dimension(:,:,:), array of desired size
!
! VARIABLES:
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE: Evenly resample given quantity.
!
! EXAMPLE:
!
! TO DO:
!
! MODIFICATION HISTORY:
!   2014-03-23 (F.C. Geneva) Written
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O parameters ---
real, allocatable, dimension(:), &
                   intent(inout)        :: qf
real, pointer, dimension(:), &
                   intent(inout)        :: qr
!
! --- Local variables ---
!
integer                                 :: m, n, r, i
!--------------------------------------------------------------------------------------------------
m = lbound(qr,1)
n = ubound(qr,1)
r = (ubound(qf,1)-lbound(qf,1)+1) / size(qr)

do i=m,n
  !
  qr(i) = qf(1+r*(i-m))
  !
end do ! i
end subroutine resample1D
!
!----------*************---------------------------------------------------------------------------
subroutine resample3D(qf, qr)
!--------------------------------------------------------------------------------------------------
! NAME:
!   resample3D ('resample_3D')
!
! PURPOSE:
!   Resample a given quantity to fit a given array.
!
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   call resample3D(qf, qr)
!
! INPUT:
!   qf:          ('q_f'), real, dimension(:,:,:), quantity to resample
!
! INPUT/OUTPUT:
!
! OUTPUT:
!   qr:          ('q_r'), real, dimension(:,:,:), array of desired size
!
! VARIABLES:
!
! ROUTINES:
!
! MODULES:
!
! SIDE EFFECTS:
!
! RESTRICTIONS:
!
! PROCEDURE: Evenly resample given quantity.
!
! EXAMPLE:
!
! TO DO:
!
! MODIFICATION HISTORY:
!   2014-03-23 (F.C. Geneva) Written
!--------------------------------------------------------------------------------------------------
!
implicit none
!
! --- I/O parameters ---
real, allocatable, dimension(:,:,:), &
                   intent(inout)        :: qf
real, pointer, dimension(:,:,:), &
                   intent(inout)        :: qr
!
! --- Local variables ---
!
integer                                 :: m1, m2, m3, n1, n2, n3, &
                                           r1, r2, r3, i1, i2, i3
!--------------------------------------------------------------------------------------------------
m1 = lbound(qr,1)
m2 = lbound(qr,2)
m3 = lbound(qr,3)
n1 = ubound(qr,1)
n2 = ubound(qr,2)
n3 = ubound(qr,3)
r1 = (dq_n1-dq_m1+1) / size(qr, 1)
r2 = (dq_n2-dq_m2+1) / size(qr, 2)
r3 = (dq_n3-dq_m3+1) / size(qr, 3)

do i3=m3,n3
  do i2=m2,n2
    do i1=m1,n1
      !
      qr(i1,i2,i3) = qf(1+r1*(i1-m1),1+r2*(i2-m2),1+r3*(i3-m3))
      !
    end do ! i1
  end do ! i2
end do ! i3
end subroutine resample3D
!
end module rhd_dq_module
