subroutine init()
use uio_bulk_module
implicit none
call uio_init(progrm='pybold')
end subroutine init
!
subroutine pyread(uio_data_reg, modelfile, parfile, imodel)
use uio_reader_module
use rhd_dq_module
implicit none
external uio_data_reg
character*(*),              intent(in)  :: modelfile
character*(*),              intent(in)  :: parfile
integer,                    intent(in)  :: imodel
integer                                 :: imodel0, nmodel
!
record => uio_data_reg
call uio_lbl_count(modelfile, startlabel='dataset', endlabel='enddataset', &
                   nlabel=nmodel)
if (trim(parfile) /= '' .and. nmodel>0) then
        use_dq = .True.
else
        use_dq = .False.
endif
if (imodel<0) then
        imodel0=nmodel+imodel+1
else
        imodel0=imodel
endif
!print*, 'nmodel: ', nmodel
if (use_dq) then
        call rhd_dq_init(parfile, modelfile, nmodel, imodel0, 6)
endif
if (nmodel>0) then
        call uio_struct_rd(modelfile, startlabel='dataset', &
                           endlabel='enddataset', nlabel=imodel0, &
                           substartlabel='box', subendlabel='endbox', &
                           headread=.true.)
else
        call uio_struct_rd(modelfile, endlabel='enddataset')
endif
if (use_dq) then
        call rhd_dq_end_init()
endif
end subroutine pyread
!
subroutine open_model(modelfile, parfile)
use uio_reader_module
implicit none
character*(*),              intent(in)  :: modelfile
character*(*),              intent(in)  :: parfile
!
r_modelfile = modelfile
r_parfile = parfile
call uio_openrd(r_nc, r_modelfile)
call uio_lbl_count(r_modelfile, r_nc, &
                   startlabel='dataset', endlabel='enddataset', &
                   nlabel=r_nmodel, maxlabel=1)
if (trim(r_parfile) /= '' .and. r_nmodel>0) then
        use_dq = .True.
else
        use_dq = .False.
endif
call uio_chpos(r_nc, 'rewind')
r_imodel = 0
if (r_nmodel>0) r_nmodel = -1
end subroutine open_model
!
subroutine count_model()
use uio_reader_module
implicit none
integer                                 :: pos
!
call uio_chpos(r_nc, 'get', posout=pos)
call uio_chpos(r_nc, 'rewind')
call uio_lbl_count(r_modelfile, r_nc, &
                   startlabel='dataset', endlabel='enddataset', nlabel=r_nmodel)
call uio_chpos(r_nc, 'goto', pos)
end subroutine count_model
!
subroutine close_model()
use uio_reader_module
implicit none
!
call uio_closrd(r_nc)
end subroutine close_model
!
subroutine read_header(uio_data_reg)
use uio_reader_module
implicit none
external uio_data_reg
integer                                 :: pos
!
record => uio_data_reg
if (r_nmodel/=0) then
        call uio_chpos(r_nc, 'get', posout=pos)
        call uio_chpos(r_nc, 'rewind')
        call uio_struct_rd(r_modelfile, channel=r_nc, headread=.true.)
        call uio_chpos(r_nc, 'goto', pos)
endif
end subroutine read_header
!
subroutine read_next(uio_data_reg)
use uio_reader_module
implicit none
external uio_data_reg
!
record => uio_data_reg
if (r_imodel < r_nmodel .or. r_nmodel < 0) then
        r_imodel = r_imodel+1
        if (use_dq) then
                call rhd_dq_init(r_parfile, r_modelfile, r_nmodel, r_imodel, 6)
        endif
        if (r_nmodel/=0) then
                call uio_struct_rd(r_modelfile, channel=r_nc, &
                                   startlabel='dataset', &
                                   endlabel='enddataset', &
                                   substartlabel='box', &
                                   subendlabel='endbox')
        else
                !print*, 'Inside else !!!'
                call uio_struct_rd(r_modelfile, channel=r_nc, &
                                   endlabel='enddataset')
        endif
        if (use_dq .and. r_imodel >=0) then
                call rhd_dq_end_init()
        endif
endif
end subroutine read_next
!
subroutine get_dq(uio_data_reg, dq)
use uio_reader_module
use rhd_dq_module
implicit none
external uio_data_reg
character*(*), intent(in)               :: dq
real, pointer, dimension(:,:,:)         :: tmp3D
!
record => uio_data_reg
select case (dq)
  case ('T')
          call rhd_dq_vars(T=tmp3D)
          label='T'
          call record('link')
  case ('P')
          call rhd_dq_vars(P=tmp3D)
          label='P'
          call record('link')
  case ('s')
          call rhd_dq_vars(s=tmp3D)
          label='s'
          call record('link')
  case ('j1')
          call rhd_dq_vars(j1=tmp3D)
          label='j1'
          call record('link')
  case ('j2')
          call rhd_dq_vars(j2=tmp3D)
          label='j2'
          call record('link')
  case ('j3')
          call rhd_dq_vars(j3=tmp3D)
          label='j3'
          call record('link')
  case ('jabs')
          call rhd_dq_vars(jabs=tmp3D)
          label='jabs'
          call record('link')
  case ('kappa')
          call rhd_dq_vars(kappa=tmp3D)
          label='kappa'
          call record('link')
  case ('Bc1')
          call rhd_dq_vars(Bc1=tmp3D)
          label='Bc1'
          call record('link')
  case ('Bc2')
          call rhd_dq_vars(Bc2=tmp3D)
          label='Bc2'
          call record('link')
  case ('Bc3')
          call rhd_dq_vars(Bc3=tmp3D)
          label='Bc3'
          call record('link')
  case ('divB')
          call rhd_dq_vars(divB=tmp3D)
          label='divB'
          call record('link')
  case ('vabs')
          call rhd_dq_vars(vabs=tmp3D)
          label='vabs'
          call record('link')
  case ('vh')
          call rhd_dq_vars(vh=tmp3D)
          label='vh'
          call record('link')
  case ('ekin')
          call rhd_dq_vars(ekin=tmp3D)
          label='ekin'
          call record('link')
  case ('plin')
          call rhd_dq_vars(plin=tmp3D)
          label='plin'
          call record('link')
  case ('vmflux')
          call rhd_dq_vars(vmflux=tmp3D)
          label='vmflux'
          call record('link')
  case ('g1')
          call rhd_dq_vars(g1=tmp3D)
          label='g1'
          call record('link')
  case ('g3')
          call rhd_dq_vars(g3=tmp3D)
          label='g3'
          call record('link')
  case ('cs')
          call rhd_dq_vars(cs=tmp3D)
          label='cs'
          call record('link')
  case ('mach')
          call rhd_dq_vars(mach=tmp3D)
          label='mach'
          call record('link')
  case ('mmu')
          call rhd_dq_vars(mmu=tmp3D)
          label='mmu'
          call record('link')
  case ('tau')
          call rhd_dq_vars(tau=tmp3D)
          label='tau'
          call record('link')
  case ('Babs')
          call rhd_dq_vars(Babs=tmp3D)
          label='Babs'
          call record('link')
  case ('Bh')
          call rhd_dq_vars(Bh=tmp3D)
          label='Bh'
          call record('link')
  case ('B2')
          call rhd_dq_vars(B2=tmp3D)
          label='B2'
          call record('link')
  case ('emag')
          call rhd_dq_vars(emag=tmp3D)
          label='emag'
          call record('link')
  case ('cA')
          call rhd_dq_vars(cA=tmp3D)
          label='cA'
          call record('link')
  case ('beta')
          call rhd_dq_vars(beta=tmp3D)
          label='beta'
          call record('link')
  case ('csca')
          call rhd_dq_vars(csca=tmp3D)
          label='csca'
          call record('link')
end select
end subroutine get_dq
!
subroutine write_model(modelfile, xb1, xb2, xb3, xc1, xc2, xc3, &
                       v1, v2, v3, rho, ei, Bb_flag, Bb1, Bb2, Bb3, &
                       B1_unit, B2_unit, B3_unit, dtime, itime, time, &
                       time_db, description, history, version, &
                       time_out_mean_last, time_out_full_last, &
                       m1, n1, m2, n2, m3, n3)
!use const_module
use rhd_gl_module
use rhd_box_module
use rhd_sub_module
!use rhd_io_module
!use tabinter_coeff_module
use uio_base_module
use rhd_prop_module
implicit none
!
! --- I/O parameters ---
character*(*),              intent(in)  :: modelfile
integer,                    intent(in)  :: itime, m1, n1, m2, n2, m3, n3
real, target, dimension(m1:n1+1), &
                         intent(inout)  :: xb1
real, target, dimension(m2:n2+1), &
                         intent(inout)  :: xb2
real, target, dimension(m3:n3+1), &
                         intent(inout)  :: xb3
real, target, dimension(m1:n1), &
                         intent(inout)  :: xc1
real, target, dimension(m2:n2), &
                         intent(inout)  :: xc2
real, target, dimension(m3:n3), &
                         intent(inout)  :: xc3
real, target, &
       dimension(m1:n1,m2:n2,m3:n3), &
                         intent(inout)  :: v1, v2, v3, rho, ei
real, target, &
     dimension(m1:n1+1,m2:n2,m3:n3), &
                         intent(inout)  :: Bb1
real, target, &
     dimension(m1:n1,m2:n2+1,m3:n3), &
                         intent(inout)  :: Bb2
real, target, &
     dimension(m1:n1,m2:n2,m3:n3+1), &
                         intent(inout)  :: Bb3
real,                       intent(in)  :: dtime, time, &
                                           time_out_mean_last,time_out_full_last
real(kind=kind(0.0D+00)),   intent(in)  :: time_db
character, dimension(:,:), &
                            intent(in)  :: history
character(len=80),          intent(in)  :: version, description, &
                                           B1_unit, B2_unit, B3_unit
logical,                    intent(in)  :: Bb_flag
!
! --- Local variables ---
type(box_type)                          :: box
character(len=80)                       :: date_end
character(len=80), dimension(nhismax)   :: history_in, history_end
character(len=80), dimension(1)         :: description_end
integer                                 :: ncbox, nhis_in, nhis_end, nchars, &
                                           i, j
!
nhis_in=size(history,1)
nchars=size(history,2)
do i=1,nhismax
  history_in(i)=repeat(' ',80)
  history_end(i)=repeat(' ',80)
end do
do i=1,nhis_in
  do j=1,nchars
    history_in(i)(j:j)=history(i,j)
  end do
end do
!print*, modelfile, xb1, xb2, xb3, xc1, xc2, xc3, &
!                       v1, v2, v3, rho, ei, Bb_flag, Bb1, Bb2, Bb3, &
!                       itime, time, time_db, history_in, &
!                       time_out_mean_last, time_out_full_last, &
!                       m1, n1, m2, n2, m3, n3, shape(history)
!
call rhd_box_Init(box)
!
box%m1=m1
box%n1=n1
box%m2=m2
box%n2=n2
box%m3=m3
box%n3=n3
!
box%itime=itime
box%time=time
box%time_db=time_db
!
box%xb1=>xb1(m1:n1+1)
box%xb2=>xb2(m2:n2+1)
box%xb3=>xb3(m3:n3+1)
box%xc1=>xc1(m1:n1)
box%xc2=>xc2(m2:n2)
box%xc3=>xc3(m3:n3)
!
box%rho=>rho(m1:n1,m2:n2,m3:n3)
box%ei=>ei(m1:n1,m2:n2,m3:n3)
box%v1=>v1(m1:n1,m2:n2,m3:n3)
box%v2=>v2(m1:n1,m2:n2,m3:n3)
box%v3=>v3(m1:n1,m2:n2,m3:n3)
!
prop%Bb1_IsInFile_flag=.false.
prop%Bb2_IsInFile_flag=.false.
prop%Bb3_IsInFile_flag=.false.
if(Bb_flag) then
  box%Bb1=>Bb1(m1:n1+1,m2:n2,m3:n3)
  box%Bb2=>Bb2(m1:n1,m2:n2+1,m3:n3)
  box%Bb3=>Bb3(m1:n1,m2:n2,m3:n3+1)
  prop%Bb1_IsInFile_flag=.true.
  prop%Bb2_IsInFile_flag=.true.
  prop%Bb3_IsInFile_flag=.true.
  prop%B1_unit=B1_unit
  prop%B2_unit=B2_unit
  prop%B3_unit=B3_unit
endif
!
! --- History ---
call uio_dattim(date_end)
history_end(1:nhis_in)=history_in(1:nhis_in)
if (nhis_in > (nhismax-1)) then
  history_end(nhismax/2-1)='...'
  history_end(nhismax/2:nhismax-3)=history_in(nhis_in-(nhismax-3)+nhismax/2:nhis_in)
  nhis_end=nhismax-1
else
  nhis_end=nhis_in
endif
nhis_end=nhis_end+1
history_end(nhis_end)='Box saved with pybold: ' // trim(date_end)
!
if (trim(description) == '') then
  description_end=(/'Pybold-generated RHD-simulation model snapshot'/)
else
  description_end(1) = description
endif
!print*, 'open,model,box,endmo,close', ncbox, dtime, &
!                  time_out_full_last, &
!                  time_out_mean_last, modelfile, &
!                  'unformatted', 'ieee_4', 1, &
!                  description, nhis_end, history_end, &
!                  trim(version)
call rhd_box_Write('open,model,box,endmo,close', ncbox, box=box, dtime=dtime, &
                  time_out_full_last=time_out_full_last, &
                  time_out_mean_last=time_out_mean_last, file=modelfile, &
                  form='unformatted', conv='ieee_4', ndes=1, &
                  description=description_end, nhis=nhis_end, &
                  history=history_end, version=trim(version))
!
! --- Open file for writing in uio-form ---
!call rhd_box_WrData('open', gl%nc_p, ncbox, file=modelfile, &
!                   form='unformatted', conv='ieee_4', ndes=1, &
!                   description=description, nhis=nhis_end, &
!                   history=history_end, version=trim(version))
! --- Write model header ---
!call rhd_box_WrData('model', gl%nc_p, ncbox, time_unit=time_unit, &
!                   itime=itime, time=time, time_db=time_db, dtime=dtime, &
!                   time_out_full_last=time_out_full_last, &
!                   time_out_mean_last=time_out_mean_last)
! --- Write modelbox ---
!call rhd_box_WrData('box', gl%nc_p, ncbox, &
!                   time_unit=time_unit, &
!                   x1_unit=x1_unit, x2_unit=x2_unit, x3_unit=x3_unit, &
!                   rho_unit=rho_unit, v_unit=v_unit, e_unit=e_unit, &
!                   box_id=box_id, &
!                   m1=m1, m2=m2, m3=m3, n1=n1, n2=n2, n3=n3, &
!                   itime=itime, time=time, time_db=time_db, dtime=dtime, &
!                   xc1=xc1, xc2=xc2, xc3=xc3, &
!                   xb1=xb1, xb2=xb2, xb3=xb3, &
!                   rho=rho, ei=ei, v1=v1, v2=v2, v3=v3, &
!                   Bb1=Bb1, Bb2=Bb2, Bb3=Bb3)
! --- Close file ---
!call rhd_box_WrData('endmo', gl%nc_p, ncbox)
! --- Close file ---
!call rhd_box_WrData('close', gl%nc_p, ncbox)
end subroutine write_model
