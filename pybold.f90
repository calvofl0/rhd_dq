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
