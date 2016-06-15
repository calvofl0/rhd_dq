subroutine get_dq(dq)
use uio_reader_module
use rhd_dq_module
implicit none
integer, intent(in)                     :: dq
real, pointer, dimension(:,:,:)         :: tmp3D
!
if (iand(dq,z'00000001').eq.z'00000001') call rhd_dq_vars(T=tmp3D)
if (iand(dq,z'00000002').eq.z'00000002') call rhd_dq_vars(P=tmp3D)
if (iand(dq,z'00000004').eq.z'00000004') call rhd_dq_vars(s=tmp3D)
if (iand(dq,z'00000008').eq.z'00000008') call rhd_dq_vars(j1=tmp3D)
if (iand(dq,z'00000010').eq.z'00000010') call rhd_dq_vars(j2=tmp3D)
if (iand(dq,z'00000020').eq.z'00000020') call rhd_dq_vars(j3=tmp3D)
if (iand(dq,z'00000040').eq.z'00000040') call rhd_dq_vars(jabs=tmp3D)
if (iand(dq,z'00000080').eq.z'00000080') call rhd_dq_vars(kappa=tmp3D)
if (iand(dq,z'00000100').eq.z'00000100') call rhd_dq_vars(Bc1=tmp3D)
if (iand(dq,z'00000200').eq.z'00000200') call rhd_dq_vars(Bc2=tmp3D)
if (iand(dq,z'00000400').eq.z'00000400') call rhd_dq_vars(Bc3=tmp3D)
if (iand(dq,z'00000800').eq.z'00000800') call rhd_dq_vars(divB=tmp3D)
if (iand(dq,z'00001000').eq.z'00001000') call rhd_dq_vars(vabs=tmp3D)
if (iand(dq,z'00002000').eq.z'00002000') call rhd_dq_vars(vh=tmp3D)
if (iand(dq,z'00004000').eq.z'00004000') call rhd_dq_vars(ekin=tmp3D)
if (iand(dq,z'00008000').eq.z'00008000') call rhd_dq_vars(plin=tmp3D)
if (iand(dq,z'00010000').eq.z'00010000') call rhd_dq_vars(vmflux=tmp3D)
if (iand(dq,z'00020000').eq.z'00020000') call rhd_dq_vars(g1=tmp3D)
if (iand(dq,z'00040000').eq.z'00040000') call rhd_dq_vars(g3=tmp3D)
if (iand(dq,z'00080000').eq.z'00080000') call rhd_dq_vars(cs=tmp3D)
if (iand(dq,z'00100000').eq.z'00100000') call rhd_dq_vars(mach=tmp3D)
if (iand(dq,z'00200000').eq.z'00200000') call rhd_dq_vars(mmu=tmp3D)
if (iand(dq,z'00400000').eq.z'00400000') call rhd_dq_vars(tau=tmp3D)
if (iand(dq,z'00800000').eq.z'00800000') call rhd_dq_vars(Babs=tmp3D)
if (iand(dq,z'01000000').eq.z'01000000') call rhd_dq_vars(Bh=tmp3D)
if (iand(dq,z'02000000').eq.z'02000000') call rhd_dq_vars(B2=tmp3D)
if (iand(dq,z'04000000').eq.z'04000000') call rhd_dq_vars(emag=tmp3D)
if (iand(dq,z'08000000').eq.z'08000000') call rhd_dq_vars(cA=tmp3D)
if (iand(dq,z'10000000').eq.z'10000000') call rhd_dq_vars(beta=tmp3D)
if (iand(dq,z'20000000').eq.z'20000000') call rhd_dq_vars(csca=tmp3D)
end subroutine get_dq
