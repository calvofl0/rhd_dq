!**************************************************************************************************
!
!   #####   #    #  #####          #####    ####
!   #    #  #    #  #    #         #    #  #    #
!   #    #  ######  #    #         #    #  #    #
!   #####   #    #  #    #         #    #  #    #
!   #   #   #    #  #    #         #    #  #    #
!   #    #  #    #  #####  ######  #####    ###### _box
!
!   Modules for Radiation-Hydrodynamics Code
!
!**************************************************************************************************
!   Fortran90
!   Flavio Calvo:          Geneva, Locarno, Freiburg
!   2014-03-26
!**************************************************************************************************
! Modules:
!   rhd_dq_box:  Derived quantities box definition
!**************************************************************************************************

!------*************-------------------------------------------------------------------------------
#ifndef PYBOLD
module rhd_dq_box
#else
module box
#endif

!--------------------------------------------------------------------------------------------------
! NAME:
!   rhd_dq_box ('rhd_derived_quantities_box')
!
! PURPOSE:
!   Derived quantities box definition
!
! CATEGORY:
!   Input/Output
!
! CALLING SEQUENCE:
!   use rhd_dq_box
!
! VARIABLES:
!   None
!
! SUBROUTINES: (contained)
!
! MODIFICATION HISTORY:
!   2014-03-26 (F.C. Geneva) Written
!--------------------------------------------------------------------------------------------------
implicit none
!
! --- Global parameters ---
!
! Files
character(len=256)                     :: dq_parfile, dq_file
integer                                :: dq_nc_p
! Fundamental quantities
real, target, allocatable, &
                   dimension(:,:,:)    :: dq_rho, dq_ei, dq_v1, dq_v2, dq_v3, &
                                          dq_Bb1, dq_Bb2, dq_Bb3
real, target, allocatable, &
                   dimension(:)        :: dq_xb1, dq_xb2, dq_xb3, &
                                          dq_xc1, dq_xc2, dq_xc3
! Derived quantities
real, target, allocatable, &
                   dimension(:,:,:)    :: dq_T, dq_P, dq_dPdrho, dq_dPdei, &
                                          dq_s, dq_j1, dq_j2, dq_j3, dq_jabs, &
                                          dq_kappa, dq_Bc1, dq_Bc2, dq_Bc3, &
                                          dq_divB
! Additional quantities
real, target, allocatable, &
                   dimension(:,:,:)    :: dq_vabs, dq_vh, dq_ekin, dq_plin, &
                                          dq_vmflux, dq_g1, dq_g3, dq_cs, &
                                          dq_mach, dq_mmu, dq_tau, dq_Bh, &
                                          dq_Babs, dq_B2, dq_emag, dq_cA, &
                                          dq_beta, dq_csca
! Global variables
real                                   :: dq_C_radHtautop, dq_grav, dq_time, &
                                          dq_deltax, dq_deltay, dq_deltaz
integer                                :: dq_itime, dq_imodel=-1, dq_nmodel, &
                                          dq_status=0, dq_nx, dq_ny, dq_nz, &
                                          dq_n1, dq_n2, dq_n3, &
                                          dq_m1, dq_m2, dq_m3
logical                                :: dq_Bb_flag
#ifndef PYBOLD
end module rhd_dq_box
#else
end module box
#endif
