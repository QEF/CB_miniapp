   program cb_cg_gamma_main

! global variables
   USE cb_module
#if defined(__MPI)
   use mp_global,            ONLY : mp_startup, mp_global_end
   use mp_world,             ONLY : world_comm
   use mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif

   implicit none
! local variables (used in the call to cegterg )
   logical, parameter :: gamma_only=.true. ! general k-point version
   complex(DP), allocatable :: evc(:,:)
   real(dp), allocatable :: eig(:)
   integer, parameter :: npol=1
   integer :: notconv, nhpsi, ig
   real(dp):: cg_iter
   integer :: max_cg_iter = 100 !default of QE 
   real(dp), allocatable :: h_diag(:,:) !in case of CG  - the preconditioner
   logical :: overlap = .false. , lrot =.false. , lscf = .true. ! lscf is true for scf calc
! additional local variables
   real(dp) :: ref=0.d0
#if defined(__MPI)
! local paralelization variables
   integer :: ndiag     ! input value of processors in the diagonalization group
   logical :: do_distr_diag_in_band_group = .false. ! whether or not the parallel diagonalization is performed inside the
                                                    ! band group or at the parallelization level above.
#endif
!------------------------------------------------------------------------
   external cb_hs_1psi, cb_s_1psi, cb_h_psi, cb_s_psi
!  subroutine cb_hs_1psi(npwx,npw,psi,hpsi,spsi)  computes H*psi and S*psi
!  subroutine cb_s_1psi (npwx,npw,psi,spsi)  computes S*psi (if needed)
!  subroutine cb_g_psi(npwx,npw,nvec,psi,eig)   computes G*psi -> psi

#if defined(__MPI)
! this call creates the parallel communicators in the MAIN code 
call mp_startup ( ndiag, diag_in_band_group = do_distr_diag_in_band_group )   
!--- THIS PART IS RELEVANT FOR THE PARALLEL USE OF THE ROUTINE IN KS_Solvers/Davidson -------------------------!
! this set the mpi communicators used internally by the routines in the Davidson library
! it passes 1) the top level communicator 
!           2) the sub-communicator of the band group 
!           3) the communicator used across band groups
!           4) whether the distributed diagonalization is performed inside the band group or at the top level
call set_mpi_comm_4_solvers( world_comm, intra_bgrp_comm, inter_bgrp_comm )
!--------------------------------------------------------------------------------------------------------------!
#endif

   call init_clocks(.true.)
   call input(gamma_only)
   call ggen(gamma_only); call export_gstart_2_solvers(gstart)
   call set_cb_potential

   if (use_overlap) write(*,*) '** TEST:  CB hamiltonian modified so as to need an overlap matrix **'
   overlap = use_overlap

   !preconditioninig is still needed in CG
   if (.not. allocated(h_diag)) allocate ( h_diag( npwx, npol ))
   do current_k=1,nks
     call init_k
     allocate( evc(npwx,nbnd), eig(nbnd) )
     call init_random_wfcs(npw,npwx,nbnd,evc)

     call start_clock('cg')
!--- THESE ARE THE RELEVANT CALLS TO THE ROUTINE IN KS_Solvers/CG ------------------------------------------!
     write(stdout,*) ' subspace rotation first'
#if defined(__MPI)
     write (6,*) 'ndiag', ndiag
     if ( ndiag == 1 ) then
#endif
        call rotate_wfc_gamma( cb_h_psi, cb_s_psi, overlap, &
                               npwx, npw, nbnd, nbnd, evc, evc, eig )
#if defined(__MPI)
     else
        call protate_wfc_gamma( cb_h_psi, cb_s_psi, overlap, &
                                npwx, npw, nbnd, nbnd, evc, evc, eig )
     endif
#endif
     nhpsi = nbnd
     write(stdout,*) ' then cg diagonalization '

     h_diag = 1.D0
     FORALL( ig = 1 : npwx )   
        h_diag(ig,:) = 1.D0 + ekin(ig) + SQRT( 1.D0 + ( ekin(ig) - 1.D0)**2)
     END FORALL 
     !in QE iter and ntry would be used here to determine the lrot, which would call rotate_wfc.
     CALL rcgdiagg( cb_hs_1psi, cb_s_1psi, h_diag, &
                    npwx, npw, nbnd, evc, eig, btype, & 
                    ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )

     nhpsi = nhpsi + cg_iter*nbnd
!--------------------------------------------------------------------------------------------------------------!

     call stop_clock('cg')
     deallocate( evc )

     if (energy_shift .and. current_k==1) ref=eig(4*ncell**3)
     call write_bands(eig,ref)
     write (6,*) 'cg_iter, nhpsi, notconv, ethr ',int(cg_iter), nhpsi, notconv, ethr
     deallocate( eig )
   end do

   call print_clock('cg')
   CALL print_clock('rotwfcg')      ; CALL print_clock('protwfcg')
   CALL print_clock('rotwfcg:hpsi') ; CALL print_clock('protwfcg:hpsi')
   CALL print_clock('rotwfcg:hc')   ; CALL print_clock('protwfcg:hc')
   CALL print_clock('rotwfcg:diag') ; CALL print_clock('pprotwfcg:diag')
   CALL print_clock('rotwfcg:evc')  ; CALL print_clock('protwfcg:evc')
   call print_clock('rcgdiagg')
   call print_clock('cg:ortho')
   call print_clock('h_psi')
   call print_clock('s_psi')
   call print_clock('hs_1psi')
   call print_clock('s_1psi')
! 
  write (6,*)
  write (6,*) ' general FFT  routines'
  call print_clock('fftw')
  call print_clock('ffts')

#if defined(__MPI)
   call mp_global_end( )
#endif

   end program cb_cg_gamma_main
