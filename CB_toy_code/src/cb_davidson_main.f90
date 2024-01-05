   program cb_davidson_main

! global variables
   USE cb_module
#if defined(__MPI)
   use mp_global,            ONLY : mp_startup, mp_global_end
   use mp_world,             ONLY : world_comm
   use mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif

   implicit none
   !
   include 'laxlib.fh'
   !
! local variables (used in the call to cegterg )
   logical, parameter :: gamma_only = .false. ! general k-point version
   complex(DP), allocatable :: evc(:,:)
   real(dp), allocatable :: eig(:)
   integer, parameter :: npol=1
   integer :: notcnv, dav_iter, nhpsi
   logical :: overlap = .false. , lrot =.false.
! additional local variables
   real(dp) :: ref=0.d0
#if defined(__MPI)
! local paralelization variables
   integer :: ndiag     ! input value of processors in the diagonalization group
   logical :: do_distr_diag_in_band_group = .false. ! whether or not the parallel diagonalization is performed inside the
                                                    ! band group or at the parallelization level above.
#endif
!------------------------------------------------------------------------
   external cb_h_psi, cb_s_psi, cb_g_psi
!  subroutine cb_h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
!  subroutine cb_s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
!  subroutine cb_g_psi(npwx,npw,nvec,psi,eig)   computes G*psi -> psi

#if defined(__MPI)
! this call creates the parallel communicators in the MAIN code 
  call mp_startup ( ndiag, diag_in_band_group = do_distr_diag_in_band_group )   
!--- THIS PART IS RELEVANT FOR THE PARALLEL USE OF THE ROUTINE IN KS_Solvers/Davidson -------------------------!
! this set the mpi communicators used internally by the routines in the Davidson library
! it passes 1) the top parent level communicator (could be different from world_comm)
!           2) the sub-communicator of the band group 
!           3) the communicator used across band groups
!           4) whether the distributed diagonalization is performed inside the band group or at the top level
 call set_mpi_comm_4_solvers( world_comm, intra_bgrp_comm, inter_bgrp_comm )

!--------------------------------------------------------------------------------------------------------------!
#endif

   call init_clocks(.true.)
   call input(gamma_only)
   call ggen(gamma_only)
   call set_cb_potential

   if (use_overlap) write(*,*) '** TEST:  CB hamiltonian modified so as to need an overlap matrix **'
   overlap = use_overlap

   do current_k=1,nks
     call init_k
     allocate( evc(npwx,nbnd), eig(nbnd) )
     call init_random_wfcs(npw,npwx,nbnd,evc)

     call start_clock('davidson')
!--- THIS IS THE RELEVANT CALL TO THE ROUTINE IN KS_Solvers/Davidson ------------------------------------------!
#if defined(__MPI)
     write (stdout,*) 'ndiag', ndiag
     if (ndiag == 1) then
#endif
        call cegterg( cb_h_psi, cb_s_psi, overlap, cb_g_psi, &
                      npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                      eig, btype, notcnv, lrot, dav_iter, nhpsi )
#if defined(__MPI)
     else
        call pcegterg(cb_h_psi, cb_s_psi, overlap, cb_g_psi, &
                      npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                      eig, btype, notcnv, lrot, dav_iter, nhpsi )
     end if
#endif
!--------------------------------------------------------------------------------------------------------------!

     call stop_clock('davidson')
     deallocate( evc )

     if (energy_shift .and. current_k==1) ref=eig(4*ncell**3)
     call write_bands(eig,ref)
     write (stdout,*) 'dav_iter, nhpsi, notcnv, ethr ', dav_iter, nhpsi, notcnv, ethr

     deallocate( eig )
   end do

   call print_clock('davidson')

   call print_clock( 'cegterg' )
   call print_clock( 'cegterg:init' )
   call print_clock( 'cegterg:diag' )
   call print_clock( 'cegterg:update' )
   call print_clock( 'cegterg:overlap' )
   call print_clock( 'cegterg:last' )

   call print_clock('h_psi')
   call print_clock('s_psi')
   call print_clock('g_psi')
! 
  write (6,*)
  write (6,*) ' general FFT  routines'
  call print_clock('fftw')
  call print_clock('ffts')

#if defined(__MPI)
   call mp_global_end( )
#endif

   end program cb_davidson_main
