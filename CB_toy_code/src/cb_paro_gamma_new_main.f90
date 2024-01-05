PROGRAM cb_paro_gamma_main

  ! global variables
  USE cb_module
#if defined(__MPI)
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE mp_world,             ONLY : world_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif

  IMPLICIT NONE
  ! local variables (used in the call to paro )
  LOGICAL, PARAMETER :: gamma_only=.TRUE.  ! Gamma-point version
  COMPLEX(DP), ALLOCATABLE :: evc(:,:)
  REAL(dp), ALLOCATABLE :: eig(:)
  INTEGER, PARAMETER :: npol=1
  INTEGER :: notconv, avg_iter, ig
  INTEGER :: max_paro_iter = 100 !default of QE
  LOGICAL :: overlap = .FALSE. , lrot =.FALSE. , lscf = .TRUE. ! lscf is true for scf calc
  ! additional local variables
  REAL(dp) :: ref=0.d0
#if defined(_OPENMP)
  INTEGER :: omp_get_max_threads
  EXTERNAL :: omp_get_max_threads
#endif
#if defined(__MPI)
  ! local paralelization variables
  INTEGER :: ndiag     ! input value of processors in the diagonalization group
  LOGICAL :: do_distr_diag_in_band_group = .FALSE. ! whether or not the parallel diagonalization is performed inside the
  ! band group or at the parallelization level above.
#endif
  !------------------------------------------------------------------------
  EXTERNAL cb_h_psi, cb_s_psi, cb_hs_psi, cb_g_1psi
  ! subroutine cb_h_psi(npwx,npw,nvec,evc,hpsi)  computes H*evc  using band parallelization
  ! subroutine cb_s_psi(npwx,npw,nvec,evc,spsi)  computes S*evc  using band parallelization
  ! subroutine cb_hs_psi(npwx,npw,nvec,evc,hpsi,spsi) computes H*evc and S*evc using band parallelization
  ! subroutine cb_g_1psi(npwx,npw,psi,eig) computes G*psi -> psi for a single band
  !
  ! ... local variables
  !
  INTEGER :: nhpsi
  !
  ! ... init local variables
  !
  avg_iter = 0.d0

#if defined(__MPI)
  ! this call creates the parallel communicators in the MAIN code
  CALL mp_startup ( ndiag, diag_in_band_group = do_distr_diag_in_band_group )
  !--- THIS PART IS RELEVANT FOR THE PARALLEL USE OF THE ROUTINE IN KS_Solvers/Davidson -------------------------!
  ! this set the mpi communicators used internally by the routines in the Davidson library
  ! it passes 1) the top level communicator
  !           2) the sub-communicator of the band group
  !           3) the communicator used across band groups
  !           4) whether the distributed diagonalization is performed inside the band group or at the top level
  CALL set_mpi_comm_4_solvers( world_comm, intra_bgrp_comm, inter_bgrp_comm )
  !--------------------------------------------------------------------------------------------------------------!
#endif
#if defined(_OPENMP)
  write (6,*) 'OPENMP is active with max_threads =', omp_get_max_threads()
#endif

  CALL init_clocks(.TRUE.)
  CALL input(gamma_only)
  CALL ggen(gamma_only); CALL export_gstart_2_solvers(gstart)
  CALL set_cb_potential

  if (use_overlap) write(*,*) '** TEST:  CB hamiltonian modified so as to need an overlap matrix **'
  overlap = use_overlap

  !preconditioninig is still needed in CG
  DO current_k=1,nks
     CALL init_k
     ALLOCATE( evc(npwx,nbnd), eig(nbnd) )
     CALL init_random_wfcs(npw,npwx,nbnd,evc)

     CALL start_clock('paro')
     !--- THIS IS THE RELEVANT CALL TO THE ROUTINE IN KS_Solvers/ParO --------------------------------------!
     CALL paro_gamma_new(cb_h_psi, cb_s_psi, cb_hs_psi, cb_g_1psi, overlap, &
               npwx, npw, nbnd, evc, eig, btype, ethr, notconv, nhpsi)

     avg_iter = avg_iter + nhpsi/float(nbnd)
     CALL stop_clock('paro')

     DEALLOCATE ( evc )

     IF (energy_shift .AND. current_k==1) ref=eig(4*ncell**3)
     CALL write_bands(eig,ref)
     WRITE (6,*) 'avg_iter, nhpsi, notconv, ethr ', avg_iter, nhpsi, notconv, ethr
     DEALLOCATE( eig )
  END DO

  CALL print_clock('paro')
  CALL print_clock('rotHSwg')          ; CALL print_clock('protHSwg')
  CALL print_clock('rotHSwg:hc')       ; CALL print_clock('protHSwg:hc')
  CALL print_clock('rotHSwg:diag')     ; CALL print_clock('protHSwg:diag')
  CALL print_clock('rotHSwg:evc')      ; CALL print_clock('protHSwg:evc')
  CALL print_clock('pcg')
  CALL print_clock('pcg:hs_1psi')
  CALL print_clock('pcg:ortho')
  CALL print_clock('s_psi')
! 
  write (6,*) 
  write (6,*) ' general FFT  routines'
  call print_clock('fftw')
  call print_clock('ffts')

#if defined(__MPI)
  CALL mp_global_end( )
#endif

END PROGRAM cb_paro_gamma_main
