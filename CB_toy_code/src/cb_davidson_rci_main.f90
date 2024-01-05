program cb_davidson_rci_main

! global variables
  USE cb_module
  use david_rci_m
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
  logical,parameter :: gamma_only=.false.  ! general k-point version.
  complex(DP), allocatable :: evc(:,:)
  real(dp), allocatable :: eig(:)
  real(dp) :: ref=0.d0
  integer, parameter :: npol=1
  integer :: notcnv, dav_iter
  logical :: overlap = .false. , lrot =.false.
  integer :: task
  type(david_rci_work_t) :: work
#if defined(__MPI)
! local paralelization variables  
  integer :: ndiag     ! input value of processors in the diagonalization group
  logical :: do_distr_diag_in_band_group = .true. ! whether or not the parallel diagonalization is performed inside the
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
!---

  call input(gamma_only)
  call ggen(gamma_only)
  call set_cb_potential

  do current_k=1,nks
    call init_k
    call david_rci_work_alloc(npw, npwx, nbnd, nbndx, npol, overlap, work)
    allocate( evc(npwx,nbnd), eig(nbnd) )
    call init_random_wfcs(npw,npwx,nbnd,evc)
    task = -1
    !---
    do while (task /= ESL_TASK_EXIT)

      call david_rci_run(npw, npwx, nbnd, nbndx, npol, evc, ethr, overlap, eig, btype, notcnv, lrot, dav_iter, work, task)

      if (iand(task, ESL_TASK_HPSI) /= 0) then
        call cb_h_psi(npwx, npw, work%ivec_size, work%psi(1,1,work%ivec_start), work%hpsi(1,1,work%ivec_start))
      end if

      if (iand(task, ESL_TASK_SPSI) /= 0) then
        call cb_s_psi(npwx, npw, work%ivec_size, work%psi(1,1,work%ivec_start), work%spsi(1,1,work%ivec_start))
      end if

      if (iand(task, ESL_TASK_GPSI) /= 0) then
        call cb_g_psi(npwx, npw, work%ivec_size, npol, work%psi(1,1,work%ivec_start), work%ew(work%ivec_start))
      end if
      
    end do
    !---
    call david_rci_work_free(work)
    deallocate( evc )

    if (energy_shift .and. current_k==1) ref=eig(4*ncell**3)
    call write_bands(eig,ref)
    write (stdout,*) 'dav_iter, notcnv, ethr ',dav_iter,notcnv, ethr

    deallocate( eig )
  end do

#if defined(__MPI)
  call mp_global_end( )
#endif

end program cb_davidson_rci_main
