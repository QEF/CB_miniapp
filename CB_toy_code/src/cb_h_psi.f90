
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE cb_h_psi( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  ! --- Wrapper routine: performs bgrp parallelization on non-distributed bands
  ! --- if suitable and required, calls old H\psi routine as cb_h_psi_
  !
  USE cb_module,        ONLY : DP
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda,m)
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda,m)
  !
  INTEGER     :: m_start, m_end
  !
  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. there is more than one band, otherwise there is nothing to parallelize
  !
  call start_clock('h_psi')
  IF (use_bgrp_in_hpsi .AND. m > 1) THEN
     ! use band parallelization here
     CALL divide(inter_bgrp_comm,m,m_start,m_end)
     hpsi(:,:) = (0.d0,0.d0) ! array must be initializad to zero to prevent adding garbage to others work
     ! call the routine if there at least one band in this band group
     IF (m_end >= m_start) CALL cb_h_psi_( lda, n, m_end-m_start+1, psi(1,m_start), hpsi(1,m_start) )
     CALL mp_sum(hpsi,inter_bgrp_comm) ! collect the result across the band group partners
  ELSE
     ! don't use band parallelization here
     CALL cb_h_psi_( lda, n, m, psi, hpsi )
  END IF
  call stop_clock('h_psi')

  RETURN
  !
END SUBROUTINE cb_h_psi
!
 subroutine cb_h_psi_(npwx,npw,nvec,psi,hpsi)
! for each input wfc
! apply the kinetic energy to each wave function
! bring the wavefunction to real space
! compute the potential energy contribution
! back to reciprocal space and add to the kinetic term

! global variables
  USE cb_module, only : DP, gstart, gamma_only_save
  USE cb_module, only : ekin, aux, fft_array, vloc, igk, dfft
  USE cb_module, only : use_overlap
  USE fft_interfaces, only : fwfft, invfft
  implicit none
! input variables
  integer, INTENT(IN) :: npwx, npw, nvec
  complex(DP), INTENT(IN) :: psi(npwx,nvec)
  complex(DP), INTENT(OUT):: hpsi(npwx,nvec)
! local variables
  integer :: ivec, ig

  do ivec=1,nvec
!$omp parallel 
!$omp workshare
    aux(:) = psi(:,ivec)
!$omp end workshare nowait
    if (use_overlap) then
!$omp workshare
       aux(1:npw) = (1.d0 + exp(-ekin(1:npw)))*psi(1:npw,ivec)
!$omp end workshare nowait
   end if

! initialize fft_array 
!$omp workshare
    fft_array(:) = CMPLX(0.d0,0.d0)
!$omp end workshare nowait

! apply the kinetic energy to each wave function
!$omp do
    do ig=1,npw
       hpsi(ig,ivec) = ekin(ig) * aux(ig)
    end do
!$omp end do nowait
! bring the wavefunction to real space
!$omp do
    do ig = 1, npw
       fft_array(dfft%nl(igk(ig))) = aux(ig)
    end do
!$omp end do nowait
    if (gamma_only_save) then
!$omp do
       do ig = gstart, npw
          fft_array(dfft%nlm(igk(ig))) = CONJG(aux(ig))
       end do
!$omp end do nowait
!$omp single
       if (gstart==2) fft_array(dfft%nl(igk(1))) = CMPLX(REAL(aux(1)),0.d0,kind=DP)
!$omp end single
    end if
!$omp end parallel
    call invfft('Wave',fft_array,dfft) !  >>> FFT G->R <<<
! compute the potential energy contribution
    fft_array(:) = vloc(:) * fft_array(:)
! back to reciprocal space and add to the kinetic term
    call fwfft('Wave',fft_array,dfft)  !  >>> FFT R->G <<<
!$omp parallel 
!$omp do
    do ig = 1, npw
       hpsi(ig,ivec) = hpsi(ig,ivec) + fft_array(dfft%nl(igk(ig)))
    end do
!$omp end do
    if (use_overlap) then
!$omp workshare
       hpsi(1:npw,ivec) = (1.d0 + exp(-ekin(1:npw)))*hpsi(1:npw,ivec)
!$omp end workshare nowait
    end if
!$omp end parallel 
  end do
  if (gamma_only_save .and. gstart==2) hpsi(1,1:nvec) = CMPLX(REAL(hpsi(1,1:nvec)),0.d0,kind=DP)

 end subroutine cb_h_psi_
!----------------------------------------------------------------------------
SUBROUTINE cb_hs_1psi( lda, n, psi, hpsi, spsi )
!----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian and the S matrix
  ! ... to a vector psi and puts the result in hpsi and spsi
  ! ... Wrapper routine - calls h_psi and s_psi
!
  USE cb_module, only : DP
  implicit none
! input variables
  INTEGER, INTENT(IN)      :: lda, n
  complex(DP), INTENT(IN) :: psi(lda)
  complex(DP), INTENT(OUT):: hpsi(lda)
  complex(DP), INTENT(OUT):: spsi(lda)
! local variables
  integer,parameter :: nvec=1
  call start_clock('hs_1psi')
  call cb_h_psi_(lda,n,nvec,psi,hpsi)
  call cb_s_psi_(lda,n,nvec,psi,spsi)
  call stop_clock('hs_1psi')

END SUBROUTINE cb_hs_1psi

!----------------------------------------------------------------------------
SUBROUTINE cb_hs_psi( lda, n, m, psi, hpsi, spsi )
!----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian and the S matrix
  ! ... to a vector psi and puts the result in hpsi and spsi
  ! ... Wrapper routine - calls h_psi and s_psi
!
  USE cb_module, only : DP
  implicit none
! input variables
  INTEGER, INTENT(IN)      :: lda, n, m
  complex(DP), INTENT(IN) :: psi(lda,m)
  complex(DP), INTENT(OUT):: hpsi(lda,m)
  complex(DP), INTENT(OUT):: spsi(lda,m)
! local variables
  call start_clock('hs_psi')
  call cb_h_psi_(lda,n,m,psi,hpsi)
  call cb_s_psi_(lda,n,m,psi,spsi)
  call stop_clock('hs_psi')

END SUBROUTINE cb_hs_psi




SUBROUTINE cb_h_1psi( lda, n, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian 
  ! ... to a vector psi and puts the result in hpsi
  ! ... Wrapper routine - calls h_psi 
  !
  ! ... No bgrp parallelization here !
  !
  USE cb_module,  ONLY: DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n
  COMPLEX (DP), INTENT(IN) :: psi(lda)
  COMPLEX (DP), INTENT(OUT) ::  hpsi(n)!
  !
  CALL start_clock( 'h_1psi' )
  CALL cb_h_psi( lda, n, 1, psi, hpsi ) ! apply H to a single wfc (no bgrp parallelization here)
  CALL stop_clock( 'h_1psi' )
  !
  RETURN
  !
END SUBROUTINE cb_h_1psi
