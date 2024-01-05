
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE cb_s_psi( lda, n, m, psi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the S
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    spsi  S*psi
  !
  ! --- Wrapper routine: performs bgrp parallelization on non-distributed bands
  ! --- if suitable and required, calls old S\psi routine as cb_s_psi_
  !
  USE cb_module,        ONLY : DP
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda,m)
  COMPLEX(DP), INTENT(OUT) :: spsi(lda,m)
  !
  INTEGER     :: m_start, m_end
  !
  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. there is more than one band, otherwise there is nothing to parallelize
  !
  call start_clock('s_psi')
  IF (use_bgrp_in_hpsi .AND. m > 1) THEN
     ! use band parallelization here
     CALL divide(inter_bgrp_comm,m,m_start,m_end)
     spsi(:,:) = (0.d0,0.d0) ! array must be initializad to zero to prevent adding garbage to others work
     ! call the routine if there at least one band in this band group
     IF (m_end >= m_start) CALL cb_s_psi_( lda, n, m_end-m_start+1, psi(1,m_start), spsi(1,m_start) )
     CALL mp_sum(spsi,inter_bgrp_comm) ! collect the result across the band group partners
  ELSE
     ! don't use band parallelization here
     CALL cb_s_psi_( lda, n, m, psi, spsi )
  END IF
  call stop_clock('s_psi')

  RETURN
  !
END SUBROUTINE cb_s_psi
!----------------------------------------------------------------------------
 subroutine cb_s_psi_(npwx,npw,nvec,psi,spsi)
!----------------------------------------------------------------------------
! for each input wfc spsi is just psi

! global variables
  USE cb_module, only : DP
  USE cb_module, only : use_overlap, ekin
  implicit none
! input variables
  integer, intent(IN) :: npwx, npw, nvec
  complex(DP),intent(IN) :: psi(npwx,nvec)
  complex(DP),intent(OUT) :: spsi(npwx,nvec)
! local variables
  integer :: ivec

! for each input wfc spsi is just psi
  if (.not.use_overlap) then
     spsi(:,:) = psi(:,:)
  else
     do ivec=1,nvec
        spsi(1:npw,ivec) = (1.d0 + exp(-ekin(1:npw)))**2 * psi(1:npw,ivec)
     end do
  end if

 end subroutine cb_s_psi_

!----------------------------------------------------------------------------
 subroutine cb_s_1psi( lda, n, psi, spsi )
!----------------------------------------------------------------------------
! global variables
  USE cb_module, only : DP
  implicit none
! input variables
  INTEGER, INTENT(IN)      :: lda, n
  complex(DP), INTENT(IN) :: psi(lda)
  complex(DP), INTENT(OUT):: spsi(lda)
! local variables
  integer,parameter :: nvec=1
  call start_clock('s_1psi')
  call cb_s_psi_(lda,n,nvec,psi,spsi)
  call stop_clock('s_1psi')

 end subroutine cb_s_1psi

