 subroutine cb_g_psi(npwx, npw, nvec, npol, psi, eig)

! global variables
  USE cb_module, ONLY : DP
  USE cb_module, ONLY : ekin
  implicit none
! input variables
  integer, intent(IN) :: npwx, npw, nvec, npol
  complex(DP), intent(INOUT) :: psi(npwx,nvec)
  real(DP), intent(IN) :: eig(nvec)
! local variables
  integer :: ivec, ig
  real(DP) :: x, denm
  
  call start_clock('g_psi')
  do ivec = 1, nvec
     do ig = 1, npw
        x = (ekin(ig) - eig(ivec))
        denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))
        psi (ig, ivec) = psi (ig, ivec) / denm
     enddo
  enddo
  call stop_clock('g_psi')

 end subroutine cb_g_psi

 subroutine cb_g_1psi(npwx, npw, psi, eig)

! global variables
  USE cb_module, ONLY : DP
  USE cb_module, ONLY : ekin
  implicit none
! input variables
  integer, intent(IN) :: npwx, npw
  complex(DP), intent(INOUT) :: psi(npwx)
  real(DP), intent(IN) :: eig
! local variables
  integer :: ivec, ig
  real(DP) :: x, denm
  
  call start_clock('g_1psi')
  do ig = 1, npw
     x = (ekin(ig) - eig)
     denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))
     psi (ig) = psi (ig) / denm
  enddo
  call stop_clock('g_1psi')

 end subroutine cb_g_1psi
