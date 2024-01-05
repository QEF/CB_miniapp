 subroutine init_k
! compute the norm of the k+g vectors and fill a kinetic energy array
! sort it dragging around the corresponding igk index

! global variables
  USE cb_module
  implicit none
! local variables
  real(DP) :: kpg(3), kpg2
  integer :: ig

! compute the norm of the k+g vectors and fill a kinetic energy array
  npw = 0
  do ig = 1, ngm
     kpg(:) = xk(:,current_k) + g(:,ig)    
     kpg2 = kpg(1)*kpg(1) + kpg(2)*kpg(2) + kpg(3)*kpg(3)
     if (kpg2 < gcutwfc) then
        npw = npw + 1 ; if (npw>npwx) stop 'something wrong in init_igk'
        igk(npw) = ig
        ekin(npw) = kpg2 * tpiba2
     end if
  end do
! sort it dragging around the corresponding igk index
  call hpsort_eps( npw, ekin, igk, eps8 ) 

 end subroutine init_k
