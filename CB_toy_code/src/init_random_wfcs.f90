 subroutine init_random_wfcs(npw,npwx,nbnd,evc)

! global variables
  use cb_module, only: DP
  use cb_module, only: ekin, gstart
  use random_numbers, only: randy
  implicit none
! input variables
  integer, intent(in) :: npw, npwx, nbnd
  complex(DP),intent(out) :: evc(npwx,nbnd)
! local variables
  integer :: ibnd, ig
  real(DP) :: alpha = 1.0d0

  evc = (0.d0,0.d0)
  do ibnd =1, nbnd
     do ig = 1, npw
        evc(ig,ibnd) = exp(-alpha*ekin(ig))*CMPLX(randy(),randy())
     end do
  end do
  if (gstart==2) evc(1,1:nbnd) = CMPLX ( REAL ( evc(1,1:nbnd) ), 0.d0, kind=DP)

 end subroutine init_random_wfcs
