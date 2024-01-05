 subroutine write_bands(eig,ref)
! write the kpoint coordinates and the number of plane waves
! write the band energies

! global variables
  USE cb_module
  implicit none
! input variables
  real(DP), intent(IN) :: eig(nbnd), ref
! local variables
  real(DP), parameter :: rytoev = 13.6056d0
  integer :: i, ibnd

! write the kpoint coordinates and the number of plane waves
  write( 6, 9021 ) ( xk(i,current_k), i=1,3 ), npw
9021 FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   bands (ev):'/ )
! write the band energies
  write( 6, 9030 ) ( (eig(ibnd)-ref) * rytoev, ibnd = 1, nbnd )
!  write( 6, 9030 ) ( eig(ibnd), ibnd = 1, nbnd )
9030 FORMAT( '  ',8F9.4 )

 end subroutine write_bands
 
