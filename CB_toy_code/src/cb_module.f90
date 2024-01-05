module cb_module

 USE fft_types,  ONLY : fft_type_descriptor
 implicit none 
 save

! ... kind definitions
 integer, parameter :: DP = selected_real_kind(14,200)
! ... numerical parameters
 real(DP) :: eps8 = 1.d-8
 real(DP) :: tpi  = 2.d0 * 3.14159265358979323846_DP
!...  I/O definitions
 integer :: stdout = 6    ! unit connected to standard output

 integer, parameter :: ncrystal = 15 
 character(len=4) :: crystal(ncrystal)
 real(DP) :: cb_alat(ncrystal)
 real(DP) :: vs3(ncrystal), vs8(ncrystal),vs11(ncrystal)
 real(DP) :: va3(ncrystal), va4(ncrystal),va11(ncrystal)

 data crystal / 'free', 'Si',  'Ge',  'Sn', 'GaP','GaAs','AlSb', 'InP','GaSb','InAs','InSb', 'ZnS','ZnSe','ZnTe','CdTe' /
 data cb_alat /  5.00,  5.43,  5.66,  6.49,  5.44,  5.64,  6.13,  5.86,  6.12,  6.04,  6.48,  5.41,  5.65,  6.07,  6.41 /
 data vs3     /  0.00, -0.21, -0.23, -0.20, -0.22, -0.23, -0.21, -0.23, -0.22, -0.22, -0.20, -0.22, -0.23, -0.22, -0.20 /
 data vs8     /  0.00, +0.04, +0.01,  0.00, +0.03, +0.01, +0.02, +0.01,  0.00,  0.00,  0.00, +0.03, +0.01,  0.00,  0.00 /
 data vs11    /  0.00, +0.08, +0.06, +0.04, +0.07, +0.06, +0.06, +0.06, +0.05, +0.05, +0.04, +0.07, +0.06, +0.05, +0.04 /
 data va3     /  0.00,  0.00,  0.00,  0.00, +0.12, +0.07, +0.06, +0.07, +0.06, +0.08, +0.06, +0.24, +0.18, +0.13, +0.15 /
 data va4     /  0.00,  0.00,  0.00,  0.00, +0.07, +0.05, +0.04, +0.05, +0.05, +0.05, +0.05, +0.14, +0.12, +0.10, +0.09 /
 data va11    /  0.00,  0.00,  0.00,  0.00, +0.02, +0.02, +0.02, +0.01, +0.01, +0.03, +0.01, +0.04, +0.03, +0.01, +0.04 /

 INTEGER  :: icrystal  ! index of the considered crystal structure
 INTEGER  :: ncell     ! number of units in each direction of the crystal supercell 
 logical  :: energy_shift ! whether in the printout VB top is shifted to zero (assuming the first k-point contains it)
 REAL(DP) :: ecutwfc   ! kinetic energy cutoff (Ry units)
 REAL(DP) :: gcutwfc   ! kinetic energy cutoff (tpiba2 units)
 INTEGER  :: current_k ! index of the current k-point 
 INTEGER  :: nks       ! number of k-points
 REAL(DP) :: alat      ! lattice parameter (in a.u.) of the computed cell
 REAL(DP) :: omega     ! volume of the crystal unit cell
 REAL(DP) :: tpiba     !  2 pi / alat
 REAL(DP) :: tpiba2    ! (2 pi / alat )**2
 REAL(DP) :: at(3,3)   ! direct lattice basis vectors
 REAL(DP) :: bg(3,3)   ! reciprocal lattice basis vectors
 REAL(DP) :: gkcut     ! maximum value of k-point length
 TYPE(fft_type_descriptor) :: dfft ! fft grid descriptor
 INTEGER  :: nnr            ! fft global array dimension
 INTEGER  :: ngm       ! number of G-vectors within the specified cutoff
 INTEGER  :: gstart=1  ! index of the first non-zero G vector (may be 2 if G=0 is locally present). used for gamma_only
 LOGICAL  :: gamma_only_save ! .T. if it's a Gamma-point  calculation 
 REAl(DP), ALLOCATABLE :: g(:,:) ! reciprocal lattice vector list
 REAl(DP), ALLOCATABLE :: gg(:)  ! reciprocal lattice vector square moduli
 INTEGER, ALLOCATABLE  :: igk(:) ! index from the |k+g| ordered list to the |g|-ordered list
 REAl(DP), ALLOCATABLE :: xk(:,:)! k-point vector list
 REAl(DP), ALLOCATABLE :: ekin(:)! kinetic energy of the |k+g| ordered list
 REAl(DP), ALLOCATABLE :: vloc(:)! effective potential on the fft grid
 COMPLEX(DP), ALLOCATABLE :: fft_array(:) ! auxiliary array for the fft operations
 COMPLEX(DP), ALLOCATABLE :: aux(:) ! auxiliary array to be used in h_psi
  
! diagonalization related parameters
 INTEGER  :: david     ! davidson inflation limit
 real(DP) :: ethr      ! eigenvalue convergence threshold
 INTEGER  :: npw       ! number of plane waves for the current k-point
 INTEGER  :: npwx      ! maximum allowed number of plane waves 
 INTEGER  :: nbnd      ! number of bands to be calculated
 INTEGER  :: nbndx     ! maximum allowed number of bands to be calculated
 INTEGER, ALLOCATABLE  :: btype(:)  ! index specifying occupied vs empty bands

 LOGICAL  :: use_overlap ! .T. if the hamiltonian is modified so as to need an overlap matrix. For testing purposes.

end module cb_module
