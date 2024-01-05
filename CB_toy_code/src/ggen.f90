 subroutine ggen(gamma_only)
! define the lattice parameter alat and then tpiba, tpiba2 and gcutwfc
! define the direct lattice vectors at(icar,ivec)
! define the reciprocal lattice vectors bg(icar,ivec)
! define fft dimensions nr1,nr2,nr3 from gcutwfc
! generate g(3,ngm) vectors and their square lengths gg(ngm)
! sort the g-list w.r.t. increasing length
! define the index nl(ig) from the g-list to the fft grid
! compute npwx and allocate igk and ekin arrays

! global variables
  USE cb_module
  USE stick_base, only : sticks_map
  USE fft_types, only : fft_type_init
  USE fft_ggen, only : fft_set_nl
#if defined(__MPI)
  use mp_bands,             ONLY : intra_bgrp_comm,nyfft
#endif
  implicit none
! input variables
  logical :: gamma_only
! local variables
  real(DP), parameter :: AA_to_au = 1.d0/0.529177d0
  real(DP), allocatable :: g_aux(:,:)
  integer, allocatable :: ind(:)
  real(DP) :: g_(3), gg_, gcutrho
  integer :: ig, i, j, k, ik
  integer ::  nr1, nr2, nr3, nr1x, nr2x, nr3x ! to be clean up
  integer ::  m1, m2, mc
  TYPE(sticks_map) :: smap
  INTEGER :: nyfft_dummy = 1
#if defined(__MPI)
  logical :: lpara =.true.  ! whether the diagonalization is done in parallel
#else
  integer :: intra_bgrp_comm = 0 ! fake communicator required input to fft_type_init (to be fixed)
  logical :: lpara=.false.
#endif

  ! store input value of gamma_only in cb_module variable to know whether it's a Gamma-point calculation
  gamma_only_save = gamma_only
! define the lattice parameter alat and omega and then tpiba, tpiba2 and gcutwfc
  alat = cb_alat(icrystal) * AA_to_au * ncell ; omega = 0.25d0 * alat**3
  tpiba = tpi / alat
  tpiba2 = tpiba*tpiba
  gcutwfc = ecutwfc / tpiba2
  gcutrho = 4.d0*gcutwfc
! define the direct lattice vectors at(icar,ivec)
  at(1,1) = 0.0d0; at(2,1) = 0.5d0; at(3,1) = 0.5d0
  at(1,2) = 0.5d0; at(2,2) = 0.0d0; at(3,2) = 0.5d0
  at(1,3) = 0.5d0; at(2,3) = 0.5d0; at(3,3) = 0.0d0
! define the reciprocal lattice vectors bg(icar,ivec)
  bg(1,1) =-1.0d0; bg(2,1) = 1.0d0; bg(3,1) = 1.0d0
  bg(1,2) = 1.0d0; bg(2,2) =-1.0d0; bg(3,2) = 1.0d0
  bg(1,3) = 1.0d0; bg(2,3) = 1.0d0; bg(3,3) =-1.0d0
! define fft dimensions nr1,nr2,nr3 from gcutwfc
  nr1 = 2 * int ( sqrt(gcutrho) * sqrt(at(1,1)**2+ at(2,1)**2+at(3,1)**2) ) + 1 ; nr1x = nr1
  nr2 = 2 * int ( sqrt(gcutrho) * sqrt(at(1,2)**2+ at(2,2)**2+at(3,2)**2) ) + 1 ; nr2x = nr2
  nr3 = 2 * int ( sqrt(gcutrho) * sqrt(at(1,3)**2+ at(2,3)**2+at(3,3)**2) ) + 1 ; nr3x = nr3
  gkcut = gcutwfc
  do ik = 1, nks
     gkcut = max(gkcut, (sqrt(gcutwfc)+sqrt(sum(xk(:,ik)**2)))**2)
  end do
!emine
!"nyfft" groups (to push FFT parallelization beyond the nz-planes limit)
! number of y-fft groups. By default =1, i.e. y-ffts are done by a single proc
#if defined(__MPI)
  nyfft_dummy=nyfft
#else
  nyfft_dummy=1
#endif


!emine
!set the data structure for the fft arrays
  CALL fft_type_init( dfft, smap, "wave", gamma_only, lpara, intra_bgrp_comm, at, bg, gkcut, gcutrho/gkcut , nyfft=nyfft_dummy)

  write (stdout,*) dfft%nr1, dfft%nr2,dfft%nr3, dfft%nnr
  nr1  = dfft%nr1  ; nr2  = dfft%nr2  ; nr3  = dfft%nr3  ; nnr = dfft%nnr 
  nr1x = dfft%nr1x ; nr2x = dfft%nr2x ; nr3x = dfft%nr3x

  ! define the clock labels ( this enables the corresponding fft too ! )
  dfft%rho_clock_label='ffts' ; dfft%wave_clock_label='fftw'




! generate g(3,ngm) vectors and their square lengths gg(ngm)
 !> count how many vectors are inside the 4*ecutwfc sphere and set ngm
  ngm = 0
  do i = -nr1/2,nr1/2
     if (gamma_only .and. i<0 ) cycle 
     do j = -nr2/2,nr2/2
        if (gamma_only .and. i==0 .and. j<0 ) cycle 
#if defined (__MPI)
 !> in the parallel case only take into account the z-columns assigne to this processor
        m1 = mod (i, nr1) + 1 ; IF (m1 < 1) m1 = m1 + nr1
        m2 = mod (j, nr2) + 1 ; IF (m2 < 1) m2 = m2 + nr2
        mc = m1 + (m2 - 1) * nr1x ; IF ( dfft%isind ( mc ) == 0) CYCLE 
#endif
        do k = -nr3/2,nr3/2
           if (gamma_only .and. i==0 .and. j==0 .and. k<0 ) cycle 
           g_(:)= bg(:,1)*i+bg(:,2)*j+bg(:,3)*k
           gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
           if (gg_ < 4.d0*gcutwfc) ngm = ngm + 1
        end do
     end do
  end do
  write (stdout,*) 'ngm =', ngm
 !> allocate g(3,ngm) and gg(ngm)
  allocate (g(3,ngm), gg(ngm))
 !> fill the arrays as needed
  ig = 0
  do i = -nr1/2,nr1/2
     if (gamma_only .and. i<0 ) cycle 
     do j = -nr2/2,nr2/2
        if (gamma_only .and. i==0 .and. j<0 ) cycle 
#if defined (__MPI)
 !> in the parallel case only take into account the z-columns assigne to this processor
        m1 = mod (i, nr1) + 1 ; IF (m1 < 1) m1 = m1 + nr1
        m2 = mod (j, nr2) + 1 ; IF (m2 < 1) m2 = m2 + nr2
        mc = m1 + (m2 - 1) * nr1x ; IF ( dfft%isind ( mc ) == 0) CYCLE 
#endif
        do k = -nr3/2,nr3/2
           if (gamma_only .and. i==0 .and. j==0 .and. k<0 ) cycle 
           g_(:)= bg(:,1)*i+bg(:,2)*j+bg(:,3)*k
           gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
           if (gg_ < 4.d0*gcutwfc) then
              ig = ig + 1; g(:,ig) = g_(:) ; gg(ig) = gg_
           end if
        end do
     end do
  end do
  if (ig.ne.ngm) stop ' something wrong in ggen'
! sort the g-list w.r.t. increasing length
  allocate ( ind(ngm), g_aux(3,ngm) ) ; ind(:) = 0 
  call hpsort_eps( ngm, gg, ind, eps8 ) ; g_aux = g
  do ig=1,ngm
     g(:,ig) = g_aux(:,ind(ig))
  end do
  if ( gg(1) < eps8 ) gstart = 2
  deallocate (ind, g_aux)

! define the index nl(ig) from the g-list to the fft grid
  CALL fft_set_nl( dfft, at, g )

! compute npwx and allocate igk and ekin arrays
  npwx = 0
  do current_k = 1, nks
     npw = 0
     do ig = 1, ngm
        g_(:) = xk(:,current_k) + g(:,ig)
        gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
        if (gg_ < gcutwfc) npw = npw + 1
     end do
     npwx = max (npw, npwx)
  end do
  write (stdout,*) 'npwx =', npwx
  allocate ( igk(npwx), ekin(npwx), aux(npwx) )

 end subroutine ggen
