 subroutine set_cb_potential
! allocate the fft grid for the potential and initialize to 0
! fill the non-vanishing fourier components according to the CB definition
! fourier transform the potential to real space

! global variables
  USE cb_module
  USE fft_interfaces, only : invfft
  implicit none
! local variables
  real(DP) :: arg
  integer :: ig, nc2
  nc2 = ncell**2

! allocate the fft grid for the potential and initialize to 0
allocate ( vloc(dfft%nnr), fft_array(dfft%nnr) ) ; fft_array(:) =CMPLX(0.d0,0.d0)
! fill the non-vanishing fourier components according to the CB definition
  do ig=1,ngm 
     if (ncell.ne.1) then
        arg = ( at(1,1)*g(1,ig) + at(2,1)*g(2,ig) + at(3,1)*g(3,ig) )/ncell       
        if (abs(arg-nint(arg))>eps8) cycle
        arg = ( at(1,2)*g(1,ig) + at(2,2)*g(2,ig) + at(3,2)*g(3,ig) )/ncell  
        if (abs(arg-nint(arg))>eps8) cycle
        arg = ( at(1,3)*g(1,ig) + at(2,3)*g(2,ig) + at(3,3)*g(3,ig) )/ncell 
        if (abs(arg-nint(arg))>eps8) cycle
     end if
     arg = -tpi * 0.125d0 * (g(1,ig) + g(2,ig) + g(3,ig)) /  ncell
     if (abs (gg(ig)- 3.d0*nc2) < eps8) fft_array(dfft%nl(ig)) = CMPLX(vs3(icrystal)*cos(arg),va3(icrystal)*sin(arg))
     if (abs (gg(ig)- 4.d0*nc2) < eps8) fft_array(dfft%nl(ig)) = CMPLX(0.d0,va4(icrystal)*sin(arg))
     if (abs (gg(ig)- 8.d0*nc2) < eps8) fft_array(dfft%nl(ig)) = CMPLX(vs8(icrystal)*cos(arg),0.d0)
     if (abs (gg(ig)-11.d0*nc2) < eps8) fft_array(dfft%nl(ig)) = CMPLX(vs11(icrystal)*cos(arg),va11(icrystal)*sin(arg))
     if (gamma_only_save) fft_array(dfft%nlm(ig)) = CONJG(fft_array(dfft%nl(ig)))
     if (gg(ig) > 11.d0*nc2 + eps8 ) exit
  end do
! fourier transform the potential to real space >>>> FFT G->R <<<< 

!  call invfft ('Smooth',fft_array,dfft)
!emine: now when you send sth to invfft it goes to invfft_y
!which goes to fft_fwinv.f90 new interface
!where we have rho or wave rahter than smooth, dense, wave etc.
!Currently in setlocal.f90 of QE, local potential is brought to real space with 
!rho grid
  call invfft ('Rho',fft_array,dfft)

   vloc(:) = REAL(fft_array(:))

  end subroutine set_cb_potential
