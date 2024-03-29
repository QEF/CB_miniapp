# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_FFT], [

have_fft=0
have_fft_include=0

  AC_MSG_CHECKING([FFT])
 
# check for FFT libraries (no check for explicit openmp)
# supported vendor replacements:
#   essl on IBM linux machines
#   ASL/Mathkeisan on Nec
#   acml on amd
if test "$fft_libs" = "" && test "$use_openmp" -eq 0
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch" in
        ppc64 | ppc64-mn )
                # check for essl
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs"
                AC_SEARCH_LIBS(dcft, essl, have_fft=1 fft_libs="$LIBS")
            ;;
        ppc64-bg | ppc64-bgq )
                # check for esslbg
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs $blas_libs"
                AC_SEARCH_LIBS(dcft, esslbg, have_fft=1 fft_libs="$LIBS")
            ;;
        necsx )
                if test "$use_fft_mathkeisan" -ne 0
                then
                   #sx5-nec or sx6-nec or sx8-nec: check in (/SX)/usr/lib
                   #sx8-nec-idris: check in /SX/opt/mathkeisan/inst/lib0
                   try_libdirs="/SX/usr/lib /SX/opt/mathkeisan/inst/lib0"
                   #check for Mathkeisan (Cray simple precision )
                   #search for initialization subroutine
                   echo $ECHO_N  "Searching in Mathkeisan" $ECHO_C
                   for dir in none $try_libdirs
                   do
                        unset ac_cv_search_zftfax # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(zftfax, fft, have_fft=1
                             try_dflags="$try_dflags try_dflags_fft_mathkeisan"
                                       fft_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_zftfax" != "no"
                        then break ; fi
                   done
                fi
                if test "$use_fft_asl" -ne 0
                then
                   #check for asl in (/SX)/usr/lib
                   try_libdirs="/SX/usr/lib"
                   #search for initialization subroutine
                   echo $ECHO_N  "Searching in Asl" $ECHO_C
                   for dir in none $try_libdirs
                   do
                        unset ac_cv_search_zfc3cl # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(zfc3cl, asl, have_fft=1
                             asl_libs="$try_loption $LIBS"
                             try_dflags="$try_dflags $try_dflags_fft_asl"
                             fft_libs="$fft_libs $asl_libs")
                        if test "$ac_cv_search_zfc3cl" != "no"
                        then break ; fi
                   done
                fi
                if test "$use_fft_para" -ne 0
                then
                   try_dflags="$try_dflags $try_dflags_fft_para"
                fi
                ;;
        esac

fi

if test "$have_fft" -eq 0 
then
   # if MKL have been found, select dfti (no need to recheck)
   if test "$have_mkl" -eq 1
   then
      try_dflags="$try_dflags -D__DFTI"
      # If not set on input, MKLROOT was set when checking blas
      try_iflags="$try_iflags -I$MKLROOT/include"
      have_fft=1
   else 
      if test "$have_armpl" -eq 1; then 
      have_fft=1 
      fi 
   fi

   # Check fftw v3 (both native and MKL interfaces)

   if test "$have_fft" -eq 0
   then
	  try_libdirs="/usr/local/lib"
          try_libdirs="$libdirs $try_libdirs $ld_library_path "
          for dir in none $try_libdirs
          do
                  unset ac_cv_search_dfftw_execute_dft # clear cached value
                  if test "$dir" = "none"
                  then
                        try_loption=
                  else
                        echo $ECHO_N "in $dir: " $ECHO_C
                        try_loption="-L$dir"
                  fi

                  CFLAGS="$test_cflags"
                  CPPFLAGS="$test_cppflags"
                  LDFLAGS=" $test_ldflags $try_loption"
                  LIBS="$fft_libs"

                  if test "$use_openmp" -eq 1
                  then
                    # Try testing openmp without -lfftw3 first, if that fails then return
                    # to previous behaviour
                    AC_SEARCH_LIBS(dfftw_execute_dft, fftw3_omp, have_fft=1
                               fft_libs="$try_loption $LIBS", , -lm)
                    if test "$have_fft" -eq 0
                    then
                      AC_SEARCH_LIBS(dfftw_execute_dft, fftw3_omp, have_fft=1
                                 fft_libs="$try_loption $LIBS -lfftw3", , -lfftw3 -lm)
                    fi 
                  else
                    AC_SEARCH_LIBS(dfftw_execute_dft, fftw3, have_fft=1
                               fft_libs="$try_loption $LIBS", , -lm)
                  fi

                  if test "$have_fft" -eq 1
                  then
                        try_dflags="$try_dflags -D__FFTW3"
                        try_incdir="$FFTW_INCLUDE $FFTW_INC $INCLUDE_PATH $CPATH $FPATH"
                        for inc in $try_incdir
                        do
                           AC_COMPILE_IFELSE([include "fftw3.f03"],have_fft_include=1,)
                           if test "$have_fft_include" -eq 1
                           then
                             try_iflags="$try_iflags -I$inc"
                             break
                           fi
                        done
                        break
                  fi

          done
   fi
fi

  AC_MSG_RESULT(${fft_libs})

fft_line="FFT_LIBS=$fft_libs"

# if no valid FFT library was found, use the local copy
if test "$have_fft" -eq 0
then
        case "$arch" in
        ppc64-bg | ppc64-bgq )
            try_dflags="$try_dflags -D__LINUX_ESSL"
            ;;
        * )
	    try_dflags="$try_dflags -D__FFTW"
            ;;
        esac
fi

AC_SUBST(fft_libs)
AC_SUBST(fft_line)
  
  ]
)
