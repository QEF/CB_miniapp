#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

if test $# = 0
then
    dirs=" LAXlib FFTXlib UtilXlib clib \
           KS_Solvers/Davidson KS_Solvers/Davidson_RCI KS_Solvers/CG KS_Solvers/PPCG \
           KS_Solvers/ParO  KS_Solvers/DENSE  \
	   CB_toy_code/src" 
          
elif
    test $1 = "-addson" 
then
    echo "The script for adding new dependencies is running"
    echo "Usage: $0 -addson DIR DEPENDENCY_DIRS"
    echo "$0 assumes that the new dependencies are in $TOPDIR/../"
#    ninput=$#
#    echo "number of input arguments: $ninput"
    dirs=$2
    shift
    shift
    add_deps=$*
    echo "dependencies in $add_deps will be searched for $dirs"
else
    dirs=$*
fi

for dir in $dirs; do

    # the following command removes a trailing slash
    DIR=`echo ${dir%/}`

    # the following would also work
    #DIR=`echo $dir | sed "s,/$,,"`

    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    LEVEL1=..
    LEVEL2=../..
    # default
    DEPENDS="$LEVEL1/include" 
    # for convenience, used later
    DEPEND1="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/FFTXlib $LEVEL1/LAXlib $LEVEL1/UtilXlib"
    DEPEND2="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/FFTXlib $LEVEL2/LAXlib $LEVEL2/UtilXlib" 
    DEPEND3="$LEVEL2/include $LEVEL2/FFTXlib $LEVEL2/LAXlib $LEVEL2/UtilXlib"
    case $DIR in 
        LAXlib )
             DEPENDS="$LEVEL1/UtilXlib " ;;
	CB_toy_code/src )
	     DEPENDS="$DEPEND2 ../../KS_Solvers/Davidson ../../KS_Solvers/CG ../../KS_Solvers/PPCG ../../KS_Solvers/ParO ../../KS_Solvers/DENSE " ;;
	KS_Solvers/Davidson | KS_Solvers/Davidson_RCI | KS_Solvers/CG | KS_Solvers/PPCG | KS_Solvers/ParO | KS_Solvers/DENSE )
	     DEPENDS="$DEPEND3" ;;
    *)
# if addson needs a make.depend file
	DEPENDS="$DEPENDS $add_deps"

    esac

    # generate dependencies file (only for directories that are present)
    if test -d $TOPDIR/../$DIR
    then
	cd $TOPDIR/../$DIR
       
	$TOPDIR/moduledep.sh $DEPENDS > make.depend
	$TOPDIR/includedep.sh $DEPENDS >> make.depend

        # handle special cases: modules for C-fortran binding,
        #   	                hdf5, MPI, FoX, libxc
        sed '/@iso_c_binding@/d' make.depend > make.depend.tmp
        sed '/@hdf5@/d;/@mpi@/d' make.depend.tmp > make.depend
        sed '/@fox_dom@/d;/@fox_wxml@/d'  make.depend > make.depend.tmp 
        sed '/@m_common_io@/d;/@xc_f03_lib_m@/d' make.depend.tmp > make.depend

        if test "$DIR" = "FFTXlib"
        then
            sed '/@mkl_dfti/d' make.depend > make.depend.tmp
            sed '/@fftw3.f/d;s/@fftw.c@/fftw.c/' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "LAXlib"
        then
            sed '/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@cudafor@/d' make.depend.tmp > make.depend
            sed '/@gbuffers@/d' make.depend > make.depend.tmp
            sed '/@cusolverdn@/d' make.depend.tmp > make.depend
            sed '/@zhegvdx_gpu@/d' make.depend > make.depend.tmp
            sed '/@dsygvdx_gpu@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "UtilXlib"
        then
            sed '/@ifcore@/d' make.depend > make.depend.tmp
            sed '/@cudafor@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "CB_toy_code/src" 
        then
            sed '/@environ_/d;/@solvent_tddfpt@/d' make.depend > make.depend.tmp
            sed '/@david_rci_m@/d' make.depend.tmp > make.depend
        fi

        rm -f make.depend.tmp

        # check for missing dependencies 
        if grep @ make.depend
        then
	   notfound=1
	   echo WARNING: dependencies not found in directory $DIR
       else
           echo directory $DIR : ok
       fi
    else
       echo directory $DIR : not present in $TOPDIR 
    fi
done
if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
