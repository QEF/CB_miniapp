!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_global
  !----------------------------------------------------------------------------
  !
  ! ... The routine mp_startup initializing MPI, plus the routine mp_global_end stopping MPI.
  ! ... Do not use this module to reference variables (e.g. communicators)
  ! ... belonging to each of the various parallelization levels:
  ! ... use the specific modules instead
  !
  USE mp_world, ONLY: mp_world_start, mp_world_end, nproc
  USE mp_bands
  !
  IMPLICIT NONE 
  !
  include 'laxlib.fh'
  !
  SAVE
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_startup ( ndiag, diag_in_band_group )
    !-----------------------------------------------------------------------
    ! ... This wrapper subroutine initializes all parallelization levels.
    ! ... 
    ! ... The np mpi processes (defined in my_world_comm) are subdivided into band
    ! ... groups (default nbgrp=1). Within each band group R & G data are distributed.
    ! ... Communication within the band group uses INTRA_BGRP_COMM 
    ! ... Communication across band groups uses INTER_BGRP_COMM
    ! ... An independent level of parallelization is used for distributed diagonalization 
    ! ... using ortho_comm (ndiag_) defined as 
    !
    USE command_line_options, ONLY : get_command_line, ndiag_, nband_
    USE parallel_include
    !
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ndiag
    LOGICAL, INTENT(IN), OPTIONAL :: diag_in_band_group
    LOGICAL :: do_distr_diag_inside_bgrp 
    INTEGER :: my_comm, ortho_parent_comm
    !
    ! start mpi
    !
    my_comm = MPI_COMM_WORLD 
    CALL mp_world_start( my_comm )
    !
    ! get input arguments
    !
    CALL get_command_line ( )
    write (6,*) 'parallel MPI defs : np=',nproc,'nb=', nband_,'  nd=', ndiag_
    !
    ! band group inizialization
    !
    CALL mp_start_bands ( nband_, my_comm )
    !
    ! distributed diag/ortho inizialization
    !
    do_distr_diag_inside_bgrp = .FALSE. ;  IF ( PRESENT(diag_in_band_group) ) do_distr_diag_inside_bgrp = diag_in_band_group
    !
    if (do_distr_diag_inside_bgrp) then
       ! used to be one diag group per bgrp with strict hierarchy: K-POINT > BAND > DIAG
       ortho_parent_comm = intra_bgrp_comm
    ELSE
       ! one diag group per individual k-point, with band group and
       ! diag group both being children of the current top level comm
       ortho_parent_comm = my_comm
    END IF

    CALL laxlib_start  ( ndiag_, my_comm, ortho_parent_comm, &
                        do_distr_diag_inside_bgrp_ = do_distr_diag_inside_bgrp )

    ndiag = ndiag_ ! copy input value to output value for latr use
    !
    RETURN
    !
  END SUBROUTINE mp_startup
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_end ( )
    !-----------------------------------------------------------------------
    !
    USE mp, ONLY : mp_comm_free
    !
    CALL mp_comm_free ( intra_bgrp_comm )
    CALL mp_comm_free ( inter_bgrp_comm )
    CALL mp_world_end( )
    !
    RETURN
    !
  END SUBROUTINE mp_global_end
  !
END MODULE mp_global
