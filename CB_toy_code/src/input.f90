 subroutine input(gamma_only)
! read the 'crystal' type (Si, Ge, GaAs, ...)
! read the linear dimension of the crystal supercell 
! read the number of bands to be computed
! read the kinetic energy cutoff (in Rydberg) 
! read the number and coordinates of the k-point set

! global variables
  USE cb_module
#if defined(__MPI)
  use mp_world,             ONLY : world_comm, root, mpime
  use mp,                   ONLY : mp_bcast
#endif
  implicit none
! input variables
  logical, intent(IN) :: gamma_only
! local variables
  character(4) :: crystal_name
  integer :: icar

  !>  system namelist
  NAMELIST / system / crystal_name, ncell, nbnd, ecutwfc, david, ethr, energy_shift, use_overlap

  !>  default values for the system namelist
  crystal_name = 'Si' ; ncell = 1 ; nbnd  = 0 ; ecutwfc = 4.d0 ; david=4 ; ethr=1.d-6; energy_shift=.true.
  use_overlap = .false.

! read the 'crystal' type (Si, Ge, GaAs, ...)
! read the linear dimension of the crystal supercell 
! read the number of bands to be computed
! read the kinetic energy cutoff (in Rydberg) 
#if defined(__MPI)
  if (mpime == root) then
#endif
     CALL input_from_file ( ) 
     read (5,system)
#if defined(__MPI)
  end if

! broadcast input
  call mp_bcast ( crystal_name, root, world_comm )
  call mp_bcast ( ncell       , root, world_comm )
  call mp_bcast ( nbnd        , root, world_comm )
  call mp_bcast ( ecutwfc     , root, world_comm )
  call mp_bcast ( david       , root, world_comm )
  call mp_bcast ( ethr        , root, world_comm )
  call mp_bcast ( energy_shift, root, world_comm )
  call mp_bcast ( use_overlap,  root, world_comm )
#endif

  do icrystal =1, ncrystal
     if (crystal_name==crystal(icrystal)) exit
  end do
  if (icrystal>ncrystal) stop ' wrong crystal name'
  if (ncell <= 0) stop ' ncell must be positive'
  if (nbnd < 0) stop ' nbnd must be positive'
  if (nbnd == 0) nbnd = 8 * ncell**3
  allocate (btype(nbnd)) ; btype = 1   ! btype=0; btype(1:min(4*ncell**3,nbnd)) = 1
  if (david < 2) stop ' david must be at least 2'
  nbndx = david * nbnd
  if (ecutwfc <= 0.d0) stop ' ecutwfc must be positive'

  if (gamma_only) then
     write (stdout,*) ' Gamma-only calculation, assuming nks=1, xk=(0,0,0) '
     nks = 1;  allocate ( xk(3,nks) ); xk = 0.d0
  else
! generic k-point calculation: read the number and coordinates of the k-point set
     write(6,*) 'nks >'
#if defined(__MPI)
     if (mpime == root) &
#endif
        read (5,*) nks
#if defined(__MPI)
     call mp_bcast ( nks, root, world_comm )
#endif

     allocate ( xk(3,nks) )

     do current_k =1,nks
#if defined(__MPI)
        if (mpime == root) &
#endif
           read (5,*) (xk(icar,current_k),icar=1,3)
     end do
#if defined(__MPI)
     call mp_bcast ( xk, root, world_comm )
#endif
  end if

 end subroutine input

!----------------------------------------------------------------------------
SUBROUTINE input_from_file( )
  !
  !! This subroutine checks command-line arguments for -i[nput] "file name"
  !! if "file name" is present, attach input unit 5 to the specified file
  !
  IMPLICIT NONE
  !
  INTEGER             :: stdin = 5, stderr = 6, ierr = 0
  CHARACTER (LEN=256) :: input_file
  LOGICAL             :: found
  !
  INTEGER :: iiarg, nargs
  !
  nargs = command_argument_count()
  found = .FALSE.
  input_file = ' '
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL get_command_argument( iiarg, input_file )
     !
     IF ( TRIM( input_file ) == '-i'     .OR. &
          TRIM( input_file ) == '-in'    .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-input' ) THEN
        !
        CALL get_command_argument( ( iiarg + 1 ) , input_file )
        found =.TRUE.
        EXIT
        !
     END IF
     !
  END DO
  !
  IF ( found ) THEN
     !
     OPEN ( UNIT = stdin, FILE = input_file, FORM = 'FORMATTED', &
            STATUS = 'OLD', IOSTAT = ierr )
     !
     ! TODO: return error code ierr (-1 no file, 0 file opened, > 1 error)
     ! do not call "errore" here: it may hang in parallel execution
     ! if this routine is called by a single processor
     !
     IF ( ierr > 0 ) WRITE (stderr, &
            '(" *** input file ",A," not found ***")' ) TRIM( input_file )
     !
  ELSE
     ierr = -1
  END IF
  !
  RETURN 
  !
END SUBROUTINE input_from_file

