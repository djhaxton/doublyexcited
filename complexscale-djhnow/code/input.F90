#include "definitions.INC"

MODULE INPUT
IMPLICIT NONE

REAL *8 :: PI = 3.14159265358979323846264
LOGICAL :: STRETCH = .FALSE.
LOGICAL :: MODBASIS = .FALSE.
LOGICAL :: VQUAD = .FALSE.
LOGICAL :: WRITEVALS = .TRUE.
LOGICAL :: WRITEVECS = .FALSE., READVECS=.FALSE.
CHARACTER(LEN=300) :: INVECFILE = "eigvecs"   ! in Bin subdirectory (hardwired)
CHARACTER(LEN=300) :: OUTVECFILE = "eigvecs"  !
CHARACTER(LEN=300) :: CORRFILE = "corr_test.txt"
CHARACTER(LEN=300) :: ENERGY_CURVE_FILE = "ecurve_test.txt"

! djh - deprecated - these should be different programs. 
! ACTION DEFINES WHAT THE EXECUTABLE PROGRAM DOES
! ACTION = 1 => DIAGONALIZE HAMILTONIAN
! ACTION = 2 => CORRUGATION CURVE
! ACTION = 3 => DIATOMIC ENERGY CURVE
! ACTION = 4 => DIATOMIC ROTATION CORRUGATION
! ACTION = 5 => COMPUTE DOUBLY EXCITED STATES WITH GPLMR

! INTEGER :: ACTION = 5

! GRID PARAMETERS

INTEGER :: N = 2
REAL *8 :: DELTA = 0.5
INTEGER :: NSMALL = 22

integer :: TFACTOR = 4

!! HMM NBIG IS USED TO ALLOCATE TINV... NOT INTENT OF ORIGINAL "NBIG" VARIABLE.
!!   in this code NBIG should equal 

!!! NO NO NO !!!  INTEGER :: NBIG = 72

INTEGER :: NCOR = 20
REAL *8 :: RADIUS = 7.0
LOGICAL :: AUTODELTA  = .TRUE.

! MOLECULE PARAMETERS
!! INTEGER :: NELEC = 1      ! NOT INPUT, PROGRAM VARIABLE, CHANGES
INTEGER :: NUCLEI = 1
INTEGER :: CHARGE = 2
INTEGER :: NUMPOINTS = 10
LOGICAL :: NUC_ON_GRID = .TRUE.

! COMPLEX SCALING ANGLE
REAL *8 :: THETA = 0.0d0

! LANCZOS PARAMETERS
! 3 is the size of the largest irreducible representation
!   so lanblocksize=3 default makes sense
INTEGER :: LANBLOCKSIZE=3  
INTEGER :: KDIM = 2001
INTEGER :: CHECK = 10
REAL *8 :: LANTOL = 1D-6
!LOGICAL :: SAVE_VECTORS = .FALSE.  ! NOT INPUT, PROGRAM VARIABLE, CHANGES

! GMRES SOLVER PARAMETERS
integer :: gm_usepetscflag=1
INTEGER :: GM_MAX_INNER_ITER = 5
INTEGER :: GM_MAXRST = 0
INTEGER, parameter :: GM_PCOPT = 0
LOGICAL :: GM_PRINT_OUTPUT = .FALSE.
REAL *8 :: GM_TOL = 1D-3
NUMBER :: GM_SHIFT

! HELIUM ORBTIAL SPECIFICATION
INTEGER :: NUMORBS = 14
LOGICAL :: SINGLET = .FALSE.
INTEGER :: NUMSKIP=0

! GPLMR PARAMETERS
INTEGER :: outmodulus=5
LOGICAL :: TAKEREALPART=.TRUE.
INTEGER :: GPLMR_NEV = 1
INTEGER :: GPLMR_NS = 5
INTEGER :: GPLMR_MAXIT = 100
REAL *8 :: GPLMR_TOL = 1D-12
NUMBER :: GPLMR_SHIFT = -999
! It would be nice to use these defaults (both false) --
!   C-norm for the eigenproblem, then Herm norm for GPLMR...
!   I assume GPLMR would prefer a shift based on the Herm norm.
!   If these variables disagree then the eigenvalues and
!   expectation values of the subspace diagonalization will
!   be the same.
LOGICAL :: HERMOPT=.FALSE.   ! define subspace eigvector problem w/ herm norm
LOGICAL :: CSHIFTOPT=.FALSE.  ! define inital shift through c-norm

! djh deprecated - these are different programs
! ACTIONS: ENERGIES VIA LANCZOS
!          H2+ CURVES
!          H2 CURVE (GROUND STATE ONLY)
!          TRANSLATE NUCLEUS
!          ROTATE DIATOMIC
!          DOUBLY-EXCITED STATE ENERGIES OF HELIUM

! RELEVANT PARAMETER DESCRIPTIONS
! ----------------------------------
! N - DETERMINES SIZE OF GRID: 2*N+1 POINTS PER SPATIAL DIMENSION
! RADIUS = PHYSICAL EXTENT OF GRID (SHOULD BE BETWEEN 3 AND 10)
! NELEC = NUMBER OF ELECTRONS, BETWEN 1 AND 3
! THETA = COMPLEX SCALING FACTOR (BETWEEN 0 AND PI/4)
! GM_SHIFT = TARGETED SHIFT FOR GMRES SOLVER
! GM_ITMAX = MAXIMUM NUMBER OF OUTER GMRES ITERATIONS
! GM_RST = MAXIMUM NUMBER OF INNER ITERATIONS PER OUTER ITERATION
! GM_TOL = RELATIVE ERROR TOLERANCE FOR GMRES
! PCOPT = PRE-CONDITIONING OPTION
!       0) NO PRE-CONDITIONING
!       1) LINEARIZED INVERSE (CHEAP BUT FAIRLY INEFFECTIVE)
!       2) SPECTRAL DECOMPOSITION OF (T-E0)
!       3) FIRST-LEVEL BLOCK INVERSION
!       4) SECOND-LEVEL BLOCK INVERSION
! There is no need for preconditioning with this algorithm.
!  we are not trying to solve (H-E)X=B with any accuracy at all.
!  Unless preconditioner is fast, it is not useful.
!  pcopt=0 hardwire.  --djh


END MODULE INPUT


subroutine checkiostat(myiostat,message,me,okay)
  implicit none
  logical,intent(out) :: okay
  integer,intent(in) :: myiostat,me
  character*(*),intent(in) :: message

  okay=.true.
  if (myiostat.ne.0) then
     okay=.false.
     if (ME==0) then
        WRITE(*,"(A8,I5,A8,I5,A3,A)") "RANK ",ME,"IOSTAT=",myiostat ," ",message
     endif
  endif

end subroutine checkiostat


subroutine get_input(inputfile,me)
  use input
  implicit none
  integer :: myiostat
  integer, intent(in) :: me
  character*(*), intent(in) :: inputfile
  logical :: ilog
  namelist/main/ &
       STRETCH ,&
       MODBASIS,&
       VQUAD ,&
       WRITEVALS,&
       WRITEVECS,&
       READVECS,&
       INVECFILE,&
       OUTVECFILE,&
       CORRFILE,&
       ENERGY_CURVE_FILE,&
       THETA
  namelist/grid/&
       N ,&
       DELTA,&
       NSMALL,&
!! NO NO NO   THIS IS USED TO ALLOCATE T^1!
!!  NOT INTENT OF ORIGINAL NBIG variable in DJH, JRJ notes     NBIG,&
!!  NBIG is an INTERNAL VARIABLE.
       TFACTOR,&
       NCOR,&
       RADIUS,&
       AUTODELTA
  namelist/molecule/&
! NOT INPUT       NELEC,&
       NUCLEI,&
       CHARGE,&
       NUMPOINTS,&
       NUC_ON_GRID
  namelist/lanczos/&
       KDIM,&
       LANBLOCKSIZE,&
       CHECK,&
       LANTOL
  namelist/gmres/&
       gm_usepetscflag, &
       GM_MAX_INNER_ITER,&
       GM_MAXRST,&
!!$       GM_PCOPT,&
       GM_PRINT_OUTPUT,&
       GM_TOL

  namelist/orbitals/&
       NUMORBS,&
       SINGLET, &
       NUMSKIP 
  namelist/gplmr/ &
       outmodulus,&
       TAKEREALPART,&
       HERMOPT,&
       CSHIFTOPT,&
       GPLMR_NEV,&
       GPLMR_NS,&
       GPLMR_MAXIT,&
       GPLMR_TOL, &
       GPLMR_SHIFT

  open(1016,file=inputfile,iostat=myiostat,status="old")
  call checkiostat(myiostat,"on open",me,ilog)
  if (.not.ilog) then
     return
  endif
  close(1016)


  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=main,iostat=myiostat)
  call checkiostat(myiostat,"namelist main",me,ilog)
  close(1016)  !! apparently I need to do this...rewind does not help


  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=grid,iostat=myiostat)
  call checkiostat(myiostat,"namelist grid",me,ilog)
  close(1016)

  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=molecule,iostat=myiostat)
  call checkiostat(myiostat,"namelist molecule",me,ilog)
  close(1016)

  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=lanczos,iostat=myiostat)
  call checkiostat(myiostat,"namelist lanczos",me,ilog)
  close(1016)

  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=gmres,iostat=myiostat)
  call checkiostat(myiostat,"namelist gmres",me,ilog)
  close(1016)

  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=orbitals,iostat=myiostat)
  call checkiostat(myiostat,"namelist orbitals",me,ilog)
  close(1016)

  open(1016,file=inputfile,iostat=myiostat,status="old")
  read(1016,nml=gplmr,iostat=myiostat)
  call checkiostat(myiostat,"namelist gplmr",me,ilog)
  close(1016)



end subroutine get_input


subroutine get_inputfile(inpfile)
  implicit none
  integer :: iarg,nargs,mylen
  character*(*), intent(out) :: inpfile
  character(len=200) :: buffer
#ifdef PGFFLAG
  integer :: myiargc
  nargs=myiargc()
#else
  nargs=iargc()
#endif

  do iarg=1,nargs
     call getarg(iarg,buffer)
     mylen=LEN_TRIM(buffer)
     if (buffer(1:4).eq.'Inp=') then
        inpfile(1:mylen-4)=buffer(5:mylen)
     endif
  end do
end subroutine get_inputfile

subroutine get_args(me)
  use input
  implicit none
  integer, intent(in) :: me
  integer :: iarg,nargs,mylen
  character(len=200) :: buffer
#ifdef PGFFLAG
  integer :: myiargc
  nargs=myiargc()
#else
  nargs=iargc()
#endif
  do iarg=1,nargs
     call getarg(iarg,buffer)
     mylen=LEN_TRIM(buffer)
     if (buffer(1:2).eq.'N=') then
        read(buffer(3:mylen),*) N
     endif
     if (buffer(1:6).eq.'Shift=') then
        read(buffer(7:mylen),*) gplmr_shift
     endif
     if (buffer(1:8).eq.'Outfile=') then
        read(buffer(9:mylen),*) outvecfile
        if (me==0) print *, "OUTFILE SET TO Bin/", TRIM(outvecfile)
     endif
     if (buffer(1:7).eq.'Infile=') then
        read(buffer(8:mylen),*) invecfile
     endif

     if (buffer(1:7).eq.'NoPetsc') then
        gm_usepetscflag=0
     endif

     if (buffer(1:5).eq.'Petsc') then
        gm_usepetscflag=1
     endif

     if (buffer(1:6).eq.'Theta=') then
        read(buffer(7:mylen),*) theta
     endif
     if (buffer(1:2).eq.'R=') then
        AUTODELTA=.true.
        read(buffer(3:mylen),*) radius
     endif
     if (buffer(1:6).eq.'Delta=') then
        AUTODELTA=.FALSE.
        read(buffer(7:mylen),*) delta
     endif

     if (buffer(1:4).eq.'Nev=') then
        read(buffer(5:mylen),*) gplmr_nev
     endif
     if (buffer(1:6).eq.'Norbs=') then
        read(buffer(7:mylen),*) numorbs
     endif
     if (buffer(1:6).eq.'Gptol=') then
        read(buffer(7:mylen),*) gplmr_tol
     endif
     if (buffer(1:5).eq.'Skip=') then
        read(buffer(6:mylen),*) numskip
     endif
     if (buffer(1:6).eq.'Gmtol=') then
        read(buffer(7:mylen),*) gm_tol
     endif
     if (buffer(1:5).eq.'Iter=') then
        read(buffer(6:mylen),*) gm_max_inner_iter
     endif
     if (buffer(1:4).eq.'Read') then
        numskip=0
        readvecs=.true.
     endif
     if (buffer(1:8).eq.'Takereal') then
        takerealpart=.true.
     endif
     if (buffer(1:10).eq.'Notakereal') then
        takerealpart=.false.
     endif
     if (buffer(1:4).eq.'Sing') then
        singlet=.true.
     endif
     if (buffer(1:4).eq.'Trip') then
        singlet=.false.
     endif
  end do
end subroutine get_args


SUBROUTINE INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM,SAVE_VECTORS,NELEC)
USE STRUCTS
USE INPUT
USE INITMOD
implicit none  ! why not????  might move this subroutine.

! Suborutine to initialize all structures for a given problem
logical,intent(IN):: SAVE_VECTORS
!character(len=200) :: inpfile="Input.Inp"

INTEGER, INTENT(IN) :: ME,NP,NELEC
TYPE(PARAMETERS), INTENT(OUT) :: PARAMS
TYPE(MPI_DATA), INTENT(OUT) :: MPIDATA
TYPE(LANPARAMS), INTENT(OUT) :: LANPAR
TYPE(MOLECULE), INTENT(OUT) :: MOL
TYPE(HAMILTONIAN), INTENT(OUT) :: HAM

! Check for valid input
IF(NELEC == 1 .AND. NP > 1 .AND. .NOT. NUC_ON_GRID) THEN
   IF(ME == 0) WRITE(*,"(A)") "ERROR: MPI not supported for one-electron for NUC_ON_GRID == .FALSE."
   stop !!lame GOTO 100
END IF

! Automatically compute resoltion based on specified box size
IF(AUTODELTA) DELTA = RADIUS/N


! Initialize molecule
CALL INIT_MOLECULE(MOL,NUCLEI)

! Define nuclear charges
MOL % CHARGE = CHARGE
! Define nuclear coordinates
MOL % RNUC = 0

! Initialize PARAMETERS structure
CALL INIT_PARAMETERS(PARAMS,N,NSMALL,TFACTOR,NCOR,NELEC,&
     DELTA,THETA,VQUAD,NUC_ON_GRID)

! Initialize MPI_DATA structure
CALL INIT_MPI_DATA(MPIDATA,PARAMS,ME,NP)

! Initialize HAMILTONIAN structure
CALL INIT_HAMILTONIAN(HAM,PARAMS,MOL,MPIDATA)

! Initialize LANPARAMS structure USED FOR ONE-E ONLY??? DJH
CALL INIT_LANPARAMS(LANPAR,NUMORBS,KDIM,CHECK, &
     LANTOL,MPIDATA % MM,.TRUE.,SAVE_VECTORS)

! Print values of structures to screen
CALL PRINTINPUT(PARAMS,MPIDATA,LANPAR)

END SUBROUTINE INIT_STRUCTS



!--------------------------------------------------------------------------------------

SUBROUTINE PRINTINPUT(PARAMS,MPIDATA,LANPAR)
USE STRUCTS
use input
implicit none

! Prints selected values from the structures to the screen

TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(LANPARAMS), INTENT(IN) :: LANPAR

IF(MPIDATA % ME == 0) THEN

   WRITE(*,"(A)") "STARTING GRID BASIS PROGRAM."
   WRITE(*,"(A)") "----------------------------------------------------------------------"
#if ISREAL
   WRITE(*,"(A)") "DATATYPE: REAL DOUBLE PRECISION"
#else
   WRITE(*,"(A)") "DATATYPE: COMPLEX DOUBLE PRECISION"
#endif

   WRITE(*,"(A,I8)") "GRID POINTS PER DIMENSION: ",PARAMS % NN
   WRITE(*,"(A,F8.2)") "GRID SPACING: ",PARAMS % DELTA
   WRITE(*,"(A,I8)") "NUMBER OF 3D GRID POINTS: ",PARAMS % M
   WRITE(*,"(A,I12)") "TOTAL BASIS DIMENSION: ",(PARAMS % M)**(PARAMS % NELEC)
   WRITE(*,"(A,I8)") "NUMBER OF ELECTRONS: ",PARAMS % NELEC
!!   WRITE(*,"(A,I8,I8)") "NBIG/NSMALL: ",PARAMS%NBIG,PARAMS%NSMALL

   WRITE(*,"(A,I8,I8)") '"NBIG":',PARAMS%NBIG
   WRITE(*,"(A,2I8)") "TFACTOR*NSMALL/NSMALL: ",PARAMS%NSMALL*PARAMS%TFACTOR,PARAMS%NSMALL
   WRITE(*,"(A,I8)") "NUMBER OF DESIRED EIGENVALUES: ",LANPAR % NUMEIGS
   WRITE(*,"(A,I8)") "KRYLOV DIMENSION: ",LANPAR % KDIM
   WRITE(*,"(A,I8)") "EXTRACTION FREQUENCY: ",LANPAR % CHECK
   WRITE(*,"(A,I8)") "NUMBER OF MPI PROCESSES: ",MPIDATA % NP
   WRITE(*,*)
   WRITE(*,"(A,I8)") "gm_usepetscflag ",gm_usepetscflag
   WRITE(*,"(A,I8)") "gm_max_inner_iter ",gm_max_inner_iter

WRITE(*,"(A)") "----------------------------------------------------------------------"

END IF

END SUBROUTINE PRINTINPUT


