#include "definitions.INC"

! --------------------------------------------------------------------------------------

PROGRAM ATOM
USE INITMOD
USE STRUCTS
USE INPUT
USE LANMOD
USE SOLVER
USE MPI
USE GPLMR
#if defined(USE_PETSC)
USE PETSCWRAP
#endif

IMPLICIT NONE

! --------------------------------------------------------------------------------------

! Declare variables and structures

INTEGER :: I, J, ME, NP, MPIERR, M, MB, FIRST, LAST
REAL *8 :: TIME, ANG_EXP
TYPE(PARAMETERS) :: PARAMS
TYPE(MPI_DATA) :: MPIDATA
TYPE(LANPARAMS), TARGET :: LANPAR
TYPE(MOLECULE) :: MOL
TYPE(HAMILTONIAN) :: HAM
TYPE(GMRES_SOLVER) :: GMSOLVER
INTEGER :: ACTION=5

! Initialize MPI variables
#if defined(USE_PETSC)
CALL PETSCINITIALIZE(PETSC_NULL_CHARACTER,MPIERR)
#else
CALL MPI_INIT(MPIERR)
#endif
CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,MPIERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,MPIERR)

! Start program timer
TIME = MPI_WTIME()

! Check for valid input
IF(NELEC == 1 .AND. NP > 1 .AND. .NOT. NUC_ON_GRID) THEN
   IF(ME == 0) WRITE(*,"(A)") "ERROR: MPI not supported for one-electron for NUC_ON_GRID == .FALSE."
   GOTO 100
END IF

! Automatically compute resoltion based on specified box size
IF(AUTODELTA) DELTA = RADIUS/N

! Initialize structures based on input parameters
CALL INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM)

#if defined(USE_PETSC)
GMSOLVER%IUNIT=GM_PRINT_OUTPUT
GMSOLVER%MAX_INNER_ITER=GM_MAX_INNER_ITER
GMSOLVER%PCOPT=GM_PCOPT
GMSOLVER%TOL=GM_TOL
#endif

! Perform the specified action
IF(ACTION == 1) THEN

   ! ACTION == 1 ==> diagonalize Hamiltonian
   CALL DIAGONALIZE(PARAMS,MPIDATA,LANPAR,HAM)

   IF(WRITEVALS) THEN
      IF(MPIDATA % ME == 0) THEN
         OPEN(4,FILE=TRIM(VALFILE))
         WRITE(4,*) LANPAR % EIGENVALUES
         CLOSE(4)
      END IF
   END IF

   IF(WRITEVECS) THEN
      CALL WRITE_EIGENVECTORS(PARAMS,LANPAR,MPIDATA)
   END IF

ELSE IF(ACTION == 2) THEN

   ! ACTION == 2 ==> compute corrugation curve
   NUC_ON_GRID = .FALSE.

   IF(NUCLEI > 1) THEN
      IF(ME == 0) WRITE(*,"(A)") "ERROR: NUCLEI must be equal to 1 to compute corrugation curve."
      GOTO 100
   ELSE
      CALL CORRUGATION(HAM,PARAMS,MPIDATA,LANPAR)
   END IF

ELSE IF(ACTION == 3) THEN

   ! ACTION == 3 ==> compute diatomic energy vs. internuclear distance
   NUC_ON_GRID = .FALSE.

   IF(NUCLEI /= 2) THEN
      IF(ME == 0) WRITE(*,"(A)") "ERROR: NUCLEI must be equal to 2 to compute diatomic energy curve."
      GOTO 100
   ELSE
      CALL ENERGYCURVE(HAM,PARAMS,MPIDATA,LANPAR)
   END IF

ELSE IF(ACTION == 4) THEN

   ! ACTION == 4 ==> compute diatomic energy for a fixed radius but varying the 
   ! angle the molecular axis makes with the +x axis
   
   NUC_ON_GRID = .FALSE.

   IF(NUCLEI /= 2) THEN
      IF(ME == 0) WRITE(*,"(A)") "ERROR: NUCLEI must be equal to 2 to compute diatomic corrugation curve."
      GOTO 100
   ELSE
      CALL ROTATE_DIATOMIC(HAM,PARAMS,MPIDATA,LANPAR)
   END IF

ELSEIF(ACTION == 5) THEN

   ! ACTION == 5 ==> compute doubly-excited states of helium
   
#if ISREAL
   IF(ME == 0) WRITE(*,"(A)") "ERROR: COMPLEX datatype needed for doubly-excited helium calculation."
   GOTO 100
#endif

   CALL DOUBLYEXCITED(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM,GMSOLVER)

ELSE

   WRITE(*,"(A)") "INVALID VALUE FOR ACTION. EXITING."
   GOTO 100

END IF

! Finalize MPI and print elapsed time of program

100 CONTINUE

TIME = MPI_WTIME() - TIME
IF(ME == 0) WRITE(*,"(A,F10.2,A)") "TOTAL TIME ELAPSED: ",TIME," SECONDS."
CALL MPI_FINALIZE(MPIERR)

!--------------------------------------------------------------------------------------

CONTAINS

!--------------------------------------------------------------------------------------

#if ISREAL == 0

! Subroutine for computing doubly-excited states with GPLMR

SUBROUTINE DOUBLYEXCITED(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM,GMSOLVER)
USE INPUT
IMPLICIT NONE

INTEGER, INTENT(IN) :: ME,NP
TYPE(PARAMETERS), INTENT(INOUT) :: PARAMS
TYPE(MPI_DATA), INTENT(INOUT) :: MPIDATA
TYPE(LANPARAMS), INTENT(INOUT), TARGET :: LANPAR
TYPE(MOLECULE), INTENT(INOUT) :: MOL
TYPE(HAMILTONIAN), INTENT(INOUT) :: HAM

INTEGER :: M,MB,FIRST,LAST
NUMBER, ALLOCATABLE :: PSI(:,:,:),E0(:), EIGS(:)
NUMBER, POINTER :: EV(:,:)
TYPE(GMRES_SOLVER) :: GMSOLVER

! EXCITED STATE EXACT ANSWERS
!------------------------------
! 2S2P SINGLET: (-0.693135,-6.865E-4)
! 2S2P TRIPLET: (-0.760492,-1.495E-4)
! 2P2P SINGLET: (-0.701946,-1.200E-3)
! 2P2P TRIPLET: -0.710500


! Set up and diagonalize 1 electron He+ Hamiltonian

NELEC = 1
SAVE_VECTORS = .TRUE.
GSORTH = .TRUE.
NUMEIGS = 10

CALL INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM)
CALL DIAGONALIZE(PARAMS,MPIDATA,LANPAR,HAM)

! Compute approximate eigenvectors for helium based
! on the eigenvectors for He+

M = PARAMS % M
MB = MPIDATA % MBLK
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

IF(ME == 0) WRITE(*,"(A)") "COMPUTING HELIUM ORBITAL APPROXIMATIONS..."

EV => LANPAR % EIGENVECTORS
ALLOCATE(PSI(1:M,FIRST:LAST,1:NUMORBS),E0(1:NUMORBS),EIGS(NUMORBS))

IF(NUMORBS == 1) THEN

   GPLMR_NEV = 1

   IF(SINGLET) THEN
      CALL HELIUM_ORBITAL(EV(:,ORB1),EV(:,ORB2),PSI(:,:,1),M,FIRST,LAST,MPIDATA,1)
   ELSE
      CALL HELIUM_ORBITAL(EV(:,ORB1),EV(:,ORB2),PSI(:,:,1),M,FIRST,LAST,MPIDATA,-1)
   ENDIF

ELSE

   GPLMR_NEV = 5

   ! SINGLETS: 2S2S, 2S2P, 2P2P
   CALL HELIUM_ORBITAL(EV(:,2),EV(:,2),PSI(:,:,1),M,FIRST,LAST,MPIDATA,1)
   CALL HELIUM_ORBITAL(EV(:,2),EV(:,3),PSI(:,:,2),M,FIRST,LAST,MPIDATA,1)
   CALL HELIUM_ORBITAL(EV(:,3),EV(:,4),PSI(:,:,3),M,FIRST,LAST,MPIDATA,1)

   ! TRIPLETS: 2S2P, 2P2P (NO 2S2S TRIPLET)
   CALL HELIUM_ORBITAL(EV(:,2),EV(:,3),PSI(:,:,4),M,FIRST,LAST,MPIDATA,-1)
   CALL HELIUM_ORBITAL(EV(:,3),EV(:,4),PSI(:,:,5),M,FIRST,LAST,MPIDATA,-1)

END IF

IF(ME == 0) WRITE(*,"(A)") "INITIALIZING HELIUM HAMILTONIAN..."

! Set up helium Hamiltonian

NELEC = 2
SAVE_VECTORS = .FALSE.
CALL INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM)

! Compute expectations of initial eigenvector guesses to use
! for the shift parameter in GPLMR

IF(ME == 0) WRITE(*,"(A)") "COMPUTING INITIAL ENERGY GUESSES..."

CALL DIAGSUBSPACE(HAM,PARAMS,MPIDATA,M,MB,NUMORBS,E0,PSI)

IF(ME == 0) THEN
   WRITE(*,"(A)") "IMPROVED INITIAL ENERGY GUESSES ARE: "
   CALL PRINTVECTOR(E0,NUMORBS)
END IF

! Compute and print the average energy of the He orbitals
GPLMR_SHIFT = SUM(E0)/NUMORBS

IF(ME == 0) THEN
   WRITE(*,"(A,F20.10,F20.10)") "AVERAGE SHIFT: ",DREAL(GPLMR_SHIFT),DIMAG(GPLMR_SHIFT)
END IF

! Call the GPLMR solver
CALL GPLMR_EIG(HAM,PARAMS,MPIDATA,GMSOLVER,PSI,EIGS, &
               PARAMS % M,MPIDATA % MBLK,GPLMR_SHIFT, &
               GPLMR_NEV,GPLMR_NS,GPLMR_MAXIT,GPLMR_TOL,MSOLVECPLX2)

DEALLOCATE(PSI,E0,EIGS)

END SUBROUTINE DOUBLYEXCITED

#endif

!--------------------------------------------------------------------------------------

SUBROUTINE INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM)
USE INPUT

! Suborutine to initialize all structures for a given problem

character(len=200) :: inpfile="Input.Inp"
INTEGER, INTENT(IN) :: ME,NP
TYPE(PARAMETERS), INTENT(OUT) :: PARAMS
TYPE(MPI_DATA), INTENT(OUT) :: MPIDATA
TYPE(LANPARAMS), INTENT(OUT) :: LANPAR
TYPE(MOLECULE), INTENT(OUT) :: MOL
TYPE(HAMILTONIAN), INTENT(OUT) :: HAM

IF(ME == 0) WRITE(*,"(A)") "INITIALIZING STRUCTURES..."


call get_inputfile(inpfile)

! Initialize variables (user input)
call GET_INPUT(inpfile,ME)

call get_args(me)


! Initialize molecule
CALL INIT_MOLECULE(MOL,NUCLEI,NELEC,DELTA)

! Define nuclear charges
MOL % CHARGE = CHARGE
! Define nuclear coordinates
MOL % RNUC = 0

! Initialize PARAMETERS structure
CALL INIT_PARAMETERS(PARAMS,N,NSMALL,NBIG,NCOR,NELEC,&
     DELTA,THETA,VQUAD,NUC_ON_GRID)

! Initialize MPI_DATA structure
CALL INIT_MPI_DATA(MPIDATA,PARAMS,ME,NP)

! Initialize HAMILTONIAN structure
CALL INIT_HAMILTONIAN(HAM,PARAMS,MOL,MPIDATA)

! Initialize LANPARAMS structure
CALL INIT_LANPARAMS(LANPAR,NUMEIGS,KDIM,CHECK, &
     LANTOL,MPIDATA % MM,GSORTH,SAVE_VECTORS)

! Print values of structures to screen
CALL PRINTINPUT(PARAMS,MPIDATA,LANPAR)

END SUBROUTINE INIT_STRUCTS

!--------------------------------------------------------------------------------------

SUBROUTINE DIAGONALIZE(PARAMS,MPIDATA,LANPAR,HAM)
USE INPUT

TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR
TYPE(HAMILTONIAN), INTENT(IN) :: HAM

! Call Lanczos solver
IF(MPIDATA % ME == 0) WRITE(*,"(A)") "CALLING LANCZOS SOLVER..."
CALL LANCZOS(PARAMS,LANPAR,HAM,MPIDATA)

! Print out computed eigenvalues
IF(MPIDATA % ME == 0) THEN
   WRITE(*,"(A)") "COMPUTED EIGENVALUES: "
   CALL PRINTVECTOR(LANPAR % EIGENVALUES, LANPAR % NUMEIGS)
END IF


END SUBROUTINE DIAGONALIZE

!--------------------------------------------------------------------------------------

SUBROUTINE ENERGYCURVE(HAM,PARAMS,MPIDATA,LANPAR)
USE INPUT, ONLY:NUMPOINTS,ENERGY_CURVE_FILE

TYPE(HAMILTONIAN), INTENT(INOUT) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR

INTEGER :: I,J,P
REAL *8 :: DELTA, DP, R
REAL *8, ALLOCATABLE :: ENERGIES(:,:), DISTANCES(:)

DELTA = PARAMS % DELTA
HAM % MOL % RNUC = 0

! Number of internuclear distance points to use for curve
P = NUMPOINTS

! Spacing between internuclear distance points
DP = 2.5D0/(1D0*P)

! Allocate arrays
ALLOCATE(ENERGIES(P,LANPAR % NUMEIGS))
ALLOCATE(DISTANCES(P))

! Loop over points and compute energies
DO I = 1,P

   IF(MPIDATA % ME == 0) WRITE(*,"(A,I5,A,I5)") "CURRENTLY AT POINT ",I," OUT OF ",P
   HAM % MOL % RNUC(1,1) = -I*DP
   HAM % MOL % RNUC(1,2) = I*DP
   DISTANCES(I) = 2*I*DP

	CALL LANCZOS(PARAMS,LANPAR,HAM,MPIDATA)
   ENERGIES(I,:) = LANPAR % EIGENVALUES
	IF(MPIDATA % ME == 0) WRITE(*,"(A,F20.10)") "GS ENERGY =  ",ENERGIES(I,1)

END DO

! Write distance/energies to ENERGY_CURVE_FILE
! The format for each colums is [DISTANCE,ENERGY_1,ENERGY_2,...,ENERGY_NUMEIGS]

IF(MPIDATA % ME == 0) THEN
   OPEN(4,FILE = ENERGY_CURVE_FILE)
   
   DO I = 1,P
   	WRITE(4,"(F20.10)",ADVANCE="NO") DISTANCES(I)
   	DO J = 1,LANPAR%NUMEIGS
   		WRITE(4,"(F20.10)",ADVANCE="NO") ENERGIES(I,J)
   	END DO
   	WRITE(4,*)
   END DO
   CLOSE(4)
   
END IF

END SUBROUTINE

!--------------------------------------------------------------------------------------

SUBROUTINE CORRUGATION(HAM,PARAMS,MPIDATA,LANPAR)
USE INPUT, ONLY:CORRFILE,NUMPOINTS

TYPE(HAMILTONIAN), INTENT(INOUT) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR

INTEGER :: I,P,FIRST,LAST,NP,ME,MBLK,MAXBLK,MPIERR
INTEGER, ALLOCATABLE :: BLOCKS(:), DISPLS(:)

REAL *8 :: DELTA, DP, R
REAL *8, ALLOCATABLE :: ENERGIES(:,:), ALLENERGIES(:,:)
REAL *8 :: VEC(3)

DELTA = PARAMS % DELTA
HAM % MOL % RNUC = 0
NP = MPIDATA % NP
ME = MPIDATA % ME

! Number of points to use for the plot
P = NUMPOINTS

! Spacing between points
R = 2D0
DP = R/P

! Define unit vector in the direction in which the nucleus is moved
VEC(1) = 1
VEC(2) = 0
VEC(3) = 0
VEC = VEC/SQRT(SUM(VEC**2))

! Allocate array
ALLOCATE(ENERGIES(P,LANPAR % NUMEIGS))

! Loop over points and compute energies
DO I = 1,P

	IF(ME == 0) WRITE(*,"(A,I5,A,I5)") "CURRENTLY AT POINT ",I," OUT OF ",P
   HAM % MOL % RNUC(:,1) = VEC*(I-1)*DP
   CALL LANCZOS(PARAMS,LANPAR,HAM,MPIDATA)
   ENERGIES(I,:) = LANPAR % EIGENVALUES
	IF(ME == 0) WRITE(*,"(A,F20.10)") "GROUND STATE = ",ENERGIES(I,1) 

END DO

! Write energy data to CORRFILE

IF(ME == 0) THEN
   OPEN(4,FILE = TRIM(CORRFILE))
   WRITE(4,"(F20.10)") ENERGIES
   CLOSE(4)
END IF

END SUBROUTINE CORRUGATION

!--------------------------------------------------------------------------------------

SUBROUTINE ROTATE_DIATOMIC(HAM,PARAMS,MPIDATA,LANPAR)
USE INPUT, ONLY:CORRFILE,NUMPOINTS

TYPE(HAMILTONIAN), INTENT(INOUT) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR

INTEGER :: I,P
REAL *8,ALLOCATABLE :: ENERGIES(:,:)
REAL *8 :: DELTA, X, Y, R, DT

! Number of points to use for the plot
P = NUMPOINTS

! Fixed distance between nuclei
R = 1D0

DELTA = PARAMS % DELTA
HAM % MOL % RNUC = 0

! Spacing between angles
DT = 2D0*(PARAMS % PI) / P

! Allocate arrray
ALLOCATE(ENERGIES(P,LANPAR % NUMEIGS))

! Lop over angles and compute energies
DO I = 1,P

	THETA = (I-1)*DT
	X = R*COS(THETA)
	Y = R*SIN(THETA)
   HAM % MOL % RNUC(1,1) = X
   HAM % MOL % RNUC(2,1) = Y
   HAM % MOL % RNUC(1,2) = -X
   HAM % MOL % RNUC(2,2) = -Y

	CALL LANCZOS(PARAMS,LANPAR,HAM,MPIDATA)
   ENERGIES(I,:) = LANPAR % EIGENVALUES
	WRITE(*,*) ENERGIES(I,1)
	
END DO

! Write energy data to CORRFILE

IF(MPIDATA % ME == 0) THEN
   OPEN(4,FILE = CORRFILE)
   WRITE(4,"(F20.10)") ENERGIES
   CLOSE(4)
END IF

DEALLOCATE(ENERGIES)

END SUBROUTINE ROTATE_DIATOMIC

!--------------------------------------------------------------------------------------

SUBROUTINE PRINTVECTOR(X,N)

! Subroutine to print a real or complex vector to the screen

INTEGER, INTENT(IN) :: N
NUMBER, INTENT(IN) :: X(1:N)
INTEGER :: I

#if ISREAL
   WRITE(*,"(F10.5)") X
#else
   DO I = 1,N
      WRITE(*,"(F20.15,F20.15)") DREAL(X(I)), DIMAG(X(I))
   END DO
#endif

END SUBROUTINE PRINTVECTOR

!--------------------------------------------------------------------------------------

SUBROUTINE PRINTINPUT(PARAMS,MPIDATA,LANPAR)

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
   WRITE(*,"(A,I8,I8)") "NBIG/NSMALL: ",PARAMS%NBIG,PARAMS%NSMALL
   WRITE(*,"(A,I8)") "NUMBER OF DESIRED EIGENVALUES: ",LANPAR % NUMEIGS
   WRITE(*,"(A,I8)") "KRYLOV DIMENSION: ",LANPAR % KDIM
   WRITE(*,"(A,I8)") "EXTRACTION FREQUENCY: ",LANPAR % CHECK
   WRITE(*,"(A,I8)") "NUMBER OF MPI PROCESSES: ",MPIDATA % NP
   WRITE(*,"(A,I8)") "MPI BLOCKS SIZE: ",MPIDATA % MM

WRITE(*,"(A)") "----------------------------------------------------------------------"

END IF

END SUBROUTINE PRINTINPUT

!--------------------------------------------------------------------------------------

SUBROUTINE WRITE_EIGENVECTORS(PARAMS,LANPAR,MPIDATA)

! Prints eigenvectors distributed over MPI processes to the file 
! specified by VECFILE

TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(LANPARAMS), TARGET, INTENT(IN) :: LANPAR
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA

INTEGER :: I, P, ME, NP, MPIERR, M
NUMBER, POINTER :: U(:)

M = PARAMS % M
ME = MPIDATA % ME
NP = MPIDATA % NP

IF(ME == 0) OPEN(4,FILE = TRIM(VECFILE))

DO I = 1,LANPAR % NUMEIGS

   U => LANPAR % EIGENVECTORS(:,I)

   IF(ME == 0) THEN
      
      WRITE(4,"(F20.15)") U

      DO P = 1,NP-1
         CALL MPI_RECV(U,M/NP,MPIDATA % MPISIZE,P,1, &
         MPIDATA % COMM,MPIDATA % STATUS,MPIERR)
         WRITE(4,"(F20.15)") U(1:(M/NP))
      END DO

   ELSE

      CALL MPI_SEND(U,M/NP,MPIDATA % MPISIZE,0,1,MPIDATA % COMM,MPIERR)

   END IF

END DO

IF(ME == 0) CLOSE(4)


END SUBROUTINE WRITE_EIGENVECTORS

!--------------------------------------------------------------------------------------

END PROGRAM ATOM
