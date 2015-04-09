#if defined(USE_PETSC)

#include "definitions.INC"
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"

MODULE PETSCWRAP

USE PETSCSYS
USE PETSCVEC
USE PETSCMAT
USE PETSCKSP
USE PETSCPC
USE SNEAK
USE HAMMOD
USE SOLVER
IMPLICIT NONE

CONTAINS

SUBROUTINE  WrapMatVec(A,X,Y,IERR)
  IMPLICIT NONE
  Mat :: A
  Vec :: X, Y
  INTEGER :: IERR

  INTEGER :: ME, LOCN, GLOBN, I, GFIRST
  INTEGER, ALLOCATABLE :: GIND(:)
  COMPLEX *16, ALLOCATABLE :: XX(:),YY(:)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,IERR)

  CALL VecGetSize(X,GLOBN,IERR)
  CALL VecGetLocalSize(X,LOCN,IERR)
  ALLOCATE(GIND(LOCN))
  ALLOCATE(XX(LOCN))
  ALLOCATE(YY(LOCN))

  ! Global indices of variables owned by ME
  IF(ME.EQ.0) THEN
    GFIRST=0
  ELSE
    GFIRST=MPIDATA_GL%MAXBLOCK*PARAMS_GL%M**(PARAMS_GL%NELEC-1)+(ME-1)*MPIDATA_GL%MBLK*PARAMS_GL%M**(PARAMS_GL%NELEC-1)
  ENDIF
  DO I=1,LOCN
    GIND(I)=GFIRST+I-1 ! 0-based
  END DO

  ! PETSc to "regular" structure
  CALL VecGetValues(X,LOCN,GIND,XX,IERR)
  CALL VecGetValues(Y,LOCN,GIND,YY,IERR)

  ! Call our matvec
  ! 1/ A*x
  CALL HAMMULT(HAM_GL,PARAMS_GL,MPIDATA_GL,XX,YY,LOCN)
  ! 2/ -shift*x
  YY=YY-SHIFT_GL*XX 

  ! Back to PETSc structure
  CALL VecSetValues(X,LOCN,GIND,XX,INSERT_VALUES,IERR)
  CALL VecAssemblyBegin(X,IERR)
  CALL VecAssemblyEnd(X,IERR)
  CALL VecSetValues(Y,LOCN,GIND,YY,INSERT_VALUES,IERR)
  CALL VecAssemblyBegin(Y,IERR)
  CALL VecAssemblyEnd(Y,IERR)

  ! Cleanup
  DEALLOCATE(GIND)
  DEALLOCATE(XX,YY)

END SUBROUTINE

SUBROUTINE  WrapPrec(P,X,Y,IERR)
  IMPLICIT NONE
  PC :: P
  Vec :: X, Y
  INTEGER :: IERR

  INTEGER :: ME, LOCN, GLOBN, I, GFIRST
  INTEGER, ALLOCATABLE :: GIND(:)
  COMPLEX *16, ALLOCATABLE :: XX(:),YY(:)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,IERR)

  CALL VecGetSize(X,GLOBN,IERR)
  CALL VecGetLocalSize(X,LOCN,IERR)
  ALLOCATE(GIND(LOCN))
  ALLOCATE(XX(LOCN))
  ALLOCATE(YY(LOCN))

  ! Global indices of variables owned by ME
  IF(ME.EQ.0) THEN
    GFIRST=0
  ELSE
    GFIRST=MPIDATA_GL%MAXBLOCK*PARAMS_GL%M**(PARAMS_GL%NELEC-1)+(ME-1)*MPIDATA_GL%MBLK*PARAMS_GL%M**(PARAMS_GL%NELEC-1)
  ENDIF
  DO I=1,LOCN
    GIND(I)=GFIRST+I-1 ! 0-based
  END DO

  ! PETSc to "regular" structure
  CALL VecGetValues(X,LOCN,GIND,XX,IERR)
  CALL VecGetValues(Y,LOCN,GIND,YY,IERR)

  ! Call one of our preconditioners
  SELECT CASE (PC_GL)
     CASE (1)
       CALL MSOLVECPLX12(XX,YY,LOCN) ! Linear inversion
     CASE (2)
       CALL MSOLVECPLX22(XX,YY,LOCN) ! Spectral decomposition
#if defined(USE_SUPERLU)
     CASE (3)
       CALL MSOLVECPLX32(XX,YY,LOCN) ! 1st-level block diagonal
     CASE (4)
       CALL MSOLVECPLX42(XX,YY,LOCN) ! 2nd-level block diagonal
     CASE (5)
       CALL MSOLVECPLX52(XX,YY,LOCN) ! 3rt-level block diagonal
#endif
     CASE DEFAULT
       ! No preconditioner
       YY=XX
  END SELECT


  ! Back to PETSc structure
  CALL VecSetValues(X,LOCN,GIND,XX,INSERT_VALUES,IERR)
  CALL VecAssemblyBegin(X,IERR)
  CALL VecAssemblyEnd(X,IERR)
  CALL VecSetValues(Y,LOCN,GIND,YY,INSERT_VALUES,IERR)
  CALL VecAssemblyBegin(Y,IERR)
  CALL VecAssemblyEnd(Y,IERR)

  ! Cleanup
  DEALLOCATE(GIND)
  DEALLOCATE(XX,YY)

END SUBROUTINE

SUBROUTINE PETSCSOLVE(B,X,N,HAM,PARAMS,MPIDATA,GMSOLVER,OUTMODULUS)
  ! Wrapper to the KPC package in PETSC.
  ! We use the matrix-free feature for the matrix and
  ! the preconditioner. We use GMRES but virtually
  ! other solver (iterative, or even direct!), could
  ! be used.

  TYPE(GMRES_SOLVER), TARGET, INTENT(INOUT) :: GMSOLVER
  TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
  TYPE(PARAMETERS), TARGET, INTENT(IN) :: PARAMS
  TYPE(MPI_DATA), TARGET, INTENT(IN) :: MPIDATA
  INTEGER, INTENT(IN) :: N , OUTMODULUS
  NUMBER, INTENT(IN) :: B(1:N)
  NUMBER, INTENT(INOUT) :: X(1:N)

  INTEGER :: ME,IERR,I,GFIRST,ITER
  INTEGER :: GIND(N)
  Vec :: VB, VX
  KSP :: solver
  Mat :: A
  PC  :: P
  DOUBLE PRECISION :: rnorm
  real*8 :: ATIME,BTIME
  INTEGER, SAVE ::  CALLEDFLAG=0
  REAL*8, SAVE :: TIMINGS(100)=0
  
  CALLEDFLAG=CALLEDFLAG+1

  ATIME=MPI_WTIME()

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,IERR)

  ! Set global stuff
  HAM_GL => HAM
  PARAMS_GL => PARAMS
  MPIDATA_GL => MPIDATA
  SHIFT_GL => GMSOLVER%SHIFT
  PC_GL => GMSOLVER%PCOPT

  ! Global indices of variables owned by ME
  IF(ME.EQ.0) THEN
    GFIRST=0
  ELSE
    GFIRST=MPIDATA%MAXBLOCK*PARAMS%M**(PARAMS%NELEC-1)+(ME-1)*MPIDATA%MBLK*PARAMS%M**(PARAMS%NELEC-1)
  ENDIF
  DO I=1,N
    GIND(I)=GFIRST+I-1 ! 0-based
  END DO

  BTIME=MPI_WTIME()  ;  TIMINGS(1)=BTIME-ATIME;   ATIME=BTIME
  
  ! Copy RHS B into a PETSC vector VB
  CALL VecCreateMPI(MPI_COMM_WORLD,N,PETSC_DETERMINE,VB,IERR)
  CALL VecSetValues(VB,N,GIND,B,INSERT_VALUES,IERR)
  CALL VecAssemblyBegin(VB,IERR)
  CALL VecAssemblyEnd(VB,IERR)
  !CALL VecView(VB,PETSC_VIEWER_STDOUT_WORLD,IERR)

  ! Copy solution vector X into a PETSC vector VX
  CALL VecCreateMPI(MPI_COMM_WORLD,N,PETSC_DETERMINE,VX,IERR)
  CALL VecSetValues(VX,N,GIND,X,INSERT_VALUES,IERR)
  CALL VecAssemblyBegin(VX,IERR)
  CALL VecAssemblyEnd(VX,IERR)
  !CALL VecView(VX,PETSC_VIEWER_STDOUT_WORLD,IERR)
  ! Initialize the linear solver
  ! It is matrix-free (we provide the matvec), and we provide the preconditioner



  BTIME=MPI_WTIME()  ;  TIMINGS(2)=BTIME-ATIME;   ATIME=BTIME

  CALL KSPCreate(MPI_COMM_WORLD,solver,IERR)
  CALL MatCreateShell(MPI_COMM_WORLD,N,N,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_NULL_INTEGER,A,IERR)
  CALL MatShellSetOperation(A,MATOP_MULT,WrapMatVec,IERR)
  !CALL KSPSetOperators(solver,A,A,DIFFERENT_NONZERO_PATTERN,IERR) ! PETSc       <3.5
  CALL KSPSetOperators(solver,A,A,IERR)

  IF(PC_GL.NE.0) THEN
     CALL KSPGetPC(solver,P,IERR)
     CALL PCSetType(P,PCSHELL,IERR)
     CALL PCShellSetApply(P,WrapPrec,IERR)
  ENDIF


  BTIME=MPI_WTIME()  ;  TIMINGS(3)=BTIME-ATIME;   ATIME=BTIME

  
  ! Tune the linear solver
  ! We use GMRES and we set some parameters
  CALL KSPSetType(solver,KSPGMRES,IERR)

!! 
#ifdef PETSC_LOW_VERFLAG
  CALL KSPSetTolerances(solver,GMSOLVER%TOL,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION, &
                        GMSOLVER%MAX_INNER_ITER,IERR) ! PETSc <3.5
#else
  CALL KSPSetTolerances(solver,GMSOLVER%TOL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,GMSOLVER%MAX_INNER_ITER,IERR)
#endif


  IF(GMSOLVER%IUNIT.NE.0) THEN
    CALL KSPMonitorSet(solver,KSPMonitorDefault,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,IERR)
  ENDIF


  BTIME=MPI_WTIME()  ;  TIMINGS(4)=BTIME-ATIME;   ATIME=BTIME


  ! Solve
  CALL KSPSolve(solver,VB,VX,IERR)

  BTIME=MPI_WTIME()  ;  TIMINGS(5)=BTIME-ATIME;   ATIME=BTIME

  CALL KSPGetIterationNumber(solver,iter,IERR)
  CALL KSPGetResidualNorm(solver,rnorm,IERR)
  ! Get solution from solver and copy back into X
  CALL KSPGetSolution(solver,VX,IERR)
  CALL VecGetValues(VX,N,GIND,X,IERR)

  BTIME=MPI_WTIME()  ;  TIMINGS(6)=BTIME-ATIME;   ATIME=BTIME
  
  ! Release stuff
  CALL VecDestroy(VB,IERR)
  CALL VecDestroy(VX,IERR)
  CALL MatDestroy(A,IERR)
  CALL KSPDestroy(solver,IERR)

  BTIME=MPI_WTIME()  ;  TIMINGS(7)=BTIME-ATIME;   ATIME=BTIME

!  IF (ME.EQ.0.AND.MOD(CALLEDFLAG,10).EQ.OUTMODULUS) THEN
!     PRINT *, "PETSCSOLVE TIMINGS"
!     WRITE(*,'(100A15)') "INIT", "VEC" ,"MAT", "SET", "SOLVE","GET","DESTROY"
!     WRITE(*,'(100F10.5)') TIMINGS(1:7)
!  ENDIF



END SUBROUTINE PETSCSOLVE

END MODULE PETSCWRAP

#endif
