#include "definitions.INC"

!--------------------------------------------------------------------------------------

MODULE SNEAK
USE STRUCTS
IMPLICIT NONE

TYPE(HAMILTONIAN), POINTER :: HAM_GL
TYPE(PARAMETERS), POINTER :: PARAMS_GL
TYPE(MPI_DATA), POINTER :: MPIDATA_GL
NUMBER, POINTER :: SHIFT_GL
INTEGER, POINTER :: PC_GL

END MODULE SNEAK

!--------------------------------------------------------------------------------------

MODULE SOLVER
USE HAMMOD
USE STRUCTS

CONTAINS


!--------------------------------------------------------------------------------------

#if ISREAL == 0
SUBROUTINE PCWRAPPER(HAM,PARAMS,MPIDATA,SHIFT,M,MB,X,PX,PRECONDITIONER)
USE SNEAK

TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
TYPE(PARAMETERS), TARGET, INTENT(IN) :: PARAMS
TYPE(MPI_DATA), TARGET, INTENT(IN) :: MPIDATA
NUMBER, TARGET, INTENT(IN) :: SHIFT
INTEGER, INTENT(IN) :: M,MB
NUMBER, INTENT(IN) :: X(1:M*MB)
NUMBER, INTENT(OUT) :: PX(1:M*MB)
EXTERNAL :: PRECONDITIONER

REAL *8 :: XR(1:M*MB,1:2), PXR(1:M*MB,1:2)
INTEGER :: NELT,IA,JA,ISYM,IWORK
REAL *8 :: A,RWORK

HAM_GL => HAM
PARAMS_GL => PARAMS
MPIDATA_GL => MPIDATA
SHIFT_GL => SHIFT

XR(:,1) = DREAL(X)
XR(:,2) = DIMAG(X)
CALL PRECONDITIONER(2*M*MB,XR,PXR,NELT,IA,JA,A,ISYM,RWORK,IWORK)
PX = DCMPLX(PXR(:,1),PXR(:,2))

END SUBROUTINE PCWRAPPER
#endif

!--------------------------------------------------------------------------------------

SUBROUTINE MYGMRES(B,X,N,HAM,PARAMS,MPIDATA,GMSOLVER)
USE SNEAK

! SOLVES (H-SHIFT)*X = B USING GMRES

TYPE(GMRES_SOLVER), TARGET, INTENT(INOUT) :: GMSOLVER
TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
TYPE(PARAMETERS), TARGET, INTENT(IN) :: PARAMS
TYPE(MPI_DATA), TARGET, INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: N 
NUMBER, INTENT(IN) :: B(1:N)
NUMBER, INTENT(INOUT) :: X(1:N)

INTEGER :: NN,I,MPIERR
INTEGER :: IUNIT,LRGW,IGWK(7),LIGW,IWORK,NELT,IA,JA,ISYM,ITOL
INTEGER, POINTER :: MAXRST,PCFLAG,ITER,ERRFLAG,PCOPT
INTEGER, POINTER :: MAX_INNER_ITER, MAX_OUTER_ITER
REAL *8, POINTER :: TOL,ERR
REAL *8 :: A,SX,SY
REAL *8, ALLOCATABLE :: RGWK(:), RWORK(:), BREAL(:,:), XREAL(:,:)
NUMBER, POINTER :: SHIFT
PROCEDURE(), POINTER :: MSOLVE=>NULL(), MVEC=>NULL()

MAX_INNER_ITER => GMSOLVER % MAX_INNER_ITER
MAX_OUTER_ITER => GMSOLVER % MAX_OUTER_ITER
MAXRST => GMSOLVER % MAXRST
PCFLAG => GMSOLVER % PCFLAG
ITER => GMSOLVER % ITER
ERRFLAG => GMSOLVER % ERRFLAG
TOL => GMSOLVER % TOL
ERR => GMSOLVER % ERR
SHIFT => GMSOLVER % SHIFT
PCOPT => GMSOLVER % PCOPT

! DEFINE ARRAY SIZE BASED ON DATATYPE
#if ISREAL == 1
   NN = N
   MVEC => MATVECREAL
   MSOLVE => MSOLVEREAL1
   IF(PCOPT == 1) THEN
      MSOLVE => MSOLVEREAL1
   ELSE
      MSOLVE => MSOLVEREAL2
   END IF
#else
   NN = 2*N
   MVEC => MATVECCPLX
   MSOLVE => MSOLVECPLX1
   IF(PCOPT == 1) THEN
      MSOLVE => MSOLVECPLX1
   ELSEIF(PCOPT == 2) THEN
      MSOLVE => MSOLVECPLX2
   ELSEIF(PCOPT == 3) THEN
      MSOLVE => MSOLVECPLX3
#if defined(USE_SUPERLU)
   ELSEIF(PCOPT == 4) THEN
      MSOLVE => MSOLVECPLX4
   ELSEIF(PCOPT == 5) THEN
      MSOLVE => MSOLVECPLX5
#endif
   END IF
#endif

! INITIALIZE GLOBAL POINTERS FOR MATVEC AND MSOLVE
HAM_GL => HAM
PARAMS_GL => PARAMS
MPIDATA_GL => MPIDATA
SHIFT_GL => SHIFT

! DEFINE POINTERS TO MATVEC AND MSOLVE ROUTINES




! I/O UNIT FOR OUTPUT:
! SET TO 0 TO SUPRESS OUTPUT
! SET TO 6 FOR STANDARD OUTPUT

IF(MPIDATA % ME == 0) THEN
   IUNIT = GMSOLVER % IUNIT
ELSE
   IUNIT = 0
END IF

NELT = 0
ITOL = 0
LIGW = 7

IGWK = 0
IGWK(1) = MAX_INNER_ITER ! HOW OFTEN TO RESTART
IGWK(2) = MAX_INNER_ITER ! NUMBER OF KRYLOV VECTOS TO ORTHOGONALIZE ON RESTART
IGWK(3) = 0
IGWK(4) = PCFLAG ! PRECONDITIONING FLAG
IGWK(5) = MAXRST ! MAXIMUM NUMBER OF ALLOWED RESTARTS

! SIZE OF WORK ARRAY
LRGW = 1 + NN*(IGWK(1)+6) + IGWK(1)*(IGWK(1)+3)

ALLOCATE(RGWK(LRGW))

#if ISREAL == 0

   ALLOCATE(BREAL(1:N,1:2),XREAL(1:N,1:2))
   BREAL(:,1) = DREAL(B)
   BREAL(:,2) = DIMAG(B)

   CALL DGMRES(NN,BREAL,XREAL,NELT,IA,JA,A,ISYM,MVEC, &
   MSOLVE,ITOL,TOL,MAX_OUTER_ITER,ITER,ERR,ERRFLAG,IUNIT,SX,SY, &
   RGWK,LRGW,IGWK,LIGW,RWORK,IWORK)

   X = DCMPLX(XREAL(:,1),XREAL(:,2))
   DEALLOCATE(BREAL,XREAL)

#else

   CALL DGMRES(NN,B,X,NELT,IA,JA,A,ISYM,MVEC, &
   MSOLVE,ITOL,TOL,MAX_OUTER_ITER,ITER,ERR,ERRFLAG,IUNIT,SX,SY, &
   RGWK,LRGW,IGWK,LIGW,RWORK,IWORK)

#endif

DEALLOCATE(RGWK)

END SUBROUTINE MYGMRES

!--------------------------------------------------------------------------------------

#if ISREAL

SUBROUTINE MATVECREAL(N,X,Y,NELT,IA,JA,A,ISYM)
USE SNEAK

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM
REAL *8, INTENT(IN) :: X(1:N), A(*)
REAL *8, INTENT(OUT) :: Y(1:N)
INTEGER :: I

CALL HAMMULT(HAM_GL,PARAMS_GL,MPIDATA_GL,X,Y,N)
Y = Y - SHIFT_GL*X

END SUBROUTINE MATVECREAL

!--------------------------------------------------------------------------------------

SUBROUTINE MSOLVEREAL1(N,R,Z,NELT,IA,JA,A,ISYM,W,IWORK)
USE SNEAK

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
NUMBER, INTENT(IN) :: R(1:N), W(*), A(*)
NUMBER, INTENT(OUT) :: Z(1:N)
REAL *8, ALLOCATABLE, SAVE :: TI(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
INTEGER :: NN, LWORK

NN = PARAMS_GL % NN
LWORK = 10*NN

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(TI(NN,NN))
   CALL KINETIC_INVERSE(HAM_GL % T,TI,NN,LWORK,SHIFT_GL)
END IF

IF(PARAMS_GL % NELEC == 1) THEN
   CALL KEMULT1(TI,R,Z,NN,MPIDATA_GL%FIRST,MPIDATA_GL%LAST,MPIDATA_GL)
ELSEIF(PARAMS_GL % NELEC == 2) THEN
   CALL KEMULT2(TI,R,Z,NN,PARAMS_GL % M, MPIDATA_GL % MBLK,MPIDATA_GL)
ELSEIF(PARAMS_GL % NELEC == 3) THEN
   CALL KEMULT3(TI,R,Z,NN,PARAMS_GL % M, MPIDATA_GL % MBLK,MPIDATA_GL)
END IF

END SUBROUTINE MSOLVEREAL1

!--------------------------------------------------------------------------------------

SUBROUTINE MSOLVEREAL2(N,R,Z,NELT,IA,JA,A,ISYM,RWORK,IWORK)
USE SNEAK

! CAUTION: WORKS ON ONE PROCESS ONLY

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: RWORK(*), A(*)
REAL *8, TARGET, INTENT(IN) :: R(1:N)
REAL *8, TARGET, INTENT(OUT) :: Z(1:N)
REAL *8, ALLOCATABLE, SAVE :: W(:,:), WT(:,:), ALLEIGS(:), TEMP(:)
REAL *8, ALLOCATABLE, SAVE :: U(:,:), V(:,:), U1(:), BIGEIGS(:,:)
REAL *8, POINTER :: RP(:,:), ZP(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
INTEGER :: NN, M, MM, LWORK, MB, I, J, NELEC, FIRST, LAST

NN = PARAMS_GL % NN
M = PARAMS_GL % M
MM = MPIDATA_GL % MM
FIRST = MPIDATA_GL % FIRST
LAST = MPIDATA_GL % LAST
MB = MPIDATA_GL % MBLK
NELEC = PARAMS_GL % NELEC

LWORK = 10*NN

IF(FIRSTTIME) THEN

   FIRSTTIME = .FALSE.
   ALLOCATE(TEMP(1:MM),W(1:NN,1:NN),WT(1:NN,1:NN))
   
   IF(NELEC == 1) THEN
   
   	ALLOCATE(ALLEIGS(FIRST:LAST))
      CALL KINETIC_EIGS(HAM_GL % T,ALLEIGS,W,WT,PARAMS_GL % THETA,NN,LWORK,FIRST,LAST)
      
   ELSEIF(NELEC == 2) THEN
   
      ALLOCATE(ALLEIGS(1:M))
      ALLOCATE(U(1:M,1:MB),V(1:M,1:MB),U1(1:M))
      ALLOCATE(BIGEIGS(1:M,FIRST:LAST))
      CALL KINETIC_EIGS(HAM_GL % T,ALLEIGS,W,WT,PARAMS_GL % THETA,NN,LWORK,1,M)
      
      DO J = FIRST,LAST
		DO I = 1,M
			BIGEIGS(I,J) = ALLEIGS(I) + ALLEIGS(J)
		END DO
      END DO
      
      !DEALLOCATE(ALLEIGS)
      
   END IF
      
END IF

IF(NELEC == 1) THEN

   CALL EIGMULT(R,TEMP,WT,NN)
   CALL KRONMULT(R,TEMP,WT,WT,WT,NN,MB,FIRST,LAST,MPIDATA_GL)
   TEMP = TEMP/(ALLEIGS - SHIFT_GL)
   CALL KRONMULT(TEMP,Z,W,W,W,NN,MB,FIRST,LAST,MPIDATA_GL)
   CALL EIGMULT(TEMP,Z,W,NN)

ELSEIF(NELEC == 2) THEN

   RP(1:M,1:MB) => R
   ZP(1:M,1:MB) => Z

   !CALL BIGEIGMULT(R,U1,WT,NN)
   !U1 = U1/(BIGEIGS - SHIFT_GL)
   !CALL BIGEIGMULT(U1,Z,W,NN)

   DO J = 1,MB
      CALL EIGMULT(RP(:,J),U1,WT,NN)
      U1 = U1/(ALLEIGS - SHIFT_GL)
      CALL EIGMULT(U1,ZP(:,J),W,NN)
   END DO

   CALL MPITRANSPOSE(RP,U,M,MB,MPIDATA_GL)

   DO J = 1,MB
      CALL EIGMULT(U(:,J),U1,WT,NN)
      U1 = U1/(ALLEIGS - SHIFT_GL)
      CALL EIGMULT(U1,V(:,J),W,NN)
   END DO

   CALL MPITRANSPOSE(V,U,M,MB,MPIDATA_GL)

   ZP = ZP + U

END IF

END SUBROUTINE MSOLVEREAL2

!--------------------------------------------------------------------------------------

#else

!--------------------------------------------------------------------------------------

SUBROUTINE MATVECCPLX(N,X,Y,NELT,IA,JA,A,ISYM)
USE SNEAK

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM
REAL *8, INTENT(IN) :: A(*)
REAL *8, INTENT(IN) :: X(1:N/2,1:2)
REAL *8, INTENT(OUT) :: Y(1:N/2,1:2)
COMPLEX *16 :: XX(1:N/2), YY(1:N/2)
INTEGER :: I, M

M = N/2
XX = DCMPLX(X(:,1),X(:,2))
CALL HAMMULT(HAM_GL,PARAMS_GL,MPIDATA_GL,XX,YY,M)

YY = YY - SHIFT_GL*XX
Y(:,1) = DREAL(YY)
Y(:,2) = DIMAG(YY)

END SUBROUTINE MATVECCPLX

!--------------------------------------------------------------------------------------

SUBROUTINE MSOLVECPLX1(N,R,Z,NELT,IA,JA,A,ISYM,W,IWORK)
USE SNEAK

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: R(1:N/2,1:2), W(*), A(*)
REAL *8, INTENT(OUT) :: Z(1:N/2,1:2)
COMPLEX *16, ALLOCATABLE, SAVE :: TI(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
INTEGER :: NN, LWORK
COMPLEX *16 :: RR(1:N/2), ZZ(1:N/2)

NN = PARAMS_GL % NN
LWORK = 10*NN

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(TI(NN,NN))
   CALL KINETIC_INVERSE(HAM_GL % T,TI,NN,LWORK,SHIFT_GL)
END IF


RR = DCMPLX(R(:,1),R(:,2))

IF(PARAMS_GL % NELEC == 1) THEN
   CALL KEMULT1(TI,RR,ZZ,NN,MPIDATA_GL%FIRST,MPIDATA_GL%LAST,MPIDATA_GL)
ELSEIF(PARAMS_GL % NELEC == 2) THEN
   CALL KEMULT2(TI,RR,ZZ,NN,PARAMS_GL % M, MPIDATA_GL % MBLK,MPIDATA_GL)
ELSEIF(PARAMS_GL % NELEC == 3) THEN
   CALL KEMULT3(TI,RR,ZZ,NN,PARAMS_GL % M, MPIDATA_GL % MBLK,MPIDATA_GL)
END IF

Z(:,1) = DREAL(ZZ)
Z(:,2) = DIMAG(ZZ)

END SUBROUTINE MSOLVECPLX1

SUBROUTINE MSOLVECPLX12(RR,ZZ,N)
USE SNEAK
! Simplified interface for PETSc

INTEGER, INTENT(IN) :: N
COMPLEX*16, INTENT(IN) :: RR(N)
COMPLEX*16, INTENT(OUT) :: ZZ(N)
COMPLEX *16, ALLOCATABLE, SAVE :: TI(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
INTEGER :: NN, LWORK

NN = PARAMS_GL % NN
LWORK = 10*NN

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(TI(NN,NN))
   CALL KINETIC_INVERSE(HAM_GL % T,TI,NN,LWORK,SHIFT_GL)
END IF


IF(PARAMS_GL % NELEC == 1) THEN
   CALL KEMULT1(TI,RR,ZZ,NN,MPIDATA_GL%FIRST,MPIDATA_GL%LAST,MPIDATA_GL)
ELSEIF(PARAMS_GL % NELEC == 2) THEN
   CALL KEMULT2(TI,RR,ZZ,NN,PARAMS_GL % M, MPIDATA_GL % MBLK,MPIDATA_GL)
ELSEIF(PARAMS_GL % NELEC == 3) THEN
   CALL KEMULT3(TI,RR,ZZ,NN,PARAMS_GL % M, MPIDATA_GL % MBLK,MPIDATA_GL)
END IF

END SUBROUTINE MSOLVECPLX12

!--------------------------------------------------------------------------------------

SUBROUTINE MSOLVECPLX2(N,R,Z,NELT,IA,JA,A,ISYM,RWORK,IWORK)
USE SNEAK

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: R(1:N/2,1:2), RWORK(*), A(*)
REAL *8, INTENT(OUT) :: Z(1:N/2,1:2)
REAL *8, ALLOCATABLE, SAVE :: W(:,:), WT(:,:)
COMPLEX *16, ALLOCATABLE, SAVE :: TEMP(:), ALLEIGS(:)
COMPLEX *16, ALLOCATABLE, SAVE :: U(:,:), V(:,:), U1(:)
COMPLEX *16, POINTER :: RP(:,:), ZP(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
INTEGER :: NN, M, MM, LWORK, MB, NELEC, J
COMPLEX *16, TARGET :: RR(1:N/2), ZZ(1:N/2)
COMPLEX *16 :: FAC
REAL *8 :: PCTHETA

NELEC = PARAMS_GL % NELEC
NN = PARAMS_GL % NN
M = PARAMS_GL % M
MB = MPIDATA_GL % MBLK
MM = MPIDATA_GL % MM
LWORK = 10*NN

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(ALLEIGS(1:M),TEMP(1:MM))
   ALLOCATE(W(1:NN,1:NN),WT(1:NN,1:NN))
   CALL KINETIC_EIGS(HAM_GL % T,ALLEIGS,W,WT,PARAMS_GL % THETA,NN,LWORK,1,M)

   IF(NELEC == 2) ALLOCATE(U(1:M,1:MB),V(1:M,1:MB),U1(1:M))

END IF

RR = DCMPLX(R(:,1),R(:,2))

IF(NELEC == 1) THEN

   CALL EIGMULT(RR,TEMP,WT,NN)
   TEMP = TEMP/(ALLEIGS - SHIFT_GL)
   CALL EIGMULT(TEMP,ZZ,W,NN)

ELSEIF(NELEC == 2) THEN

   RP(1:M,1:MB) => RR
   ZP(1:M,1:MB) => ZZ

   DO J = 1,MB
      CALL EIGMULT(RP(:,J),U1,WT,NN)
      U1 = U1/(ALLEIGS - SHIFT_GL)
      CALL EIGMULT(U1,ZP(:,J),W,NN)
   END DO

   CALL MPITRANSPOSE(RP,U,M,MB,MPIDATA_GL)

   DO J = 1,MB
      CALL EIGMULT(U(:,J),U1,WT,NN)
      U1 = U1/(ALLEIGS - SHIFT_GL)
      CALL EIGMULT(U1,V(:,J),W,NN)
   END DO

   CALL MPITRANSPOSE(V,U,M,MB,MPIDATA_GL)

   ZP = ZP + U

END IF

Z(:,1) = DREAL(ZZ)
Z(:,2) = DIMAG(ZZ)

END SUBROUTINE MSOLVECPLX2

SUBROUTINE MSOLVECPLX22(RR,ZZ,N)
USE SNEAK

! Simplified interface for PETSc

INTEGER, INTENT(IN) :: N
REAL *8, ALLOCATABLE, SAVE :: W(:,:), WT(:,:)
COMPLEX *16, ALLOCATABLE, SAVE :: TEMP(:), ALLEIGS(:)
COMPLEX *16, ALLOCATABLE, SAVE :: U(:,:), V(:,:), U1(:)
COMPLEX *16, POINTER :: RP(:,:), ZP(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
INTEGER :: NN, M, MM, LWORK, MB, NELEC, J
COMPLEX *16, TARGET :: RR(N), ZZ(N)
COMPLEX *16 :: FAC
REAL *8 :: PCTHETA

NELEC = PARAMS_GL % NELEC
NN = PARAMS_GL % NN
M = PARAMS_GL % M
MB = MPIDATA_GL % MBLK
MM = MPIDATA_GL % MM
LWORK = 10*NN

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(ALLEIGS(1:M),TEMP(1:MM))
   ALLOCATE(W(1:NN,1:NN),WT(1:NN,1:NN))
   CALL KINETIC_EIGS(HAM_GL % T,ALLEIGS,W,WT,PARAMS_GL % THETA,NN,LWORK,1,M)

   IF(NELEC == 2) ALLOCATE(U(1:M,1:MB),V(1:M,1:MB),U1(1:M))

END IF

IF(NELEC == 1) THEN

   CALL EIGMULT(RR,TEMP,WT,NN)
   TEMP = TEMP/(ALLEIGS - SHIFT_GL)
   CALL EIGMULT(TEMP,ZZ,W,NN)

ELSEIF(NELEC == 2) THEN

   RP(1:M,1:MB) => RR
   ZP(1:M,1:MB) => ZZ

   DO J = 1,MB
      CALL EIGMULT(RP(:,J),U1,WT,NN)
      U1 = U1/(ALLEIGS - SHIFT_GL)
      CALL EIGMULT(U1,ZP(:,J),W,NN)
   END DO

   CALL MPITRANSPOSE(RP,U,M,MB,MPIDATA_GL)

   DO J = 1,MB
      CALL EIGMULT(U(:,J),U1,WT,NN)
      U1 = U1/(ALLEIGS - SHIFT_GL)
      CALL EIGMULT(U1,V(:,J),W,NN)
   END DO

   CALL MPITRANSPOSE(V,U,M,MB,MPIDATA_GL)

   ZP = ZP + U

END IF

END SUBROUTINE MSOLVECPLX22

!--------------------------------------------------------------------------------------

SUBROUTINE MSOLVECPLX3old(N,R,Z,NELT,IA,JA,A,ISYM,RWORK,IWORK)
USE SNEAK
! Block-diagonal preconditioner based on the NNxNN diagonal blocks
! of the whole problem.
!
! This version uses the "all-in-one" ZGESV routine to solve for each
! NNxNN diagonal block. This is cheap memory-wise but there are
! many redundant operations (factorization of the same matrix at each
! call/iteration).
!

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: R(1:N/2,1:2)
REAL *8, INTENT(IN):: A(*), RWORK(*)
REAL *8, INTENT(OUT) :: Z(1:N/2,1:2)

INTEGER :: I,K,I1,I2,INFO,NN,MM,NELEC,NBLOCKS
COMPLEX *16 :: RR(1:N/2), ZZ(1:N/2)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE, SAVE :: DD(:,:)
INTEGER, ALLOCATABLE, SAVE :: IPIV(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
MM = MPIDATA_GL % MM

T => HAM_GL % T
V => HAM_GL % V
RR = DCMPLX(R(:,1),R(:,2))

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(DD(NN,NN),IPIV(NN))
END IF

NBLOCKS=NN**(3*NELEC-1)

DO K = 1,NBLOCKS

   I1 = (K-1)*NN + 1
   I2 = K*NN
   V => HAM_GL % V(I1:I2)

   DD = T
   DO I = 1,NN
      DD(I,I) = 3D0*NELEC*T(I,I) + V(I) - SHIFT_GL
   END DO

   ZZ(I1:I2) = RR(I1:I2)
   CALL ZGESV(NN,1,DD,NN,IPIV,ZZ(I1:I2),NN,INFO)
   Z(I1:I2,1) = DREAL(ZZ(I1:I2))
   Z(I1:I2,2) = DIMAG(ZZ(I1:I2))

END DO

END SUBROUTINE MSOLVECPLX3old

!--------------------------------------------------------------------------------------

SUBROUTINE MSOLVECPLX3(N,R,Z,NELT,IA,JA,A,ISYM,RWORK,IWORK)
USE SNEAK
USE MPI
! Block-diagonal preconditioner based on the NNxNN diagonal blocks
! of the whole problem.
!
! This version uses ZGETRF (factorization) + ZGETRS (triangular solve)
! instead of the "all-in-one" ZGESV routine. The factors are stored
! and reused instead of being recomputed at each iteration.
! This should be faster:  (NN^(NELEC-1)*NN^3+#iter*NN^(NELEC-1)N^2 flops
! intead of #iter*NN^(NELEC-1)*NN^3) but it uses much more memory:
! NN^(NELEC-1)**NN^2 instead of NN^2.
!

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: R(1:N/2,1:2)
REAL *8, INTENT(IN):: A(*), RWORK(*)
REAL *8, INTENT(OUT) :: Z(1:N/2,1:2)

INTEGER :: I,K,I1,I2,INFO,NN,MM,NELEC,IERR,NBLOCKS,NDIM,NP,ME
COMPLEX *16 :: RR(1:N/2), ZZ(1:N/2)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE, SAVE :: DD(:,:)
INTEGER, ALLOCATABLE, SAVE :: IPIV(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
MM = MPIDATA_GL % MM
NP = MPIDATA_GL % NP
ME = MPIDATA_GL % ME

T => HAM_GL % T
V => HAM_GL % V
RR = DCMPLX(R(:,1),R(:,2))

NBLOCKS=NN**(3*NELEC-1)/NP
NDIM=NBLOCKS*NN

IF(FIRSTTIME) THEN
  IF(ME.EQ.0) WRITE(*,*)'1st-level preconditioner'
  ALLOCATE(DD(NDIM,NN),IPIV(NDIM))
END IF

DO K = 1,NBLOCKS
  I1 = (K-1)*NN + 1
  I2 = K*NN

  IF(FIRSTTIME) THEN
    V => HAM_GL % V(I1:I2)
    CALL ZLACPY('N',NN,NN,T,NN,DD(I1,1),NDIM)
    DO I = 1,NN
      DD(I1+I-1,I) = 3D0*NELEC*T(I,I) + V(I) - SHIFT_GL
    END DO
    CALL ZGETRF(NN,NN,DD(I1,1),NDIM,IPIV(I1),IERR)
  END IF

  ZZ(I1:I2) = RR(I1:I2)
  CALL ZGETRS('N',NN,1,DD(I1,1),NDIM,IPIV(I1),ZZ(I1:I2),NN,IERR)
  Z(I1:I2,1) = DREAL(ZZ(I1:I2))
  Z(I1:I2,2) = DIMAG(ZZ(I1:I2))

END DO

! Subsequent iterations dont recompute factors    
FIRSTTIME = .FALSE.

END SUBROUTINE MSOLVECPLX3

SUBROUTINE MSOLVECPLX32(XX,ZZ,N)
USE SNEAK
USE MPI
! Block-diagonal preconditioner based on the NNxNN diagonal blocks
! of the whole problem.
!
! This version uses ZGETRF (factorization) + ZGETRS (triangular solve)
! instead of the "all-in-one" ZGESV routine. The factors are stored
! and reused instead of being recomputed at each iteration.
! This should be faster:  (NN^(NELEC-1)*NN^3+#iter*NN^(NELEC-1)N^2 flops
! intead of #iter*NN^(NELEC-1)*NN^3) but it uses much more memory:
! NN^(NELEC-1)**NN^2 instead of NN^2.
!
! Simplified interface for PETSc
!

INTEGER, INTENT(IN) :: N

INTEGER :: I,K,I1,I2,INFO,NN,MM,NELEC,IERR,NBLOCKS,NDIM,NP,ME
COMPLEX *16 :: XX(N),ZZ(N)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE, SAVE :: DD(:,:)
INTEGER, ALLOCATABLE, SAVE :: IPIV(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
MM = MPIDATA_GL % MM
NP = MPIDATA_GL % NP
ME = MPIDATA_GL % ME

T => HAM_GL % T
V => HAM_GL % V

NBLOCKS=NN**(3*NELEC-1)/NP
NDIM=NBLOCKS*NN

IF(FIRSTTIME) THEN
  IF(ME.EQ.0) WRITE(*,*)'1st-level preconditioner'
  ALLOCATE(DD(NDIM,NN),IPIV(NDIM))
END IF

ZZ=XX

DO K = 1,NBLOCKS
  I1 = (K-1)*NN + 1
  I2 = K*NN

  IF(FIRSTTIME) THEN
    V => HAM_GL % V(I1:I2)
    CALL ZLACPY('N',NN,NN,T,NN,DD(I1,1),NDIM)
    DO I = 1,NN
      DD(I1+I-1,I) = 3D0*NELEC*T(I,I) + V(I) - SHIFT_GL
    END DO
    CALL ZGETRF(NN,NN,DD(I1,1),NDIM,IPIV(I1),IERR)
  END IF

  CALL ZGETRS('N',NN,1,DD(I1,1),NDIM,IPIV(I1),ZZ(I1:I2),NN,IERR)

END DO

! Subsequent iterations dont recompute factors    
FIRSTTIME = .FALSE.

END SUBROUTINE MSOLVECPLX32

!--------------------------------------------------------------------------------------

#if defined(USE_SUPERLU)

SUBROUTINE MSOLVECPLX4(N,R,Z,NELT,IA,JA,A,ISYM,RWORK,IWORK)
USE SNEAK
USE SUPERLU_MOD
! Block-diagonal preconditioner based on the NN^2xNN^2 diagonal blocks
! of the whole problem. Each block is sparse and the solution for a
! block is computed with SuperLU.
!
! This version performs factorization+solve at each call; it doesnt
! reuse any information between iterations. This is cheap memory-wise
! but performs a lot of flops.
! 

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: R(1:N/2,1:2), RWORK(*), A(*)
REAL *8, INTENT(OUT) :: Z(1:N/2,1:2)

INTEGER :: I,IP,J,JP,K,I1,I2,B1,B2,INFO,NN,NELEC,NBLOCKS,NNZ,NP,ME
INTEGER, ALLOCATABLE :: II(:), JJ(:)
COMPLEX *16 :: ZZ(1:N/2)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE :: TIJ(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
COMPLEX *16 :: T1, T2

INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: GRID(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: OPTIONS(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SCALEPERMSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: LUSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SOLVESTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: ASLU(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: STAT(:)

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
NP = MPIDATA_GL % NP
ME = MPIDATA_GL % ME
T => HAM_GL % T
ZZ = DCMPLX(R(:,1),R(:,2))

NBLOCKS=NN**(3*NELEC-2)/NP
IF(FIRSTTIME) THEN
  ALLOCATE(GRID(NBLOCKS))
  ALLOCATE(OPTIONS(NBLOCKS))
  ALLOCATE(SCALEPERMSTRUCT(NBLOCKS))
  ALLOCATE(LUSTRUCT(NBLOCKS))
  ALLOCATE(SOLVESTRUCT(NBLOCKS))
  ALLOCATE(ASLU(NBLOCKS))
  ALLOCATE(STAT(NBLOCKS))
END IF

! Allocate data for the sparse matrix
! Each block has the same size and number of elements
NNZ=NN**2*(2*NN-1)
IF(FIRSTTIME) THEN
  IF(ME.EQ.0) WRITE(*,*)'2nd-level preconditioner'
  ALLOCATE(II(NNZ*NBLOCKS))
  ALLOCATE(JJ(NNZ*NBLOCKS))
  ALLOCATE(TIJ(NNZ*NBLOCKS))
END IF

DO K = 1,NBLOCKS
  !IF((K+1)*10/NBLOCKS>K*10/NBLOCKS) write(*,'(I4,A)',ADVANCE='NO')(K+1)*100/NBLOCKS,'%'

  I1 = (K-1)*NN**2 + 1
  I2 = K*NN**2
  B1 = (K-1)*NNZ + 1
  B2 = K*NNZ

  ! Compute the values of the sparse matrix
  IF(FIRSTTIME) THEN
    CALL N2BLOCKS(I1,I2,NNZ,II(B1:B2),JJ(B1:B2),TIJ(B1:B2))
  END IF
  
  ! Solve with SuperLU_DIST
  CALL SUPERLUSOLVE(II(B1:B2),JJ(B1:B2),TIJ(B1:B2),NNZ,ZZ(I1:I2),NN**2, &
       GRID(K),OPTIONS(K),SCALEPERMSTRUCT(K),LUSTRUCT(K),SOLVESTRUCT(K),&
       ASLU(K),STAT(K),FIRSTTIME)

  ! Put back in solution vector  
  Z(I1:I2,1) = DREAL(ZZ(I1:I2))
  Z(I1:I2,2) = DIMAG(ZZ(I1:I2))

END DO

Z(:,1) = DREAL(ZZ)
Z(:,2) = DIMAG(ZZ)

! Subsequent iterations dont recompute factors    
FIRSTTIME = .FALSE.

END SUBROUTINE MSOLVECPLX4

SUBROUTINE MSOLVECPLX42(XX,ZZ,N)
USE SNEAK
USE SUPERLU_MOD
! Block-diagonal preconditioner based on the NN^2xNN^2 diagonal blocks
! of the whole problem. Each block is sparse and the solution for a
! block is computed with SuperLU.
!
! This version performs factorization+solve at each call; it doesnt
! reuse any information between iterations. This is cheap memory-wise
! but performs a lot of flops.
!
! Simplified interface for PETSc
! 

INTEGER, INTENT(IN) :: N

INTEGER :: I,IP,J,JP,K,I1,I2,B1,B2,INFO,NN,NELEC,NBLOCKS,NNZ,NP,ME
INTEGER, ALLOCATABLE :: II(:), JJ(:)
COMPLEX *16 :: XX(N), ZZ(N)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE :: TIJ(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
COMPLEX *16 :: T1, T2

INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: GRID(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: OPTIONS(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SCALEPERMSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: LUSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SOLVESTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: ASLU(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: STAT(:)

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
NP = MPIDATA_GL % NP
ME = MPIDATA_GL % ME
T => HAM_GL % T

NBLOCKS=NN**(3*NELEC-2)/NP
IF(FIRSTTIME) THEN
  ALLOCATE(GRID(NBLOCKS))
  ALLOCATE(OPTIONS(NBLOCKS))
  ALLOCATE(SCALEPERMSTRUCT(NBLOCKS))
  ALLOCATE(LUSTRUCT(NBLOCKS))
  ALLOCATE(SOLVESTRUCT(NBLOCKS))
  ALLOCATE(ASLU(NBLOCKS))
  ALLOCATE(STAT(NBLOCKS))
END IF

ZZ=XX

! Allocate data for the sparse matrix
! Each block has the same size and number of elements
NNZ=NN**2*(2*NN-1)
IF(FIRSTTIME) THEN
  IF(ME.EQ.0) WRITE(*,*)'2nd-level preconditioner'
  ALLOCATE(II(NNZ*NBLOCKS))
  ALLOCATE(JJ(NNZ*NBLOCKS))
  ALLOCATE(TIJ(NNZ*NBLOCKS))
END IF

DO K = 1,NBLOCKS
  !IF((K+1)*10/NBLOCKS>K*10/NBLOCKS) write(*,'(I4,A)',ADVANCE='NO')(K+1)*100/NBLOCKS,'%'

  I1 = (K-1)*NN**2 + 1
  I2 = K*NN**2
  B1 = (K-1)*NNZ + 1
  B2 = K*NNZ

  ! Compute the values of the sparse matrix
  IF(FIRSTTIME) THEN
    CALL N2BLOCKS(I1,I2,NNZ,II(B1:B2),JJ(B1:B2),TIJ(B1:B2))
  END IF
  
  ! Solve with SuperLU_DIST
  CALL SUPERLUSOLVE(II(B1:B2),JJ(B1:B2),TIJ(B1:B2),NNZ,ZZ(I1:I2),NN**2, &
       GRID(K),OPTIONS(K),SCALEPERMSTRUCT(K),LUSTRUCT(K),SOLVESTRUCT(K),&
       ASLU(K),STAT(K),FIRSTTIME)

END DO

! Subsequent iterations dont recompute factors    
FIRSTTIME = .FALSE.

END SUBROUTINE MSOLVECPLX42

SUBROUTINE N2BLOCKS(I1,I2,NNZ,II,JJ,TIJ)
! Returns the sparse matrix corresponding to the
! NN^2xNN^2 diagonal blocks of the whole matrix,
! i.e., I kron T + T kron I + (3NELEC-2) diag(T) + V(I1:I2) - SHIFT_GL 

! NNZ = NUMBER OF NONZERO ELEMENTS = NN^2*(2*NN-1)
! II = ARRAY OF I INDICES
! JJ = ARRAY OF J INDICES
! TIJ = MATRIX ELEMENT OF (I,J)
USE SNEAK

COMPLEX *16, POINTER :: T(:,:), V(:)
INTEGER, INTENT(IN) :: I1, I2, NNZ
INTEGER, INTENT(INOUT) :: II(1:NNZ), JJ(1:NNZ)
COMPLEX *16, INTENT(INOUT) :: TIJ(1:NNZ)

INTEGER :: I,J,IP,JP,K,NN,NELEC
COMPLEX *16 :: T1, T2

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
T => HAM_GL % T
V => HAM_GL % V(I1:I2)

K = 1

DO JP = 1,NN
DO J = 1,NN
  DO IP = 1,NN
  DO I = 1,NN

    T1 = 0
    T2 = 0
    IF(I == IP) T1 = T(J,JP) ! T kron I
    IF(J == JP) T2 = T(I,IP) ! I kron T

    IF(T1 /= 0 .OR. T2 /= 0) THEN
      II(K) = I + (J-1)*NN
      JJ(K) = IP + (JP-1)*NN
      TIJ(K) = T1 + T2
      IF(II(K).EQ.JJ(K)) THEN
        ! Diagonal values already contains 2 diag(T)
        TIJ(K)=1.5D0*NELEC*TIJ(K)+V(II(K))-SHIFT_GL
      END IF
      K = K + 1
    END IF

  END DO
  END DO
END DO
END DO

END SUBROUTINE N2BLOCKS

#endif

!--------------------------------------------------------------------------------------

#if defined(USE_SUPERLU)

SUBROUTINE MSOLVECPLX5(N,R,Z,NELT,IA,JA,A,ISYM,RWORK,IWORK)
USE SNEAK
USE SUPERLU_MOD
! Block-diagonal preconditioner based on the NN^3xNN^3 diagonal blocks
! of the whole problem. Each block is sparse and the solution for a
! block is computed with SuperLU.
!
! This version performs factorization+solve at each call; it doesnt
! reuse any information between iterations. This is cheap memory-wise
! but performs a lot of flops.
! 

INTEGER, INTENT(IN) :: N, NELT, IA(*), JA(*), ISYM, IWORK(*)
REAL *8, INTENT(IN) :: R(1:N/2,1:2), RWORK(*), A(*)
REAL *8, INTENT(OUT) :: Z(1:N/2,1:2)

INTEGER :: I,IP,J,JP,K,I1,I2,B1,B2,INFO,NN,NELEC,NBLOCKS,NNZ,NP,ME
INTEGER, ALLOCATABLE :: II(:), JJ(:)
COMPLEX *16 :: ZZ(1:N/2)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE :: TIJ(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
COMPLEX *16 :: T1, T2

INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: GRID(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: OPTIONS(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SCALEPERMSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: LUSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SOLVESTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: ASLU(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: STAT(:)

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
NP = MPIDATA_GL % NP
ME = MPIDATA_GL % ME
T => HAM_GL % T
ZZ = DCMPLX(R(:,1),R(:,2))

NBLOCKS=NN**(3*NELEC-3)/NP ! NP should divide N^3
IF(FIRSTTIME) THEN
  ALLOCATE(GRID(NBLOCKS))
  ALLOCATE(OPTIONS(NBLOCKS))
  ALLOCATE(SCALEPERMSTRUCT(NBLOCKS))
  ALLOCATE(LUSTRUCT(NBLOCKS))
  ALLOCATE(SOLVESTRUCT(NBLOCKS))
  ALLOCATE(ASLU(NBLOCKS))
  ALLOCATE(STAT(NBLOCKS))
END IF

! Allocate data for the sparse matrix
! Each block has the same size and number of elements
NNZ=NN**3*(3*NN-2)
IF(FIRSTTIME) THEN
  IF(ME.EQ.0) WRITE(*,*)'3rd-level preconditioner'
  ALLOCATE(II(NNZ*NBLOCKS))
  ALLOCATE(JJ(NNZ*NBLOCKS))
  ALLOCATE(TIJ(NNZ*NBLOCKS))
END IF


DO K = 1,NBLOCKS
  IF(ME.EQ.0.AND.FIRSTTIME.AND.(K+1)*10/NBLOCKS>K*10/NBLOCKS) write(*,'(I4,A)',ADVANCE='NO')(K+1)*100/NBLOCKS,'%'

  I1 = (K-1)*NN**3 + 1
  I2 = K*NN**3
  B1 = (K-1)*NNZ + 1
  B2 = K*NNZ

  ! Compute the values of the sparse matrix
  IF(FIRSTTIME) THEN
    CALL N3BLOCKS(I1,I2,NNZ,II(B1:B2),JJ(B1:B2),TIJ(B1:B2))
  END IF  

  ! Solve with SuperLU_DIST
  CALL SUPERLUSOLVE(II(B1:B2),JJ(B1:B2),TIJ(B1:B2),NNZ,ZZ(I1:I2),NN**3, &
       GRID(K),OPTIONS(K),SCALEPERMSTRUCT(K),LUSTRUCT(K),SOLVESTRUCT(K),&
       ASLU(K),STAT(K),FIRSTTIME)

  ! Put back in solution vector  
  Z(I1:I2,1) = DREAL(ZZ(I1:I2))
  Z(I1:I2,2) = DIMAG(ZZ(I1:I2))

END DO

Z(:,1) = DREAL(ZZ)
Z(:,2) = DIMAG(ZZ)

! Subsequent iterations dont recompute factors    
FIRSTTIME = .FALSE.

END SUBROUTINE MSOLVECPLX5

SUBROUTINE MSOLVECPLX52(XX,ZZ,N)
USE SNEAK
USE SUPERLU_MOD
! Block-diagonal preconditioner based on the NN^3xNN^3 diagonal blocks
! of the whole problem. Each block is sparse and the solution for a
! block is computed with SuperLU.
!
! This version performs factorization+solve at each call; it doesnt
! reuse any information between iterations. This is cheap memory-wise
! but performs a lot of flops.
! 
! Simplified interface for PETSc

INTEGER, INTENT(IN) :: N

INTEGER :: I,IP,J,JP,K,I1,I2,B1,B2,INFO,NN,NELEC,NBLOCKS,NNZ,NP,ME
INTEGER, ALLOCATABLE :: II(:), JJ(:)
COMPLEX *16 :: XX(N),ZZ(N)
COMPLEX *16, POINTER :: T(:,:), V(:)
COMPLEX *16, ALLOCATABLE :: TIJ(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.
COMPLEX *16 :: T1, T2

INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: GRID(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: OPTIONS(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SCALEPERMSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: LUSTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: SOLVESTRUCT(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: ASLU(:)
INTEGER(SUPERLU_PTR),ALLOCATABLE,SAVE :: STAT(:)

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
NP = MPIDATA_GL % NP
ME = MPIDATA_GL % ME
T => HAM_GL % T

NBLOCKS=NN**(3*NELEC-3)/NP ! NP should divide N^3
IF(FIRSTTIME) THEN
  ALLOCATE(GRID(NBLOCKS))
  ALLOCATE(OPTIONS(NBLOCKS))
  ALLOCATE(SCALEPERMSTRUCT(NBLOCKS))
  ALLOCATE(LUSTRUCT(NBLOCKS))
  ALLOCATE(SOLVESTRUCT(NBLOCKS))
  ALLOCATE(ASLU(NBLOCKS))
  ALLOCATE(STAT(NBLOCKS))
END IF

! Allocate data for the sparse matrix
! Each block has the same size and number of elements
NNZ=NN**3*(3*NN-2)
IF(FIRSTTIME) THEN
  IF(ME.EQ.0) WRITE(*,*)'3rd-level preconditioner'
  ALLOCATE(II(NNZ*NBLOCKS))
  ALLOCATE(JJ(NNZ*NBLOCKS))
  ALLOCATE(TIJ(NNZ*NBLOCKS))
END IF

ZZ=XX

DO K = 1,NBLOCKS
  IF(ME.EQ.0.AND.FIRSTTIME.AND.(K+1)*10/NBLOCKS>K*10/NBLOCKS) write(*,'(I4,A)',ADVANCE='NO')(K+1)*100/NBLOCKS,'%'

  I1 = (K-1)*NN**3 + 1
  I2 = K*NN**3
  B1 = (K-1)*NNZ + 1
  B2 = K*NNZ

  ! Compute the values of the sparse matrix
  IF(FIRSTTIME) THEN
    CALL N3BLOCKS(I1,I2,NNZ,II(B1:B2),JJ(B1:B2),TIJ(B1:B2))
  END IF  

  ! Solve with SuperLU_DIST
  CALL SUPERLUSOLVE(II(B1:B2),JJ(B1:B2),TIJ(B1:B2),NNZ,ZZ(I1:I2),NN**3, &
       GRID(K),OPTIONS(K),SCALEPERMSTRUCT(K),LUSTRUCT(K),SOLVESTRUCT(K),&
       ASLU(K),STAT(K),FIRSTTIME)

END DO

! Subsequent iterations dont recompute factors    
FIRSTTIME = .FALSE.

END SUBROUTINE MSOLVECPLX52

SUBROUTINE N3BLOCKS(I1,I2,NNZ,II,JJ,TIJ)
! Returns the sparse matrix corresponding to the
! NN^3xNN^3 diagonal blocks of the whole matrix,
! i.e., I kron T + (3NELEC-1) diag(T) + V(I1:I2) - SHIFT_GL 

! NNZ = NUMBER OF NONZERO ELEMENTS = NN^2*(2*NN-1)
! II = ARRAY OF I INDICES
! JJ = ARRAY OF J INDICES
! TIJ = MATRIX ELEMENT OF (I,J)
USE SNEAK

COMPLEX *16, POINTER :: T(:,:), V(:)
INTEGER, INTENT(IN) :: I1, I2, NNZ
INTEGER, INTENT(INOUT) :: II(1:NNZ), JJ(1:NNZ)
COMPLEX *16, INTENT(INOUT) :: TIJ(1:NNZ)

INTEGER :: I,J,IP,JP,K,KP,Z,NN,NELEC
COMPLEX *16 :: T1, T2, T3

NN = PARAMS_GL % NN
NELEC = PARAMS_GL % NELEC
T => HAM_GL % T
V => HAM_GL % V(I1:I2)

Z = 1

DO KP = 1,NN
DO K = 1,NN
  DO JP = 1,NN
  DO J = 1,NN
    DO IP = 1,NN
    DO I = 1,NN
  
      T1 = 0
      T2 = 0
      T3 = 0
      IF(I == IP .AND. K == KP) T1 = T(J,JP)
      IF(J == JP .AND. K == KP) T2 = T(I,IP)
      IF(J == JP .AND. I == IP) T3 = T(K,KP)
  
      IF(T1 /= 0 .OR. T2 /= 0 .OR. T3/=0) THEN
        II(Z) = I + (J-1)*NN + (K-1)*NN**2
        JJ(Z) = IP + (JP-1)*NN + (KP-1)*NN**2
        TIJ(Z) = T1 + T2 + T3
        IF(II(Z).EQ.JJ(Z)) THEN
          ! Diagonal values already contains 3 diag(T)
          TIJ(Z)=NELEC*TIJ(Z)+V(II(Z))-SHIFT_GL
        END IF
        Z = Z + 1
      END IF
  
    END DO
    END DO
  END DO
  END DO
END DO
END DO

END SUBROUTINE N3BLOCKS

#endif

!--------------------------------------------------------------------------------------

#endif

!--------------------------------------------------------------------------------------

SUBROUTINE KINETIC_EIGS(T,ALLEIGS,W,WT,THETA,NN,LWORK,FIRST,LAST)
USE SNEAK

INTEGER, INTENT(IN) :: NN, LWORK, FIRST, LAST
NUMBER, INTENT(IN) :: T(1:NN,1:NN)
NUMBER, INTENT(OUT) :: ALLEIGS(FIRST:LAST)
REAL *8, INTENT(IN) :: THETA
REAL *8, INTENT(OUT) :: W(1:NN,1:NN), WT(1:NN,1:NN)
REAL *8 :: EIGS(1:NN), WORK(1:LWORK), VR(1:NN,1:NN)
INTEGER :: I,J,K,L,INFO
COMPLEX *16 :: SCALE

#if ISREAL
   VR = T
#else
	VR = DREAL(T)/(COS(2D0*THETA))
#endif

CALL DSYEV('V','U',NN,VR,NN,EIGS,WORK,LWORK,INFO)

W = VR
WT = TRANSPOSE(W)

DO L = FIRST,LAST
   K = (L-1)/NN**2 + 1
   J = (L-(K-1)*NN**2-1)/NN + 1
   I = L-(K-1)*NN**2-(J-1)*NN
   ALLEIGS(L) = EIGS(I) + EIGS(J) + EIGS(K)
END DO

SCALE = DCMPLX(COS(THETA),SIN(THETA))
ALLEIGS = ALLEIGS / SCALE**2

END SUBROUTINE KINETIC_EIGS

!--------------------------------------------------------------------------------------

SUBROUTINE KINETIC_INVERSE(T,TI,NN,LWORK,SHIFT)

INTEGER, INTENT(IN) :: NN, LWORK
NUMBER, INTENT(IN) :: T(1:NN,1:NN), SHIFT
NUMBER, INTENT(OUT) :: TI(1:NN,1:NN)

NUMBER :: TEMP(1:NN,1:NN), WORK(1:LWORK)
INTEGER :: I, INFO, IPIV(1:NN)

TEMP = T
TI = 0

DO I = 1,NN
   TI(I,I) = 1
   TEMP(I,I) = TEMP(I,I) - SHIFT
END DO

#if ISREAL
   CALL DSYSV('U',NN,NN,TEMP,NN,IPIV,TI,NN,WORK,LWORK,INFO)
#else
   CALL ZGESV(NN,NN,TEMP,NN,IPIV,TI,NN,INFO)
#endif

END SUBROUTINE KINETIC_INVERSE

!--------------------------------------------------------------------------------------

SUBROUTINE INVERSEITERATION(HAM,PARAMS,MPIDATA,GMSOLVER,PSI,E0,MM)

TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(HAMILTONIAN), INTENT(IN) :: HAM
TYPE(GMRES_SOLVER), INTENT(INOUT) :: GMSOLVER
INTEGER, INTENT(IN) :: MM
NUMBER, INTENT(INOUT) :: PSI(1:MM)
NUMBER, INTENT(IN) :: E0
NUMBER :: PSI2(1:MM), HPSI(1:MM), E
REAL *8 :: MYERR, MYTOL
INTEGER :: I,J,ERRFLAG

! INITIALIZE VAIRBALES

I = 0
E = E0
MYTOL = 1D-4
MYERR = 1D0

! START INVERSE ITERATION

IF(MPIDATA % ME == 0) WRITE(*,"(A)") "STARTING INVERSE ITERATION ..."

DO WHILE(MYERR > MYTOL)

   I = I + 1
   PSI = PSI/SQRT(MPIDOT(PSI,PSI,MM,MPIDATA))
   PSI2 = PSI
   CALL MYGMRES(PSI2,PSI,MM,HAM,PARAMS,MPIDATA,GMSOLVER)
   E = E0 + 1D0/MPIDOT(PSI,PSI2,MM,MPIDATA)
   PSI = PSI/SQRT(MPIDOT(PSI,PSI,MM,MPIDATA))
   CALL HAMMULT(HAM,PARAMS,MPIDATA,PSI,HPSI,MM)
   HPSI = HPSI - E*PSI

#if ISREAL
   MYERR = SQRT(MPIDOT(HPSI,HPSI,MM,MPIDATA))
#else
   MYERR = SQRT(ABS(MPIDOT(CONJG(HPSI),HPSI,MM,MPIDATA)))
#endif

   IF(MPIDATA % ME == 0) WRITE(*,"(A,I5,A,E12.4E2)") " ITERATION ",I,", ERR = ",MYERR

END DO

!IF(MPIDATA % ME == 0) WRITE(*,"(A,F10.5,F10.5)") "COMPUTED ENERGY: ",DREAL(E),DIMAG(E)
IF(MPIDATA % ME == 0) WRITE(*,"(A,F10.5)") "COMPUTED ENERGY: ",E

END SUBROUTINE INVERSEITERATION

!--------------------------------------------------------------------------------------

END MODULE SOLVER
