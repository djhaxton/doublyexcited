#include "definitions.INC"

MODULE LANMOD
USE STRUCTS
USE MPI
USE HAMMOD
IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------


!! want to calculate 14 orbitals, full shells n=1,2,3.
!!  That way I can get the full spectrum in my initial guess
!!  diagonalization.... got a better way?

SUBROUTINE BLOCKLANCZOS(PARAMS,LANPAR,HAM,MPIDATA,BB)

TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN), TARGET :: MPIDATA
TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR
TYPE(HAMILTONIAN), INTENT(IN), TARGET :: HAM
INTEGER, INTENT(IN) :: BB !! blocksize

INTEGER :: NN,M,MM,KDIM,CHECK,NUMEIGS,ME,FIRST,LAST,MBLK,III,IB,JB,QQ
INTEGER :: LWORK,I,J,K,INFO,IND(LANPAR % NUMEIGS)
REAL *8 :: MAXERR,TOL
real*8, ALLOCATABLE :: GUESSV(:,:) 
NUMBER, ALLOCATABLE :: V(:,:)  !! current block
NUMBER, ALLOCATABLE :: W(:,:)
NUMBER, ALLOCATABLE :: A(:,:,:), B(:,:,:), U(:), VNEW(:,:),VOLD(:,:), A1(:), B1(:), RWORK(:)
NUMBER, ALLOCATABLE :: VMAT(:,:,:), Z(:,:,:), WORK(:), TDM(:,:,:,:), TDM1(:,:,:,:), VL(:,:,:)
NUMBER, POINTER :: T(:,:), VX(:)


T => HAM % T
VX => HAM % V

NN = PARAMS % NN
M = PARAMS % M
MM = MPIDATA % MM
ME = MPIDATA % ME
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST
MBLK = MPIDATA % MBLK

NUMEIGS = LANPAR % NUMEIGS
KDIM = LANPAR % KDIM
CHECK = LANPAR % CHECK
TOL = LANPAR % TOL
!LWORK = 6*KDIM*BB
LWORK = 12*KDIM*BB

ALLOCATE(V(MM,BB),W(MM,BB),U(MM),A(BB,BB,KDIM),B(BB,BB,KDIM),VOLD(MM,BB),A1(KDIM*BB),B1(KDIM*BB),VNEW(MM,BB),GUESSV(MM,BB))

ALLOCATE(Z(BB,KDIM,BB*KDIM),WORK(LWORK))

IF(LANPAR % GSORTH) ALLOCATE(VMAT(MM,BB,KDIM))

#if ISREAL == 1
   DO I = 1,NUMEIGS
      IND(I) = I
   END DO
#else
   ALLOCATE(TDM(BB,KDIM,BB,KDIM),TDM1(BB,KDIM,BB,KDIM),VL(BB,KDIM,KDIM*BB),RWORK(LWORK))
   TDM = 0
   TDM1 = 0
#endif

!V=1d0
!do IB=1,BB
!   V(IB,IB)=V(IB,IB)+IB
!enddo

call RANDOM_NUMBER(GUESSV)
V=GUESSV
call RANDOM_NUMBER(GUESSV)
V=V+GUESSV*(0d0,1d-4)

do QQ=1,2
do IB=1,BB
   CALL GRAMSCHMIDT(V,V(:,IB),IB,MPIDATA)
   V(:,IB) = V(:,IB)/SQRT(MPIDOT(V(:,IB),V(:,IB),MM,MPIDATA))
enddo
enddo


B(:,:,1) = 0

DO J = 1,KDIM-1

   IF(LANPAR % GSORTH) VMAT(:,:,J) = V(:,:)
   do IB=1,BB
      CALL HAMMULT(HAM,PARAMS,MPIDATA,V(:,IB),W(:,IB),MM)
      do JB=1,BB
         A(JB,IB,J) = MPIDOT(V(:,JB),W(:,IB),MM,MPIDATA)
      enddo
   enddo

!   print *, "ALPHABLOCK", A(:,:,J)

   VOLD(:,:) = V(:,:)
   VNEW(:,:) = W(:,:)

   DO QQ=1,2
      DO IB=1,BB
         IF(LANPAR % GSORTH) THEN
            CALL GRAMSCHMIDT(RESHAPE(VMAT(:,:,:),(/MM,KDIM*BB/)),VNEW(:,IB),J*BB+1,MPIDATA)
         ENDIF
         
         CALL GRAMSCHMIDT(VOLD,VNEW(:,IB),BB+1,MPIDATA)
         CALL GRAMSCHMIDT(VNEW,VNEW(:,IB),IB,MPIDATA)
         VNEW(:,IB) = VNEW(:,IB)/SQRT(MPIDOT(VNEW(:,IB),VNEW(:,IB),MM,MPIDATA))
      enddo
   ENDDO


   V(:,:)=VNEW(:,:)

   DO IB=1,BB
      DO JB=1,BB
         B(IB,JB,J+1) = MPIDOT(V(:,IB),W(:,JB),MM,MPIDATA)
      ENDDO
   ENDDO

#if ISREAL == 0
   TDM(:,J,:,J) = A(:,:,J)
   TDM(:,J,:,J+1) = TRANSPOSE(B(:,:,J+1))
   IF(J > 1) TDM(:,J,:,J-1) = B(:,:,J)
#endif
   
   IF(MOD(J,CHECK) == 0) THEN
      
#if ISREAL == 1

!      A1(1:J) = A(1:J)
!      B1(1:J) = B(1:J)
!      CALL DSTEV('V',J,A1(1:J),B1(2:J),Z(1:J,1:J),J,WORK,INFO)
      print *, "NOT ALLOWED"; stop
#else

      TDM1(:,:,:,:) = TDM(:,:,:,:)

      CALL ZGEEV('V','V',J*BB,TDM1,KDIM*BB,A1,VL, &
           KDIM*BB,Z,KDIM*BB,WORK,LWORK,RWORK,INFO)
      if (info.ne.0) then
         print *, "INFO ZGEEV",info;         stop
      endif
      CALL EIGSORT(DREAL(A1(1:J*BB)),IND,J*BB,NUMEIGS)

#endif

      IF(.not.LANPAR % GSORTH) THEN
         print *, "OOGA NEED ORTH"; stop
      endif

      DO K = 1,NUMEIGS
         
         U = 0
!            DO I = 1,J*BB
         DO III=1,J
            DO IB=1,BB
               U = U + Z(IB,III,IND(K))*VMAT(:,IB,III)
            ENDDO
         END DO

!!  herm normalize them for reliability. 

         U = U/SQRT(MPIHERMDOT(U,U,MM,MPIDATA))

         IF(LANPAR % SAVE_VECTORS) LANPAR % EIGENVECTORS(:,K) = U
         
         CALL HAMMULT(HAM,PARAMS,MPIDATA,U,W(:,1),MM)
         U = W(:,1) - A1(IND(K))*U

!! Residual certainly should be herm norm!  Zero c-norm residual
!! is not meaningful.         LANPAR % RESIDUAL(K) = SQRT(MPIDOT(U,U,MM,MPIDATA))
         LANPAR % RESIDUAL(K) = ABS(SQRT(MPIHERMDOT(U,U,MM,MPIDATA)))
      END DO

      MAXERR = MAXVAL(LANPAR % RESIDUAL)

      IF(ME == 0) THEN
         WRITE(*,"(A,I5,A,E12.4E2)") "J = ",J,",   ERROR = ",MAXERR
!         WRITE(*,'(100E10.2)') LANPAR%RESIDUAL(1:NUMEIGS)
!         WRITE(*,'(T10,2F20.10)') (A1(IND(K)),K=1,NUMEIGS)
      END IF
      
      IF(MAXERR < LANPAR % TOL) THEN

         IF(MPIDATA % ME == 0) THEN
            WRITE(*,"(A,I4)") "LANCZOS CONVERGED AT ITERATION ",J
         END IF

         LANPAR % EIGENVALUES = A1(IND)

         EXIT

      END IF
   END IF

   IF(J == KDIM-1) THEN
      IF(MPIDATA % ME == 0) THEN
         WRITE(*,"(A)") "LANCZOS FAILED TO CONVERGE."
      ENDIF
      stop
   ENDIF

END DO

DEALLOCATE(V,W,U,A,B,VOLD,A1,B1,Z,WORK,VNEW,GUESSV)
IF(LANPAR % GSORTH) DEALLOCATE(VMAT)

#if ISREAL == 0
   DEALLOCATE(TDM,TDM1,VL,RWORK)
#endif

END SUBROUTINE BLOCKLANCZOS      

!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------

!! DOES NOT NORMALIZE, GRAM-SCHMIDT ORTHOGONALIZATION ONLY

SUBROUTINE GRAMSCHMIDT(U,V,N,MPIDATA)

TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: N
NUMBER, INTENT(IN) :: U(:,:)
NUMBER, INTENT(INOUT) :: V(:)
INTEGER :: K
NUMBER :: D

DO K = 1,N-1
   D = MPIDOT(U(:,K),V,SIZE(V),MPIDATA)
   V = V - D*U(:,K)
END DO

END SUBROUTINE GRAMSCHMIDT

!--------------------------------------------------------------------------------------

SUBROUTINE INIT_RANDOM_SEED()

! INITIALIZES RANDOM SEED BASED ON SYSTEM CLOCK
!   hmm not good IMHO... want results to be reproducible... and independent
!   of number of processors also -djh

INTEGER :: I_SEED !!, DT_SEED(8)
INTEGER, ALLOCATABLE :: A_SEED(:)
REAL :: R

CALL RANDOM_SEED(SIZE = I_SEED)
ALLOCATE(A_SEED(I_SEED)); A_SEED=0
CALL RANDOM_SEED(GET = A_SEED)
!!CALL DATE_AND_TIME(VALUES = DT_SEED)
A_SEED(I_SEED) = 449               !!DT_SEED(8)
A_SEED(1) = 124                    !!DT_SEED(8)*DT_SEED(7)*DT_SEED(6)
CALL RANDOM_SEED(PUT = A_SEED)
DEALLOCATE(A_SEED)

END SUBROUTINE INIT_RANDOM_SEED

!--------------------------------------------------------------------------------------

SUBROUTINE EIGSORT(X,IND,M,N)

! RETURNS THE INDICES OF THE N SMALLEST VALUES OF X

INTEGER, INTENT(IN) :: M,N
REAL *8, INTENT(IN) :: X(1:M)
INTEGER, INTENT(INOUT) :: IND(1:N)
REAL *8 :: VAL, Y(1:M), MX
INTEGER :: I,J,K,L
LOGICAL :: FOUND

Y = X
MX = MAXVAL(X)

DO K = 1,N
   I = MINLOC(Y,1,.TRUE.)
   IND(K) = I
   Y(I) = MX+1
END DO

END SUBROUTINE EIGSORT

!--------------------------------------------------------------------------------------

END MODULE LANMOD







!!$ #if ISREAL == 0
!!$ 
!!$ SUBROUTINE ARNOLDI(PARAMS,LANPAR,HAM,MPIDATA)
!!$ 
!!$ TYPE(PARAMETERS), INTENT(IN) :: PARAMS
!!$ TYPE(MPI_DATA), INTENT(IN), TARGET :: MPIDATA
!!$ TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR
!!$ TYPE(HAMILTONIAN), INTENT(INOUT), TARGET :: HAM
!!$ 
!!$ INTEGER :: NN,M,MM,KDIM,CHECK,NUMEIGS,ME,FIRST,LAST
!!$ INTEGER :: LWORK,I,J,K,INFO,IND(LANPAR % NUMEIGS)
!!$ REAL *8 :: MAXERR,TOL
!!$ 
!!$ NUMBER, ALLOCATABLE :: V(:), W(:), HESS(:,:), TEMP(:,:), EIGS(:), U(:)
!!$ NUMBER, ALLOCATABLE :: VMAT(:,:), Z(:,:), WORK(:), RWORK(:), VL(:,:)
!!$ 
!!$ 
!!$ NN = PARAMS % NN
!!$ M = PARAMS % M
!!$ MM = MPIDATA % MM
!!$ ME = MPIDATA % ME
!!$ FIRST = MPIDATA % FIRST
!!$ LAST = MPIDATA % LAST
!!$ 
!!$ NUMEIGS = LANPAR % NUMEIGS
!!$ KDIM = LANPAR % KDIM
!!$ CHECK = LANPAR % CHECK
!!$ TOL = LANPAR % TOL
!!$ LWORK = 6*KDIM
!!$ 
!!$ ALLOCATE(V(MM),W(MM),HESS(KDIM,KDIM),TEMP(KDIM,KDIM),EIGS(KDIM),U(MM))
!!$ ALLOCATE(VMAT(MM,KDIM),Z(KDIM,KDIM),WORK(LWORK),RWORK(2*KDIM),VL(KDIM,KDIM))
!!$ 
!!$ V = 1
!!$ V = V/SQRT(MPIDOT(V,V,MM,MPIDATA))
!!$ 
!!$ DO J = 1,KDIM-1
!!$ 
!!$    VMAT(:,J) = V
!!$    CALL HAMMULT(HAM,PARAMS,MPIDATA,V,W,MM)
!!$    
!!$    DO K = 1,J
!!$       HESS(K,J) = MPIDOT(VMAT(:,K),W,MM,MPIDATA)
!!$       W = W - HESS(K,J)*VMAT(:,K)
!!$    END DO
!!$    
!!$    HESS(J+1,J) = SQRT(MPIDOT(W,W,MM,MPIDATA))
!!$    V = W/HESS(J+1,J)
!!$    
!!$    IF(MOD(J,CHECK) == 0) THEN
!!$    
!!$       TEMP = HESS
!!$       CALL ZGEEV('N','V',J,TEMP(1:J,1:J),J,EIGS(1:J),VL(1:J,1:J), &
!!$                   J,Z(1:J,1:J),J,WORK,LWORK,RWORK(1:2*J),INFO)
!!$       CALL EIGSORT(DREAL(EIGS(1:J)),IND,J,NUMEIGS)
!!$ 
!!$       DO K = 1,NUMEIGS
!!$          U = 0
!!$          DO I = 1,J
!!$             U = U + Z(I,IND(K))*VMAT(:,I)
!!$          END DO
!!$ 
!!$          IF(LANPAR % SAVE_VECTORS) LANPAR % EIGENVECTORS(:,K) = U
!!$ 
!!$          CALL HAMMULT(HAM,PARAMS,MPIDATA,U,W,MM)
!!$          U = W - EIGS(IND(K))*U
!!$          LANPAR % RESIDUAL(K) = ABS(SQRT(MPIDOT(U,U,MM,MPIDATA)))
!!$       END DO
!!$ 		
!!$       !LANPAR % RESIDUAL = ABS(Z(J,1:NUMEIGS))
!!$       MAXERR = MAXVAL(LANPAR % RESIDUAL)
!!$ 
!!$       IF(ME == 0) THEN
!!$          WRITE(*,"(A,I5,A,E12.4E2)") "J = ",J,",   ERROR = ",MAXERR
!!$       END IF
!!$       
!!$       IF(MAXERR < LANPAR % TOL) THEN
!!$ 
!!$          IF(MPIDATA % ME == 0) THEN
!!$             WRITE(*,"(A,I4)") "ARNOLDI CONVERGED AT ITERATION ",J
!!$          END IF
!!$ 
!!$          LANPAR % EIGENVALUES = EIGS(IND)
!!$ 
!!$          EXIT
!!$ 
!!$       END IF
!!$ END IF
!!$ 
!!$ IF(J == KDIM-1 .AND. MPIDATA % ME == 0) WRITE(*,"(A)") "ARNOLDI FAILED TO CONVERGE."
!!$ END DO
!!$ 
!!$ DEALLOCATE(V,W,HESS,TEMP,EIGS,VMAT,Z,WORK,RWORK,VL)
!!$ 
!!$ END SUBROUTINE ARNOLDI    
!!$ #endif
!!$ 


!!$ 
!!$ 
!!$ SUBROUTINE LANCZOS(PARAMS,LANPAR,HAM,MPIDATA)
!!$ 
!!$ TYPE(PARAMETERS), INTENT(IN) :: PARAMS
!!$ TYPE(MPI_DATA), INTENT(IN), TARGET :: MPIDATA
!!$ TYPE(LANPARAMS), INTENT(INOUT) :: LANPAR
!!$ TYPE(HAMILTONIAN), INTENT(IN), TARGET :: HAM
!!$ 
!!$ INTEGER :: NN,M,MM,KDIM,CHECK,NUMEIGS,ME,FIRST,LAST,MBLK
!!$ INTEGER :: LWORK,I,J,K,INFO,IND(LANPAR % NUMEIGS)
!!$ REAL *8 :: MAXERR,TOL
!!$ NUMBER, ALLOCATABLE :: V(:), W(:)
!!$ NUMBER, ALLOCATABLE :: A(:), B(:), U(:), VOLD(:), A1(:), B1(:), RWORK(:), VT(:)
!!$ NUMBER, ALLOCATABLE :: VMAT(:,:), Z(:,:), WORK(:), TDM(:,:), TDM1(:,:), VL(:,:)
!!$ NUMBER, POINTER :: T(:,:), VX(:)
!!$ 
!!$ 
!!$ T => HAM % T
!!$ VX => HAM % V
!!$ 
!!$ NN = PARAMS % NN
!!$ M = PARAMS % M
!!$ MM = MPIDATA % MM
!!$ ME = MPIDATA % ME
!!$ FIRST = MPIDATA % FIRST
!!$ LAST = MPIDATA % LAST
!!$ MBLK = MPIDATA % MBLK
!!$ 
!!$ NUMEIGS = LANPAR % NUMEIGS
!!$ KDIM = LANPAR % KDIM
!!$ CHECK = LANPAR % CHECK
!!$ TOL = LANPAR % TOL
!!$ LWORK = 6*KDIM
!!$ 
!!$ ALLOCATE(V(MM),W(MM),U(MM),A(KDIM),B(KDIM),VOLD(MM),A1(KDIM),B1(KDIM),VT(MM))
!!$ ALLOCATE(Z(KDIM,KDIM),WORK(LWORK))
!!$ IF(LANPAR % GSORTH) ALLOCATE(VMAT(MM,KDIM))
!!$ 
!!$ #if ISREAL == 1
!!$    DO I = 1,NUMEIGS
!!$       IND(I) = I
!!$    END DO
!!$ #else
!!$    ALLOCATE(TDM(KDIM,KDIM),TDM1(KDIM,KDIM),VL(KDIM,KDIM),RWORK(LWORK))
!!$    TDM = 0
!!$    TDM1 = 0
!!$ #endif
!!$ 
!!$ V = 1D0
!!$ V = V/SQRT(MPIDOT(V,V,MM,MPIDATA))
!!$ B(1) = 0
!!$ 
!!$ DO J = 1,KDIM-1
!!$ 
!!$    IF(LANPAR % GSORTH) VMAT(:,J) = V
!!$    CALL HAMMULT(HAM,PARAMS,MPIDATA,V,W,MM)
!!$   A(J) = MPIDOT(V,W,MM,MPIDATA)
!!$ 
!!$  print *, "ALPHA", A(J)
!!$ 
!!$    W = W - A(J)*V - B(J)*VOLD
!!$    VOLD = V
!!$    
!!$    B(J+1) = SQRT(MPIDOT(W,W,MM,MPIDATA))
!!$    V = W/B(J+1)
!!$    IF(LANPAR % GSORTH) CALL GRAMSCHMIDT(VMAT,V,J+1,MPIDATA)
!!$    
!!$ #if ISREAL == 0
!!$    TDM(J,J) = A(J)
!!$    TDM(J,J+1) = B(J+1)
!!$    IF(J > 1) TDM(J,J-1) = B(J)
!!$ #endif
!!$ 
!!$ IF(MOD(J,CHECK) == 0) THEN
!!$ 
!!$ #if ISREAL == 1
!!$       A1(1:J) = A(1:J)
!!$       B1(1:J) = B(1:J)
!!$       CALL DSTEV('V',J,A1(1:J),B1(2:J),Z(1:J,1:J),J,WORK,INFO)
!!$ #else
!!$       TDM1 = TDM
!!$ 
!!$       print *, "MAT00"
!!$       do i=1,J
!!$          write(*,'(100F10.5)') TDM1(1:J,i)
!!$       enddo
!!$       print *, "TEMPSTOP00"; stop
!!$ 
!!$       CALL ZGEEV('N','V',J,TDM1(1:J,1:J),J,A1(1:J),VL(1:J,1:J), &
!!$       J,Z(1:J,1:J),J,WORK,LWORK,RWORK(1:2*J),INFO)
!!$       if (info.ne.0) then
!!$          print *, "INFO ZGEEV",info;         stop
!!$       endif
!!$ 
!!$ 
!!$ !print *, "EIGS"
!!$ !write(*,'(2F20.10)') A1(1:J)
!!$ !print *, "TEMPSTOP8800FIRST"
!!$ !stop
!!$ 
!!$ 
!!$       CALL EIGSORT(DREAL(A1(1:J)),IND,J,NUMEIGS)
!!$ #endif
!!$ 
!!$       IF(LANPAR % GSORTH) THEN
!!$ 
!!$          DO K = 1,NUMEIGS
!!$ 
!!$             U = 0
!!$             DO I = 1,J
!!$                U = U + Z(I,IND(K))*VMAT(:,I)
!!$             END DO
!!$ 
!!$             IF(LANPAR % SAVE_VECTORS) LANPAR % EIGENVECTORS(:,K) = U
!!$ 
!!$             CALL HAMMULT(HAM,PARAMS,MPIDATA,U,W,MM)
!!$             U = W - A1(IND(K))*U
!!$             print *, "FSDLKJ"; stop
!!$ ! residual should be from hermitian norm            
!!$             LANPAR % RESIDUAL(K) = SQRT(MPIDOT(U,U,MM,MPIDATA))
!!$          END DO
!!$       ELSE
!!$          LANPAR % RESIDUAL = ABS(Z(J,IND))
!!$       END IF
!!$ 
!!$       MAXERR = MAXVAL(LANPAR % RESIDUAL)
!!$ 
!!$       IF(ME == 0) THEN
!!$          WRITE(*,"(A,I5,A,E12.4E2)") "J = ",J,",   ERROR = ",MAXERR
!!$       END IF
!!$       
!!$       IF(MAXERR < LANPAR % TOL) THEN
!!$ 
!!$          IF(MPIDATA % ME == 0) THEN
!!$             WRITE(*,"(A,I4)") "LANCZOS CONVERGED AT ITERATION ",J
!!$          END IF
!!$ 
!!$          LANPAR % EIGENVALUES = A1(IND)
!!$ 
!!$          EXIT
!!$ 
!!$       END IF
!!$ END IF
!!$ 
!!$ IF(J == KDIM-1 .AND. MPIDATA % ME == 0) WRITE(*,"(A)") "LANCZOS FAILED TO CONVERGE."
!!$ END DO
!!$ 
!!$ DEALLOCATE(V,W,U,A,B,VOLD,A1,B1,Z,WORK)
!!$ IF(LANPAR % GSORTH) DEALLOCATE(VMAT)
!!$ 
!!$ #if ISREAL == 0
!!$    DEALLOCATE(TDM,TDM1,VL,RWORK)
!!$ #endif
!!$ 
!!$ END SUBROUTINE LANCZOS      
!!$ 
!!$ 