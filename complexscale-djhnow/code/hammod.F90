
! This module contains procedures related to the 
! basis and Hamiltonian, including all the procedures needed
! for computing the matrix-vector product Y = H*X.

#include "definitions.INC"

MODULE HAMMOD
USE STRUCTS
USE MPI
USE QUAD
IMPLICIT NONE
CONTAINS

! ----------------------------------------------------------------------------

! HAMILTONIAN MATVEC OPERATION, SUPPORTED FOR NELEC = 1,2 OR 3

SUBROUTINE HAMMULT(HAM,PARAMS,MPIDATA,X,Y,MM)

TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(HAMILTONIAN), INTENT(IN) :: HAM

INTEGER, INTENT(IN) :: MM
NUMBER, INTENT(IN) :: X(1:MM)
NUMBER, INTENT(OUT) :: Y(1:MM)
NUMBER :: Z(1:MM)
INTEGER :: NN, M, MB, FIRST, LAST

INTEGER, SAVE :: ICALLED=0

ICALLED=ICALLED+1

   
NN = PARAMS % NN
M = PARAMS % M
MB = MPIDATA % MBLK
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

IF (MPIDATA%ME==0.AND.MOD(ICALLED,100).eq.0) THEN
   PRINT *, "ICALLED", icalled
ENDIF

! Perform the appropriate kinetic energy matvec based on 
! the number of electrons

IF(PARAMS % NELEC == 1) THEN
  CALL KEMULT1(HAM % T,X,Y,NN,FIRST,LAST,MPIDATA)
ELSE IF(PARAMS % NELEC == 2) THEN
  CALL KEMULT2(HAM % T,X,Y,NN,M,MB,MPIDATA)
ELSE IF(PARAMS % NELEC == 3) THEN
  CALL KEMULT3(HAM % T,X,Y,NN,M,MB,MPIDATA)
END IF

! Special case where all nuclei are at grid points
IF(PARAMS % NUC_ON_GRID) THEN
   Y = Y + (HAM % V) * X
! General case when nuclei are not at grid points
ELSE
   CALL PEMULT(HAM,PARAMS,MPIDATA,X,Z,MM)
   Y = Y + Z
END IF


END SUBROUTINE HAMMULT

! ----------------------------------------------------------------------------

SUBROUTINE KEMULT1(T,X,Y,NN,FIRST,LAST,MPIDATA)

! Performs the parallel kinetic energy matvec Y = T*X for 1 electron

TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: NN,FIRST,LAST
NUMBER, INTENT(IN) :: T(1:NN,1:NN)
NUMBER, INTENT(IN) ::  X(FIRST:LAST)
NUMBER, INTENT(OUT) :: Y(FIRST:LAST)
INTEGER :: I,J,K,L,II,JJ,MPIERR
NUMBER :: YY(1:NN**3), S
INTEGER, ALLOCATABLE, SAVE :: IINDEX(:,:,:), JINDEX(:,:,:), KINDEX(:,:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(IINDEX(NN,NN,2),JINDEX(NN,NN,2),KINDEX(NN,NN,2))
   CALL GETINDICES(IINDEX,JINDEX,KINDEX,NN,FIRST,LAST)
END IF

DO K = 1,NN
DO J = 1,NN
DO I = 1,NN

   JJ = I + (J-1)*NN + (K-1)*NN**2
   S = 0

#ifdef PGFFLAG
   IF (IINDEX(J,K,1).LE.IINDEX(J,K,2)) THEN
#endif
   DO L = IINDEX(J,K,1),IINDEX(J,K,2)
      II = L + (J-1)*NN + (K-1)*NN**2
      S = S + T(L,I)*X(II)
   END DO
#ifdef PGFFLAG
   ENDIF
   IF (JINDEX(I,K,1).LE.JINDEX(I,K,2)) THEN
#endif
   DO L = JINDEX(I,K,1),JINDEX(I,K,2)
      II = I + (L-1)*NN + (K-1)*NN**2
      S = S + T(L,J)*X(II)
   END DO
#ifdef PGFFLAG
   ENDIF
   IF (KINDEX(I,J,1).LE.KINDEX(I,J,2)) THEN
#endif
   DO L = KINDEX(I,J,1),KINDEX(I,J,2) 
      II = I + (J-1)*NN + (L-1)*NN**2
      S = S + T(L,K)*X(II)
   END DO
#ifdef PGFFLAG
   ENDIF
#endif
   YY(JJ) = S

END DO
END DO
END DO


CALL MPI_REDUCE_SCATTER(YY,Y,MPIDATA % BLOCKS, &
MPIDATA % MPISIZE,MPI_SUM,MPIDATA % COMM,MPIERR)

END SUBROUTINE KEMULT1

! ----------------------------------------------------------------------------

SUBROUTINE KEMULT2(T,X,Y,NN,M,MB,MPIDATA)

! Performs the parallel kinetic energy matvec for 2 electrons
! The two-electron vector X is reshaped into a matrix
! The 1-electron version of T is applied to each column of X
! The result is then transposed and T is applied again to each column
! The operation is summarized as Y = T*X + (T*X')' where ' denotes transposition


TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: NN, M, MB
NUMBER, INTENT(IN) :: T(1:NN,1:NN), X(1:M,1:MB)
NUMBER, INTENT(OUT) :: Y(1:M,1:MB)

INTEGER :: I,J
NUMBER, ALLOCATABLE, SAVE :: W(:,:), Z(:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(W(M,MB),Z(M,MB))
END IF


! Y = T*X

!$OMP PARALLEL DO SHARED(T,X,Y,NN,MB) PRIVATE(J)
DO J = 1,MB
   CALL KEMULTS(T,X(:,J),Y(:,J),NN)
END DO
!$OMP END PARALLEL DO

! W = TRANSPOSE(X)
CALL MPITRANSPOSE(X,W,M,MB,MPIDATA)

! Z = T*W

!$OMP PARALLEL DO SHARED(T,W,Z,NN,MB) PRIVATE(J)
DO J = 1,MB
   CALL KEMULTS(T,W(:,J),Z(:,J),NN)
END DO
!$OMP END PARALLEL DO

! W = TRANSPOSE(Z)
CALL MPITRANSPOSE(Z,W,M,MB,MPIDATA)

! Y = T*X + TRANSPOSE(T*TRANSPOSE(X))
Y = Y + W

END SUBROUTINE KEMULT2

!--------------------------------------------------------------------------------------

SUBROUTINE KEMULT3(T,X,Y,NN,M,MB,MPIDATA)

! Performs the kinetic energy matvec for 3 electron
! Analogous to KEMULT2

TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: NN,M,MB
NUMBER, INTENT(IN) :: T(1:NN,1:NN), X(1:M,1:M,1:MB)
NUMBER, INTENT(OUT) :: Y(1:M,1:M,1:MB)
NUMBER :: W(1:M,1:M,1:MB), Z(1:M,1:M,1:MB)
INTEGER :: I,J,K

!$OMP PARALLEL DO SHARED(T,X,Y,NN,M,MB) PRIVATE(J,K)
DO K = 1,MB
   DO J = 1,M
      CALL KEMULTS(T,X(:,J,K),Y(:,J,K),NN)
   END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SHARED(T,X,W,NN,M,MB) PRIVATE(I,K)
DO K = 1,MB
   DO I = 1,M
      CALL KEMULTS(T,X(I,:,K),W(I,:,K),NN)
   END DO
END DO
!$OMP END PARALLEL DO

Y = Y + W

DO I = 1,M
   CALL MPITRANSPOSE(X(I,:,:),W(I,:,:),M,MB,MPIDATA)
   DO J = 1,MB
      CALL KEMULTS(T,W(I,:,J),Z(I,:,J),NN)
   END DO
   CALL MPITRANSPOSE(Z(I,:,:),W(I,:,:),M,MB,MPIDATA)
   Y(I,:,:) = Y(I,:,:) + W(I,:,:)
END DO


END SUBROUTINE KEMULT3

!--------------------------------------------------------------------------------------

SUBROUTINE KEMATRIX(T,NN,DELTA)

! Defines the 1D kinetic energy matrix elements
! T is NN x NN where NN = 2*N+1

INTEGER, INTENT(IN) :: NN
REAL *8, INTENT(IN) :: DELTA
NUMBER, INTENT(OUT) :: T(1:NN,1:NN)
INTEGER :: I,J
REAL *8 :: PI

PI = 4D0*ATAN(1D0)

DO J = 1,NN
DO I = 1,NN

   IF(I == J) THEN
      T(I,J) = (PI/DELTA)**2/6D0
   ELSE
      T(I,J) = (-1D0)**(I-J)/(DELTA**2*(I-J)**2)
   END IF
   
END DO
END DO

END SUBROUTINE KEMATRIX

!--------------------------------------------------------------------------------------

SUBROUTINE KEMULTS(T,X,Y,NN)

! Performs the non-parallel kinetic energy matvec for 1 electron

INTEGER, INTENT(IN) :: NN
NUMBER, INTENT(IN) :: T(1:NN,1:NN), X(1:NN,1:NN,1:NN)
NUMBER,INTENT(OUT) :: Y(1:NN,1:NN,1:NN)

INTEGER :: I,J,K,L

DO K = 1,NN
DO J = 1,NN
DO I = 1,NN

   Y(I,J,K) = SUM(T(:,I)*X(:,J,K)) &
            + SUM(T(:,J)*X(I,:,K)) &
            + SUM(T(:,K)*X(I,J,:))

END DO
END DO
END DO

END SUBROUTINE KEMULTS

!--------------------------------------------------------------------------------------

SUBROUTINE PEMULT(HAM,PARAMS,MPIDATA,X,Y,MM)

! General subroutine for the potentiel matvec, which includes all terms:
! nuclei-nuclei repulsion, electron-nuclei repulsion and electron-electron repulsion
 
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
INTEGER, INTENT(IN) :: MM
NUMBER, INTENT(IN) :: X(1:MM)
NUMBER, INTENT(OUT) :: Y(1:MM)
NUMBER :: Z(1:MM)

INTEGER :: I,J,NN, M, MB, FIRST, LAST
INTEGER :: N,NBIG,NCOR,NUCLEI
REAL *8 :: DELTA,R,S
REAL *8, ALLOCATABLE :: B(:,:,:,:)
INTEGER, POINTER :: CHARGE(:)
REAL *8, POINTER :: RNUC(:,:)

  print *, "fixme real data types. 1113"; stop
!!!  b is real but basis is complex.


N = PARAMS % N
NBIG = PARAMS % NBIG
NCOR = PARAMS % NCOR
DELTA = PARAMS % DELTA
NN = PARAMS % NN
M = PARAMS % M
MB = MPIDATA % MBLK
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

NUCLEI = HAM % MOL % NUCLEI
CHARGE => HAM % MOL % CHARGE
RNUC => HAM % MOL % RNUC

S = 0

DO I = 1,NUCLEI-1
DO J = I+1,NUCLEI

   R = SQRT(ABS(RNUC(1,I)-RNUC(1,J))**2 &
          + ABS(RNUC(2,I)-RNUC(2,J))**2 &
          + ABS(RNUC(3,I)-RNUC(3,J))**2)

   S = S + (1D0*CHARGE(I)*CHARGE(J))/R

END DO
END DO

ALLOCATE(B(-NCOR:NCOR,-N:N,1:3,1:NUCLEI))
CALL BASIS_MATRIX(N,NCOR,NUCLEI,HAM%MOL%RNUC,DELTA,B)

IF(PARAMS % NELEC == 1) THEN
   CALL PEMULT_CORR(HAM,PARAMS,MPIDATA,B,X,Y,N,NCOR)
ELSE IF(PARAMS % NELEC == 2) THEN
   CALL PEMULT_CORR_2(HAM,PARAMS,MPIDATA,B,X,Y,M,MB)
   CALL TWO_BODY_MULT(HAM,PARAMS,MPIDATA,X,Z,M,MB)
   Y = Y + Z
END IF

Y = Y + S*X

DEALLOCATE(B)

END SUBROUTINE PEMULT

!--------------------------------------------------------------------------------------

NUMBER FUNCTION MPIDOT(X,Y,N,MPIDATA)

! Computes the dot product of 2 vectors that are identically 
! distributed over MPI processes

INTEGER, INTENT(IN) :: N
NUMBER, INTENT(IN) :: X(1:N), Y(1:N)
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA

NUMBER :: S, SS
INTEGER :: MPIERR

! compute local dot product
SS = SUM(X*Y)

! sum the local dot products over all processes
CALL MPI_ALLREDUCE(SS,S,1,MPIDATA%MPISIZE,MPI_SUM,MPIDATA % COMM,MPIERR)
MPIDOT = S

END FUNCTION MPIDOT


NUMBER FUNCTION MPIHERMDOT(X,Y,N,MPIDATA)
IMPLICIT NONE !!!
! Computes the dot product of 2 vectors that are identically 
! distributed over MPI processes

INTEGER, INTENT(IN) :: N
NUMBER, INTENT(IN) :: X(1:N), Y(1:N)
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA

NUMBER :: S, SS,XC(1:N)
INTEGER :: MPIERR

! compute local dot product
XC=CONJG(X)
SS = SUM(XC*Y)

! sum the local dot products over all processes
CALL MPI_ALLREDUCE(SS,S,1,MPIDATA%MPISIZE,MPI_SUM,MPIDATA % COMM,MPIERR)
MPIHERMDOT = S

END FUNCTION MPIHERMDOT

!--------------------------------------------------------------------------------------

NUMBER FUNCTION NOMPIDOT(X,Y,N)
IMPLICIT NONE !!!
INTEGER, INTENT(IN) :: N
NUMBER, INTENT(IN) :: X(1:N), Y(1:N)

! compute local dot product
NOMPIDOT = SUM(X*Y)

END FUNCTION NOMPIDOT


NUMBER FUNCTION NOMPIHERMDOT(X,Y,N)
IMPLICIT NONE !!!
INTEGER, INTENT(IN) :: N
NUMBER, INTENT(IN) :: X(1:N), Y(1:N)
NUMBER :: XC(1:N)

! compute local dot product
XC=CONJG(X)
NOMPIHERMDOT = SUM(XC*Y)

END FUNCTION NOMPIHERMDOT

!--------------------------------------------------------------------------------------

SUBROUTINE POTENTIAL1(HAM,PARAMS,MOL,MPIDATA)

! Defines the total potentential energy for a 1-electron system 
! Only called if NUC_ON_GRID == .TRUE.

TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MOLECULE), TARGET, INTENT(IN) :: MOL
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA

INTEGER, POINTER :: CHARGE(:)
REAL *8, POINTER :: RNUC(:,:)
NUMBER, POINTER :: V(:)
INTEGER :: I,J,K,P,NI,IP,JP,KP,N,NN,FIRST,LAST
REAL *8 :: A, S, R, DELTA
LOGICAL :: ONGRID

DELTA = PARAMS % DELTA
N = PARAMS % N
NN = PARAMS % NN
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

V(FIRST:LAST) => HAM % V
RNUC => HAM % MOL % RNUC
CHARGE => HAM % MOL % CHARGE

! Compute nuclei repulsion energy

V = 0
   
DO I = 1,HAM%MOL%NUCLEI-1
DO J = I+1,HAM%MOL%NUCLEI
   
   R = SQRT(ABS(RNUC(1,I)-RNUC(1,J))**2 &
          + ABS(RNUC(2,I)-RNUC(2,J))**2 &
          + ABS(RNUC(3,I)-RNUC(3,J))**2)
             
   V = V + (1D0*CHARGE(I)*CHARGE(J))/R
   
END DO
END DO

! Compute nuclei-electron attraction energy

DO NI = 1,MOL%NUCLEI

   A = -CHARGE(NI)/DELTA**3

   IP = FLOOR(ABS(RNUC(1,NI))/DELTA)
   JP = FLOOR(ABS(RNUC(2,NI))/DELTA)
   KP = FLOOR(ABS(RNUC(3,NI))/DELTA)
   
   DO P = FIRST,LAST

      K = (P-1)/NN**2 + 1
      J = (P-(K-1)*NN**2-1)/NN + 1
      I = P-(K-1)*NN**2-(J-1)*NN
      I = I-N-1;  J = J-N-1;  K = K-N-1
      S = HAM % TINV(I-IP,J-JP,K-KP)
      V(P) = V(P) + A*S

   END DO

END DO

END SUBROUTINE POTENTIAL1

!--------------------------------------------------------------------------------------

SUBROUTINE POTENTIAL2(HAM,PARAMS,MOL,MPIDATA)

! Defines the total potentential energy for a 2-electron system
! Only called if NUC_ON_GRID == .TRUE.

TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MOLECULE), TARGET, INTENT(IN) :: MOL
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA

INTEGER, POINTER :: CHARGE(:)
REAL *8, POINTER :: RNUC(:,:), TINV(:,:,:)
NUMBER, POINTER :: V(:,:)
INTEGER :: I,J,K,NI,IP,JP,KP,N,NN,M,NBIG,FIRST,LAST
INTEGER :: I1,J1,K1,I2,J2,K2,P,Q
REAL *8 :: A, S, R, DELTA
LOGICAL :: ONGRID

DELTA = PARAMS % DELTA
N = PARAMS % N
NN = PARAMS % NN
M = PARAMS % M
NBIG = PARAMS % NBIG
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG) => HAM % TINV
V(1:M,FIRST:LAST) => HAM % V
RNUC => HAM % MOL % RNUC
CHARGE => HAM % MOL % CHARGE

! Compute nuclei-repulsion energy

V = 0
   
DO I = 1,HAM%MOL%NUCLEI-1
DO J = I+1,HAM%MOL%NUCLEI
   
   R = SQRT(ABS(RNUC(1,I)-RNUC(1,J))**2 &
          + ABS(RNUC(2,I)-RNUC(2,J))**2 &
          + ABS(RNUC(3,I)-RNUC(3,J))**2)
             
   V = V + (1D0*CHARGE(I)*CHARGE(J))/R
   
END DO
END DO

! Compute nuclei-electron attraction and
! electron-electron repulsion energies

DO NI = 1,MOL % NUCLEI

   A = -CHARGE(NI)/DELTA**3
   IP = 0; JP = 0; KP = 0;

   DO Q = FIRST,LAST
 
   K1 = (Q-1)/NN**2 + 1
   J1 = (Q-(K1-1)*NN**2-1)/NN + 1
   I1 = Q-(K1-1)*NN**2-(J1-1)*NN
   I1 = I1-N-1;  J1 = J1-N-1;  K1 = K1-N-1
         
      DO P = 1,M
 
         K2 = (P-1)/NN**2 + 1
         J2 = (P-(K2-1)*NN**2-1)/NN + 1
         I2 = P-(K2-1)*NN**2-(J2-1)*NN
         I2 = I2-N-1;  J2 = J2-N-1;  K2 = K2-N-1
         
         V(P,Q) = V(P,Q) + A*TINV(I1-IP,J1-JP,K1-KP) &
                         + A*TINV(I2-IP,J2-JP,K2-KP) &
                         + TINV(I1-I2,J1-J2,K1-K2)/DELTA**3

      END DO

   END DO


END DO

END SUBROUTINE POTENTIAL2

!--------------------------------------------------------------------------------------

SUBROUTINE POTENTIAL3(HAM,PARAMS,MOL,MPIDATA)

! Defines the total potentential energy for a 3-electron system
! Only called if NUC_ON_GRID == .TRUE.

TYPE(HAMILTONIAN), TARGET, INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MOLECULE), TARGET, INTENT(IN) :: MOL
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA

INTEGER, POINTER :: CHARGE(:)
REAL *8, POINTER :: RNUC(:,:), TINV(:,:,:)
NUMBER, POINTER :: V(:,:,:)
INTEGER :: I,J,K,NI,IP,JP,KP,N,NN,M,NBIG,FIRST,LAST
INTEGER :: I1,J1,K1,I2,J2,K2,I3,J3,K3,II,JJ,KK
REAL *8 :: A, S, R, DELTA
LOGICAL :: ONGRID

DELTA = PARAMS % DELTA
N = PARAMS % N
NN = PARAMS % NN
M = PARAMS % M
NBIG = PARAMS % NBIG
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

TINV => HAM % TINV
V(1:M,1:M,FIRST:LAST) => HAM % V
RNUC => HAM % MOL % RNUC
CHARGE => HAM % MOL % CHARGE

! Compute nuclei-repulsion energy

V = 0
   
DO I = 1,HAM%MOL%NUCLEI-1
DO J = I+1,HAM%MOL%NUCLEI
   
   R = SQRT(ABS(RNUC(1,I)-RNUC(1,J))**2 &
          + ABS(RNUC(2,I)-RNUC(2,J))**2 &
          + ABS(RNUC(3,I)-RNUC(3,J))**2)
             
   V = V + (1D0*CHARGE(I)*CHARGE(J))/R
   
END DO
END DO

! Compute nuclei-electron attraction and
! electron-electron repulsion energies

DO NI = 1,MOL % NUCLEI

   A = -CHARGE(NI)/DELTA**3
   IP = 0; JP = 0; KP = 0;

   DO KK = FIRST,LAST

      K3 = (KK-1)/NN**2 + 1
      J3 = (KK-(K3-1)*NN**2-1)/NN + 1
      I3 = KK-(K3-1)*NN**2-(J3-1)*NN

      I3 = I3-N-1;  J3 = J3-N-1;  K3 = K3-N-1

      DO JJ = 1,M

         K2 = (JJ-1)/NN**2 + 1
         J2 = (JJ-(K2-1)*NN**2-1)/NN + 1
         I2 = JJ-(K2-1)*NN**2-(J2-1)*NN

         I2 = I2-N-1;  J2 = J2-N-1;  K2 = K2-N-1

         DO II = 1,M

            K1 = (II-1)/NN**2 + 1
            J1 = (II-(K1-1)*NN**2-1)/NN + 1
            I1 = II-(K1-1)*NN**2-(J1-1)*NN

            I1 = I1-N-1;  J1 = J1-N-1;  K1 = K1-N-1

            V(II,JJ,KK) = V(II,JJ,KK) &
                       +  A*TINV(I1-IP,J1-JP,K1-KP) &
                       +  A*TINV(I2-IP,J2-JP,K2-KP) &
                       +  A*TINV(I3-IP,J3-JP,K3-KP) &
                       + TINV(I1-I2,J1-J2,K1-K2)/DELTA**3 &
                       + TINV(I1-I3,J1-J3,K1-K3)/DELTA**3 &
                       + TINV(I2-I3,J2-J3,K2-K3)/DELTA**3

         END DO
      END DO
   END DO

END DO

END SUBROUTINE POTENTIAL3

!--------------------------------------------------------------------------------------

ELEMENTAL NUMBER FUNCTION BASIS(X,I,DELTA)

! Evaluates the i'th 1D basis function at the point(s) X
! This is an elemental function that can be called if X is an array

INTEGER, INTENT(IN) :: I
REAL *8, INTENT(IN) :: X, DELTA
REAL *8 :: R,PI


PI = 4D0*ATAN(1D0)
R = PI*(X/DELTA - I)

IF(ABS(R) < 1D-14) THEN
   BASIS = 1D0
ELSE
   BASIS = SIN(R)/R
END IF

END FUNCTION BASIS

!--------------------------------------------------------------------------------------

SUBROUTINE GETINDICES(IINDEX,JINDEX,KINDEX,NN,FIRST,LAST)

! Defines all of the boundary indices of a 3D array 
! that is split arbitrarily over MPI processes

INTEGER, INTENT(IN) :: NN,FIRST,LAST
INTEGER, INTENT(OUT) :: IINDEX(:,:,:),JINDEX(:,:,:),KINDEX(:,:,:)
INTEGER :: I,J,K,L,II,P,TEMP(1:NN)

DO J = 1,NN
   DO K = 1,NN
      P = 0
      DO L = 1,NN
         II = L + (J-1)*NN + (K-1)*NN**2
         IF(II >= FIRST .AND. II <=LAST) THEN
            P = P+1
            TEMP(P) = L
         END IF
      END DO
      IINDEX(J,K,1) = MINVAL(TEMP(1:P))
      IINDEX(J,K,2) = MAXVAL(TEMP(1:P))
   END DO
END DO

DO I = 1,NN
   DO K = 1,NN
      P = 0
      DO L = 1,NN
         II = I + (L-1)*NN + (K-1)*NN**2
         IF(II >= FIRST .AND. II <=LAST) THEN
            P = P+1
            TEMP(P) = L
         END IF
      END DO
      JINDEX(I,K,1) = MINVAL(TEMP(1:P))
      JINDEX(I,K,2) = MAXVAL(TEMP(1:P))
   END DO
END DO

DO I = 1,NN
   DO J = 1,NN
      P = 0
      DO L = 1,NN
         II = I + (J-1)*NN + (L-1)*NN**2
         IF(II >= FIRST .AND. II <=LAST) THEN
         P = P+1
         TEMP(P) = L
         END IF
      END DO
      KINDEX(I,J,1) = MINVAL(TEMP(1:P))
      KINDEX(I,J,2) = MAXVAL(TEMP(1:P))
   END DO
END DO

END SUBROUTINE GETINDICES

!--------------------------------------------------------------------------------------

SUBROUTINE EIGMULT(X,Y,W,NN)

! Computes a matrix-vector product of the form
! Y = KRON(W,KRON(W,W))*X
! where W is an NN x NN matrix for a 1D system

! Warning: does not work if X and Y are distributed over processes
! Before calling, X must be gathered into a collective vector
! After calling, Y must be scattered to each process

INTEGER, INTENT(IN) :: NN
REAL *8, INTENT(IN) :: W(1:NN,1:NN)
NUMBER, INTENT(IN) :: X(1:NN,1:NN,1:NN)
NUMBER, INTENT(OUT) :: Y(1:NN,1:NN,1:NN)
NUMBER :: X1(1:NN,1:NN,1:NN), X2(1:NN,1:NN,1:NN)
INTEGER :: I,J,K

DO J = 1,NN
DO I = 1,NN
   X1(I,J,:) = MATMUL(W,X(I,J,:))
END DO
END DO

DO K = 1,NN
DO I = 1,NN
   X2(I,:,K) = MATMUL(W,X1(I,:,K))
END DO
END DO

DO K = 1,NN
DO J = 1,NN
   Y(:,J,K) = MATMUL(W,X2(:,J,K))
END DO
END DO

END SUBROUTINE EIGMULT

!--------------------------------------------------------------------------------------

SUBROUTINE KRONMULT(X,Y,W1,W2,W3,NN,MB,FIRST,LAST,MPIDATA)

TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: NN,MB,FIRST,LAST
REAL *8, INTENT(IN) :: W1(1:NN,1:NN), W2(1:NN,1:NN), W3(1:NN,1:NN)
NUMBER, INTENT(IN) :: X(FIRST:LAST)
NUMBER, INTENT(OUT) :: Y(FIRST:LAST)
NUMBER :: X1(1:NN,1:NN,1:NN), X2(1:NN,1:NN,1:NN), S
INTEGER :: I,J,K,L,II,MPIERR
INTEGER, ALLOCATABLE, SAVE :: IINDEX(:,:,:), JINDEX(:,:,:), KINDEX(:,:,:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

IF(FIRSTTIME) THEN
   FIRSTTIME = .FALSE.
   ALLOCATE(IINDEX(NN,NN,2),JINDEX(NN,NN,2),KINDEX(NN,NN,2))
   CALL GETINDICES(IINDEX,JINDEX,KINDEX,NN,FIRST,LAST)
END IF

DO J = 1,NN
DO I = 1,NN

   DO K = 1,NN
      S = 0 
      DO L = KINDEX(I,J,1),KINDEX(I,J,2)
         II = I + (J-1)*NN + (L-1)*NN**2
         S = S + W3(K,L)*X(II)
      END DO
      X1(I,J,K) = S
   END DO
   
END DO
END DO

CALL MPI_ALLREDUCE(MPI_IN_PLACE,X1,NN**3,MPIDATA%MPISIZE,MPI_SUM,MPIDATA % COMM,MPIERR)

IF(MPIDATA % ME == 0) THEN

DO K = 1,NN
DO I = 1,NN
   X2(I,:,K) = MATMUL(W2,X1(I,:,K))
END DO
END DO

DO K = 1,NN
DO J = 1,NN
   X1(:,J,K) = MATMUL(W1,X2(:,J,K))
END DO
END DO

END IF

CALL MPI_SCATTERV(X1,MPIDATA%BLOCKS,MPIDATA%DISPLS,MPIDATA%MPISIZE, &
                  Y,MB,MPIDATA%MPISIZE,0,MPIDATA%COMM,MPIERR)
                  
END SUBROUTINE KRONMULT

!--------------------------------------------------------------------------------------

SUBROUTINE MPITRANSPOSE(X,Y,M,MB,MPIDATA)
IMPLICIT NONE !!!
! Transposes a matrix X split over MPI processes by column chunks 
! Y = TRANSPOSE(X) is split the same way as X
! Uses the MPI_ALLTOALL subroutine

INTEGER, INTENT(IN) :: M,MB
NUMBER, INTENT(IN) :: X(1:M,1:MB)
TYPE(MPI_DATA), TARGET, INTENT(IN) :: MPIDATA
NUMBER, INTENT(OUT) :: Y(1:M,1:MB)

NUMBER, ALLOCATABLE, SAVE :: XXT(:,:,:), YY(:,:,:), XT(:,:)
INTEGER :: I, MPIERR, NP, MAXBLOCK
INTEGER, POINTER :: BLOCKS(:)
INTEGER, ALLOCATABLE, SAVE :: B(:)
LOGICAL, SAVE :: FIRSTTIME = .TRUE.

NP = MPIDATA % NP
MAXBLOCK = MPIDATA % MAXBLOCK
BLOCKS => MPIDATA % BLOCKS

IF(FIRSTTIME) THEN

   FIRSTTIME = .FALSE.
   ALLOCATE(B(NP),XT(MAXBLOCK,M))
   ALLOCATE(XXT(MAXBLOCK,MAXBLOCK,NP))
   ALLOCATE(YY(MAXBLOCK,MAXBLOCK,NP))

   B(1) = 1
   IF (NP.GT.1)   B(2) = MAXBLOCK + 1

   DO I = 3,NP
      B(I) = B(I-1) + BLOCKS(I)
   END DO

END IF

IF(NP == 1) THEN
   Y = TRANSPOSE(X)
ELSE

   XT(1:MB,1:M) = TRANSPOSE(X)

   DO I = 1,NP
      XXT(:,1:BLOCKS(I),I) = XT(:,B(I):B(I)+BLOCKS(I)-1)
   ENDDO

   CALL MPI_ALLTOALL(XXT,MAXBLOCK**2,MPIDATA % MPISIZE,YY, &
   MAXBLOCK**2,MPIDATA % MPISIZE,MPIDATA % COMM,MPIERR)

   DO I = 1,NP
   Y(B(I):B(I)+BLOCKS(I)-1,:) = YY(1:BLOCKS(I),1:MB,I)
   ENDDO

END IF

END SUBROUTINE MPITRANSPOSE

!--------------------------------------------------------------------------------------

SUBROUTINE DERIVATIVE(F,NN,DELTA)

! Computes the first derivative matrix elements in the sinc basis

INTEGER, INTENT(IN) :: NN
REAL *8, INTENT(IN) :: DELTA
REAL *8, INTENT(OUT) :: F(1:NN,1:NN)
INTEGER :: I,J

F = 0

DO J = 1,NN
DO I = 1,NN

   IF(I /= J) F(I,J) = (-1D0)**(I-J)/(DELTA*(I-J))

END DO
END DO

END SUBROUTINE DERIVATIVE

!--------------------------------------------------------------------------------------

SUBROUTINE GAUSSIAN(V,Z,NN)

! Computes a model gaussian potential that has a known resonance
! Use for testing complex scaling methods
! Z is an array of complex points where the potential is evaluates
! V is the potential at those points
 
INTEGER, INTENT(IN) :: NN
COMPLEX *16, INTENT(IN) :: Z(1:NN)
COMPLEX *16, INTENT(OUT) :: V(1:NN,1:NN,1:NN)
INTEGER :: I,J,K,L
COMPLEX *16 :: V1,V2,XX,YY,ZZ
REAL *8 :: LAMBDA, GAMMA, A(3)

LAMBDA = 0.35
GAMMA = 0.13
A = 0
A(3) = 2

DO K = 1,NN
   DO J = 1,NN
      DO I = 1,NN

         XX = Z(I); YY = Z(J); ZZ = Z(K)
   
         V1 = EXP(-GAMMA*( (XX-A(1))**2 + (YY-A(2))**2 + (ZZ-A(3))**2))
         V2 = EXP(-((XX+A(1))**2 + (YY+A(2))**2 + (ZZ+A(3))**2))
         V(I,J,K) = LAMBDA*(XX**2 + YY**2 + ZZ**2)*(V1 + V2)

      END DO
   END DO
END DO


END SUBROUTINE GAUSSIAN

!--------------------------------------------------------------------------------------

SUBROUTINE HARMONIC(V,Z,NN)

! Computes the 3D harmonic potential (used as a test case for debugging)

INTEGER, INTENT(IN) :: NN
COMPLEX *16, INTENT(IN) :: Z(1:NN)
COMPLEX *16, INTENT(OUT) :: V(1:NN,1:NN,1:NN)
INTEGER :: I,J,K

DO K = 1,NN
   DO J = 1,NN
      DO I = 1,NN

         V(I,J,K) = 0.5*(Z(I)**2 + Z(J)**2 + Z(K)**2)

      END DO
   END DO
END DO

END SUBROUTINE HARMONIC

!--------------------------------------------------------------------------------------

SUBROUTINE EIGENVALUES(H,E,ND,NUMEIGS)
IMPLICIT NONE

! A wrapper to LAPACK to compute eigenvalues of an ND X ND matrix H
! E is the output vector of NUMEIGS eigenvalues
! Works for either real symmetric or general complex matrices

INTEGER, INTENT(IN) :: ND, NUMEIGS
NUMBER, INTENT(IN) :: H(1:ND,1:ND)
NUMBER, INTENT(OUT) :: E(1:NUMEIGS)
INTEGER :: LWORK,INFO
REAL *8, ALLOCATABLE :: RWORK(:)
NUMBER :: A(1:ND,1:ND)
NUMBER, ALLOCATABLE :: WORK(:), VL(:,:), VR(:,:), EE(:)

LWORK = 6*ND
ALLOCATE(WORK(LWORK),RWORK(2*ND),EE(ND))
A = H

! real symmetric solver
#if ISREAL == 1
   CALL DSYEV('N','U',ND,A,ND,EE,WORK,LWORK,INFO)
! general comlex solver
#else
   CALL ZGEEV('N','N',ND,A,ND,EE,VL,ND,VR,ND,WORK,LWORK,RWORK,INFO)
#endif

E = EE(1:NUMEIGS)
DEALLOCATE(WORK,RWORK,EE)

END SUBROUTINE EIGENVALUES

!--------------------------------------------------------------------------------------

! SIGN specifies whether a singlet (+1) or a triplet (-1) wave function is desired

SUBROUTINE HELIUM_ORBITAL(XLIST,IX,JX,TOTX,PSI,M,FIRST,LAST,SIGN,MPIDATA,ME)
  IMPLICIT NONE !!!
  TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
  INTEGER, INTENT(IN) :: M, FIRST, LAST, SIGN,ME,IX,JX,TOTX
  NUMBER, INTENT(OUT) :: PSI(1:M,FIRST:LAST)
  INTEGER :: I,J,MPIERR
  NUMBER,INTENT(IN) :: XLIST(1:M,TOTX)
  NUMBER :: X1(1:M),X2(1:M)

  IF (IX.LT.1.OR.JX.LT.1.OR.IX.GT.TOTX.OR.JX.GT.TOTX) then
     PRINT *, "BOOGA ERROR 45.";     STOP
  ENDIF

  X1(:)=XLIST(:,IX)
  X2(:)=XLIST(:,JX)
  
  psi(:,:)=0d0
  
  DO J = FIRST,LAST
     DO I = 1,M
        PSI(I,J) = X1(I)*X2(J) + SIGN*X2(I)*X1(J)
     END DO
  END DO

  IF (IX.EQ.JX) THEN
     PSI(:,:)=PSI(:,:)/2
  ELSE
     PSI(:,:)=PSI(:,:)/SQRT(2D0)
  ENDIF

!  NORM_PSI = SQRT(MPIDOT(PSI,PSI,SIZE(PSI),MPIDATA))
!  PSI = PSI / NORM_PSI

END SUBROUTINE HELIUM_ORBITAL



SUBROUTINE HELIUM_ORBITALGATHER(HOWMANY,V1,X1,M,FIRST,LAST,MPIDATA,ME,hermopt)
  IMPLICIT NONE !!!  in case it is moved somewhere else
  logical, intent(IN) :: hermopt
  TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
  INTEGER, INTENT(IN) :: M, FIRST, LAST, ME,HOWMANY
  NUMBER, INTENT(IN) :: V1(FIRST:LAST,HOWMANY) 
  NUMBER, INTENT(OUT) :: X1(1:M,HOWMANY)
  INTEGER ::MPIERR,ii,jj

  do ii=1,howmany
     CALL MPI_GATHERV(V1(:,ii),MPIDATA % MM, MPIDATA % MPISIZE, X1(:,ii), MPIDATA % BLOCKS, &
          MPIDATA % DISPLS, MPIDATA % MPISIZE, 0, MPIDATA % COMM, MPIERR)
     
     CALL MPI_BCAST(X1(:,ii),M,MPIDATA % MPISIZE, 0, MPIDATA % COMM, MPIERR)
  enddo

!! could go ahead and gram-schmidt orthogonalize now, in a manner depending on hermopt,
!!  to avoid doing generalized eigenproblem 
!!     ..ok do that

  if (HERMOPT) then
     do ii=1,howmany
        do jj=1,ii-1
           x1(:,ii)=x1(:,ii)-x1(:,jj)*NOMPIHERMDOT(x1(:,jj),x1(:,ii),M)
        enddo
        x1(:,ii)=x1(:,ii)/SQRT(NOMPIHERMDOT(x1(:,ii),x1(:,ii),M))
     enddo
  else
     do ii=1,howmany
        do jj=1,ii-1
           x1(:,ii)=x1(:,ii)-x1(:,jj)*NOMPIDOT(x1(:,jj),x1(:,ii),M)
        enddo
        x1(:,ii)=x1(:,ii)/SQRT(NOMPIDOT(x1(:,ii),x1(:,ii),M))
     enddo
  endif

END SUBROUTINE HELIUM_ORBITALGATHER



! --------------------------------------------------------------------------------------

! FYI - this does NOT return a HERMITIAN-orthonormal subspace.  
!  ( If HERMOPT=0, returns C-orthogonal, hermitian normed subspace. )
!  ( If HERMOPT.ne.0, returns hermitian normed subspace, no orthogonality for theta/=0. )
!   -djh

SUBROUTINE DIAGSUBSPACE(HAM,PARAMS,MPIDATA,M,MB,NUMOUT,OUTE0,OUTPSI,me,etarget,hermopt,numorbs,singlet,EV,FIRST,LAST,np)
  IMPLICIT NONE !!!
! Diagonalizes the totbas dimensional subspace spanned by the 
! 2 electron vectors contains in PSI... returns
  LOGICAL,intent(IN) :: SINGLET
  INTEGER, INTENT(IN) :: M,MB,me,numout,numorbs,first,last,np
  TYPE(HAMILTONIAN), INTENT(IN) :: HAM
  TYPE(PARAMETERS), INTENT(IN) :: PARAMS
  TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
  TYPE(MPI_DATA) :: OLDMPIDATA
  NUMBER, INTENT(IN) :: etarget,EV(1:M,numorbs)
  NUMBER :: TEMPPSI(1:M,1:MB),  HPSI(1:M,1:MB), &
       TEMPPSI2(1:M,1:MB)
  NUMBER, INTENT(OUT) :: OUTE0(1:NUMOUT),OUTPSI(1:M,1:MB,NUMOUT)
  NUMBER, ALLOCATABLE :: WORK(:), A(:,:),OVL(:,:),&
       VL(:,:), VR(:,:),ALPHA(:),&
       BETA(:),E0(:)
  integer, allocatable :: itable(:,:)
  REAL *8, allocatable :: RWORK(:)
  INTEGER :: I,J,INFO,LWORK,totbas,kk,iorb,jorb,isinglet,itriplet
  LOGICAL :: HERMOPT

  if (ME==0) THEN
     PRINT *, "GO DIAGSUBSPACE"
  ENDIF

  HPSI(:,:)=0d0

  if (SINGLET) then
     totbas=NUMORBS*(NUMORBS+1)/2
  else
     totbas=NUMORBS*(NUMORBS-1)/2
  endif
  IF (NUMOUT.gt.totbas) then
     if (ME==0) PRINT *, "Reset NUMOUT to be less than tot orbital basis size",numout,totbas, " due to numorbs= ",numorbs
     stop
  endif

  allocate(A(1:TOTBAS,1:TOTBAS),OVL(1:TOTBAS,1:TOTBAS),&
       VL(1:TOTBAS,1:TOTBAS), VR(1:TOTBAS,1:TOTBAS),ALPHA(1:TOTBAS),&
       BETA(1:TOTBAS),itable(2,1:totbas),e0(1:TOTBAS),RWORK(1:4*TOTBAS))
  LWORK = 6*TOTBAS
  ALLOCATE(WORK(LWORK))

  itriplet=1; isinglet=-1
  if (singlet) then
     itriplet=0
     isinglet=1
  endif
  
  kk=0
  
  do iorb=1,numorbs
     do jorb=iorb+itriplet,numorbs
        kk=kk+1
        itable(:,kk)=(/iorb,jorb/)
     enddo
  enddo
  if (kk.ne.totbas) then
     print *, "OOGA"; stop
  endif
  
  DO J = 1,totbas

     if (ME==0.and.(mod(j,10).eq.0)) print *, "   ...in diagsubspace... ",J," OF ",totbas

!! (X1,X2,PSI,M,FIRST,LAST,SIGN,ME)

     CALL HELIUM_ORBITAL(EV(:,:),itable(1,j),itable(2,j),numorbs,TEMPPSI(:,:),M,FIRST,LAST,isinglet,MPIDATA,me)

     CALL HAMMULT(HAM,PARAMS,MPIDATA,TEMPPSI,HPSI,M*MB)

     DO I = 1,totbas

        CALL HELIUM_ORBITAL(EV(:,:),itable(1,i),itable(2,i),numorbs,TEMPPSI2(:,:),M,FIRST,LAST,isinglet,MPIDATA,me)

        if (.not.HERMOPT) then
           A(I,J) = MPIDOT(TEMPPSI2(:,:),HPSI,M*MB,MPIDATA)
           OVL(I,J) = MPIDOT(TEMPPSI2(:,:),TEMPPSI(:,:),M*MB,MPIDATA)
        else
           A(I,J) = MPIHERMDOT(TEMPPSI2(:,:),HPSI,M*MB,MPIDATA)
           OVL(I,J) = MPIHERMDOT(TEMPPSI2(:,:),TEMPPSI(:,:),M*MB,MPIDATA)
        endif
     END DO
  END DO
  
! always do generalized eigenproblem.  If degeneracy in lan solve for
! one electron orbitals, basis can be non c orthonormal.  Could GS orthog.
! Instead doing generalized eigenproblem regardless of HERMOPT.
!
!    .. no now did northonorm with helium_orbital regardless of HERMOPT. 
!
!  CALL ZGGEV('V','V',TOTBAS,A,TOTBAS,OVL,TOTBAS,ALPHA,BETA,VL,TOTBAS,VR,TOTBAS,WORK,LWORK,RWORK,INFO)
!  E0=ALPHA/BETA
!  if (HERMOPT) then
!     VR(:,:)=VL(:,:) !!... for some reason it's the left eigvects we want
!     !! ... it's like VL and VR are backwards.  
!  endif
  
   CALL ZGEEV('V','V',TOTBAS,A,TOTBAS,E0,VL,TOTBAS,VR,TOTBAS,WORK,LWORK,RWORK,INFO)

  if (info.ne.0) then
     if (me==0) then
        print *, "INFO zggev 77 DIAGSUBSPACE",info; stop
     endif
  endif
  
  DEALLOCATE(WORK)
  
  call mysort(VR,E0,totbas,totbas,etarget)
  OUTE0(:)=E0(1:NUMOUT)
  OUTPSI(:,:,:)=0d0
  DO i = 1,TOTBAS

     CALL HELIUM_ORBITAL(EV(:,:),itable(1,i),itable(2,i),numorbs,TEMPPSI(:,:),M,FIRST,LAST,isinglet,MPIDATA,me)
     
     CALL ZGEMM('N','N',M*MB, NUMOUT,1,(1d0,0d0),TEMPPSI,M*MB,VR(I,1),TOTBAS,(1d0,0d0),OUTPSI,M*MB)
  enddo

!  always herm normed methinks?  Was not before.  Could have been feeding huge vectors into GPLMR...?
!  Or should a full GS orthog be performed? 
!     i.e. mix eigenvectors to make a hermitian orthonormal space?  Then set sigma based on average
!     of expectation values of those vectors?  Maybe should do that... for the moment,
!     just norming; the hamiltonian remains diagonal within the space.

  do i=1,numout
     OUTPSI(:,:,I)=OUTPSI(:,:,I)/SQRT(MPIHERMDOT(OUTPSI(:,:,I),OUTPSI(:,:,I),M*MB,MPIDATA))
  END DO

  deallocate(A,ovl,vl,vr,alpha,beta,itable,e0,rwork)

END SUBROUTINE DIAGSUBSPACE




SUBROUTINE GET_EXPECTS(HAM,PARAMS,MPIDATA,M,MB,NUMVECTS,PSI,EXPECTS,cshiftopt)
IMPLICIT NONE !!!

! Diagonalizes the NUMVECTS dimensional subspace spanned by the 
! 2 electron vectors contains in PSI
LOGICAL,INTENT(IN) :: cshiftopt
TYPE(HAMILTONIAN), INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: M,MB,NUMVECTS
NUMBER :: HPSI(1:M,1:MB)
NUMBER, INTENT(IN) :: PSI(1:M,1:MB,1:NUMVECTS)
NUMBER, INTENT(OUT) :: EXPECTS(1:NUMVECTS)
INTEGER :: I

DO I = 1,NUMVECTS
   CALL HAMMULT(HAM,PARAMS,MPIDATA,PSI(:,:,I),HPSI,M*MB)

   if (cshiftopt) then
      EXPECTS(I)=MPIDOT(PSI(:,:,I),HPSI,M*MB,MPIDATA) &
           /MPIDOT(PSI(:,:,I),PSI(:,:,I),M*MB,MPIDATA)
   else

!! default.  right?!?!?!?!?
!!$  so same as eig for HERMOPT, but different from eig for .not.HERMOPT.
      EXPECTS(I)=MPIHERMDOT(PSI(:,:,I),HPSI,M*MB,MPIDATA) &
           /MPIHERMDOT(PSI(:,:,I),PSI(:,:,I),M*MB,MPIDATA)

   endif

!!$ to check: same as eig
!!$if (HERMOPT) THEN
!!$   EXPECTS(I)=MPIHERMDOT(PSI(:,:,I),HPSI,M*MB,MPIDATA) &
!!$ELSE
!!$   EXPECTS(I)=MPIDOT(PSI(:,:,I),HPSI,M*MB,MPIDATA)
!!$endif
!!        

END DO

END SUBROUTINE GET_EXPECTS

! --------------------------------------------------------------------------------------

SUBROUTINE TWO_POTENTIAL(TINV,V2,DELTA,NBIG,N,NN,M1,FIRST,LAST)

! Computes V2, the electron-electron repulsion potential using TINV

INTEGER, INTENT(IN) :: N,NBIG,NN,M1,FIRST,LAST
REAL *8, INTENT(IN) :: DELTA,TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG)
NUMBER, INTENT(OUT) :: V2(1:M1,FIRST:LAST)
INTEGER :: I1,J1,K1,I2,J2,K2,II,JJ

DO JJ = FIRST,LAST

   K1 = (JJ-1)/NN**2 + 1
   J1 = (JJ-(K1-1)*NN**2-1)/NN + 1
   I1 = JJ-(K1-1)*NN**2-(J1-1)*NN
   I1 = I1-N-1;  J1 = J1-N-1;  K1 = K1-N-1
      
   DO II = 1,M1
   
      K2 = (II-1)/NN**2 + 1
      J2 = (II-(K2-1)*NN**2-1)/NN + 1
      I2 = II-(K2-1)*NN**2-(J2-1)*NN
      I2 = I2-N-1;  J2 = J2-N-1;  K2 = K2-N-1

      V2(II,JJ) = TINV(I2-I1,J2-J1,K2-K1)/DELTA**3
      
   END DO
END DO


END SUBROUTINE TWO_POTENTIAL

!--------------------------------------------------------------------------------------

SUBROUTINE ANGULARMULT(X,Y,MM,N,NN,DELTA,OPT,MPIDATA)

! Applies the angular momentum operators to a one-electron wave function
! The value of OPT specifies which component of the angular momentum to apply

INTEGER, INTENT(IN) :: MM,N,NN,OPT
NUMBER, INTENT(IN) :: X(1:MM)
NUMBER, INTENT(OUT) :: Y(1:MM)
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
REAL *8, INTENT(IN) :: DELTA

REAL *8 :: EYE(-N:N,-N:N), D(-N:N,-N:N), DIPOLE(-N:N,-N:N)
NUMBER :: Z1(1:MM), Z2(1:MM)
INTEGER :: K,MB,FIRST,LAST

MB = MPIDATA % MBLK
FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST

CALL DERIVATIVE(D,NN,DELTA)

EYE = 0
DIPOLE = 0
DO K = -N,N
   EYE(K,K) = 1D0
   DIPOLE(K,K) = K*DELTA
END DO


IF(OPT == 1) THEN
   CALL KRONMULT(X,Z1,EYE,DIPOLE,D,NN,MB,FIRST,LAST,MPIDATA)
   CALL KRONMULT(X,Z2,EYE,D,DIPOLE,NN,MB,FIRST,LAST,MPIDATA)
ELSEIF(OPT == 2) THEN
   CALL KRONMULT(X,Z1,DIPOLE,EYE,D,NN,MB,FIRST,LAST,MPIDATA)
   CALL KRONMULT(X,Z2,D,EYE,DIPOLE,NN,MB,FIRST,LAST,MPIDATA)
ELSEIF(OPT == 3) THEN
   CALL KRONMULT(X,Z1,DIPOLE,D,EYE,NN,MB,FIRST,LAST,MPIDATA)
   CALL KRONMULT(X,Z2,D,DIPOLE,EYE,NN,MB,FIRST,LAST,MPIDATA)
END IF

Y = Z1 - Z2

END SUBROUTINE ANGULARMULT

! --------------------------------------------------------------------------------------

NUMBER FUNCTION ANGULAREXP(PSI,PARAMS,MPIDATA,MM,OPT)

! Computes the expectation of an angular momentum operator

INTEGER, INTENT(IN) :: MM,OPT
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
NUMBER, INTENT(IN) :: PSI(1:MM)
NUMBER :: LPSI(1:MM), L2PSI(1:MM), LLPSI(1:MM)

LLPSI = 0

CALL ANGULARMULT(PSI,LPSI,MM,PARAMS % N, PARAMS % NN,PARAMS % DELTA,1,MPIDATA)
CALL ANGULARMULT(LPSI,L2PSI,MM,PARAMS % N, PARAMS % NN,PARAMS % DELTA,1,MPIDATA)
LLPSI = LLPSI + L2PSI

CALL ANGULARMULT(PSI,LPSI,MM,PARAMS % N, PARAMS % NN,PARAMS % DELTA,2,MPIDATA)
CALL ANGULARMULT(LPSI,L2PSI,MM,PARAMS % N, PARAMS % NN,PARAMS % DELTA,2,MPIDATA)
LLPSI = LLPSI + L2PSI

CALL ANGULARMULT(PSI,LPSI,MM,PARAMS % N, PARAMS % NN,PARAMS % DELTA,3,MPIDATA)
CALL ANGULARMULT(LPSI,L2PSI,MM,PARAMS % N, PARAMS % NN,PARAMS % DELTA,3,MPIDATA)
LLPSI = LLPSI + L2PSI


ANGULAREXP = MPIDOT(PSI,LLPSI,MM,MPIDATA)

END FUNCTION ANGULAREXP

! --------------------------------------------------------------------------------------

SUBROUTINE EXPECTATIONS(HAM,PARAMS,MPIDATA,LANPAR,EXT,EXV,NUMEIGS)

! Computes expectations of T and  V
! Used to compute Virial rations <T> / <V>
! EXT(J) = < PSI(J) | T | PSI(J) >
! EXT(J) = < PSI(J) | V | PSI(J) >

TYPE(HAMILTONIAN), INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
TYPE(LANPARAMS), TARGET, INTENT(IN) :: LANPAR
INTEGER, INTENT(IN) :: NUMEIGS
NUMBER, INTENT(OUT) :: EXT(1:NUMEIGS), EXV(1:NUMEIGS)

INTEGER :: K,MM
NUMBER, ALLOCATABLE :: Y(:)
NUMBER, POINTER :: EV(:,:)

EV => LANPAR % EIGENVECTORS
MM = MPIDATA % MM
ALLOCATE(Y(MM))

DO K = 1,NUMEIGS

   CALL KEMULT1(HAM % T,EV(:,K),Y,PARAMS % NN,MPIDATA % FIRST,MPIDATA % LAST,MPIDATA)
   EXT(K) = MPIDOT(EV(:,K),Y,MM,MPIDATA)
   EXV(K) = MPIDOT(EV(:,K),HAM % V * EV(:,K),MM,MPIDATA)

END DO

DEALLOCATE(Y)

END SUBROUTINE EXPECTATIONS

! --------------------------------------------------------------------------------------

SUBROUTINE BASIS_MATRIX(N,NCOR,NUCLEI,RNUC,DELTA,B)

! Computes an array of coefficients in terms of the basis functions
! used for the nuclear potential matvec when nuclei are not on grid points

INTEGER, INTENT(IN) :: N,NCOR,NUCLEI
REAL *8, INTENT(IN) :: DELTA, RNUC(1:3,1:NUCLEI)
REAL *8, INTENT(OUT) :: B(-NCOR:NCOR,-N:N,1:3,1:NUCLEI)

INTEGER :: I,J,K,L
REAL *8 :: R

!! no probably is fine print *, "fixme real data types. 1112"; stop

! B = NCOR x N x NUCELI x 3

DO L = 1,NUCLEI
   DO K = 1,3
      DO J = -N,N
         DO I = -NCOR,NCOR
            R = I*DELTA - RNUC(K,L)
            B(I,J,K,L) = BASIS(R,J,DELTA)
         END DO
      END DO
   END DO
END DO

END SUBROUTINE BASIS_MATRIX

!--------------------------------------------------------------------------------------

SUBROUTINE PEMULT_CORR(HAM,PARAMS,MPIDATA,B,X,Y,N,NCOR)

! Computes the nuclear potential matvec Y = V*X for a one-electron
! system when nuclei are not on grid points

INTEGER, INTENT(IN) :: N,NCOR
TYPE(HAMILTONIAN), TARGET,  INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
REAL *8, INTENT(IN) :: B(:,:,:,:)
NUMBER, INTENT(IN) :: X(-N:N,-N:N,-N:N)
NUMBER, INTENT(OUT) :: Y(-N:N,-N:N,-N:N)

INTEGER :: I,J,K,NI,NUCLEI
INTEGER, POINTER :: CHARGE(:)
REAL *8, POINTER :: TINV(:,:,:)

REAL *8 :: DELTA,BT(-N:N,-NCOR:NCOR)
NUMBER :: X1(-NCOR:NCOR,-N:N,-N:N)
NUMBER :: X2(-NCOR:NCOR,-NCOR:NCOR,-N:N)
NUMBER :: X3(-NCOR:NCOR,-NCOR:NCOR,-NCOR:NCOR)
NUMBER :: Y2(-N:N,-N:N,-NCOR:NCOR)
NUMBER :: Y1(-N:N,-NCOR:NCOR,-NCOR:NCOR)

!!! probably fine print *, "fixme real data types. 1111"; stop

NUCLEI = HAM % MOL % NUCLEI
DELTA = PARAMS % DELTA
CHARGE => HAM % MOL % CHARGE
TINV => HAM % TINV(-NCOR:NCOR,-NCOR:NCOR,-NCOR:NCOR)

Y = 0

DO NI = 1,NUCLEI

   DO K = -N,N
   DO J = -N,N
      X1(:,J,K) = MATMUL(B(:,:,1,NI),X(:,J,K))
   END DO
   END DO

   DO K = -N,N
   DO I = -NCOR,NCOR
      X2(I,:,K) = MATMUL(B(:,:,2,NI),X1(I,:,K))
   END DO
   END DO

   DO J = -NCOR,NCOR
   DO I = -NCOR,NCOR
      X3(I,J,:) = MATMUL(B(:,:,3,NI),X2(I,J,:))
   END DO
   END DO

   X3 = -CHARGE(NI)*X3*TINV/DELTA**3

   BT = TRANSPOSE(B(:,:,1,NI))

   DO K = -NCOR,NCOR
   DO J = -NCOR,NCOR
      Y1(:,J,K) = MATMUL(BT,X3(:,J,K))
   END DO
   END DO

   BT = TRANSPOSE(B(:,:,2,NI))
   DO K = -NCOR,NCOR
   DO I = -N,N
      Y2(I,:,K) = MATMUL(BT,Y1(I,:,K))
   END DO
   END DO

   BT = TRANSPOSE(B(:,:,3,NI))
   DO J = -N,N
   DO I = -N,N
      Y(I,J,:) = Y(I,J,:) + MATMUL(BT,Y2(I,J,:))
   END DO
   END DO

END DO

END SUBROUTINE PEMULT_CORR

!--------------------------------------------------------------------------------------

SUBROUTINE PEMULT_CORR_2(HAM,PARAMS,MPIDATA,B,X,Y,M,MB)

! Computes the nuclear potential matvec Y = V*X for a two-electron
! system when nuclei are not on grid points
! Analagous to KEMULT2

INTEGER, INTENT(IN) :: M,MB
TYPE(HAMILTONIAN), TARGET,  INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
REAL *8, INTENT(IN) :: B(:,:,:,:)
NUMBER, INTENT(IN) :: X(1:M,1:MB)
NUMBER, INTENT(OUT) :: Y(1:M,1:MB)
NUMBER :: W(1:M,1:MB), Z(1:M,1:MB)
INTEGER :: J,N,NCOR

N = PARAMS % N
NCOR = PARAMS % NCOR

!$OMP PARALLEL DO SHARED(HAM,PARAMS,MPIDATA,B,X,Y,N,NCOR) PRIVATE(J)
DO J = 1,MB
   CALL PEMULT_CORR(HAM,PARAMS,MPIDATA,B,X(:,J),Y(:,J),N,NCOR)
END DO
!$OMP END PARALLEL DO

CALL MPITRANSPOSE(X,W,M,MB,MPIDATA)

!$OMP PARALLEL DO SHARED(HAM,PARAMS,MPIDATA,B,W,Z,N,NCOR) PRIVATE(J)
DO J = 1,MB
  CALL PEMULT_CORR(HAM,PARAMS,MPIDATA,B,W(:,J),Z(:,J),N,NCOR)
END DO
!$OMP END PARALLEL DO

CALL MPITRANSPOSE(Z,W,M,MB,MPIDATA)

Y = Y + W

END SUBROUTINE PEMULT_CORR_2

!--------------------------------------------------------------------------------------

SUBROUTINE TWO_BODY_MULT(HAM,PARAMS,MPIDATA,X,Y,M,MB)

! Computes the electron-electron potential matvec when nuclei are not at grid points

TYPE(HAMILTONIAN), TARGET,  INTENT(IN) :: HAM
TYPE(PARAMETERS), INTENT(IN) :: PARAMS
TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
INTEGER, INTENT(IN) :: M,MB
NUMBER, INTENT(IN) :: X(1:M,MPIDATA % FIRST : MPIDATA % LAST)
NUMBER, INTENT(OUT) :: Y(1:M,MPIDATA % FIRST : MPIDATA % LAST)
INTEGER :: I1,J1,K1,I2,J2,K2,P,Q,FIRST,LAST,N,NN,NBIG
REAL *8 :: DELTA
REAL *8, POINTER :: TINV(:,:,:)

FIRST = MPIDATA % FIRST
LAST = MPIDATA % LAST
N = PARAMS % N
NN = PARAMS % NN
NBIG = PARAMS % NBIG
DELTA = PARAMS % DELTA
TINV(-NBIG:NBIG,-NBIG:NBIG,-NBIG:NBIG) => HAM % TINV

!$OMP PARALLEL DO SHARED(N,NN,X,Y,TINV,DELTA) PRIVATE(P,Q,I1,J1,K1,I2,J2,K2)
DO Q = FIRST,LAST

   K1 = (Q-1)/NN**2 + 1
   J1 = (Q-(K1-1)*NN**2-1)/NN + 1
   I1 = Q-(K1-1)*NN**2-(J1-1)*NN
   I1 = I1-N-1;  J1 = J1-N-1;  K1 = K1-N-1

   DO P = 1,M

      K2 = (P-1)/NN**2 + 1
      J2 = (P-(K2-1)*NN**2-1)/NN + 1
      I2 = P-(K2-1)*NN**2-(J2-1)*NN
      I2 = I2-N-1;  J2 = J2-N-1;  K2 = K2-N-1

      Y(P,Q) = X(P,Q)*TINV(I1-I2,J1-J2,K1-K2)/DELTA**3

   END DO

END DO
!$OMP END PARALLEL DO

END SUBROUTINE TWO_BODY_MULT

!--------------------------------------------------------------------------------------





subroutine mysort(in, values,n,lda,etarget)
  implicit none
  NUMBER :: in(lda,n), out(lda,n), values(lda), newvals(lda),etarget
  integer :: lda, n,i,j,whichlowest, flag
  real*8 :: lowval 
  integer :: taken(lda), order(lda)

  taken=0;  order=-1
  do j=1,n
     whichlowest=-1; flag=0;     lowval=1.d+30  !! is not used (see flag)
     do i=1,n
        if ( taken(i) .eq. 0 ) then
           if ((flag.eq.0) .or.(abs(values(i)-etarget) .le. lowval)) then
              flag=1;              lowval=abs(values(i)-etarget); whichlowest=i
           endif
        endif
     enddo
     if ((whichlowest.gt.n).or.(whichlowest.lt.1)) then
        write(*,*) taken;        write(*,*) 
         write(*,*) "WHICHLOWEST ERROR, J=",j," WHICHLOWEST=", whichlowest
         stop
     endif
     if (taken(whichlowest).ne.0) then
        write(*,*) "TAKEN ERROR.";         stop
     endif
     taken(whichlowest)=1;     order(j)=whichlowest
     out(:,j)=in(:,order(j));     newvals(j)=values(order(j))
  enddo

  values(1:n)=newvals(1:n)
  in(:,:)=out(:,:)

end subroutine mysort


END MODULE HAMMOD
