MODULE FOURIER
IMPLICIT NONE

CONTAINS

SUBROUTINE FOURIER_INVERSE(TINV,NB)

INTEGER, INTENT(IN) :: NB
REAL *8, INTENT(OUT) :: TINV(-NB:NB,-NB:NB,-NB:NB)

INTEGER :: I,J,K,L
REAL *8 :: PI
COMPLEX *16 :: X(-NB:NB,-NB:NB,-NB:NB), Y(-NB:NB,-NB:NB,-NB:NB)
COMPLEX *16 :: T1,T2,T3,T(-NB:NB)

PI = 4D0*ATAN(1D0)

DO J = -NB,NB
	T(J) = ((-1D0)**J)/J**2
	IF(J == 0) T(J) = PI**2/6D0
END DO

DO K = -NB,NB
DO J = -NB,NB
DO I = -NB,NB
	T1 = 0; T2 = 0; T3 = 0
	IF(J == 0 .AND. K == 0) T1 = T(I)
	IF(I == 0 .AND. K == 0) T2 = T(J)
	IF(I == 0 .AND. J == 0) T3 = T(K)
	X(I,J,K) = T1+T2+T3
END DO
END DO
END DO


CALL DFT3(X,Y,NB)
CALL IDFT3(1D0/Y,X,NB)
TINV = 2D0*PI*DREAL(X)


!WRITE(*,"(F20.10)") DREAL(X(0,0,0))
!OPEN(4,FILE="C:\Users\Jeremiah\Postdoc\MATLAB\tinv.txt")
!WRITE(4,"(F20.10)") DREAL(X)
!CLOSE(4)

END SUBROUTINE FOURIER_INVERSE

! ---------------------------------------------------------------

SUBROUTINE DFT1(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N)
INTEGER :: J,K,NN
COMPLEX *16 :: A
REAL *8 :: PI

NN = 2*N+1
PI = 4D0*ATAN(1D0)
A = -2D0*PI*DCMPLX(0D0,1D0)/NN

DO K = -N,N
	Y(K) = 0
	DO J = -N,N
		Y(K) = Y(K) + X(J)*EXP(A*J*K)
   END DO
END DO

END SUBROUTINE DFT1

! ---------------------------------------------------------------

SUBROUTINE IDFT1(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N)
INTEGER :: J,K,NN
COMPLEX *16 :: A
REAL *8 :: PI

NN = 2*N+1
PI = 4D0*ATAN(1D0)
A = 2D0*PI*DCMPLX(0D0,1D0)/NN

DO K = -N,N
	Y(K) = 0
	DO J = -N,N
		Y(K) = Y(K) + X(J)*EXP(A*J*K)/NN
   END DO
END DO

END SUBROUTINE IDFT1

! ---------------------------------------------------------------

SUBROUTINE DFT3(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N,-N:N,-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N,-N:N,-N:N)

COMPLEX *16 :: Y1(-N:N,-N:N,-N:N)
INTEGER :: I,J

DO J = -N,N
DO I = -N,N
	CALL DFT1(X(:,I,J),Y(:,I,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
	CALL DFT1(Y(I,:,J),Y1(I,:,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
	CALL DFT1(Y1(I,J,:),Y(I,J,:),N)
END DO
END DO

END SUBROUTINE DFT3

! ---------------------------------------------------------------

SUBROUTINE IDFT3(X,Y,N)

INTEGER, INTENT(IN) :: N
COMPLEX *16, INTENT(IN) :: X(-N:N,-N:N,-N:N)
COMPLEX *16, INTENT(OUT) :: Y(-N:N,-N:N,-N:N)

COMPLEX *16 :: Y1(-N:N,-N:N,-N:N)
INTEGER :: I,J

DO J = -N,N
DO I = -N,N
	CALL IDFT1(X(:,I,J),Y(:,I,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
	CALL IDFT1(Y(I,:,J),Y1(I,:,J),N)
END DO
END DO

DO J = -N,N
DO I = -N,N
	CALL IDFT1(Y1(I,J,:),Y(I,J,:),N)
END DO
END DO

END SUBROUTINE IDFT3

! ---------------------------------------------------------------

END MODULE FOURIER