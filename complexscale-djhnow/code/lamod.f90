! parallel dense linear algebra
MODULE LAMOD
     !
     USE MPI
     !
     IMPLICIT NONE 
     !
     PRIVATE
     !
     PUBLIC ::                        &
        normalize,                    &
        crossmult_diag,               &
        crossmult,                    &
        colnorms,                     &
        cholQR                                    

 
CONTAINS     


   SUBROUTINE normalize(X, M, MB, k)
      ! Normalize all columns of X to have unit length 
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN)    :: M, MB, k
      COMPLEX*16, INTENT(INOUT) :: X(M*MB, k)
      !
      ! ... local variables
      !
      REAL*8  :: loc_scalprod(k)     ! local scalar products 
      REAL*8  :: colnorms(k)         ! column norms squared 
      INTEGER :: n, j, mpierr = 0
      !
      REAL*8 , EXTERNAL :: ddot
      !
      n = 2*M*MB
      DO j= 1, k
         loc_scalprod(j) = ddot(n, X(:,j), 1, X(:,j), 1)
      END DO
      !    
      CALL MPI_Allreduce(loc_scalprod, colnorms, k, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) 
      !
      DO j= 1, k
         X(:,j) = X(:,j)/SQRT(colnorms(j))
      END DO
      !
   END SUBROUTINE normalize


   SUBROUTINE crossmult_diag(X, Y, D, M, MB, k)
      ! Compute D = diag(X'Y) 
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN)    ::  M, MB, k
      COMPLEX*16, INTENT(IN)    ::  X(M*MB, k), Y(M*MB, k)
      COMPLEX*16, INTENT(OUT)    ::  D(k) 
      !
      ! ... local variables
      !
      COMPLEX*16  :: loc_scalprod(k)     ! local scalar products 
      INTEGER     :: n, j, mpierr = 0
      !
      COMPLEX*16, EXTERNAL :: zdotc
      !
      n = M*MB
      DO j= 1, k
        !loc_scalprod(j) = zdotc(SIZE(X,1), X(:,j), 1, Y(:,j), 1)
        loc_scalprod(j) = sum(conjg(X(:,j))*Y(:,j))
      END DO
      !    
      CALL MPI_Allreduce(loc_scalprod, D, k, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierr) 
      !
  END SUBROUTINE crossmult_diag

  
   SUBROUTINE colnorms(X, nrm, M, MB, k)
      ! Compute nrm = sqrt(diag(X'X)) 
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN)    ::  M, MB, k
      COMPLEX*16, INTENT(IN)    ::  X(M*MB, k)
      REAL*8, INTENT(OUT)       ::  nrm(k) 
      !
      ! ... local variables
      !
      REAL*8  :: loc_scalprod(k)     ! local scalar products 
      INTEGER :: n, j, mpierr = 0
      !
      REAL*8 , EXTERNAL :: ddot
      !
      n = 2*M*MB
      DO j= 1, k
         loc_scalprod(j) = ddot(n, X(:,j), 1, X(:,j), 1)
      END DO
      !    
      CALL MPI_Allreduce(loc_scalprod, nrm, k, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) 
      nrm = SQRT( nrm )
      !
  END SUBROUTINE colnorms


!  SUBROUTINE crossmult(X, Y, Z, M, MB, k)
  SUBROUTINE crossmult(X, k1, Y, k2, Z, ldz, M, MB)

      IMPLICIT NONE
      INTEGER,    INTENT(IN)       ::  M, MB, k1, k2, ldz ! (ldz >= k1)
      COMPLEX*16, INTENT(INOUT)    ::  X(M*MB, k1), Y(M*MB, k2)
      COMPLEX*16, INTENT(OUT)      ::  Z(ldz,*)
      COMPLEX*16  :: XY_loc(k1,k2)     ! local X'Y 
      INTEGER :: n, j, mpierr = 0
      COMPLEX*16  ::  recv_buf(k1,k2)
      
      n  = M*MB
      
      CALL zgemm('C','N', k1, k2, n, dcmplx(1.0), X, n, Y, n, dcmplx(0.0), XY_loc, k1)
      CALL MPI_Allreduce(XY_loc, recv_buf, k1*k2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierr) 

      DO j = 1 , k2
         Z(1:k1,j) = recv_buf(:,j)  
      END DO 
        
  END SUBROUTINE crossmult


  SUBROUTINE cholQR(X, k, R, ldr, M, MB)
      ! Compute Cholesky based QR decomposition of X 
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN)       ::  M, MB, k, ldr
      COMPLEX*16, INTENT(INOUT)    ::  X(M*MB, k)
!      COMPLEX*16, INTENT(OUT)      ::  R(ldr,*) 

      COMPLEX*16, INTENT(OUT)      ::  R(ldr,k) ! upper right triangule  is the R matrix in the QR decomposition

      !
      ! ... local variables
      !
      COMPLEX*16  :: XTX_loc(k,k)     ! local X'X 
      real*8 :: rsum
      INTEGER :: n,  mpierr = 0, info = 0,ii,jj
      integer, save :: fff=0
      !
      n  = M*MB
      !
      ! R = X'X
      !
!      CALL zgemm('C','N', k, k, n, dcmplx(1.0), X, n, X, n, dcmplx(0.0), XTX_loc, k)
!      CALL MPI_Allreduce(XTX_loc, R, k**2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierr) 
      CALL crossmult(X, k, X, k, R, ldr, M, MB)
      !    
      ! ... Run Cholesky

      rsum=0
      do ii=1,k
      do jj=1,k
         rsum=rsum+abs(R(jj,ii))
      enddo
      enddo
      if (rsum.eq.0d0.and.fff.eq.0) then
         print *, 'Warning rsum=0 in CholQ'
         fff=1
      else
         CALL zpotrf('U', k, R, ldr, info)
      endif
      IF (info /= 0) THEN
         WRITE(*,*) 'Cholesky failed!!!', info
         do ii=1,min(5,k)
            write(*,'(1000E10.2)') R(1:min(5,k),ii)
         enddo
        STOP
      END IF
      !
      ! ... Perform triangular solve
      CALL ztrsm('R','U', 'N', 'N', n, k, dcmplx(1.0), R, ldr, X, n)  
      ! 
  END SUBROUTINE cholQR


END MODULE LAMOD



