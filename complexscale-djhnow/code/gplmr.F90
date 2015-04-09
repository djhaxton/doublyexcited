#include "definitions.INC"

MODULE GPLMR

CONTAINS

#if ISREAL==0

SUBROUTINE GPLMR_EIG(HAM,PARAMS,MPIDATA,GMSOLVER,X,E,M,MB,SIGMAIN,NEV,NS,MAXIT,TOL,PRECONDITIONER,OUTMODULUS)
   USE HAMMOD
   USE STRUCTS
   USE SOLVER
   USE LAMOD   
   USE MPI  
#if defined(USE_PETSC)
   USE PETSCWRAP
#endif
   IMPLICIT NONE 
   INTEGER, INTENT(IN)      ::  M, MB, NEV, NS, MAXIT ,OUTMODULUS
   TYPE(HAMILTONIAN), INTENT(IN) :: HAM
   TYPE(PARAMETERS), INTENT(IN) :: PARAMS
   TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
   TYPE(GMRES_SOLVER), INTENT(INOUT) :: GMSOLVER
   REAL *8, INTENT(IN) :: TOL
   COMPLEX*16, INTENT(IN)   ::  SIGMAIN
   COMPLEX*16 :: SIGMA
   COMPLEX*16, INTENT(INOUT)      ::  X(M, MB, NEV) 
   COMPLEX*16, INTENT(OUT)        ::  E(NEV)      
   EXTERNAL :: PRECONDITIONER

   COMPLEX*16      ::  W(M, MB, NEV), S(M, MB, NEV*NS), P(M, MB, NEV)
   COMPLEX*16      ::  AX(M, MB, NEV), AW(M, MB, NEV), AS(M, MB, NEV*NS), AP(M, MB, NEV)   
   COMPLEX*16      ::  ASIGMAX(M, MB, NEV), ASIGMAW(M, MB, NEV), ASIGMAS(M, MB, NEV*NS), ASIGMAP(M, MB, NEV)   
   COMPLEX*16      ::  KP((3+NS)*NEV, (3+NS)*NEV), MP((3+NS)*NEV, (3+NS)*NEV), BUF(M,MB,NEV)  
   COMPLEX*16      ::  KPLOC((3+NS)*NEV, (3+NS)*NEV), MPLOC((3+NS)*NEV, (3+NS)*NEV) 
   COMPLEX*16      ::  R(NS*NEV, NEV), COORDX(NEV,NEV), COORDW(NEV,NEV), COORDS(NEV*NS,NEV), COORDP(NEV,NEV) 
   REAL*8          ::  RESNRM(NEV), MAXRES
   INTEGER         ::  ITER, J, J1, J2, TMP, JPREV, CUR_S_SIZE, L, MMB, NPROJ, ELEM_COUNT, LDR, MAXBLOCK, N
   INTEGER         ::  ME, IERROR
   INTEGER, PARAMETER :: ROOT = 0 !! ??? why hardwire ??? huge dimension! too much memory. MAX_OUTER_ITER = 200

   REAL*8            :: TGEMM, TGEMM0, TGEMM1,      &
                      TSPMV, TSPMV0, TSPMV1,      &
                      TPREC, TPREC0, TPREC1,      &
                      TOT_TIME, TOT_TIME0, TOT_TIME1

   SIGMA=SIGMAIN

   TGEMM = 0.0 
   TSPMV = 0.0 
   TPREC = 0.0 
   TOT_TIME = 0.0

   MAXBLOCK = MPIDATA % MAXBLOCK
   N = PARAMS % NN

   TOT_TIME0 = MPI_WTIME() 

   LDR  = NS*NEV          
   ITER = 1  
   MMB = M*MB
   NPROJ = (3+NS)*NEV      
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, IERROR) 

   IF (ME == ROOT) WRITE(*,*) "GPLMR_EIG STARTED ... "
   IF (ME == ROOT) WRITE(*,*) "THE TARGETED SHIFT IS "
   IF (ME == ROOT) PRINT*, SIGMA
  
   CALL NORMALIZE(X, M, MB, NEV)

   TSPMV0 = MPI_WTIME()
   CALL MATBLOCK(HAM,PARAMS,MPIDATA,X,AX,M,MB,NEV)
   TSPMV1 = MPI_WTIME()
   TSPMV  = TSPMV + (TSPMV1 - TSPMV0)

   CALL CROSSMULT_DIAG(X, AX, E, M, MB, NEV)

   DO J = 1, NEV 
      W(:, :, J)  = AX(: , :, J) - X(:, :, J)*E(J)  
   END DO 

   CALL COLNORMS(W, RESNRM, M, MB, NEV)  

   MAXRES = MAXVAL(RESNRM)
   
   DO WHILE (ITER <= MAXIT .AND. MAXRES > TOL)

      CALL CHOLQR(X, NEV, R, LDR, M, MB)
      
      CALL ZTRSM('R','U', 'N', 'N', MMB, NEV, DCMPLX(1.0), R, LDR, AX, MMB)  
      
      IF (GMSOLVER%MAX_INNER_ITER.GT.0) THEN
         TPREC0 = MPI_WTIME()
         CALL PRECBLOCK(HAM,PARAMS,MPIDATA,GMSOLVER,W,BUF,SIGMA,&
              M,MB,NEV,GMSOLVER%MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
         TPREC1 = MPI_WTIME()
         TPREC  = TPREC + (TPREC1 - TPREC0)
         W = BUF
      ENDIF

      CALL CROSSMULT(X, NEV, W, NEV, R, LDR, M, MB) 
      
      CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), X, MMB, R, LDR, DCMPLX(1.0), W, MMB) 
      
      CALL CHOLQR(W, NEV, R, LDR, M, MB)       
      TSPMV0 = MPI_WTIME()
      
      CALL MATBLOCK(HAM,PARAMS,MPIDATA,W,AW,M,MB,NEV)
      TSPMV1 = MPI_WTIME()
      TSPMV  = TSPMV + (TSPMV1 - TSPMV0)
      
      IF ( NS > 0 ) THEN
         
         DO J = 1, NEV
            S(:, :, J) = AW(:, :, J) - W(:, :, J)*E(J)
         END DO
         
         IF (GMSOLVER%MAX_INNER_ITER.GT.0) THEN
            TPREC0 = MPI_WTIME()
            CALL PRECBLOCK(HAM,PARAMS,MPIDATA,GMSOLVER,S,BUF,SIGMA, &
                 M,MB,NEV,GMSOLVER%MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
            TPREC1 = MPI_WTIME() 
            TPREC  = TPREC + (TPREC1 - TPREC0)
            S(:, :, 1:NEV) = BUF
         ENDIF

         CALL CROSSMULT(X, NEV, S, NEV, R, LDR, M, MB)                                           ! R = X'*S
         CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), X, MMB, R, LDR, DCMPLX(1.0), S, MMB) ! S = S - XR
         CALL CROSSMULT(W, NEV, S, NEV, R, LDR, M, MB)                                                  ! R = W'*S
         CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), W, MMB, R, LDR, DCMPLX(1.0), S, MMB) ! S = S - WR
         CALL CHOLQR(S, NEV, R, LDR, M, MB)       
         TSPMV0 = MPI_WTIME()

         CALL MATBLOCK(HAM,PARAMS,MPIDATA,S,AS,M,MB,NEV)
         TSPMV1 = MPI_WTIME()
         TSPMV  = TSPMV + (TSPMV1 - TSPMV0)       

         IF ( NS > 1 ) THEN
             ! ...  CONSTRUCT MORE
             DO L = 2, NS
                J1 = (L-1)*NEV+1
                J2 = L*NEV
                CUR_S_SIZE = J1 - 1   ! CURRENT SIZE OF S BLOCK 
                TMP = 1
                DO J = J1, J2
                   JPREV = J - NEV
                   S(:, :, J) = AS(:, :, JPREV) - S(:, :, JPREV)*E(TMP)
                   TMP = TMP + 1
                END DO

                IF (GMSOLVER%MAX_INNER_ITER.GT.0) THEN
                   TPREC0 = MPI_WTIME()
                   CALL PRECBLOCK(HAM,PARAMS,MPIDATA,GMSOLVER,S(:,:,J1),BUF,SIGMA, &
                        M,MB,NEV,GMSOLVER%MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
                   TPREC1 = MPI_WTIME() 
                   TPREC  = TPREC + (TPREC1 - TPREC0)
                   S(:, :, J1:J2) = BUF
                ENDIF

                CALL CROSSMULT(X, NEV, S(:,:,J1), NEV, R, LDR, M, MB)                                                   ! R = X'*S
                CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), X, MMB, R, LDR, DCMPLX(1.0), S(:,:,J1), MMB)          ! S = S - XR
                CALL CROSSMULT(W, NEV, S(:,:,J1), NEV, R, LDR, M, MB)                                                        ! R = W'*S
                CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), W, MMB, R, LDR, DCMPLX(1.0), S(:,:,J1), MMB)     ! S = S - WR
                CALL CROSSMULT(S, CUR_S_SIZE, S(:,:,J1), NEV, R, LDR, M, MB)                                       ! R = S_PREV'*S
                CALL ZGEMM('N', 'N', MMB, NEV, CUR_S_SIZE, DCMPLX(-1.0), S, MMB, R, LDR, DCMPLX(1.0), S(:,:,J1), MMB)     ! S = S - SPREVR

                CALL CHOLQR(S(:,:,J1), NEV, R, LDR, M, MB)
                TSPMV0 = MPI_WTIME() 
                CALL MATBLOCK(HAM,PARAMS,MPIDATA,S(:,:,J1),AS(:,:,J1),M,MB,NEV)
                TSPMV1 = MPI_WTIME() 
                TSPMV  = TSPMV + (TSPMV1 - TSPMV0)
             END DO
        END IF 
        
     END IF

     IF (ITER > 1 ) THEN
        !
        ! ... ORTHOGONOLIZE P AGAINST X AND W
        CALL CROSSMULT(X, NEV, P, NEV, R, LDR, M, MB)                                                      ! R = X'*P
        CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), X, MMB, R, LDR, DCMPLX(1.0), P, MMB)     ! P = P - XR
        CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), AX, MMB, R, LDR, DCMPLX(1.0), AP, MMB)   ! AP = AP - AXR
        CALL CROSSMULT(W, NEV, P, NEV, R, LDR, M, MB)                                                      ! R = W'*P
        CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), W, MMB, R, LDR, DCMPLX(1.0), P, MMB)     ! P = P - WR
        CALL ZGEMM('N', 'N', MMB, NEV, NEV, DCMPLX(-1.0), AW, MMB, R, LDR, DCMPLX(1.0), AP, MMB)   ! AP = AP - AWR
        !
        ! ... ORTHOGONOLIZE P AGAINST ALL S BLOCKS  
        CALL CROSSMULT(S, LDR, P, NEV, R, LDR, M, MB)                                                      ! R = S'*P
        CALL ZGEMM('N', 'N', MMB, NEV, LDR, DCMPLX(-1.0), S, MMB, R, LDR, DCMPLX(1.0), P, MMB)     ! P = P - SR
        CALL ZGEMM('N', 'N', MMB, NEV, LDR, DCMPLX(-1.0), AS, MMB, R, LDR, DCMPLX(1.0), AP, MMB)   ! AP = AP - ASR
        !
        ! ORTHONORMALIZE COLUMNS OF P AND UPDATE AP ACCORDINGLY
        CALL CHOLQR(P, NEV, R, LDR, M, MB)       
        CALL ZTRSM('R','U', 'N', 'N', MMB, NEV, DCMPLX(1.0), R, LDR, AP, MMB)  
        !
     END IF
     !
     ! ... CONSTRUCT HARMONIC PROBLEM
     !

     ASIGMAX = AX - SIGMA*X 
     ASIGMAW = AW - SIGMA*W 
     ASIGMAS = AS - SIGMA*S
     !
     ! ... CONSTRUCT DIAGONAL BLOCKS OF K AND M
     !
     ! K(1:NEV, 1:NEV) = ASIGMAX'*ASIGMAX, M(1:NEV, 1:NEV) = ASIGMAX'*X     
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAX, MMB, ASIGMAX, MMB, DCMPLX(0.0), KPLOC, NPROJ)
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAX, MMB, X, MMB, DCMPLX(0.0), MPLOC, NPROJ)
     ! K(NEV+1:2*NEV, NEV+1:2*NEV) = ASIGMAW'*ASIGMAW, M(NEV+1:2*NEV, NEV+1:2*NEV) = ASIGMAW'*W    
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAW, MMB, ASIGMAW, MMB, DCMPLX(0.0), KPLOC(NEV+1, NEV+1), NPROJ)
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAW, MMB, W, MMB, DCMPLX(0.0), MPLOC(NEV+1, NEV+1), NPROJ)
     ! K(2*NEV+1:3*NEV, 2*NEV+1:3*NEV) = ASIGMAS'*ASIGMAS, M(2*NEV+1:3*NEV, 2*NEV+1:3*NEV) = ASIGMAS'*S    
     CALL ZGEMM('C','N', LDR, LDR, MMB, DCMPLX(1.0), ASIGMAS, MMB, ASIGMAS, MMB, DCMPLX(0.0), KPLOC(2*NEV+1, 2*NEV+1), NPROJ)
     CALL ZGEMM('C','N', LDR, LDR, MMB, DCMPLX(1.0), ASIGMAS, MMB, S, MMB, DCMPLX(0.0), MPLOC(2*NEV+1, 2*NEV+1), NPROJ)
     !
     ! ... CONSTRUCT UPPER OFF-DIAGONAL BLOCKS OF K AND M
     !
     ! K(1:NEV, NEV+1:2*NEV) = ASIGMAX'*ASIGMAW, M(1:NEV, NEV+1:2*NEV) = ASIGMAX'*W     
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAX, MMB, ASIGMAW, MMB, DCMPLX(0.0), KPLOC(1,NEV+1), NPROJ)
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAX, MMB, W, MMB, DCMPLX(0.0), MPLOC(1,NEV+1), NPROJ)
     ! K(1:NEV, 2NEV+1:3*NEV) = ASIGMAX'*ASIGMAS, M(1:NEV, 2*NEV+1:3*NEV) = ASIGMAX'*S     
     CALL ZGEMM('C','N', NEV, LDR, MMB, DCMPLX(1.0), ASIGMAX, MMB, ASIGMAS, MMB, DCMPLX(0.0), KPLOC(1,2*NEV+1), NPROJ)
     CALL ZGEMM('C','N', NEV, LDR, MMB, DCMPLX(1.0), ASIGMAX, MMB, S, MMB, DCMPLX(0.0), MPLOC(1,2*NEV+1), NPROJ)
     ! K(NEV+1:2*NEV, 2NEV+1:3*NEV) = ASIGMAW'*ASIGMAS, M(NEV+1:2NEV, 2*NEV+1:3*NEV) = ASIGMAW'*S     
     CALL ZGEMM('C','N', NEV, LDR, MMB, DCMPLX(1.0), ASIGMAW, MMB, ASIGMAS, MMB, DCMPLX(0.0), KPLOC(NEV+1,2*NEV+1), NPROJ)
     CALL ZGEMM('C','N', NEV, LDR, MMB, DCMPLX(1.0), ASIGMAW, MMB, S, MMB, DCMPLX(0.0), MPLOC(NEV+1,2*NEV+1), NPROJ)
     !
     ! ... CONSTRUCT LOWER OFF-DIAGONAL BLOCKS OF K AND M
     !
     ! K(NEV+1:2*NEV, 1:NEV) = K(1:NEV, NEV+1:2*NEV)', M(NEV+1:2*NEV, 1:NEV+1) = W'*ASIGMAX     
     KPLOC(NEV+1:2*NEV, 1:NEV)   = CONJG(TRANSPOSE(KPLOC(1:NEV, NEV+1:2*NEV)))  
     CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAW, MMB, X, MMB, DCMPLX(0.0), MPLOC(NEV+1,1), NPROJ)
     ! K(2NEV+1:3*NEV, 1:NEV) = K(1:NEV, 2NEV+1:3*NEV)', M(2NEV+1:3NEV, 1:NEV) = S'*ASIGMAX     
     KPLOC(2*NEV+1:2*NEV+LDR, 1:NEV)   = CONJG(TRANSPOSE(KPLOC(1:NEV, 2*NEV+1:2*NEV+LDR)))  
     CALL ZGEMM('C','N', LDR, NEV, MMB, DCMPLX(1.0), ASIGMAS, MMB, X, MMB, DCMPLX(0.0), MPLOC(2*NEV+1,1), NPROJ)
     ! K(2NEV+1:3*NEV, NEV+1:2NEV) = K(NEV+1:2NEV, 2NEV+1:3*NEV)', M(2NEV+1:3NEV, NEV+1:2*NEV) = S'*ASIGMAW     
     KPLOC(2*NEV+1:2*NEV+LDR, NEV+1:2*NEV)   = CONJG(TRANSPOSE(KPLOC(NEV+1:2*NEV, 2*NEV+1:2*NEV+LDR)))  
     CALL ZGEMM('C','N', LDR, NEV, MMB, DCMPLX(1.0), ASIGMAS, MMB, W, MMB, DCMPLX(0.0), MPLOC(2*NEV+1,NEV+1), NPROJ)
     !
     ! ... FINISH CONSTRUCTION OF THE PROJECTED MATRICES
     !
     IF (ITER > 1 ) THEN
        ASIGMAP = AP - SIGMA*P 
        ! K(1:NEV, 3NEV+1:4*NEV) = ASIGMAX'*ASIGMAP, M(1:NEV, 3*NEV+1:4*NEV) = ASIGMAX'*P     
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAX, MMB, ASIGMAP, MMB, DCMPLX(0.0), KPLOC(1,2*NEV+LDR+1), NPROJ)
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAX, MMB, P, MMB, DCMPLX(0.0), MPLOC(1,2*NEV+LDR+1), NPROJ)
        ! K(NEV+1:2NEV, 3NEV+1:4*NEV) = ASIGMAW'*ASIGMAP, M(NEV+1:2NEV, 3*NEV+1:4*NEV) = ASIGMAW'*P     
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAW, MMB, ASIGMAP, MMB, DCMPLX(0.0), KPLOC(NEV+1,2*NEV+LDR+1), NPROJ)
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAW, MMB, P, MMB, DCMPLX(0.0), MPLOC(NEV+1,2*NEV+LDR+1), NPROJ)
        ! K(2NEV+1:3NEV, 3NEV+1:4*NEV) = ASIGMAS'*ASIGMAP, M(2NEV+1:3NEV, 3*NEV+1:4*NEV) = ASIGMAS'*P     
        CALL ZGEMM('C','N', LDR, NEV, MMB, DCMPLX(1.0), ASIGMAS, MMB, ASIGMAP, MMB, DCMPLX(0.0), KPLOC(2*NEV+1,2*NEV+LDR+1), NPROJ)
        CALL ZGEMM('C','N', LDR, NEV, MMB, DCMPLX(1.0), ASIGMAS, MMB, P, MMB, DCMPLX(0.0), MPLOC(2*NEV+1,2*NEV+LDR+1), NPROJ)
        ! K(3*NEV+1:4*NEV+1, 3*NEV+1:4*NEV) = ASIGMAP'*ASIGMAP, M(3*NEV+1:4*NEV+1, 3*NEV+1:4*NEV) = ASIGMAP'*P   
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAP, &
                    MMB, ASIGMAP, MMB, DCMPLX(0.0), KPLOC(2*NEV+LDR+1, 2*NEV+LDR+1), NPROJ)
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAP, MMB, P, MMB, DCMPLX(0.0), MPLOC(2*NEV+LDR+1, 2*NEV+LDR+1), NPROJ)
        !
        ! K(3NEV+1:4*NEV,1:NEV) = K(1:NEV, 3NEV+1:4*NEV)', M(3*NEV+1:4*NEV, 1:NEV) = P'*ASIGMAX    
        KPLOC(2*NEV+LDR+1:NPROJ, 1:NEV)   = CONJG(TRANSPOSE(KPLOC(1:NEV, 2*NEV+LDR+1:NPROJ)))  
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAP, MMB, X, MMB, DCMPLX(0.0), MPLOC(2*NEV+LDR+1,1), NPROJ)
        ! K(3NEV+1:4*NEV,NEV+1:2NEV) = K(NEV+1:2NEV, 3NEV+1:4*NEV)', M(3*NEV+1:4*NEV, NEV+1:2NEV) = P'*ASIGMAW    
        KPLOC(2*NEV+LDR+1:NPROJ, NEV+1:2*NEV)   = CONJG(TRANSPOSE(KPLOC(NEV+1:2*NEV, 2*NEV+LDR+1:NPROJ)))  
        CALL ZGEMM('C','N', NEV, NEV, MMB, DCMPLX(1.0), ASIGMAP, MMB, W, MMB, DCMPLX(0.0), MPLOC(2*NEV+LDR+1,NEV+1), NPROJ)
        ! K(3NEV+1:4*NEV,2NEV+1:3NEV) = K(2NEV+1:3NEV, 3NEV+1:4*NEV)', M(3*NEV+1:4*NEV, 2NEV+1:2NEV) = P'*ASIGMAS    
        KPLOC(2*NEV+LDR+1:NPROJ, 2*NEV+1:2*NEV+LDR)   = CONJG(TRANSPOSE(KPLOC(2*NEV+1:2*NEV+LDR, 2*NEV+LDR+1:NPROJ)))  
        CALL ZGEMM('C','N', NEV, LDR, MMB, DCMPLX(1.0), ASIGMAP, MMB, S, MMB, DCMPLX(0.0), MPLOC(2*NEV+LDR+1,2*NEV+1), NPROJ)
        ! 
     END IF
     !
     ! ... ALLREDUCE KPLOC AND MPLOC
     !
     ELEM_COUNT = NPROJ**2 
     !

     CALL MPI_ALLREDUCE(KPLOC, KP, ELEM_COUNT, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERROR) 
     CALL MPI_ALLREDUCE(MPLOC, MP, ELEM_COUNT, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERROR) 
     !
     ! ... SOLVE THE NON-HERMITIAN GENERALIZED PROJECTED EIGENVALUE PROBLEM
     !
!PRINT *, MP
!IF (ME == 0) PRINT*, KP
!IF (ME == 0) PRINT*, MP


     IF (ITER > 1) THEN
        CALL SOLVE_GEVP(KP, MP, NPROJ, NPROJ, NEV)
     ELSE
        CALL SOLVE_GEVP(KP, MP, NPROJ, NPROJ-NEV, NEV)        
     END IF

     COORDX = KP(1:NEV,1:NEV) 
     COORDW = KP(NEV+1:2*NEV,1:NEV) 
     COORDS = KP(2*NEV+1:2*NEV+LDR,1:NEV) 
     ! 
     ! ... FORM P 
     !
     IF (ITER > 1) THEN
        COORDP = KP(2*NEV+LDR+1 : NPROJ,1:NEV) 
        CALL ZGEMM('N','N', MMB, NEV, NEV, DCMPLX(1.0), P, MMB, COORDP, NEV, DCMPLX(0.0), BUF, MMB)
        P = BUF
        CALL ZGEMM('N','N', MMB, NEV, NEV, DCMPLX(1.0), AP, MMB, COORDP, NEV, DCMPLX(0.0), BUF, MMB) 
        AP = BUF
     ELSE
        P  = DCMPLX(0.0)          
        AP = DCMPLX(0.0)          
     END IF
     ! P = P + S*COORDS, AP = AP + AS*COORDS
     CALL ZGEMM('N','N', MMB, NEV, LDR, DCMPLX(1.0), S, MMB, COORDS, LDR, DCMPLX(1.0), P, MMB)
     CALL ZGEMM('N','N', MMB, NEV, LDR, DCMPLX(1.0), AS, MMB, COORDS, LDR, DCMPLX(1.0), AP, MMB)
     ! P = P + W*COORDW, AP = AP + AW*COORDW
     CALL ZGEMM('N','N', MMB, NEV, NEV, DCMPLX(1.0), W, MMB, COORDW, NEV, DCMPLX(1.0), P, MMB)
     CALL ZGEMM('N','N', MMB, NEV, NEV, DCMPLX(1.0), AW, MMB, COORDW, NEV, DCMPLX(1.0), AP, MMB)
     !
     ! ... UPDATE X AND AX:  X = X*COORDX + P, AX = AX*COORDX + AP
     !
     CALL ZGEMM('N','N', MMB, NEV, NEV, DCMPLX(1.0), X, MMB, COORDX, NEV, DCMPLX(0.0), BUF, MMB)
     X = BUF + P
     CALL ZGEMM('N','N', MMB, NEV, NEV, DCMPLX(1.0), AX, MMB, COORDX, NEV, DCMPLX(0.0), BUF, MMB)
     AX = BUF + AP
     !
     ! ... COMPUTE NEW RESIDUAL  
!     CALL NORMALIZE(X, M, MB, NEV)
     CALL COLNORMS(X, RESNRM, M, MB, NEV)         ! TEMPORARILY STORE NORMS OF X'S COLUMNS IN RESNRM 
     !
     ! NORMALIZE COLUMNS OF X AND UPDATE AX ACCORDINGLY    
     DO J = 1, NEV
        X(:,:,J)  = X(:,:,J)/RESNRM(J)
        AX(:,:,J) = AX(:,:,J)/RESNRM(J)
     END DO
     !
     CALL CROSSMULT_DIAG(X, AX, E, M, MB, NEV)    ! COMPUTE RQS E = DIAG(X'AX) 
     !
     DO J = 1, NEV 
        W(:, :, J)  = AX(: , :, J) - X(:, :, J)*E(J)  
     END DO 
     !
     ! REPORT RESIDUAL NORMS
     CALL COLNORMS(W, RESNRM, M, MB, NEV)   ! RESNRM = SQRT( DIAG(W'*W) )
     !
     MAXRES = MAXVAL(RESNRM)

     IF (ME == ROOT.and.(maxres.le.tol.or.mod(ITER,outmodulus).eq.0)) THEN
        WRITE(*,'( "ITERATION ", I5)') ITER
     !

        DO J = 1, NEV

           WRITE(*,"(A,I4,A,F20.12,F20.12,A,ES10.2E2)") "EIG ",J," = ",DREAL(E(J)),DIMAG(E(J)),",  ERROR = ",RESNRM(J)

        END DO
     ENDIF
        
     !
     ITER = ITER + 1

!!!  PRINT ITERATION TIMINGS 
!
     TOT_TIME1 = MPI_WTIME() 
     TOT_TIME  = TOT_TIME + (TOT_TIME1 - TOT_TIME0)

     IF (ME == ROOT.and.(maxres.le.tol.or.mod(ITER,outmodulus).eq.0)) THEN
        WRITE(*,'("TOTAL  TIME   = ", F7.1)') TOT_TIME
        WRITE(*,'("MATVEC  TIME  = ", F7.1,  " PERCENTAGE ", F4.1)')  TSPMV, 100*TSPMV/TOT_TIME
        WRITE(*,'("INVERSE  TIME = ", F7.1, " PERCENTAGE ", F4.1)') TPREC, 100*TPREC/TOT_TIME
     ENDIF

     TOT_TIME0 = MPI_WTIME() 

!!!  END PRINT ITERATION TIMINGS 
!! moving this up
!   MAXRES = MAXVAL(RESNRM)

   END DO
   ! ... END MAIN LOOP

!TOT_TIME1 = MPI_WTIME() 
!TOT_TIME  = TOT_TIME + (TOT_TIME1 - TOT_TIME0)
!
!   ! ... PRINT OUT TIMINGS
!   WRITE(*,'("TOTAL  TIME = ", F5.1)'), TOT_TIME
!   WRITE(*,'("MATVEC  TIME = ", F5.1,  " PERCENTAGE ", F4.1)'),  TSPMV, 100*TSPMV/TOT_TIME
!   WRITE(*,'("PREC  TIME = ", F5.1, " PERCENTAGE ", F4.1)'), TPREC, 100*TPREC/TOT_TIME
!!   WRITE(*,'("DGEMM  TIME = ", F5.1, " PERCENTAGE ", F4.1)'), TGEMM, 100*TGEMM/TOT_TIME

100 CONTINUE

!---------------------------------------------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------------------------------------------

    ! FIND (RIGHT) EIGENVECTORS OF (K,M) ASSOCIATED WITH NEV SMALLEST EIGENVALUES 
    SUBROUTINE SOLVE_GEVP(KP, MP, LD, K, NEV)
       !
       !
       IMPLICIT NONE 
       !
       ! ... I/O VARIABLES
       !
       ! LEADING DIMENSION OF K AND M, THEIR SIZE, AND THE DESIRED NUMBER OF SMALLEST PAIRS OF (K,M)
       INTEGER, INTENT(IN)          ::  LD, K, NEV
       COMPLEX *16, INTENT(INOUT)   ::  KP(LD, *)  ! ON ENTRY, CONTAINS K OF THE PENCIL (K,M). 
                                                  ! ON EXIT, INITIAL NEV COLUMNS ARE THE TARGETED EIGEVECTORS 
       COMPLEX *16, INTENT(IN)      ::  MP(LD, *)  ! CONTAINS M OF THE PENCIL (K,M)
       !
       ! ... LOCAL VARIABLES
       !
!       COMPLEX*16                   :: ALPHA(K), BETA(K), VR(LD, K)
!       COMPLEX*16, ALLOCATABLE      :: WORK(:), VL(:,:) ! NOT REFERENCED IN ZGGEV
       COMPLEX*16                   :: ALPHA(K), BETA(K), VR(LD, K), VL(1,1)
       COMPLEX*16, ALLOCATABLE      :: WORK(:)
       INTEGER                      :: LWORK, INFO, IDX(K)
       REAL*8                       :: RWORK(8*K), ALAMBDA(K)
       !
       LWORK = 10*MAX(1,2*K)
       ALLOCATE( WORK(LWORK) ) 
       !
       CALL ZGGEV('N', 'V', K, KP, LD, MP, LD, ALPHA, BETA, VL, LD, VR, LD, WORK, LWORK, RWORK, INFO)
       !
       IF (INFO == 0) THEN
         ! ... SORT EIGENVALUES IN ASCENDINGLY IN THEIR ABSOLUTE VALUE 
         ALAMBDA = ABS(ALPHA/BETA)                                !!! ASSUME THAT BETA ARE NOT TOO SMALL IN ABS VALUE  
         !CALL KB07AD(ALAMBDA, K, IDX)   ! SORT EIGENVALUES BY THEIR MODULUS 
         CALL EIGSORT(ABS(ALAMBDA),IDX,SIZE(ALAMBDA),K)
         KP(:,1:NEV) = VR(:,IDX(1:NEV))            
!ALAMBDA = ALPHA/BETA + SIGMA
!PRINT*, ALAMBDA(IDX)
       ELSE
          WRITE(*,*) "SOLVE OF THE PROJECTED EIGENVALUE PROBLEM FAILED"
          STOP
       END IF
                 
       
       DEALLOCATE( WORK ) 


    END SUBROUTINE SOLVE_GEVP 

!---------------------------------------------------------------------------------------------------------------

    SUBROUTINE MATBLOCK(HAM,PARAMS,MPIDATA,X,AX,M,MB,NEV)

       IMPLICIT NONE 
       TYPE(HAMILTONIAN), INTENT(IN) :: HAM
       TYPE(PARAMETERS), INTENT(IN) :: PARAMS
       TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
       INTEGER, INTENT(IN)       ::  M, MB, NEV         
       COMPLEX*16, INTENT(IN)    ::  X(M, MB, NEV)
       COMPLEX*16, INTENT(OUT)   ::  AX(M, MB, NEV)
       INTEGER                   ::  ME, NP, IERROR
       INTEGER                   ::  J
       
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, IERROR) 
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, IERROR) 
       
       DO J = 1, NEV
          CALL HAMMULT(HAM,PARAMS,MPIDATA,X(:,:,J),AX(:,:,J),M*MB)
       END DO
       
    END SUBROUTINE MATBLOCK 
 
 !---------------------------------------------------------------------------------------------------------------

    SUBROUTINE PRECBLOCK(HAM,PARAMS,MPIDATA,GMSOLVER,X,PX,SIGMA,M,MB,NEV,MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)

       IMPLICIT NONE 
       TYPE(HAMILTONIAN), INTENT(IN) :: HAM
       TYPE(PARAMETERS), INTENT(IN) :: PARAMS
       TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
       TYPE(GMRES_SOLVER), INTENT(INOUT) :: GMSOLVER
       
       INTEGER, INTENT(IN)       ::  M, MB, NEV, MAX_INNER_ITER, OUTMODULUS
       COMPLEX *16, INTENT(IN)   ::  SIGMA
       COMPLEX*16, INTENT(IN)    ::  X(M*MB, NEV)
       COMPLEX*16, INTENT(OUT)   ::  PX(M*MB, NEV) ! PRECONDITIONED BLOCK X
       INTEGER                   ::  ME, NP, IERROR
       INTEGER                   ::  PCFLAG, MAX_RST
       INTEGER                   ::  J
       REAL*8                    ::  TOL, REL_RES ! TOLERANCE, GMRES RELATIVE RESIDUAL NORM UPON COMPLETION, 
       INTEGER                   ::  GMRES_ITER   ! NUMBER OF GMRES ITERATIONS TO CONVERGE
       INTEGER                   :: IDIST, ISEED(4)
       EXTERNAL                  :: PRECONDITIONER

!       print *, "USEPETSCFLAG IS", gmsolver%usepetscflag;stop

       IF (GMSOLVER%MAX_INNER_ITER == 0) THEN
          print *, "What? just don't call me.  right?"
          stop
       else IF (GMSOLVER%MAX_INNER_ITER == 1) THEN
          call matblock(HAM,PARAMS,MPIDATA,X,PX,M,MB,NEV)
          return
       endif

       if (gmsolver%usepetscflag==0) then
          call PRECBLOCK_PETSC(HAM,PARAMS,MPIDATA,GMSOLVER,X,PX,SIGMA,M,MB,NEV,MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
       else
          call PRECBLOCK_BASIC(HAM,PARAMS,MPIDATA,GMSOLVER,X,PX,SIGMA,M,MB,NEV,MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
       endif

    end SUBROUTINE PRECBLOCK


 !---------------------------------------------------------------------------------------------------------------
 
    SUBROUTINE PRECBLOCK_PETSC(HAM,PARAMS,MPIDATA,GMSOLVER,X,PX,SIGMA,M,MB,NEV,MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)

       IMPLICIT NONE 
       TYPE(HAMILTONIAN), INTENT(IN) :: HAM
       TYPE(PARAMETERS), INTENT(IN) :: PARAMS
       TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
       TYPE(GMRES_SOLVER), INTENT(INOUT) :: GMSOLVER
       
       INTEGER, INTENT(IN)       ::  M, MB, NEV, MAX_INNER_ITER, OUTMODULUS
       COMPLEX *16, INTENT(IN)   ::  SIGMA
       COMPLEX*16, INTENT(IN)    ::  X(M*MB, NEV)
       COMPLEX*16, INTENT(OUT)   ::  PX(M*MB, NEV) ! PRECONDITIONED BLOCK X
       INTEGER                   ::  ME, NP, IERROR
       INTEGER                   ::  PCFLAG, MAX_RST
       INTEGER                   ::  J
       REAL*8                    ::  TOL, REL_RES ! TOLERANCE, GMRES RELATIVE RESIDUAL NORM UPON COMPLETION, 
       INTEGER                   ::  GMRES_ITER   ! NUMBER OF GMRES ITERATIONS TO CONVERGE
       INTEGER                   :: IDIST, ISEED(4)
       EXTERNAL                  :: PRECONDITIONER

       IF (GMSOLVER%MAX_INNER_ITER == 0) THEN
          print *, "What? just don't call me.  right?"
          stop
       ENDIF

#if defined(USE_PETSC)
       DO J = 1, NEV
          CALL PETSCSOLVE(X(:,J),PX(:,J),M*MB,HAM,PARAMS,MPIDATA,GMSOLVER,OUTMODULUS)
       END DO
#else

       PX(:,:)=0d0
       print *, "ACK! DON't CALL ME WITHOUT PETSC"; stop

#endif


!!$ OLD LOGIC DJH DOES NOT UNDERSTAND.
!!$
!!$       !IF (MAX_INNER_ITER < 1) THEN
!!$          DO J = 1, NEV
!!$#if defined(USE_PETSC)
!!$             CALL PETSCSOLVE(X(:,J),PX(:,J),M*MB,HAM,PARAMS,MPIDATA,GMSOLVER,OUTMODULUS)
!!$#else
!!$
!!$             CALL PCWRAPPER(HAM,PARAMS,MPIDATA,SIGMA,M,MB,X(:,J),PX(:,J),PRECONDITIONER)
!!$#endif
!!$             !CALL PCSOLVE(HAM,PARAMS,MPIDATA,X(:,J),PX(:,J),SIGMA,M*MB)
!!$          END DO
!!$       !ELSE
!!$          !DO J = 1, NEV
!!$          !   CALL MYGMRES(X(:,J),PX(:,J),M*MB,HAM,PARAMS,MPIDATA,GMSOLVER)
!!$          !END DO
!!$       !END IF

       
    END SUBROUTINE PRECBLOCK_PETSC
 

 !-------------------------------------------------------------------------------------------------------------

!!$ SUBROUTINE PCSOLVE(HAM,PARAMS,MPIDATA,X,PX,SIGMA,MM)
!!$ IMPLICIT NONE
!!$ 
!!$ TYPE(HAMILTONIAN), INTENT(IN) :: HAM
!!$ TYPE(PARAMETERS), INTENT(IN) :: PARAMS
!!$ TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
!!$ INTEGER, INTENT(IN) :: MM
!!$ COMPLEX *16, INTENT(IN) :: SIGMA
!!$ COMPLEX *16, TARGET, INTENT(IN) :: X(1:MM)
!!$ COMPLEX *16, TARGET, INTENT(OUT) :: PX(1:MM)
!!$ 
!!$ REAL *8, ALLOCATABLE, SAVE :: W(:,:), WT(:,:)
!!$ COMPLEX *16, ALLOCATABLE, SAVE :: TEMP(:), ALLEIGS(:)
!!$ COMPLEX *16, ALLOCATABLE, SAVE :: U(:,:), V(:,:), U1(:)
!!$ COMPLEX *16, POINTER :: XP(:,:), PXP(:,:)
!!$ LOGICAL, SAVE :: FIRSTTIME = .TRUE.
!!$ INTEGER :: NN, M, MB, LWORK, NELEC, J
!!$ 
!!$ NELEC = PARAMS % NELEC
!!$ NN = PARAMS % NN
!!$ M = PARAMS % M
!!$ MB = MPIDATA % MBLK
!!$ LWORK = 10*NN
!!$ 
!!$ IF(FIRSTTIME) THEN
!!$ 
!!$    FIRSTTIME = .FALSE.
!!$    ALLOCATE(ALLEIGS(1:M),TEMP(1:MM))
!!$    ALLOCATE(W(1:NN,1:NN),WT(1:NN,1:NN))
!!$    CALL KINETIC_EIGS(HAM % T,ALLEIGS,W,WT,PARAMS % THETA,NN,LWORK,1,M)
!!$ 
!!$    IF(NELEC == 2) ALLOCATE(U(1:M,1:MB),V(1:M,1:MB),U1(1:M))
!!$ 
!!$ END IF
!!$ 
!!$ 
!!$ IF(NELEC == 1) THEN
!!$ 
!!$    CALL EIGMULT(X,TEMP,WT,NN)
!!$    TEMP = TEMP/(ALLEIGS - SIGMA)
!!$    CALL EIGMULT(TEMP,PX,W,NN)
!!$ 
!!$ ELSEIF(NELEC == 2) THEN
!!$ 
!!$    XP(1:M,1:MB) => X
!!$    PXP(1:M,1:MB) => PX
!!$ 
!!$    DO J = 1,MB
!!$       CALL EIGMULT(XP(:,J),U1,WT,NN)
!!$       U1 = U1/(ALLEIGS - SIGMA)
!!$       CALL EIGMULT(U1,PXP(:,J),W,NN)
!!$    END DO
!!$ 
!!$    CALL MPITRANSPOSE(XP,U,M,MB,MPIDATA)
!!$ 
!!$    DO J = 1,MB
!!$       CALL EIGMULT(U(:,J),U1,WT,NN)
!!$       U1 = U1/(ALLEIGS - SIGMA)
!!$       CALL EIGMULT(U1,V(:,J),W,NN)
!!$    END DO
!!$ 
!!$    CALL MPITRANSPOSE(V,U,M,MB,MPIDATA)
!!$ 
!!$    PXP = PXP + U
!!$ 
!!$ END IF
!!$ 
!!$ END SUBROUTINE PCSOLVE


!---------------------------------------------------------------------------------------------------------------

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

!-----------------------------------------------------------------------------------------------------------------

END SUBROUTINE GPLMR_EIG

#endif

END MODULE GPLMR



#ifndef PGFFLAG

SUBROUTINE PRECBLOCK_BASIC(HAM,PARAMS,MPIDATA,GMSOLVER,X,PX,SIGMA,M,MB,NEV,MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
  use hr_basicsolvemod
  USE HAMMOD
  USE STRUCTS
  USE MPI  

  IMPLICIT NONE 
  TYPE(HAMILTONIAN), INTENT(IN) :: HAM
  TYPE(PARAMETERS), INTENT(IN) :: PARAMS
  TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
  TYPE(GMRES_SOLVER), INTENT(INOUT) :: GMSOLVER
  
  INTEGER, INTENT(IN)       ::  M, MB, NEV, MAX_INNER_ITER, OUTMODULUS
  COMPLEX *16, INTENT(IN)   ::  SIGMA
  COMPLEX*16, INTENT(IN)    ::  X(M*MB, NEV)
  COMPLEX*16, INTENT(OUT)   ::  PX(M*MB, NEV) ! PRECONDITIONED BLOCK X
  INTEGER                   ::  ME, NP
  INTEGER                   ::  PCFLAG, MAX_RST
  INTEGER                   ::  J
  REAL*8                    ::  TOL, REL_RES ! TOLERANCE, GMRES RELATIVE RESIDUAL NORM UPON COMPLETION, 
  INTEGER                   ::  GMRES_ITER   ! NUMBER OF GMRES ITERATIONS TO CONVERGE
  INTEGER                   :: IDIST, ISEED(4)
  EXTERNAL                  :: PRECONDITIONER
  integer :: printflag
  
  printflag=1
  call basicblocklansolve(nev,M*MB, max_inner_iter, M**2, x,px,printflag,&
       allhmultsub,allhmultsub,allhdots0,sigma)
contains

  subroutine allhmultsub(in,out,howmany,checksize)
    integer, intent(in) :: howmany,checksize
    integer :: i
    complex*16 :: in(checksize,howmany),out(checksize,howmany)
    
!! call from petscsolve       call hammult(ham_gl,params_gl,mpidata_gl,xx,yy,locn)

    if (checksize.ne.M*MB) then
       print *, "CHECKSIZE!!!", checksize,M*MB; stop
    endif
    call hammult(ham,params,mpidata,in,out,M*MB)
    
  end subroutine allhmultsub
  
  function qhdot(in,out,n)
    implicit none
    integer :: n,i
    complex*16 :: in(n),out(n),qhdot,csum
    csum=0
    do i=1,n
       csum=csum+ CONJG(in(i))*out(i)
    enddo
    qhdot=csum
  end function qhdot
  
  subroutine allhdots0(bravectors,ketvectors,n,num1,num2,outdots)
    use mpi
    implicit none
    integer :: id,jd,num1,num2,n,IERROR,qq
    complex*16 :: bravectors(n,num1), ketvectors(n,num2), outdots(num1,num2), tempdots(num1,num2) !, checkdots(num1,num2)

    if (n.ne.M*MB) then
       print *, "CHECKSIZE!!!", n,M*MB; stop
    endif

    do id=1,num1
       do jd=1,num2
          outdots(id,jd)= qhdot(bravectors(:,id),ketvectors(:,jd),n)

!! check
!          checkdots(id,jd)= mpihermdot(bravectors(:,id),ketvectors(:,jd),n,mpidata)

       enddo
    enddo

    call mpi_allreduce(outdots,tempdots,num1*num2,MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERROR)
    if (IERROR.ne.0) then
       print *, "IERROR REDUCE PRECBLOCK_BASIC ",ierror;stop
    endif
    outdots(:,:)=tempdots(:,:)

!    checkdots(:,:)=checkdots(:,:)-outdots(:,:)
!    do id=1,num1
!       do jd=1,num2
!          if (abs(checkdots(id,jd)).gt.1d-10) then
!             print *, "CHECKDOT ERROR"
!             do qq=1,num2
!                write(*,'(100F10.5)') checkdots(:,qq)
!             enddo
!             stop
!          endif
!       enddo
!    enddo

  end subroutine allhdots0
  
end SUBROUTINE PRECBLOCK_BASIC


#else


 NOTDONE... ASSIGNPOINTER

module submod
  use structs

  TYPE(MPI_DATA),pointer :: MPIDATA

  integer :: M,MB

contains

  subroutine allhmultsub(in,out,howmany,checksize)
    integer, intent(in) :: howmany,checksize
    integer :: i
    complex*16 :: in(checksize,howmany),out(checksize,howmany)
    
!! call from petscsolve       call hammult(ham_gl,params_gl,mpidata_gl,xx,yy,locn)

    if (checksize.ne.M*MB) then
       print *, "CHECKSIZE!!!", checksize,M*MB; stop
    endif
    call hammult(ham,params,mpidata,in,out,M*MB)
    
  end subroutine allhmultsub
  
  function qhdot(in,out,n)
    implicit none
    integer :: n,i
    complex*16 :: in(n),out(n),qhdot,csum
    csum=0
    do i=1,n
       csum=csum+ CONJG(in(i))*out(i)
    enddo
    qhdot=csum
  end function qhdot
  
  subroutine allhdots0(bravectors,ketvectors,n,num1,num2,outdots)
    use mpi
    implicit none
    integer :: id,jd,num1,num2,n,IERROR,qq
    complex*16 :: bravectors(n,num1), ketvectors(n,num2), outdots(num1,num2), tempdots(num1,num2), checkdots(num1,num2), &
         mpihermdot

    if (n.ne.M*MB) then
       print *, "CHECKSIZE!!!", n,M*MB; stop
    endif

    do id=1,num1
       do jd=1,num2
          outdots(id,jd)= qhdot(bravectors(:,id),ketvectors(:,jd),n)

!! check
          checkdots(id,jd)= mpihermdot(bravectors(:,id),ketvectors(:,jd),n,mpidata)

       enddo
    enddo

    call mpi_allreduce(outdots,tempdots,num1*num2,MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERROR)
    if (IERROR.ne.0) then
       print *, "IERROR REDUCE PRECBLOCK_BASIC ",ierror;stop
    endif
    outdots(:,:)=tempdots(:,:)

    checkdots(:,:)=checkdots(:,:)-outdots(:,:)

    do id=1,num1
       do jd=1,num2
          if (abs(checkdots(id,jd)).gt.1d-10) then
             print *, "CHECKDOT ERROR"
             do qq=1,num2
                write(*,'(100F10.5)') checkdots(:,qq)
             enddo
             stop
          endif
       enddo
    enddo

  end subroutine allhdots0

end module submod



SUBROUTINE PRECBLOCK_BASIC(HAM,PARAMS,MPIDATA,GMSOLVER,X,PX,SIGMA,M,MB,NEV,MAX_INNER_ITER,PRECONDITIONER,OUTMODULUS)
  use hr_basicsolvemod
  use submod
  USE HAMMOD
  USE STRUCTS
  USE MPI  

  IMPLICIT NONE 
  TYPE(HAMILTONIAN), INTENT(IN) :: HAM
  TYPE(PARAMETERS), INTENT(IN) :: PARAMS
  TYPE(MPI_DATA), INTENT(IN) :: MPIDATA
  TYPE(GMRES_SOLVER), INTENT(INOUT) :: GMSOLVER
  
  INTEGER, INTENT(IN)       ::  M, MB, NEV, MAX_INNER_ITER, OUTMODULUS
  COMPLEX *16, INTENT(IN)   ::  SIGMA
  COMPLEX*16, INTENT(IN)    ::  X(M*MB, NEV)
  COMPLEX*16, INTENT(OUT)   ::  PX(M*MB, NEV) ! PRECONDITIONED BLOCK X
  INTEGER                   ::  ME, NP
  INTEGER                   ::  PCFLAG, MAX_RST
  INTEGER                   ::  J
  REAL*8                    ::  TOL, REL_RES ! TOLERANCE, GMRES RELATIVE RESIDUAL NORM UPON COMPLETION, 
  INTEGER                   ::  GMRES_ITER   ! NUMBER OF GMRES ITERATIONS TO CONVERGE
  INTEGER                   :: IDIST, ISEED(4)
  EXTERNAL                  :: PRECONDITIONER
  integer :: printflag
  
  printflag=1
  call basicblocklansolve(nev,M*MB, max_inner_iter, M**2, x,px,printflag,&
       allhmultsub,allhmultsub,allhdots0,sigma)

end SUBROUTINE PRECBLOCK_BASIC


#endif
  
  
