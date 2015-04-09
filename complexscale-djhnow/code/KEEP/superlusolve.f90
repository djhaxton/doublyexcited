#if defined(USE_SUPERLU)

SUBROUTINE SUPERLUSOLVE(II,JJ,TIJ,NNZ,B,NM,GRID,OPTIONS,SCALEPERMSTRUCT,LUSTRUCT,SOLVESTRUCT,A,STAT,DOFROMSCRATCH)
! Solve the system A x =b using SuperLU_DIST
! A is a sparse matrix defined by triplets
! II, JJ, TIV (NNZ nonzero elements)
! B is the right-hand side (size LOCN)
  USE SUPERLU_MOD
  INCLUDE 'mpif.h'

  INTEGER, INTENT(IN) :: NNZ, NM
  INTEGER, INTENT(IN) :: II(NNZ), JJ(NNZ)
  COMPLEX *16, INTENT(IN) :: TIJ(NNZ)
  COMPLEX *16, INTENT(INOUT) :: B(NM)

  INTEGER                 :: ME, NP, NPROW, NPCOL, MPIERR, INFO, I, J, IND, COUNT
  INTEGER, ALLOCATABLE    :: INDEX(:), ROWS(:), COLS(:), JPTR(:) 
  COMPLEX*16, ALLOCATABLE :: VAL(:), TMP(:)
  INTEGER(SUPERLU_PTR)    :: GRID
  INTEGER(SUPERLU_PTR)    :: OPTIONS
  INTEGER(SUPERLU_PTR)    :: SCALEPERMSTRUCT
  INTEGER(SUPERLU_PTR)    :: LUSTRUCT
  INTEGER(SUPERLU_PTR)    :: SOLVESTRUCT
  INTEGER(SUPERLU_PTR)    :: A
  INTEGER(SUPERLU_PTR)    :: STAT 
  DOUBLE PRECISION        :: BERR 
  LOGICAL                 :: DOFROMSCRATCH

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,MPIERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,MPIERR)

  IF(DOFROMSCRATCH) THEN  
    CALL F_CREATE_GRIDINFO_HANDLE(GRID)
    CALL F_CREATE_OPTIONS_HANDLE(OPTIONS)
    CALL F_CREATE_SCALEPERM_HANDLE(SCALEPERMSTRUCT)
    CALL F_CREATE_LUSTRUCT_HANDLE(LUSTRUCT)
    CALL F_CREATE_SOLVESTRUCT_HANDLE(SOLVESTRUCT)
    CALL F_CREATE_SUPERMATRIX_HANDLE(A)
    CALL F_CREATE_SUPERLUSTAT_HANDLE(STAT)

    ! Initialize the 2D grid for SuperLU_DIST
    ! For now, we assume that each process as at least one
    ! instance of SuperLU and that no instance of SuperLU
    ! is shared
    NPROW = 1 !FLOOR(SQRT(REAL(NP)))
    NPCOL = 1 !NP/NPROW
    !CALL F_SUPERLU_GRIDINIT(MPI_COMM_WORLD, NPROW, NPCOL, GRID)
    CALL F_SUPERLU_GRIDMAP(MPI_COMM_WORLD, NPROW, NPCOL, ME, 1,  GRID)

    ! Some processes might be left out
    CALL GET_GRIDINFO(GRID, ME, NPROW, NPCOL)
    IF(ME >= NPROW * NPCOL) THEN 
      GOTO 100
    END IF

    ! Convert triplet to compressed column
    !  1/ Sort column indices and apply the sorting permutation
    !     to row indices and values
    ALLOCATE(COLS(NNZ),INDEX(NNZ))
    COLS = JJ
    DO I=1,NNZ
      INDEX(I)=I
    END DO
    CALL QSORT(INDEX,COLS,1,NNZ,NNZ)
    ALLOCATE(ROWS(NNZ),VAL(NNZ),TMP(NNZ))
    DO I=1,NNZ
      ROWS(I)=II(INDEX(I))
      VAL(I)=TIJ(INDEX(I))
    END DO
    !  2/ Compute JPTR. NB: this piece of code assumes that there
    !  are no empty rows/columns.
    ALLOCATE(JPTR(NM+1))
    JPTR(1)=1
    IND=2
    COUNT=0
    DO I=2,NNZ
      COUNT = COUNT + 1
      IF(COLS(I)/=COLS(I-1)) THEN
        JPTR(IND)= JPTR(IND-1) + COUNT
        IND = IND + 1
        COUNT = 0
      END IF 
    END DO
    JPTR(NM+1)=NNZ+1 
    DEALLOCATE(COLS)
    !  3/ Sort rows within each column
    DO I=1,NNZ
      INDEX(I)=I
    END DO
    DO J=1,NM
      CALL QSORT(INDEX,ROWS,JPTR(J),JPTR(J+1)-1,NNZ)
      DO I=JPTR(J),JPTR(J+1)-1
        TMP(I)=VAL(INDEX(I))
      END DO
    END DO
    VAL=TMP
    DEALLOCATE(INDEX,TMP)

    !  4/ 0-based indices for C-based SuperLU_DIST
    JPTR = JPTR - 1
    ROWS = ROWS - 1

    ! Distribute the matrix to the process gird
    CALL  F_ZCREATE_DIST_MATRIX(A, NM, NM, NNZ, VAL, ROWS, JPTR, GRID)

    ! Initialize stuff 
    CALL F_SET_DEFAULT_OPTIONS(OPTIONS)
    CALL SET_SUPERLU_OPTIONS(OPTIONS,PrintStat=NO)
    CALL GET_SUPERMATRIX(A, NROW=NM, NCOL=NM)
    CALL F_SCALEPERMSTRUCTINIT(NM, NM, SCALEPERMSTRUCT)
    CALL F_LUSTRUCTINIT(NM, NM, LUSTRUCT)
    CALL F_PSTATINIT(STAT)
  ELSE
    CALL SET_SUPERLU_OPTIONS(OPTIONS,Fact=FACTORED)
    CALL SET_SUPERLU_OPTIONS(OPTIONS,PrintStat=NO)
  END IF

  ! Call SuperLU_DIST
  CALL F_PZGSSVX(OPTIONS, A, SCALEPERMSTRUCT, B, NM, 1, &
                 GRID, LUSTRUCT, SOLVESTRUCT, BERR, STAT, INFO)

  !IF (INFO == 0) THEN
  !   WRITE(*,*)'Backward error: ', BERR
  !ELSE
  !   WRITE(*,*)'INFO from f_pdgssvx = ', INFO
  !ENDIF

  ! Deallocate the storage allocated by SuperLU_DIST
  !CALL F_PSTATFREE(STAT)
  !CALL F_DESTROY_COMPROWLOC_MAT_DIST(A)
  !CALL F_SCALEPERMSTRUCTFREE(SCALEPERMSTRUCT)
  !CALL F_DESTROY_LU(NM, GRID, LUSTRUCT)
  !CALL F_LUSTRUCTFREE(LUSTRUCT)
  !CALL F_ZSOLVEFINALIZE(OPTIONS, SOLVESTRUCT)

  ! Release the SuperLU process grid
100   CONTINUE
  !CALL F_SUPERLU_GRIDEXIT(GRID)

  ! Deallocate things
  !CALL F_DESTROY_GRIDINFO_HANDLE(GRID)
  !CALL F_DESTROY_OPTIONS_HANDLE(OPTIONS)
  !CALL F_DESTROY_SCALEPERM_HANDLE(SCALEPERMSTRUCT)
  !CALL F_DESTROY_LUSTRUCT_HANDLE(LUSTRUCT)
  !CALL F_DESTROY_SOLVESTRUCT_HANDLE(SOLVESTRUCT)
  !CALL F_DESTROY_SUPERMATRIX_HANDLE(A)
  !CALL F_DESTROY_SUPERLUSTAT_HANDLE(STAT)

  IF(DOFROMSCRATCH) THEN 
    DEALLOCATE(ROWS,JPTR,VAL)
  END IF

END SUBROUTINE SUPERLUSOLVE

RECURSIVE SUBROUTINE QSORT(INDEX, VALUES, IBEG, IEND, N)
  ! Simple in-place Quick Sort
  ! VALUES is sorted in increasing numerical order, and
  ! INDEX is sorted accordingly so that it contains the
  ! sorting permutation (assuming INDEX=1:N)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: IBEG, IEND, N
  INTEGER, INTENT(INOUT) :: INDEX(N)
  ! Here VALUES is an array of integers but it could have any
  ! type (as long as ">" works)
  INTEGER, INTENT(INOUT)    :: VALUES(N)

  INTEGER :: I, POS, IPIV
  INTEGER :: INDPIV, INDTMP
  REAL    :: VALPIV, VALTMP

  IF(IEND>IBEG) THEN
    ! Choose a pivot
    IPIV   = IBEG + (IEND-IBEG)/2
    VALPIV = VALUES(IPIV)
    INDPIV = INDEX(IPIV)

    ! Put pivot at the end
    VALTMP       = VALUES(IPIV)
    INDTMP       = INDEX(IPIV)
    VALUES(IPIV) = VALUES(IEND)
    INDEX(IPIV)  = INDEX(IEND)
    VALUES(IEND) = VALTMP
    INDEX(IEND)  = INDTMP          

    ! Push values larger than VALPIV to the end 
    POS = IBEG
    DO I = IBEG, IEND - 1
      IF(VALUES(I) <= VALPIV) THEN
        VALTMP      = VALUES(I)
        INDTMP      = INDEX(I)
        VALUES(I)   = VALUES(POS)
        INDEX(I)    = INDEX(POS)
        VALUES(POS) = VALTMP
        INDEX(POS)  = INDTMP          
        POS         = POS + 1
      END IF
    END DO

    ! Put pivot to its final place
    VALTMP       = VALUES(POS)
    INDTMP       = INDEX(POS)
    VALUES(POS)  = VALUES(IEND)
    INDEX(POS)   = INDEX(IEND)
    VALUES(IEND) = VALTMP
    INDEX(IEND)  = INDTMP 

    ! Sort the two sublists
    CALL QSORT(INDEX, VALUES, IBEG, POS-1, N)
    CALL QSORT(INDEX, VALUES, POS+1, IEND, N)
  END IF

END SUBROUTINE QSORT

#endif
