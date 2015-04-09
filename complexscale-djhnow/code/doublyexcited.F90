#include "definitions.INC"

! --------------------------------------------------------------------------------------
! ACTION == 5 ==> compute doubly-excited states of helium
! --------------------------------------------------------------------------------------
  
! EXCITED STATE EXACT ANSWERS
!------------------------------
! 2S2P SINGLET: (-0.693135,-6.865E-4)
! 2S2P TRIPLET: (-0.760492,-1.495E-4)
! 2P2P SINGLET: (-0.701946,-1.200E-3)
! 2P2P TRIPLET: -0.710500


MODULE MULTMOD
  USE STRUCTS

  TYPE(PARAMETERS) :: PARAMS
  TYPE(HAMILTONIAN) :: HAM
  TYPE(MPI_DATA) :: MPIDATA

END MODULE

SUBROUTINE MYMATMUL(IN,OUT,HOWMANY,CHECKSIZE)
  USE MULTMOD
  USE HAMMOD
  IMPLICIT NONE
  INTEGER :: HOWMANY,CHECKSIZE,II,THISSIZE
  NUMBER :: IN(CHECKSIZE,HOWMANY),OUT(CHECKSIZE,HOWMANY)

  THISSIZE = PARAMS%M * MPIDATA%MBLK

  IF (CHECKSIZE.NE.THISSIZE) THEN
     PRINT *, "SIZE ERRORRRR 66", CHECKSIZE,THISSIZE
     STOP
  ENDIF

  DO II=1,HOWMANY
    CALL HAMMULT(HAM,PARAMS,MPIDATA,IN(:,II),OUT(:,II),THISSIZE)
  ENDDO

!!  CALL MATBLOCK(HAM,PARAMS,MPIDATA,IN,OUT,PARAMS%M,MPIDATA%MBLK,HOWMANY)

END SUBROUTINE MYMATMUL


SUBROUTINE MYDOTSUB(BRAS,KETS,SIZE,NN,MM,OUTDOTS)
  USE MPI
  IMPLICIT NONE
  INTEGER :: SIZE,NN,MM,II,JJ,IERR
  NUMBER :: BRAS(SIZE,NN), KETS(SIZE,MM), OUTDOTS(NN,MM),TEMPDOTS(NN,MM)
  DO JJ=1,MM
    DO II=1,NN
      TEMPDOTS(II,JJ)=DOT_PRODUCT(BRAS(:,II),KETS(:,JJ))
    ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(TEMPDOTS,OUTDOTS,MM*NN,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)

END SUBROUTINE MYDOTSUB

PROGRAM ATOM
  USE MULTMOD
  USE STRUCTS
  USE INITMOD
  USE INPUT
  USE LANMOD
  USE SOLVER
  USE MPI
  USE GPLMR
  USE UTILMOD
#if defined(USE_PETSC)
  USE PETSCWRAP
#endif
  IMPLICIT NONE
  
  ! --------------------------------------------------------------------------------------
  ! Declare variables and structures
  
  INTEGER :: I, J, ME, NP, MPIERR, M, MB, FIRST, LAST, TOTGET,MAXGET,&
    THISSIZE,TOTSIZE
  NUMBER, ALLOCATABLE :: PSI(:,:,:),E0(:), EIGS(:),EXPECTS(:)
  NUMBER, ALLOCATABLE :: myEV(:,:)
  REAL *8 :: TIME, ANG_EXP
!!  TYPE(PARAMETERS) :: PARAMS
!!  TYPE(MPI_DATA) :: MPIDATA
  TYPE(LANPARAMS), TARGET :: LANPAR
  TYPE(MOLECULE) :: MOL
!!  TYPE(HAMILTONIAN) :: HAM
  TYPE(GMRES_SOLVER) :: GMSOLVER
  character(len=200) :: inpfile="Input.Inp"
  EXTERNAL :: MYMATMUL,MYDOTSUB

  
#if ISREAL
  IF(ME == 0) WRITE(*,"(A)") "ERROR: COMPLEX datatype needed for doubly-excited helium calculation."
  stop
#endif
  call system("mkdir Bin")
  
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
  
! Initialize variables (user input)
  call get_inputfile(inpfile);  call GET_INPUT(inpfile,ME);  call get_args(me)

!! action=5 defaults from atom.f90
! Initialize structures based on input parameters
  
  CALL INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM,.FALSE.,2)
  

  GMSOLVER%IUNIT=GM_PRINT_OUTPUT
  GMSOLVER%MAX_INNER_ITER=GM_MAX_INNER_ITER
  GMSOLVER%PCOPT=GM_PCOPT
  GMSOLVER%TOL=GM_TOL
  GMSOLVER%usepetscflag=gm_usepetscflag

! ... NELEC and SAVE_VECTORS are not input variables apparently because they
! are changed at different steps of the calculation.  So this whole routine
! must be called again. argh.

  CALL INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM,.TRUE.,1)
  M = PARAMS % M
  MB = MPIDATA % MBLK
  FIRST = MPIDATA % FIRST
  LAST = MPIDATA % LAST

  IF(.not.READVECS) THEN  
     CALL DIAGONALIZE(PARAMS,MPIDATA,LANPAR,HAM,LANBLOCKSIZE)
  endif

  MAXGET=GPLMR_NEV; TOTGET=GPLMR_NEV+NUMSKIP
  IF(.not.READVECS) THEN  
     MAXGET=NUMORBS*(NUMORBS-1)/2
     IF (SINGLET) THEN
        MAXGET=NUMORBS*(NUMORBS+1)/2
     ENDIF
     if (GPLMR_NEV+NUMSKIP.gt.MAXGET) THEN
        IF (ME==0) PRINT *, "FOR INITIAL GUESS, GPLRM_NEV TOO BIG FOR NUMORBS.",GPLMR_NEV,NUMSKIP,NUMORBS
        STOP
     ENDIF
     TOTGET=MIN(MAXGET,GPLMR_NEV+NUMSKIP+5)
  ENDIF

  IF(ME == 0) WRITE(*,*) "ALLOCATING PSI... ",totget


  ALLOCATE(PSI(1:M,FIRST:LAST,1:TOTGET),E0(1:TOTGET),EIGS(TOTGET),EXPECTS(1:TOTGET));  PSI(:,:,:)=0d0
     
  IF(ME == 0) WRITE(*,*) "     ....OK ALLOC"

! Compute approximate eigenvectors for helium based
! on the eigenvectors for He+

  if (.not.READVECS ) then
     ALLOCATE(myEV(1:M,NUMORBS))
     IF(ME == 0) WRITE(*,"(A)") "    ...ORBITALGATHER..."
     CALL HELIUM_ORBITALGATHER(numorbs,LANPAR%EIGENVECTORS(:,:),myEV(:,:),M,FIRST,LAST,MPIDATA,ME,hermopt)
  endif

  CALL INIT_STRUCTS(ME,NP,PARAMS,MPIDATA,LANPAR,MOL,HAM,.FALSE.,2)

  IF(READVECS) THEN
     CALL READ_EIGENVECTORS(PSI,TOTGET,MPIDATA,PARAMS%M,invecfile)
  ELSE

! Compute expectations of initial eigenvector guesses to use
! for the shift parameter in GPLMR

     IF(ME == 0)        WRITE(*,"(A)") "COMPUTING INITIAL GUESSES..."
     CALL DIAGSUBSPACE(HAM,PARAMS,MPIDATA,M,MB,TOTGET,E0,PSI,me,gplmr_shift,&
          hermopt,numorbs,singlet,myEV(:,:),FIRST,LAST,np)
     deallocate(myEV)
  ENDIF

  if (NUMSKIP.GT.0) THEN
     PSI(:,:,1:TOTGET-NUMSKIP)=PSI(:,:,1+NUMSKIP:TOTGET)
     IF (.NOT.READVECS) E0(1:TOTGET-NUMSKIP)=E0(1+NUMSKIP:TOTGET)
  ENDIF

  CALL GET_EXPECTS(HAM,PARAMS,MPIDATA,M,MB,TOTGET-NUMSKIP,PSI,EXPECTS,cshiftopt)

  IF(ME == 0.and..not.READVECS) THEN
     WRITE(*,"(A)") "SUBSPACE EIGENVALS:"
     WRITE(*,'(2F20.14,T50,F20.14)')(E0(i),ABS(E0(i)-EXPECTs(i)),i=1,GPLMR_NEV)
     print *, " --xx-- "
     WRITE(*,'(2F20.14,T50,F20.14)')(E0(i),ABS(E0(i)-EXPECTs(i)),i=1+GPLMR_NEV,TOTGET-NUMSKIP)
  END IF

  IF(ME == 0) THEN
     WRITE(*,"(A)") "SUBSPACE EXPECTATION VALS:"
     WRITE(*,'(2F20.14)') EXPECTS(1:GPLMR_NEV)
     print *, " --xx-- "
     WRITE(*,'(2F20.14)') EXPECTS(1+GPLMR_NEV:TOTGET-NUMSKIP)
  END IF
     
  
! set shift to equal average... no set it to equal first
!  I assume it should be set using a quantity computed with the
!  hermitian norm!!!  Now it is.  EXPECTS not E0 with MPIHERMDOT
!  defining EXPECTS.  Works better now.

!if (abs(999+GPLMR_SHIFT).eq.0d0) then

   GPLMR_SHIFT = SUM(EXPECTS(1:GPLMR_NEV))/GPLMR_NEV

!  GPLMR_SHIFT=EXPECTS(1)
!endif

   IF(TAKEREALPART) GPLMR_SHIFT = REAL(GPLMR_SHIFT,8)

   IF(ME == 0)  WRITE(*,"(A,F20.10,F20.10)") "CHOSEN SHIFT: ",DREAL(GPLMR_SHIFT),DIMAG(GPLMR_SHIFT)


! -------------------------------------------- !
! -------- CALCULATE STATE WITH GPLMR -------- !
! -------------------------------------------- !

  THISSIZE=PARAMS % M * MPIDATA % MBLK
  TOTSIZE=PARAMS % M * PARAMS % M 

!  not debugged - djh.  just for personal learning purposes.
!  IF (ME==0) PRINT *, "CALLING HARMONIC RITZ."   
!  CALL HARMONIC_RITZ(GPLMR_NEV,1,THISSIZE,6,TOTSIZE,PSI,EIGS,&
!    1,1,1d-3,  MYMATMUL,  MYDOTSUB,  GPLMR_SHIFT, 100d0, &
!    6,6,6)
!  IF (ME==0) PRINT *, "   ...CALLED HARMONIC RITZ."

  CALL GPLMR_EIG(HAM,PARAMS,MPIDATA,GMSOLVER,PSI,EIGS, &
               PARAMS % M,MPIDATA % MBLK,GPLMR_SHIFT, &
               GPLMR_NEV,GPLMR_NS,GPLMR_MAXIT,GPLMR_TOL,MSOLVECPLX2,&
               OUTMODULUS)

  if (ME==0) then
     print *, "FINAL EIGVALS:"
     write(*,'(2F20.12)') EIGS(1:GPLMR_NEV);    print *
  endif


  IF(WRITEVECS) THEN
     if (ME==0) print *, "WRITING eigvects..."
     CALL WRITE_EIGENVECTORS(PSI,GPLMR_NEV,MPIDATA,PARAMS%M,outvecfile)
     if (ME==0) print *, "   ..written."
  END IF

  DEALLOCATE(PSI,E0,EIGS,EXPECTS)

! Finalize MPI and print elapsed time of program

  TIME = MPI_WTIME() - TIME
  IF(ME == 0) WRITE(*,"(A,F10.2,A)") "TOTAL TIME ELAPSED: ",TIME," SECONDS."
  CALL MPI_FINALIZE(MPIERR)
  
  
END PROGRAM ATOM

