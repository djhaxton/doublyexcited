
!! getting lapack segfaults gnu edison, sometimes.

subroutine invsqrtmat(A,N)
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=exp(A_in)
  implicit none
  integer :: N,lwork,i,k,j
  complex*16 :: eig(N)  ,CVL(N,N),CVR(N,N),CA(N,N)
  complex*16 :: A(N,N)
  complex*16 :: VL(N,N),VR(N,N),work(80*N),rwork(40*N)

  lwork=80*N

!print *, "AAAHA";

  j=0
  call zgeev('V','V',N,A,N,eig,VL,N,VR,N,work,lwork,rwork,i)


  VL(:,:)=TRANSPOSE(VR(:,:))

!print *, "LAAAddd";
  call invmatsmooth(VL,N,N,0d0)
!print *, "OOOD";

!! apply the function
  do k=1,N
     eig(k) = 1d0/sqrt( eig(k) )
     CVR(:,k)=VR(:,k)*eig(k)
  enddo
  CVL(:,:)=VL(:,:)

!! rebuild the matrix

  call ZGEMM('N','T',N,N,N,(1d0,0d0),CVR,N,CVL,N,(0d0,0d0),CA,N)
  A(:,:)=CA(:,:)  !! OK IMP CONV

end subroutine invsqrtmat



subroutine invmatsmooth(A,N,LDA,tol)  !! inverse of ANY matrix.
!! input :
!! A - an N by N matrix
!! N - the dimension of A
!! output:
!! A - an N by N matrix that is A=(A_in)**-1
  implicit none
  integer :: N,lwork,i,j,k,LDA
  real*8 :: SV(N),tol
  complex*16 :: A(LDA,N), SAVEA(LDA,N)
  complex*16 :: U(N,N),VT(N,N),work(5*N),zwork(5*N)
  lwork=5*N

  SAVEA=A

!! do the svd

  call zgesvd('A','A',N,N,A,LDA,SV,U,N,VT,N,zwork,lwork,work,i)

  if (i.ne.0) then
     print *, "ERR SVD",i;
     do i=1,N
        print *, SAVEA(:,i)
     enddo
     stop
  endif
!! apply the function
  do k=1,N
     if (SV(k).ne.0d0) then
        if(SV(k).lt.tol*SV(1)) then    !! it is positive ,whatever
           SV(k)= 1d0 / tol / SV(1)
        else
           SV(k) = 1d0 / SV(k)
        endif
     endif

  enddo
!! rebuild the matrix
  do j=1,N
    do i=1,N
      A(i,j) = 0d0
      do k=1,N
         A(i,j) = A(i,j) + CONJG(VT(k,i)) * SV(k) * CONJG(U(j,k))
       enddo
    enddo
  enddo

end subroutine invmatsmooth




subroutine symorthogmat(A,N)  
  implicit none
  integer :: N,lwork,i,j,k
  real*8, allocatable :: SV(:)
  complex*16 :: A(N,N),Asave(N,N)
  complex*16, allocatable :: U(:,:),VT(:,:),work(:),zwork(:)
  real*8, allocatable :: rwork(:)
  do i=1,N;  do j=1,i
     if (abs((CONJG(A(i,j)))-A(j,i)).gt.1.d-7) then
        print *, "SYM ERR ORTHOG", abs((A(i,j))-A(j,i));           stop
     endif
  enddo;  enddo
  

  Asave(:,:)=A(:,:)

!! do the svd
  lwork=10*N
  allocate(U(N,N),VT(N,N),SV(N),work(lwork),zwork(lwork),rwork(lwork))
  call zgesvd('A','A',N,N,A,N,SV,U,N,VT,N,zwork,lwork,work,i)
  if (i.ne.0) then
     print *, "ERR herm SVD",i,N
     do k=1,N
        writE(*,'(10000E9.2)') Asave(:,k)
     enddo
     stop
  endif
!! apply the function
  do k=1,N
    if(abs(SV(k)/SV(1)).lt.1d-10) then
      print *, "SING symorthogmat!", sv(k)
      SV(k) = 1d0 / sqrt(SV(k))     
!      SV(k) = 0d0
    else
      SV(k) = 1d0 / sqrt(SV(k))     
    endif
  enddo
!! rebuild the matrix
  do j=1,N
     do i=1,N
        A(i,j) = 0d0
        do k=1,N

           A(i,j) = A(i,j) + U(i,k) * SV(k) * VT(k,j)                  !HC OF INVERSE

        enddo
     enddo
  enddo
  deallocate(U,VT,SV,work,zwork,rwork)

end subroutine symorthogmat


  




