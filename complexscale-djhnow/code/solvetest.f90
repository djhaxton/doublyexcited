
!! Calculates bound states of 1D and 2D harmonic oscillator using sinc DVR or finite difference.

module mymodule

!! ADJUSTABLE PARAMETERS !!

  integer, parameter :: &
       size=30, &          !! number of points (size of 1D problem)
       numvects=3        !! number of vectors to return. 

! complex*16 :: etarget=(5.85d0,-0.2d0)   !! TARGET ENERGY

 complex*16 :: etarget=0d0

  integer :: &
      order=10           !! krylov dimension
  real*8 :: length=20d0, & !! size of "box", determines spacing
       springconstant=1d0    !! spring constant.  harmonic oscillator.

!! THOSE WERE THE ADJUSTABLE PARAMETERS, ABOVE. 

  real*8 :: pi, spacing , const
  complex*16, allocatable :: potential(:),points(:),kematrix(:,:),hammatrix(:,:)

  complex*16, parameter :: eitheta=(12d0,1d0)/13d0

end module mymodule

!! matrix vector multiplication routine to be passed as argument
!!  to lanczos subroutine

module mysubroutines
contains

subroutine mymultsub(in,out,howmany,checksize)
  use mymodule
  integer, intent(in) :: howmany,checksize
  integer :: i
  complex*16 :: in(size,howmany),out(size,howmany)

  if (size.ne.checksize) then
     print *, "AAUGH!!",size,checksize
     stop
  endif

!  do i=1,howmany
!     out(:,i)=in(:,i)*potential(:)*eitheta**2
!  enddo
!  call zgemm('n','n',size,howmany,size,CONJG(eitheta**2),kematrix,size,in,size,(1d0,0d0),out,size)

  call zgemm('n','n',size,howmany,size,(1d0,0d0),hammatrix,size,in,size,(0d0,0d0),out,size)

end subroutine mymultsub



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
  implicit none
  integer :: id,jd,num1,num2,n
  complex*16 :: bravectors(n,num1), ketvectors(n,num2), outdots(num1,num2)
  do id=1,num1
     do jd=1,num2
        outdots(id,jd)= qhdot(bravectors(:,id),ketvectors(:,jd),n)
     enddo
  enddo
end subroutine allhdots0

end module


program myprogram
  use mymodule
  use mysubroutines
  use hr_basicsolvemod
  implicit none

  integer :: i
  complex*16 :: vectors(size,numvects),solvevectors(size,numvects),tempmatrix(size,size), &
       solvevectors2(size,numvects),multvectors(size,numvects),multvectors2(size,numvects)
  integer :: ipiv(size)
  real*8 :: realvectors(size,numvects)

  real*8 :: ketot(1-size:size-1),realpoints(size)
!! NOTE DIFFERENT (BETTER) DIMENSIONS HERE for ketot(1-size:size-1)
!!    versus versions of SineDVR.F90 elsewhere, ketot(-size:size)

  print *, "GO MYPROGRAM."
  call system("date")

  allocate(potential(size),points(size),kematrix(size,size),hammatrix(size,size))

  pi=4*atan(1d0)

  spacing = length/(size+1)


  call sinedvr(ketot,realpoints,size,spacing)
  points(:)=realpoints(:)
  potential(:)=0.5d0*springconstant*points(:)**2

  do  i=1,size
     kematrix(:,i)=ketot(1-i:size-i)
  enddo

  hammatrix(:,:)=kematrix(:,:)*CONJG(eitheta**2)
  do  i=1,size
     hammatrix(i,i)=hammatrix(i,i) + potential(i)*eitheta**2
  enddo


  print *, "CALL SOLVE.",etarget

  call system("date")

  call RANDOM_NUMBER(realvectors(:,:))
  vectors=realvectors

  call RANDOM_NUMBER(realvectors(:,:))
  vectors=vectors+realvectors**3*(0d0,1d0)

  print *, "GOSOLVE"
  call basicblocklansolve(numvects,size,order,size,vectors,solvevectors,1,mymultsub,mymultsub,allhdots0,etarget)

  tempmatrix(:,:)=hammatrix(:,:)
  do i=1,size
     tempmatrix(i,i)=tempmatrix(i,i)-etarget
  enddo

  solvevectors2(:,:)=vectors(:,:)
  call zgesv(size,numvects,tempmatrix,size,ipiv,solvevectors2,size,i)
  if (i.ne.0) then
     print *, "INFO ZGESV",i;stop
  endif

  tempmatrix(:,:)=hammatrix(:,:)
  do i=1,size
     tempmatrix(i,i)=tempmatrix(i,i)-etarget
  enddo
  
  call ZGEMM('N','N',size,numvects,size,(1d0,0d0),tempmatrix,size,solvevectors,size,(0d0,0d0),multvectors,size)
  call ZGEMM('N','N',size,numvects,size,(1d0,0d0),tempmatrix,size,solvevectors2,size,(0d0,0d0),multvectors2,size)



  print *, "OKVECTORS:"
  do i=1,size
!     write(*,'(100F15.8)') solvevectors(i,:)-solvevectors2(i,:)/ &
!          sqrt(abs(solvevectors(i,:)**2+solvevectors2(i,:)**2))

!     write(*,'(100F12.6)') multvectors(i,:),multvectors2(i,:),vectors(i,:)

     write(*,'(100F12.4)') vectors(i,:)/multvectors(i,:),vectors(i,:)/multvectors2(i,:)

!     write(*,'(100F12.4)') vectors(i,:)/multvectors(i,:),vectors(i,:)/multvectors2(i,:),vectors(i,:)


  enddo
  print *, "OKSTOP."



end program myprogram
