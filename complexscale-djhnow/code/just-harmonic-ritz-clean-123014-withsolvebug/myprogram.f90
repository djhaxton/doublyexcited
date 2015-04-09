
!! Calculates bound states of 1D and 2D harmonic oscillator using sinc DVR or finite difference.

module mymodule

!! ADJUSTABLE PARAMETERS !!

  integer :: fdopt=0       !! 0 = sinc DVR  1 = finite difference

  integer, parameter :: &
       size=79, &          !! number of points (size of 1D problem)
       numvects=1,&        !! number of vectors to return.   numvects>=numconv
       numconv=1           !! number of vectors to attempt to converge

!!  complex*16 :: etarget=(0.123d0,0.456d0)   !! TARGET ENERGY

 complex*16 :: etarget=(5.85d0,-0.2d0)   !! TARGET ENERGY

! complex*16 :: etarget=(6.5d0,0d0) 

! complex*16 :: etarget=(0.71d0,0.15d0)   !! TARGET ENERGY



!! ENERGY RANGE for Krylov iterations (too small will cause gram-schmidt breakdown)
  real*8 :: range=40d0
!! ENERGY RANGE: will adjust etarget if an eigenvalue this close to etarget is obtained
  real*8 :: improverange=2d0
  real*8 :: improvethresh=2d0

  integer :: &
      order=12,&           !! krylov dimension
      genorder=12, &       !! krylov dimension for generating operator ~ (H-etarget)^-1*exp(-|H-etarget|^2/range^2)
      preporder=29         !! krylov dimension for preparing operator exp(-|H-etarget|^2/range^2)
  real*8 :: lanthresh=1d-6 !! lanczos convergence criterion
  real*8 :: length=20d0, & !! size of "box", determines spacing
       springconstant=1d0    !! spring constant.  harmonic oscillator.
  complex*16 ::      facs(2)=(/&
       (1.1d0,0d0),&
       (0.9d0,0d0)/)

!! THOSE WERE THE ADJUSTABLE PARAMETERS, ABOVE. 

  real*8 :: pi, spacing , const
  complex*16,allocatable :: potential(:),points(:),kematrix(:,:)

end module mymodule

!! matrix vector multiplication routine to be passed as argument
!!  to lanczos subroutine

module mysubroutines
contains

subroutine mymultsub(in,out,howmany,checksize)
  use mymodule
  integer, intent(in) :: howmany,checksize
  integer :: i
  complex*16 :: in(size,howmany),out(size,howmany),  eitheta


  eitheta=(12d0,1d0)/13d0


  if (size.ne.checksize) then
     print *, "AAUGH!!",size,checksize
     stop
  endif
  if (fdopt==0) then
     do i=1,howmany
        out(:,i)=in(:,i)*potential(:)*eitheta**2
     enddo
!     call zgemm('n','n',size,howmany,size,(1d0,0d0),kematrix,size,in,size,(1d0,0d0),out,size)
     call zgemm('n','n',size,howmany,size,CONJG(eitheta**2),kematrix,size,in,size,(1d0,0d0),out,size)
  else
     do i=1,howmany
        out(:,i)=in(:,i)*(potential(:)+1d0/spacing**2)
        out(1:size-1,i)=out(1:size-1,i)-in(2:size,i)/spacing**2/2
        out(2:size,i)=out(2:size,i)-in(1:size-1,i)/spacing**2/2
     enddo
  endif

!  out(:,:)=out(:,:)*(12d0,1d0)/13d0

!  out(:,:)=out(:,:)*(1d0,0.11111111d0)

!  out(:,:)=out(:,:)*(1d0,0.01d0)

!  out(:,:)=out(:,:)*(12d0,1d0)/13d0


end subroutine mymultsub


subroutine my2dmultsub(in,out,howmany,checksize)
  use mymodule
  integer, intent(in) :: howmany,checksize
  integer :: i
  complex*16 :: in(size,size,howmany),out(size,size,howmany), intrans(size,size,howmany), &
       outtrans(size,size,howmany)
  if (size**2.ne.checksize) then
     print *, "2D  AAUGH!!",size,checksize
     stop
  endif
  call mymultsub(in,out,howmany*size,size)
  out=out*facs(1)
  do i=1,size
     intrans(:,i,:)=in(i,:,:)
  enddo
  call mymultsub(intrans,outtrans,howmany*size,size)
  do i=1,size
     out(:,i,:)=out(:,i,:)+outtrans(i,:,:)*facs(2)
  enddo
end subroutine my2dmultsub


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
  implicit none

  integer :: i
  complex*16 :: vectors(size,numvects),values(numvects)
  complex*16 :: twovectors(size,size,numvects),twovalues(numvects)
  real*8 :: ketot(1-size:size-1),realpoints(size)
!! NOTE DIFFERENT (BETTER) DIMENSIONS HERE for ketot(1-size:size-1)
!!    versus versions of SineDVR.F90 elsewhere, ketot(-size:size)

  print *, "GO MYPROGRAM."
  call system("date")

  allocate(potential(size),points(size),kematrix(size,size))

  pi=4*atan(1d0)

  spacing = length/(size+1)


  call sinedvr(ketot,realpoints,size,spacing)
  points(:)=realpoints(:)
  potential(:)=0.5d0*springconstant*points(:)**2

  do  i=1,size
     kematrix(:,i)=ketot(1-i:size-i)
  enddo


!  print *, "ENTER ETARGET"
!  read(*,*) etarget
!  print *, "ENTER RANGE"
!  read(*,*) range

  print *, "CALL LANCZOS.",etarget,range

  call system("date")

if (1==0) then

  call harmonic_ritz(numvects, numconv, size, order,size,vectors,values,1,0,lanthresh,&
       mymultsub,allhdots0,etarget,range,improverange,improvethresh,-1,19,29)

  call system("date")
  print *, "VALUES"
  write(*,'(2F20.10)') values(:)


else

  call harmonic_ritz(numvects, numconv, size**2, order,size**2,twovectors,twovalues,1,0,lanthresh,&
       my2dmultsub,allhdots0,etarget,range,improverange,improvethresh,genorder,genorder,preporder)

  print *, "VALUES"
  write(*,'(2F20.10)') twovalues(:)

endif

end program myprogram
