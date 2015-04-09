
!! ketot is matrix elements of -1/2 d^2/dx^2 in sinc function basis

!! NOTE DIFFERENT (BETTER) DIMENSIONS HERE for ketot(1-size:size-1)
!!    versus versions of SineDVR.F90 elsewhere, ketot(-size:size)

!! spacing and gridpoints is input; ketot and points is output
!! one column of ketot is stored; they are all the same
!!   do  kematrix(i,j)=ketot(i-j)   with kematrix(gridpoints,gridpoints)


subroutine sineDVR(ketot, points,gridpoints,spacing)
  implicit none

  integer :: i,gridpoints

  real*8 :: ketot(1-gridpoints:gridpoints-1),points(gridpoints)

  real*8 ::  mypi,spacing

  myPi   = 4d0*atan(1d0)

  do i=1,gridpoints
     points(i) = i*spacing
  enddo

  points(:)=points(:)-(points(1)+points(gridpoints))/2d0

  do i=1-gridpoints,gridpoints-1
     if (i==0) then
        ketot(i) = myPi**2/3.d0
     else
        ketot(i) = 2.d0/real(i,8)**2 * (-1)**(i)
     endif
  enddo

  ketot=ketot/2.d0/spacing**2 

end subroutine sineDVR



  





