
!! include this file in another to define what's below
!! (files blocklanczos_real.F90 and blocklanczos_complex.F90 provided for convenience)
!! or uncomment the following and compile this file as complex or real version
! 
!  uncomment
!#define DATATYPE complex*16
!#define REALFLAG 0
!     or
!#define DATATYPE real
!#define REALFLAG 1
!   ,  and
!#define OFLWR write(*,*)
!#define CFL
!#define mpifileptr *
!#define CFLST stop
!#define OFL

!
!  CALL THIS SUBROUTINE (AT END OF FILE)
!
!subroutine harmonic_ritz( &
!     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
!     inprintflag,guessflag,lanthresh,&
!     multsub,indotsub,   etarget,range, &
!     invmult_order, invexpmult_order, expmult_order)


!!! ON INPUT
!!
!!   lanblocknum        size of krylov block; number of eigenvectors calculated
!!   numvects           number of those to converge within lanthresh tolerance
!!   lansize            size of problem per processor
!!   order              Krylov space dimension
!!   maxiter            size of problem
!!   outvectors         if guessflag=1, guess vectors    outvectors(lansize,lanblocknum)
!!   inprintflag        print flag
!!   guessflag          use outvectors for initial krylov block?
!!   lanthresh          convergence criterion
!!   multsub            matrix vector multiplication subroutine.  arbitrary complex nonsymmetric.
!!                      multiples NN vectors.  must be programmed as 
!!                           multsub(in(lansize,NN),out(lansize,NN),NN,lansize)
!!   indotsub           dot product subroutine.
!!                      computes HERMITIAN dot products between a set of NN vectors and one of MM vectors.
!!                           indotsub(bravectors(lansize,NN),ketvectors(lansize,MM),lansize,NN,MM,result(NN,MM)
!!   etarget            TARGET ENERGY. 
!!   range              RANGE (in units of energy).  small range will cause breakdown of gram schmidt proceedure.
!!   invmult_order      NOT USED  (please set to -1)
!!   invexpmult_order   Krylov space dimension for subroutine that generates each member of main Krylov space
!!   expmult_order      Krylov space dimension for subroutine that generates initial Krylov vector
!!                        if guessflag=0
!!! ON OUTPUT
!!
!!   outvectors    eigenvectors
!!   outvalues     eigenvalues



! naah, just use REALFLAG.
!#if REALFLAG == 0
!#define MYGEMM ZGEMM
!#define MYZERO (0d0,0d0)
!#define MYONE (1d0,0d0)
!#elseif REALFLAG == 1
!#define MYGEMM DGEMM
!#define MYZERO 0d0
!#define MYONE 1d0
!#else
!    FAIL FAIL FAIL FAIL
!#endif



module hr_grammod
contains

!! allows real valued indefinite (possibly negative) inner product.
!!  norm of vector (1 or -1) output in mysign.
!!  vectors 1:m need to be orthonormal
!!  their norms are in prevsigns(1:m)

subroutine myhgramschmidt_many(n, howmany, m, previous, vector,inalldotsub,prevsigns,mysign,myfactors)
  implicit none
  
! n is the length of the vectors; m is how many to orthogonalize to

  integer :: n,m,mysign,prevsigns(m),howmany
  DATATYPE :: previous(n,howmany,m), vector(n,howmany),tempdots(m,howmany),myhdots(m,1),&
       tempprev(n,m,howmany)
  real*8 :: myfactors(howmany)
  integer :: i,j,k
  DATATYPE :: norm,temp(1,1)
  external :: inalldotsub

  do i=1,howmany
     tempprev(:,:,i)=previous(:,i,:)
  enddo
  do j=1,2 
     if (m.ne.0) then
        myhdots(:,:)=0d0
        do k=1,howmany
           call inalldotsub(tempprev(:,:,k),vector(:,k),n,m,1,tempdots(:,k))
           myhdots(:,1)=myhdots(:,1)+tempdots(:,k)*myfactors(k)
        enddo
     endif
     do i=1,m
       vector(:,:)=vector-previous(:,:,i)* myhdots(i,1) * prevsigns(i)
     enddo
     norm=0d0
     do k=1,howmany
        call inalldotsub(vector(:,k),vector(:,k),n,1,1,temp)
        norm=norm+temp(1,1)*myfactors(k)
     enddo
     if (abs(imag(norm)).gt.abs(real(norm,8))*1d-8) then
        OFLWR "AAUG IMAGNORM ", norm
     endif
     if (real(norm,8).ge.0d0) then
        mysign= 1
     else
        mysign= (-1)
     endif
     norm=sqrt(abs(norm))
     vector=vector/norm
     if (abs(norm).lt.1e-7) then
        OFLWR "Gram schmidt norm",norm,m; CFL
     endif
  enddo

end subroutine myhgramschmidt_many

end module hr_grammod




module hr_setmod

contains

subroutine fast_harmritz_setup( &
     lanblocknum, lansize,order,maxiter,initvectors,etarget,&
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,&
     lanham,lansimpleham,lansimpleovl,lansimplemultovl,lanvects)
  use hr_grammod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order
  external :: ingenmultsub,inh0multsub,indotsub
  integer, parameter :: manynum=2
  DATATYPE, intent(in) ::       initvectors(lansize,lanblocknum),etarget
  DATATYPE ::       lanham(lanblocknum,order,lanblocknum,order),lansimpleham(lanblocknum,order,lanblocknum,order),&
       lansimpleovl(lanblocknum,order,lanblocknum,order),lansimplemultovl(lanblocknum,order,lanblocknum,order),&
       lanvects(lansize,lanblocknum,order), &
       lansimplemultvects(lansize,lanblocknum,order), lanmanyvects(lansize,manynum,lanblocknum,order)
  real*8 :: factors(manynum)
  integer ::  iorder,i,thislanblocknum, printflag,lansigns(lanblocknum,order)

!! should be okay actually
  if (order*lanblocknum.gt.maxiter) then
     OFLWR "What in the world?  Order*lanblocknum is greater than lansize.", order,lanblocknum,lansize; CFLST
  endif

  printflag=inprintflag

  do iorder=1,order

!! DOES PARTIAL BLOCK IF AT THE END (leaving this)
     thislanblocknum=lanblocknum
     if (iorder*lanblocknum.ge.maxiter) then
        thislanblocknum=maxiter-(iorder-1)*lanblocknum
     endif
     if (iorder.eq.1) then
        lanvects(:,:,1)=initvectors(:,:)
     else
        call ingenmultsub(lanvects(:,:,iorder-1),lanvects(:,:,iorder),thislanblocknum,lansize,indotsub,inh0multsub)
     endif
     call inh0multsub(lanvects(:,:,iorder),lansimplemultvects(:,:,iorder),thislanblocknum,lansize)

     lansimplemultvects(:,:,iorder)=lansimplemultvects(:,:,iorder)-etarget*lanvects(:,:,iorder)
     lanmanyvects(:,1,:,:)=lansimplemultvects(:,:,:)
     lanmanyvects(:,2,:,:)=lanvects(:,:,:)
     factors(1)=1d0
     factors(2)=0d0
!!$     factors(2)=(-1)*abs(etarget)**2

     do i=1,thislanblocknum
        call myhgramschmidt_many(lansize, manynum, (iorder-1)*lanblocknum+i-1,  &
             lanmanyvects,lanmanyvects(:,:,i,iorder),&
             indotsub,lansigns,lansigns(i,iorder),factors)
     enddo

     lansimplemultvects(:,:,:)=lanmanyvects(:,1,:,:)
     lanvects(:,:,:)=lanmanyvects(:,2,:,:)

  enddo

  call indotsub(lansimplemultvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lanham)

  do iorder=1,order
     do i=1,lanblocknum
        lanham(:,:,i,iorder)=lanham(:,:,i,iorder)*lansigns(:,:)
     enddo
  enddo

  call indotsub(lanvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleovl)
  lansimplemultvects(:,:,:)=lansimplemultvects(:,:,:)+etarget*lanvects(:,:,:)
  call indotsub(lanvects(:,:,:),lansimplemultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleham)
  call indotsub(lansimplemultvects(:,:,:),lansimplemultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimplemultovl)

end subroutine fast_harmritz_setup


subroutine harm_orthog_oneblock(invects, &
     lanblocknum, lansize,etarget,&
     inprintflag,&
     inh0multsub,indotsub,&
     lantrans)
  use hr_grammod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,lanblocknum,inprintflag
  external :: inh0multsub,indotsub
  integer, parameter :: manynum=2
  DATATYPE, intent(in) :: etarget
  DATATYPE ::    invects(lansize,lanblocknum), tempvects(lansize,lanblocknum), &
       lantrans(lanblocknum,lanblocknum),&
       tempdots(lanblocknum,lanblocknum),&
       indots(lanblocknum,lanblocknum),&
       lansimplemultvects(lansize,lanblocknum)
  integer ::  i,printflag,j,k

  printflag=inprintflag

  call inh0multsub(invects(:,:),lansimplemultvects(:,:),lanblocknum,lansize)
  lansimplemultvects(:,:)=lansimplemultvects(:,:)-etarget*invects(:,:)
  call indotsub(lansimplemultvects(:,:),lansimplemultvects(:,:),lansize,lanblocknum,lanblocknum,indots)
  lantrans(:,:)=indots(:,:)

!! 2nd argument is C (or N, would be the same -- T does fail)

#if REALFLAG == 0
  call hermsqrtmatsmooth(lantrans,lanblocknum,lanblocknum,1d-10)
  call ZGEMM('N','C',lansize,lanblocknum,lanblocknum,(1d0,0d0),invects,lansize,lantrans,lanblocknum,&
       (0d0,0d0),tempvects,lansize)
#else
  call realsqrtmatsmooth(lantrans,lanblocknum,lanblocknum,1d-10)
  call DGEMM('N','T',lansize,lanblocknum,lanblocknum,1d0,invects,lansize,lantrans,lanblocknum,&
       0d0,tempvects,lansize)
#endif

  invects(:,:)=tempvects(:,:)

!! check

  call inh0multsub(invects(:,:),lansimplemultvects(:,:),lanblocknum,lansize)

  lansimplemultvects(:,:)=lansimplemultvects(:,:)-etarget*invects(:,:)

  call indotsub(lansimplemultvects(:,:),lansimplemultvects(:,:),lansize,lanblocknum,lanblocknum,tempdots)
  do i=1,lanblocknum
     do j=1,lanblocknum
        if (i==j) tempdots(i,j)=tempdots(i,j)-1d0
        if (abs(tempdots(i,j)).gt.1d-10) then
           print *, "AACK DOTS:"
           do k=1,lanblocknum
              write(*,'(100E20.3)') tempdots(:,k)
           enddo
           stop
        endif
     enddo
  enddo

end subroutine harm_orthog_oneblock



subroutine orthog_and_ham_oneblock(invects, &
     lanblocknum, lansize,&
     inprintflag,&
     inh0multsub,indotsub,&
     lanham)
  use hr_grammod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,lanblocknum,inprintflag
  external :: inh0multsub,indotsub
  DATATYPE :: &
       invects(lansize,lanblocknum), &
       tempvects(lansize,lanblocknum), &
       lantrans(lanblocknum,lanblocknum),&
       lanham(lanblocknum,lanblocknum),&
       indots(lanblocknum,lanblocknum),&
       lansimplemultvects(lansize,lanblocknum)
  integer ::  printflag

  printflag=inprintflag
  call indotsub(invects(:,:),invects(:,:),lansize,lanblocknum,lanblocknum,indots)
  lantrans(:,:)=indots(:,:)

!! 2nd argument is C (or N, would be the same -- T does fail)

#if REALFLAG == 0
  call hermsqrtmatsmooth(lantrans,lanblocknum,lanblocknum,1d-10)
  call ZGEMM('N','C',lansize,lanblocknum,lanblocknum,(1d0,0d0),invects,lansize,lantrans,lanblocknum,&
       (0d0,0d0),tempvects,lansize)
#else
  call realsqrtmatsmooth(lantrans,lanblocknum,lanblocknum,1d-10)
  call DGEMM('N','T',lansize,lanblocknum,lanblocknum,1d0,invects,lansize,lantrans,lanblocknum,&
       0d0,tempvects,lansize)
#endif

  invects(:,:)=tempvects(:,:)
  call inh0multsub(invects(:,:),lansimplemultvects(:,:),lanblocknum,lansize)
  call indotsub(invects(:,:),lansimplemultvects(:,:),lansize,lanblocknum,lanblocknum,lanham)

end subroutine orthog_and_ham_oneblock

end module hr_setmod





!!!!!!!!!!!!!!!!!!!!
!!!! ETARGETMOD !!!!
!!!!!!!!!!!!!!!!!!!!

module hr_etargetmod

!! TARGET ENERGY.  The numvects vectors
!!  closest to etarget are calculated.
  
  DATATYPE :: etarget = -1d3 

!! RANGE determines the operator used to generate the Krylov space.
!!   = (H-E)^-1 exp(-|H-E|^2/range^2)   or a pretty close approximation thereof

  real*8 :: range=100d0

  integer :: invmult_order=19
  integer :: invexpmult_order=19
  integer :: expmult_order=19

end module hr_etargetmod



!! fixed order, no restart, not for solving things per se, just a finite-order approximation
!!  to the inverse, error unchecked.

!! modeflag=1: usual.  calculates approximation to (h0-etarget)^-1 x invectors.
!!
!!   with QQ = exp(-(h0-etarget)^* (h0-etarget)/range^2)
!!          2:    approximation to (h0-etarget)^-1 x Q x invectors or something
!!                          range not used for modeflag=1
!!          3:    approx to Q x invectors
!!

module hr_basicsolvemod

contains

subroutine basicblocklansolve(&
     lanblocknum, lansize,order, maxiter, invectors, outvectors, &
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,     etarget, range, modeflag)
  use hr_setmod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,modeflag
  real*8, intent(in) :: range
  DATATYPE, intent(in) :: etarget,invectors(lansize,lanblocknum)
  DATATYPE, intent(out) :: outvectors(lansize,lanblocknum)
  DATATYPE ::   initvectors(lansize,lanblocknum), &
       lanham(lanblocknum*order,lanblocknum*order),lansimpleovl(lanblocknum*order,lanblocknum*order),&
       lansimplemultovl(lanblocknum,order,lanblocknum,order),&
       lansimpleham(lanblocknum,order,lanblocknum,order),&
       lanvects(lansize,lanblocknum,order),        lantrans(lanblocknum,lanblocknum), &
       solvevecs(lanblocknum,order,lanblocknum),&
       rhs(lanblocknum,order,lanblocknum)
  integer ::  thisdim,printflag
  external :: ingenmultsub,inh0multsub,indotsub

  if (order.le.0) then
     print *, "ORDERERROR", order; stop
  endif

!! should be okay actually
  if (order*lanblocknum.gt.maxiter) then
     OFLWR "What in the world?  Order*lanblocknum is greater than lansize.", order,lanblocknum,lansize; CFLST
  endif

  printflag=inprintflag
  thisdim=min(maxiter,order*lanblocknum)

  initvectors(:,:)=invectors(:,:)


  call harm_orthog_oneblock(initvectors, &
       lanblocknum, lansize,etarget,&
       inprintflag,&
       inh0multsub,indotsub,&
       lantrans)

  call fast_harmritz_setup( &
       lanblocknum, lansize,order,maxiter,initvectors,etarget,&
       inprintflag,&
       ingenmultsub,inh0multsub,indotsub,&
       lanham,lansimpleham,lansimpleovl,lansimplemultovl,lanvects)

  rhs(:,:,:)=0d0
  rhs(:,1,:)=lantrans(:,:)

  select case(modeflag)

  case(1)  !! usual, (H-E)^-1

!!no, don't need for harmritz  call complexsolve(lanham,thisdim,order*lanblocknum,rhs,lanblocknum,solvevecs)
!! instead just multiply by lanham as long as lansimplemultovl is identity which it should be

  case(2)

!!similarly to exponentiate (-1) x |(H-E)|^2, exponentiate (+1) x overlap
     lansimpleovl=lansimpleovl/range**2
#if REALFLAG == 0
     call cexpmat(lansimpleovl,thisdim,order*lanblocknum)
#else
     call rexpmat(lansimpleovl,thisdim,order*lanblocknum)
#endif

!! which way?  Haven't done the math, not sure.
     lanham=MATMUL(lanham,lansimpleovl)

  case(3)

     lanham=lansimpleovl/range**2
#if REALFLAG == 0
     call cexpmat(lanham,thisdim,order*lanblocknum)
#else
     call rexpmat(lanham,thisdim,order*lanblocknum)
#endif

  case default

     print *, "AAUGHH MODEFLAG", modeflag; stop
     
  end select

#if REALFLAG == 0
  call ZGEMM('N','N',thisdim,lanblocknum,thisdim,(1d0,0d0),lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       (0d0,0d0),solvevecs,order*lanblocknum)
  call ZGEMM('N','N',lansize,lanblocknum,thisdim,(1d0,0d0),lanvects(:,:,:),lansize,solvevecs,order*lanblocknum,&
       (0d0,0d0),outvectors,lansize)
#else
  call DGEMM('N','N',thisdim,lanblocknum,thisdim,1d0,lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       0d0,solvevecs,order*lanblocknum)
  call DGEMM('N','N',lansize,lanblocknum,thisdim,1d0,lanvects(:,:,:),lansize,solvevecs,order*lanblocknum,&
       0d0,outvectors,lansize)
#endif

end subroutine basicblocklansolve

end module hr_basicsolvemod



module hr_multmod

contains

subroutine invmultsub(in,out,howmany,checksize,allhdots0,mymultsub)
  use hr_etargetmod
  use hr_basicsolvemod
  implicit none
  integer, intent(in) :: howmany,checksize
  DATATYPE :: in(checksize,howmany),out(checksize,howmany)
  external :: allhdots0,mymultsub

  if (invmult_order*howmany.gt.checksize) then
     print *, "HMM, checkme. if doing parallel, maybe ok."; stop
  endif

  call basicblocklansolve(howmany,checksize,invmult_order,&
     invmult_order*howmany,in,out,1,mymultsub,mymultsub,allhdots0,etarget,range,1)

!  does NOT perform the same.  I lose the target vector more often.
!  do i=1,howmany
!     call basicblocklansolve(1,checksize,invmult_order,invmult_order,in(:,i),out(:,i),1,mymultsub,mymultsub,allhdots0,etarget,range,1)
!  enddo

end subroutine invmultsub



subroutine invexpmultsub(in,out,howmany,checksize,allhdots0,mymultsub)
  use hr_etargetmod
  use hr_basicsolvemod
  implicit none
  integer, intent(in) :: howmany,checksize
  DATATYPE :: in(checksize,howmany),out(checksize,howmany)
  external :: allhdots0,mymultsub

  if (invexpmult_order*howmany.gt.checksize) then
     print *, "HMM, checkme. if doing parallel, maybe ok."; stop
  endif

  call basicblocklansolve(howmany,checksize,invexpmult_order,&
     invexpmult_order*howmany,in,out,1,mymultsub,mymultsub,allhdots0,&
     etarget,range,2)

end subroutine invexpmultsub



subroutine expmultsub(in,out,howmany,checksize,allhdots0,mymultsub)
  use hr_etargetmod
  use hr_basicsolvemod
  implicit none
  integer, intent(in) :: howmany,checksize
  DATATYPE :: in(checksize,howmany),out(checksize,howmany)
  external :: allhdots0,mymultsub

  if (expmult_order*howmany.gt.checksize) then
     print *, "HMM, checkme. if doing parallel, maybe ok."; stop
  endif

  call basicblocklansolve(howmany,checksize,expmult_order,&
     expmult_order*howmany,in,out,1,mymultsub,mymultsub,allhdots0,etarget,range,3)

end subroutine expmultsub



subroutine idmultsub(in,out,howmany,checksize,allhdots0,mymultsub) !! Ok unused
  implicit none
  integer, intent(in) :: howmany,checksize
  DATATYPE :: in(checksize,howmany),out(checksize,howmany)
  external :: allhdots0,mymultsub

  out(:,:)=in(:,:)

end subroutine idmultsub




subroutine idfuncsub(values,thisdim)  !! Ok unused
end subroutine idfuncsub

subroutine invfuncsub(values,thisdim)
  use hr_etargetmod
  implicit none
  integer :: thisdim
  DATATYPE :: values(thisdim)
  values=1d0/values+etarget
!!$     values=1d0/values
end subroutine invfuncsub


end module hr_multmod




module hr_coreblocklanmod

contains

subroutine newblocklanczos0( &
     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
     inprintflag,guessflag,lanthresh,&
     inprepmultsub,ingenmultsub,inh0multsub,indotsub,infunctionsub)
  use hr_eigenmod
  use hr_etargetmod
  use hr_setmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,numvects,guessflag
  real*8, intent(in) :: lanthresh
  DATATYPE, intent(out) :: outvalues(lanblocknum)  , outvectors(lansize,lanblocknum)
  real*8 :: stopsum,  error(lanblocknum), temprealvectors(lansize,lanblocknum)
  DATATYPE ::   initvectors(lansize,lanblocknum), values(order*lanblocknum), &
       lanham(lanblocknum*order,lanblocknum*order),lansimpleovl(lanblocknum,order,lanblocknum,order),&
       lansimplemultovl(lanblocknum,order,lanblocknum,order),&
       laneigvects(lanblocknum*order,order*lanblocknum),lansimpleham(lanblocknum,order,lanblocknum,order),&
       lanvects(lansize,lanblocknum,order), tempvectors(lansize,lanblocknum), tempvects(lanblocknum*order,order*lanblocknum),&
       tempvects3(lanblocknum*order,order*lanblocknum),tempvects2(lanblocknum*order,order*lanblocknum),&
       simpledot(order*lanblocknum),simplehamdot(order*lanblocknum),simplequaddot(order*lanblocknum),&
       simpleexpect(order*lanblocknum),simplevardot(order*lanblocknum),errormeasure(order*lanblocknum), &
       lastvalues(lanblocknum), csums(1,1)
  DATATYPE :: littlelanham(lanblocknum,lanblocknum),littlelaneigvects(lanblocknum,lanblocknum),littlevalues(lanblocknum)
  integer ::  flag,j,i, thisdim,restarts,sortorder(order*lanblocknum),printflag
  external :: inprepmultsub,ingenmultsub,inh0multsub,indotsub,infunctionsub

!  print *, "CHECK:", &
!     lanblocknum, numvects, lansize,order, maxiter; stop


  if (numvects.gt.lanblocknum) then
     OFLWR "numvects must be less than or equal to lanblocknum"; CFLST
  endif

!! should be okay actually
  if (order*lanblocknum.gt.maxiter) then
     OFLWR "What in the world?  Order*lanblocknum is greater than lansize.", order,lanblocknum,lansize; CFLST
  endif

  printflag=inprintflag
  thisdim=min(maxiter,order*lanblocknum)

  values=0;   restarts=0
  initvectors=0; 
  lanham=0;  laneigvects=0; lanvects=0d0;

  lastvalues(:)=1d10;

     

  call random_seed()
  if (guessflag==0) then
!     do i=1,lanblocknum
!        call RANDOM_NUMBER(temprealvect)
!        initvectors(:,i)=temprealvect(:)*(0.1d0,1d0)**i
!        call RANDOM_NUMBER(temprealvect)
!        initvectors(:,i)=initvectors(:,i)+temprealvect(:)*(0.6d0,8d0)**(i+1)
!     enddo
     call RANDOM_NUMBER(temprealvectors)
     initvectors(:,:)=temprealvectors(:,:)
#if REALFLAG == 0
     call RANDOM_NUMBER(temprealvectors)
     initvectors(:,:)=initvectors(:,:)+(0d0,1d0)*temprealvectors(:,:)**3
#endif

  else
     initvectors(:,:)=outvectors(:,:)
  endif


  call inprepmultsub(initvectors,tempvectors,lanblocknum,lansize,indotsub,inh0multsub)
  initvectors(:,:)=tempvectors(:,:)

  flag=0
  do while (flag.eq.0)

     call fast_harmritz_setup( &
          lanblocknum, lansize,order,maxiter,initvectors,etarget,&
          inprintflag,&
          ingenmultsub,inh0multsub,indotsub,&
          lanham,lansimpleham,lansimpleovl,lansimplemultovl,lanvects)

#if REALFLAG == 0
     call get_eigen_c_plain(lanham,thisdim,order*lanblocknum,laneigvects,values)
#else
     call get_eigen_two(lanham,thisdim,order*lanblocknum,laneigvects,values)
#endif

     call infunctionsub(values,thisdim)
     
#if REALFLAG == 0
     call ZGEMM('N','N',thisdim,thisdim,thisdim,(1d0,0d0),lansimpleham,order*lanblocknum,laneigvects,order*lanblocknum,&
          (0d0,0d0),tempvects,order*lanblocknum)
     call ZGEMM('N','N',thisdim,thisdim,thisdim,(1d0,0d0),lansimpleovl,order*lanblocknum,laneigvects,order*lanblocknum,&
          (0d0,0d0),tempvects2,order*lanblocknum)
     call ZGEMM('N','N',thisdim,thisdim,thisdim,(1d0,0d0),lansimplemultovl,order*lanblocknum,laneigvects,order*lanblocknum,&
          (0d0,0d0),tempvects3,order*lanblocknum)
#else
     call DGEMM('N','N',thisdim,thisdim,thisdim,1d0,lansimpleham,order*lanblocknum,laneigvects,order*lanblocknum,&
          0d0,tempvects,order*lanblocknum)
     call DGEMM('N','N',thisdim,thisdim,thisdim,1d0,lansimpleovl,order*lanblocknum,laneigvects,order*lanblocknum,&
          0d0,tempvects2,order*lanblocknum)
     call DGEMM('N','N',thisdim,thisdim,thisdim,1d0,lansimplemultovl,order*lanblocknum,laneigvects,order*lanblocknum,&
          0d0,tempvects3,order*lanblocknum)
#endif
     
     do i=1,thisdim
        simplehamdot(i)=DOT_PRODUCT(laneigvects(1:thisdim,i),tempvects(1:thisdim,i))
        simpledot(i)=DOT_PRODUCT(laneigvects(1:thisdim,i),tempvects2(1:thisdim,i))
        simplequaddot(i)=DOT_PRODUCT(laneigvects(1:thisdim,i),tempvects3(1:thisdim,i))
        simpleexpect(i)=simplehamdot(i)/simpledot(i)
        
!!$ PREV              simplevardot(i)=simpleexpect(i)-values(i)

!!$ KEEPME-CHECK-SHOULD EQUAL ZERO FOR HARMONIC RITZ.  simplequaddot should be = 1.
!!$              simplevardot(i)=(simplequaddot(i) &
!!$                   -2*real(conjg(etarget)*simplehamdot(i))+abs(etarget)**2*simpledot(i))/conjg(simplehamdot(i)-etarget*simpledot(i)) &
!!$                   -values(i)+etarget

!!$ KEEPME-ERROR RELATIVE TO EXPECTATION VALUE <H^2>-<H>^2
!!$  USING THIS ONE

        simplevardot(i)=simplequaddot(i)/simpledot(i)-abs(simpleexpect(i))**2

!!$ ERROR RELATIVE TO EIGENVALUE
!!$              simplevardot(i)=simplequaddot(i)/simpledot(i)-2*real(conjg(simpleexpect(i))*values(i))+abs(values(i))**2

     enddo

!!$           call mycsort_getorder(simpleexpect,thisdim,etarget,sortorder)

     errormeasure(1:thisdim)=abs(simplevardot(1:thisdim)) 
!!$                *exp(abs(simpleexpect(1:thisdim)-etarget)/range) 

#if REALFLAG == 0 
     call mycsort_getorder(errormeasure,thisdim,(0d0,0d0),sortorder)
     
     call mycsort_vecs(laneigvects,thisdim,order*lanblocknum,sortorder)
     call mycsort_vals(values,thisdim,sortorder)
     call mycsort_vals(simpleexpect,thisdim,sortorder)
     call mycsort_vals(simplevardot,thisdim,sortorder)

     call ZGEMM('N','N',lansize,lanblocknum,thisdim,(1d0,0d0),lanvects(:,:,:),lansize,laneigvects,order*lanblocknum,&
          (0d0,0d0),outvectors,lansize)
#else
     call myrsort_getorder(errormeasure,thisdim,0d0,sortorder)
     
     call myrsort_vecs(laneigvects,thisdim,order*lanblocknum,sortorder)
     call myrsort_vals(values,thisdim,sortorder)
     call myrsort_vals(simpleexpect,thisdim,sortorder)
     call myrsort_vals(simplevardot,thisdim,sortorder)

     call DGEMM('N','N',lansize,lanblocknum,thisdim,1d0,lanvects(:,:,:),lansize,laneigvects,order*lanblocknum,&
          0d0,outvectors,lansize)
#endif


     do  j=1, lanblocknum
        call indotsub(outvectors(:,j:j),outvectors(:,j:j),lansize,1,1,csums)
        outvectors(:,j)=outvectors(:,j)/sqrt(csums(1,1))
     enddo

!! rediagonalize H in first block (rayleigh-ritz in first block)

!! orthogs outvectors and gets littlelanham

     call orthog_and_ham_oneblock(outvectors, &
          lanblocknum, lansize,&
          inprintflag,&
          inh0multsub,indotsub,&
          littlelanham)



#if REALFLAG == 0
     call get_eigen_c_plain(littlelanham,lanblocknum,lanblocknum,littlelaneigvects,littlevalues)
     call ZGEMM('N','N',lansize,lanblocknum,lanblocknum,(1d0,0d0),outvectors(:,:),lansize,littlelaneigvects,lanblocknum,&
          (0d0,0d0),tempvectors,lansize)
#else
     call get_eigen_two(littlelanham,lanblocknum,lanblocknum,littlelaneigvects,littlevalues)
     call DGEMM('N','N',lansize,lanblocknum,lanblocknum,1d0,outvectors(:,:),lansize,littlelaneigvects,lanblocknum,&
          0d0,tempvectors,lansize)
#endif
     outvectors(:,:)=tempvectors(:,:)

     do  j=1, lanblocknum
        call indotsub(outvectors(:,j:j),outvectors(:,j:j),lansize,1,1,csums)
        outvectors(:,j)=outvectors(:,j)/sqrt(csums(1,1))
     enddo


     if (printflag.ne.0) then
        OFL; write(mpifileptr,'(A25,1000F19.12)') "H-R EIGVALS            ", values(1:lanblocknum); CFL
        OFL; write(mpifileptr,'(A25,1000F19.12)') "R-R EXPECTATION VALS   ", simpleexpect(1:lanblocknum); CFL
        OFL; write(mpifileptr,'(A25,1000F38.31)') "ERRORS                 ", sqrt(abs(simplevardot(1:lanblocknum))); CFL
     endif
           
     call inh0multsub(outvectors(:,:),tempvectors(:,:),lanblocknum,lansize)
     
     do  j=1, lanblocknum
        call indotsub(outvectors(:,j:j),tempvectors(:,j:j),lansize,1,1,outvalues(j))
        tempvectors(:,j) = tempvectors(:,j) - outvalues(j) * outvectors(:,j)            
        call indotsub(tempvectors(:,j:j),tempvectors(:,j:j),lansize,1,1,csums)
        error(j)=sqrt(abs(csums(1,1)))
     enddo

     if (printflag.ne.0) then
        OFL; write(mpifileptr,'(A25,1000F19.12)') "EXPECTS AFTER REDIAG  ", outvalues(:); CFL
        OFL; write(mpifileptr,'(A25,1000F38.31)') "ERRORS AFTER REDIAG   ", error(:); CFL
     endif

     call myrsort_getorder(error(1:lanblocknum),lanblocknum,0d0,sortorder)

     stopsum=error(sortorder(numvects))

!!$     call mycsort_getorder(outvalues,lanblocknum,etarget,sortorder)
!!$     stopsum=0d0
!!$     do i=1,numvects
!!$        if (stopsum.lt.error(sortorder(i))) then
!!$           stopsum=error(sortorder(i))
!!$        endif
!!$     enddo

     
     if (stopsum.gt.lanthresh) then
!! If I'm at maximum order, i.e. krylov dimension = dimension of problem, I should get the right answer.
        if (thisdim.ge.maxiter) then
           OFLWR "MAX DIM REACHED, NO CONVERGENCE -- THRESHOLDS TOO HIGH? BUG?",stopsum,lanthresh; CFLST
        endif
     else
        flag=1
     endif

     if (flag.ne.0) then
        if (printflag.ne.0) then
           OFLWR "Converged. ",stopsum,lanthresh ;           CFL
           OFLWR "    restarts=",restarts
        endif
     else
        OFLWR "  Not converged. restarting.", stopsum,lanthresh,restarts; CFL
        restarts=restarts+1
     endif
     
     initvectors(:,:)=outvectors(:,:)
     
  enddo
  
end subroutine newblocklanczos0

end module hr_coreblocklanmod




subroutine harmonic_ritz( &
     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
     inprintflag,guessflag,lanthresh,&
     inh0multsub,indotsub,   inetarget,inrange, &
     in_invmult_order, in_invexpmult_order, in_expmult_order)

  use hr_coreblocklanmod
  use hr_etargetmod
  use hr_multmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,numvects,guessflag, &
       in_invmult_order, in_invexpmult_order, in_expmult_order
  DATATYPE, intent(in) :: inetarget
  real*8, intent(in) :: lanthresh,inrange
  DATATYPE, intent(out) :: outvalues(lanblocknum)  , outvectors(lansize,lanblocknum)
  external :: inh0multsub,indotsub

  etarget=inetarget
  range=inrange
  invmult_order=in_invmult_order
  invexpmult_order=in_invexpmult_order
  expmult_order=in_expmult_order

!  print *, "xCHECK:", &
!     lanblocknum, numvects, lansize,order, maxiter
!  print *, etarget,range,invmult_order,invexpmult_order,expmult_order; stop
  

  call newblocklanczos0( &
     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
     inprintflag,guessflag,lanthresh,&
     expmultsub,invexpmultsub,inh0multsub,indotsub, invfuncsub)

end subroutine harmonic_ritz














