
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


#define C1 (1d0,0d0)
#define C0 (0d0,0d0)
#define CI (0d0,1d0)




module hr_grammod
implicit none

contains

function threefactors(ereal,eimag)
  implicit none
  DATATYPE, dimension(3,3) :: threefactors
  real*8, intent(in) :: ereal,eimag

  threefactors(:,:)=RESHAPE((/ &
       C1,   -CI,      -ereal*C1,&
       CI,    C1,      -eimag*C1,&
       -ereal*C1,-eimag*C1,   (ereal**2+eimag**2)*C1 &
       /),(/3,3/))


!!$ USUAL HARMONIC RITZ (using hermitian matvecs only)
!!$  threefactors(:,:)=RESHAPE((/ &
!!$       C1,   -CI,      -ereal*C1,&
!!$       CI,    C1,      -eimag*C1,&
!!$       -ereal*C1,-eimag*C1,   (ereal**2+eimag**2)*C1 &
!!$       /),(/3,3/))


!!$ USUAL HARMONIC RITZ (ILL FOUNDED - MATRIX (1 1,-1-1) defeats it)
!!$   threefactors(:,:)=RESHAPE((/ &
!!$       C1,    C1,      -ereal*C1,&
!!$       C1,    C1,      eimag*CI,&
!!$       -ereal*C1,-eimag*CI,  (ereal**2+eimag**2)*C1 &
!!$       /),(/3,3/))

  
end function threefactors


!! allows real valued indefinite (possibly negative) inner product.
!!  norm of vector (1 or -1) output in mysign.
!!  vectors 1:m need to be orthonormal
!!  their norms are in prevsigns(1:m)

subroutine myhgramschmidt_many(n, howmany, m, previous, vector,inalldotsub,prevsigns,mysign,myfactors)
  implicit none
  
! n is the length of the vectors; m is how many to orthogonalize to

  integer :: n,m,mysign,prevsigns(m),howmany
  DATATYPE :: previous(n,howmany,m), vector(n,howmany),tempdots(m,howmany,howmany),myhdots(m,1),&
       tempprev(n,m,howmany)
  DATATYPE :: myfactors(howmany,howmany)
  integer :: i,j,k,kk
  DATATYPE :: norm,temp(1,1)
  external :: inalldotsub

  tempdots(:,:,:)=0d0

  do i=1,howmany
     tempprev(:,:,i)=previous(:,i,:)
  enddo
  do j=1,2 
     if (m.ne.0) then
        myhdots(:,:)=0d0
        do k=1,howmany
           do kk=1,howmany
              if (myfactors(kk,k).ne.0d0) then
                 call inalldotsub(tempprev(:,:,kk),vector(:,k),n,m,1,tempdots(:,kk,k))
                 myhdots(:,1)=myhdots(:,1)+tempdots(:,kk,k)*myfactors(kk,k)
              endif
           enddo
        enddo
     endif
     do i=1,m
       vector(:,:)=vector-previous(:,:,i)* myhdots(i,1) * prevsigns(i)
     enddo
     norm=0d0
     do k=1,howmany
        do kk=1,howmany
           if (myfactors(kk,k).ne.0d0) then
              call inalldotsub(vector(:,kk),vector(:,k),n,1,1,temp)
              norm=norm+temp(1,1)*myfactors(kk,k)
           endif
        enddo
     enddo
     if (abs(imag(norm)).gt.abs(real(norm,8))*1d-8) then
        OFLWR "AAUG IMAGNORM ", norm
     endif
     if (real(norm,8).ge.0d0) then
        mysign= 1
     else
        print *, "NEGSIGN",norm; stop   !! should never happen 
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





subroutine fast_harmritz_setup_new( &
     lanblocknum, lansize,order,maxiter,initvectors,in_etarget,&
     inprintflag,&
     ingenmultsub,inh0multsub,inhermmultsub,inantimultsub,&
     indotsub,&
     lanham,lansimpleham,lanquad,lansimpleovl,lanvects,whichway)
  use hr_grammod
  implicit none 

  integer, intent(in) :: whichway  !! 1=orthog wrt (H-E)^dag(H-E)  0=herm orthog

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order
  external :: ingenmultsub,inh0multsub,indotsub,inhermmultsub,inantimultsub
  DATATYPE, intent(in) ::       initvectors(lansize,lanblocknum),in_etarget
  DATATYPE ::       lanham(lanblocknum,order,lanblocknum,order),lansimpleham(lanblocknum,order,lanblocknum,order),&
       lansimpleovl(lanblocknum*order,lanblocknum*order),&
       lanquad(lanblocknum*order,lanblocknum*order),&
       lanhermmultvects(lansize,lanblocknum,order), &
       lanantimultvects(lansize,lanblocknum,order), &
       lanvects(lansize,lanblocknum,order), &
       lansimplemultvects(lansize,lanblocknum,order), lanmanyvects(lansize,3,lanblocknum,order)
  integer ::  iorder,i,thislanblocknum, printflag,lansigns(lanblocknum,order)
  real*8 :: ereal,eimag
  DATATYPE :: herefactors(3,3)
!! should be okay actually
  if (order*lanblocknum.gt.maxiter) then
     OFLWR "What in the world?  Order*lanblocknum is greater than lansize.", order,lanblocknum,lansize; CFLST
  endif

  if (whichway.eq.1) then
     herefactors(:,:)=threefactors(ereal,eimag)
  else
     herefactors(:,:)=0d0
     herefactors(3,3)=1d0
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

     call inhermmultsub(lanvects(:,:,iorder),lanhermmultvects(:,:,iorder),thislanblocknum,lansize)
     call inantimultsub(lanvects(:,:,iorder),lanantimultvects(:,:,iorder),thislanblocknum,lansize)

     lanmanyvects(:,1,:,:)=lanhermmultvects(:,:,:)
     lanmanyvects(:,2,:,:)=lanantimultvects(:,:,:)
     lanmanyvects(:,3,:,:)=lanvects(:,:,:)

     ereal=real(in_etarget,8)
     eimag=imag(in_etarget)



     do i=1,thislanblocknum
        call myhgramschmidt_many(lansize, 3, (iorder-1)*lanblocknum+i-1,  &
             lanmanyvects,lanmanyvects(:,:,i,iorder),&
             indotsub,lansigns,lansigns(i,iorder),herefactors)
     enddo

     lanhermmultvects(:,:,:)=lanmanyvects(:,1,:,:)
     lanantimultvects(:,:,:)=lanmanyvects(:,2,:,:)
     lanvects(:,:,:)=lanmanyvects(:,3,:,:)

  enddo


  lansimplemultvects(:,:,:)=lanhermmultvects(:,:,:)+lanantimultvects(:,:,:)*(0d0,1d0)



  call indotsub(lanvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleovl)
  call indotsub(lanvects(:,:,:),lansimplemultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleham)


  lansimplemultvects(:,:,:)=lansimplemultvects(:,:,:)-in_etarget*lanvects(:,:,:)

  call indotsub(lansimplemultvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lanham)

  call indotsub(lansimplemultvects(:,:,:),lansimplemultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lanquad)


  do iorder=1,order
     do i=1,lanblocknum
        if (lansigns(i,iorder).ne.1) then
           OFLWR "WTFFFF...",lansigns(i,iorder),i,iorder; CFLST
        endif
!!        lanham(:,:,i,iorder)=lanham(:,:,i,iorder)*lansigns(:,:)
     enddo
  enddo



end subroutine fast_harmritz_setup_new




subroutine fast_harmritz_setup_3( &
     lanblocknum, lansize,order,maxiter,initvectors,&
     inprintflag,&
     ingenmultsub,inh0multsub,inhermmultsub,inantimultsub,indotsub,&
     lanhermham,lanantiham,harmvalues,lanvects)
  use hr_grammod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order
  external :: ingenmultsub,inhermmultsub,inantimultsub,indotsub,inh0multsub
  DATATYPE, intent(in) ::       initvectors(lansize,lanblocknum)
  DATATYPE ::       lanantiham(lanblocknum*order,lanblocknum*order), &
       lanhermham(lanblocknum*order,lanblocknum*order), &
       lansimpleovl(lanblocknum*order,lanblocknum*order),&
       lanvects(lansize,lanblocknum,order), &
       tempvects(lansize,lanblocknum,order), &
       lanmanyvects(lansize,3,lanblocknum,order), &
       lanhermmultvects(lansize,lanblocknum,order), &
       lanantimultvects(lansize,lanblocknum,order), &
       harmvects(lanblocknum*order,lanblocknum*order)
  real*8 ::  harmvalues(lanblocknum*order)
  integer ::  iorder,i,thislanblocknum, printflag,lansigns(lanblocknum,order),thisdim


  thisdim=min(maxiter,order*lanblocknum)

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

     call inhermmultsub(lanvects(:,:,iorder),lanhermmultvects(:,:,iorder),thislanblocknum,lansize)
     call inantimultsub(lanvects(:,:,iorder),lanantimultvects(:,:,iorder),thislanblocknum,lansize)

     lanmanyvects(:,1,:,:)=lanhermmultvects(:,:,:)
     lanmanyvects(:,2,:,:)=lanantimultvects(:,:,:)
     lanmanyvects(:,3,:,:)=lanvects(:,:,:)

     do i=1,thislanblocknum
        call myhgramschmidt_many(lansize, 3, (iorder-1)*lanblocknum+i-1,  &
             lanmanyvects,lanmanyvects(:,:,i,iorder),&
             indotsub,lansigns,lansigns(i,iorder),threefactors(0d0,0d0))
     enddo

     lanhermmultvects(:,:,:)=lanmanyvects(:,1,:,:)
     lanantimultvects(:,:,:)=lanmanyvects(:,2,:,:)
     lanvects(:,:,:)=lanmanyvects(:,3,:,:)

  enddo


  call indotsub(lanvects(:,:,:),lanvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lansimpleovl)


#if REALFLAG == 0
  call get_eigen_herm(lansimpleovl,thisdim,order*lanblocknum,harmvects,harmvalues)
#else
  call get_eigen_two(lansimpleovl,thisdim,order*lanblocknum,harmvects,harmvalues)
#endif
  harmvalues(:)=1d0/harmvalues(:)

  do i=1,thisdim

     harmvects(:,i)=harmvects(:,i)*sqrt(harmvalues(i))  !! NOW NORMALIZED.

  enddo


  tempvects(:,:,:)=0d0
#if REALFLAG == 0
  call ZGEMM('N','N',lansize,thisdim,thisdim,(1d0,0d0),lanvects(:,:,:),lansize,harmvects,order*lanblocknum,&
       (0d0,0d0),tempvects,lansize)
#else
  call DGEMM('N','N',lansize,thisdim,thisdim,1d0,lanvects(:,:,:),lansize,harmvects,order*lanblocknum,&
       0d0,tempvects,lansize)
#endif
  
  lanvects(:,:,:)=tempvects(:,:,:)


  do iorder=1,order

     call inhermmultsub(lanvects(:,:,iorder),lanhermmultvects(:,:,iorder),thislanblocknum,lansize)
     call inantimultsub(lanvects(:,:,iorder),lanantimultvects(:,:,iorder),thislanblocknum,lansize)

  enddo

  call indotsub(lanvects(:,:,:),lanhermmultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lanhermham)

  call indotsub(lanvects(:,:,:),lanantimultvects(:,:,:),lansize,order*lanblocknum,order*lanblocknum,lanantiham)

!  print *, "HERMHAM.", lanhermham(1:3,1:3); stop


!! NOW, THE BASIS SIMULTANEOUSLY DIAGONALIZES (H^\dag H) and the identity operator
!!   with respect to the hermitian inner product  (hermitian orthogonal)
!!   
!! <psi_i | H^dag H | psi_j > = \delta_ij harmvalues(i)
!! <psi_i | H^dag H | psi_j > = \delta_ij 1/harmvalues(i)


!! NOW MAKE IT NORMALIZED

!  j=0
!  do iorder=1,order
!     thislanblocknum=lanblocknum
!     if (iorder*lanblocknum.ge.maxiter) then
!        thislanblocknum=maxiter-(iorder-1)*lanblocknum
!     endif
!     do i=1,thislanblocknum
!        j=j+1

!        lanvects(:,i,iorder)=lanvects(:,i,iorder)/harmvalues(j)
!     enddo
!  enddo

end subroutine fast_harmritz_setup_3





subroutine harm_orthog_oneblock_new(invects, &
     lanblocknum, lansize,in_etarget,&
     inprintflag,&
     inhermmultsub,inantimultsub,&
     indotsub,&
     lantrans,whichway)
  use hr_grammod
  use hr_eigenmod
  implicit none 

  integer, intent(in) :: whichway   !! 1=orthog wrt (H-E)^dag(H-E)  0=herm orthog

  integer,intent(in) :: lansize,lanblocknum,inprintflag
  external :: indotsub,inantimultsub,inhermmultsub

  DATATYPE, intent(in) :: in_etarget
  DATATYPE ::    invects(lansize,lanblocknum), tempvects(lansize,lanblocknum), &
       lantrans(lanblocknum,lanblocknum),&
       tempdots(lanblocknum,lanblocknum),&
       checkdots(lanblocknum,lanblocknum),&
       lanmanyvectors(lansize,lanblocknum,3),&
       lanhermmultvects(lansize,lanblocknum),&
       lanantimultvects(lansize,lanblocknum)
  DATATYPE :: herefactors(3,3)
  real*8 :: ereal,eimag
  integer ::  i,printflag,j,k,jj,kk,flag

  ereal=real(in_etarget,8)
  eimag=imag(in_etarget)

  printflag=inprintflag

  if (whichway.eq.1) then
     herefactors(:,:)=threefactors(ereal,eimag)
  else
     herefactors(:,:)=0d0
     herefactors(3,3)=1d0
  endif

  do jj=1,2

     call inhermmultsub(invects(:,:),lanhermmultvects(:,:),lanblocknum,lansize)
     call inantimultsub(invects(:,:),lanantimultvects(:,:),lanblocknum,lansize)
     
     lanmanyvectors(:,:,1)=lanhermmultvects(:,:)
     lanmanyvectors(:,:,2)=lanantimultvects(:,:)
     lanmanyvectors(:,:,3)=invects(:,:)
     
     
     checkdots(:,:)=0d0
     do k=1,3
        do kk=1,3
           if (herefactors(kk,k).ne.0d0) then
              call indotsub(lanmanyvectors(:,:,kk),lanmanyvectors(:,:,k),lansize,lanblocknum,lanblocknum,tempdots)
              checkdots(:,:)=checkdots+tempdots(:,:)*herefactors(kk,k)
           endif
        enddo
     enddo

     if (jj.eq.1) then
        lantrans(:,:)=checkdots(:,:)

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

     else    !! check

        flag=0
        do i=1,lanblocknum
           do j=1,lanblocknum
              if (i==j) checkdots(i,j)=checkdots(i,j)-1d0
              if (abs(checkdots(i,j)).gt.1d-10) then
                 flag=1
              endif
           enddo
        enddo
        
        if (flag.eq.1) then
           print *, "AACK DxxxOTS:"
           do k=1,lanblocknum
              write(*,'(100E20.3)') checkdots(:,k)
           enddo
           stop
        endif

     endif
     
  enddo

end subroutine harm_orthog_oneblock_new







end module hr_setmod








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


subroutine myblocklansolve(&
     lanblocknum, lansize,order, maxiter, invectors, outvectors, &
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,     in_etarget, range, modeflag)
  use hr_setmod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,modeflag
  real*8, intent(in) :: range
  DATATYPE, intent(in) :: in_etarget,invectors(lansize,lanblocknum)
  DATATYPE, intent(out) :: outvectors(lansize,lanblocknum)
  external :: ingenmultsub,inh0multsub,indotsub

call newblocklansolve(&
     lanblocknum, lansize,order, maxiter, invectors, outvectors, &
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,     in_etarget, range, modeflag)
end subroutine myblocklansolve


subroutine basicblocklansolve(&
     lanblocknum, lansize,order, maxiter, invectors, outvectors, &
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,     in_etarget, range, modeflag)
  use hr_setmod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,modeflag
  real*8, intent(in) :: range
  DATATYPE, intent(in) :: in_etarget,invectors(lansize,lanblocknum)
  DATATYPE, intent(out) :: outvectors(lansize,lanblocknum)
  DATATYPE ::   initvectors(lansize,lanblocknum), &
       lanham(lanblocknum*order,lanblocknum*order),lansimpleovl(lanblocknum*order,lanblocknum*order),&
       lansimpleham(lanblocknum*order,lanblocknum*order),&
       lanquad(lanblocknum*order,lanblocknum*order),&
       lanvects(lansize,lanblocknum*order),        lantrans(lanblocknum,lanblocknum), &
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


  call harm_orthog_oneblock_new(initvectors, &
       lanblocknum, lansize,in_etarget,&
       inprintflag,&
       h0multsub_herm_shell, & 
       h0multsub_antiherm_shell, & 
       indotsub,&
       lantrans,1)

  call fast_harmritz_setup_new( &
       lanblocknum, lansize,order,maxiter,initvectors,in_etarget,&
       inprintflag,&
       ingenmultsub,&
       inh0multsub,&
       h0multsub_herm_shell, & 
       h0multsub_antiherm_shell, & 
       indotsub,&
       lanham,lansimpleham,lanquad,lansimpleovl,lanvects,1)





  rhs(:,:,:)=0d0
  rhs(:,1,:)=lantrans(:,:)

  select case(modeflag)

  case(1)  !! usual, (H-E)^-1

!!no, don't need for harmritz  call complexsolve(lanham,thisdim,order*lanblocknum,rhs,lanblocknum,solvevecs)
!! instead just multiply by lanham as long as variance matrix is identity which it should be

  case(2)

     lansimpleovl=lansimpleovl/range**2
#if REALFLAG == 0
     call cexpmat(lansimpleovl,thisdim,order*lanblocknum)
#else
     call rexpmat(lansimpleovl,thisdim,order*lanblocknum)
#endif

!! which way?  Haven't done the math, not sure.
!     lanham=MATMUL(lanham,lansimpleovl)

     lanham=MATMUL(lansimpleovl,lanham)

  case(3)

!! to exponentiate (-1) x |(H-E)|^2, exponentiate (+1) x overlap

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
  call ZGEMM('N','N',lansize,lanblocknum,thisdim,(1d0,0d0),lanvects(:,:),lansize,solvevecs,order*lanblocknum,&
       (0d0,0d0),outvectors,lansize)
#else
  call DGEMM('N','N',thisdim,lanblocknum,thisdim,1d0,lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       0d0,solvevecs,order*lanblocknum)
  call DGEMM('N','N',lansize,lanblocknum,thisdim,1d0,lanvects(:,:),lansize,solvevecs,order*lanblocknum,&
       0d0,outvectors,lansize)
#endif

contains 

#include "shell.h.f90"

end subroutine basicblocklansolve






subroutine newblocklansolve(&
     lanblocknum, lansize,order, maxiter, invectors, outvectors, &
     inprintflag,&
     ingenmultsub,inh0multsub,indotsub,     in_etarget, range, modeflag)
  use hr_setmod
  use hr_eigenmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,modeflag
  real*8, intent(in) :: range
  DATATYPE, intent(in) :: in_etarget,invectors(lansize,lanblocknum)
  DATATYPE, intent(out) :: outvectors(lansize,lanblocknum)
  DATATYPE ::   initvectors(lansize,lanblocknum), &
       lanham(lanblocknum*order,lanblocknum*order),lansimpleovl(lanblocknum*order,lanblocknum*order),&
       lansimpleham(lanblocknum*order,lanblocknum*order),&
       lanquad(lanblocknum*order,lanblocknum*order),&
       lanvects(lansize,lanblocknum*order),        lantrans(lanblocknum,lanblocknum), &
       solvevecs(lanblocknum*order,lanblocknum),&
       rhs(lanblocknum*order,lanblocknum)
  integer ::  thisdim,printflag,i
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


  call harm_orthog_oneblock_new(initvectors, &
       lanblocknum, lansize,in_etarget,&
       inprintflag,&
       h0multsub_herm_shell, & 
       h0multsub_antiherm_shell, & 
       indotsub,&
       lantrans,0)   !! herm orthog.


  call fast_harmritz_setup_new( &
       lanblocknum, lansize,order,maxiter,initvectors,in_etarget,&
       inprintflag,&
       ingenmultsub,&
       inh0multsub,&
       h0multsub_herm_shell, & 
       h0multsub_antiherm_shell, & 
       indotsub,&
       lanham,lansimpleham,lanquad,lansimpleovl,lanvects,0)



  select case(modeflag)

  case(1)  !! usual, (H-E)^-1

     rhs(:,:)=0d0
     rhs(1:lanblocknum,:)=lantrans(:,:)
     
     do i=1,lanblocknum
        rhs(:,i)=MATMUL(lanham(:,:),rhs(:,i))
     enddo
        
     call complexsolve(lanquad,thisdim,order*lanblocknum,rhs,lanblocknum,solvevecs)

  case(2)

     print *, "HMM, doesn't really work"; stop

     rhs(:,:)=0d0
     rhs(1:lanblocknum,:)=lantrans(:,:)

     do i=1,lanblocknum
        rhs(:,i)=MATMUL(lanham(:,:),rhs(:,i))
     enddo

     lanquad(:,:)=lanquad(:,:)/range**2

#if REALFLAG == 0
     call cexpmat(lanquad,thisdim,order*lanblocknum)
#else
     call rexpmat(lanquad,thisdim,order*lanblocknum)
#endif
     do i=1,thisdim
        lanquad(i,i)=lanquad(i,i)-1d0
     enddo
        
     call complexsolve(lanquad,thisdim,order*lanblocknum,rhs,lanblocknum,solvevecs)

  case(3)

     rhs(:,:)=0d0
     rhs(1:lanblocknum,:)=lantrans(:,:)

!! to exponentiate (-1) x |(H-E)|^2, exponentiate (+1) x overlap

     lanham=(-1)*lanquad/range**2
#if REALFLAG == 0
     call cexpmat(lanham,thisdim,order*lanblocknum)
#else
     call rexpmat(lanham,thisdim,order*lanblocknum)
#endif

#if REALFLAG == 0
  call ZGEMM('N','N',thisdim,lanblocknum,thisdim,(1d0,0d0),lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       (0d0,0d0),solvevecs,order*lanblocknum)
#else
  call DGEMM('N','N',thisdim,lanblocknum,thisdim,1d0,lanham(:,:),order*lanblocknum,rhs,order*lanblocknum,&
       0d0,solvevecs,order*lanblocknum)
#endif

  case default

     print *, "AAUGHH MODEFLAG", modeflag; stop
     
  end select


#if REALFLAG == 0
  call ZGEMM('N','N',lansize,lanblocknum,thisdim,(1d0,0d0),lanvects(:,:),lansize,solvevecs,order*lanblocknum,&
       (0d0,0d0),outvectors,lansize)
#else
  call DGEMM('N','N',lansize,lanblocknum,thisdim,1d0,lanvects(:,:),lansize,solvevecs,order*lanblocknum,&
       0d0,outvectors,lansize)
#endif

contains 

#include "shell.h.f90"

end subroutine newblocklansolve




end module hr_basicsolvemod



!!!!!!!!!!!!!!!!!!!!
!!!! ETARGETMOD !!!!
!!!!!!!!!!!!!!!!!!!!

module hr_etargetmod

!! TARGET ENERGY.  The numvects vectors
!!  closest to etarget are calculated.
  
  DATATYPE :: etarget = -1d3 

!! RANGE determines the operator used to generate the Krylov space.
!!   = (H-E)^-1 exp(-|H-E|^2/range^2)   or a pretty close approximation thereof

!! DON'T CHANGE THESE; WILL HAVE NO EFFECT; ARE SET ON CALL TO HARMONIC_RITZ.
  real*8 :: range=-999d0
  real*8 :: improverange=-999d0
  real*8 :: improvethresh=-999d0

  integer :: invmult_order= -1
  integer :: invexpmult_order= -1
  integer :: expmult_order= -1

  DATATYPE :: etargettest

end module hr_etargetmod



module hr_multmod

contains

subroutine invmultsub3(in,out,howmany,checksize,allhdots0,mymultsub)
  use hr_etargetmod
  use hr_basicsolvemod
  implicit none
  integer, intent(in) :: howmany,checksize
  integer :: i
  DATATYPE :: in(checksize,howmany),out(checksize,howmany),temp(checksize,howmany)
  external :: allhdots0,mymultsub
  
  temp(:,:)=in(:,:)
  do i=1,3
     call invmultsub(temp,out,howmany,checksize,allhdots0,mymultsub)
     temp(:,:)=out(:,:)
  enddo
end subroutine invmultsub3


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

  call myblocklansolve(howmany,checksize,invmult_order,invmult_order*howmany,in,out,1,mymultsub,mymultsub,allhdots0,etarget,range,1)

!  does NOT perform the same.  I lose the target vector more often.
!  do i=1,howmany
!     call myblocklansolve(1,checksize,invmult_order,invmult_order,in(:,i),out(:,i),1,mymultsub,mymultsub,allhdots0,etarget,range,1)
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

  call myblocklansolve(howmany,checksize,invexpmult_order,invexpmult_order*howmany,in,out,1,mymultsub,mymultsub,allhdots0,etarget,range,2)

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

  call myblocklansolve(howmany,checksize,expmult_order,expmult_order*howmany,in,out,1,mymultsub,mymultsub,allhdots0,etarget,range,3)

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



subroutine betterharmonic( &
     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
     inprintflag,guessflag,lanthresh,&
     inprepmultsub,ingenmultsub,inh0multsub,indotsub)
  use hr_eigenmod
  use hr_etargetmod
  use hr_setmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,numvects,guessflag
  real*8, intent(in) :: lanthresh
  DATATYPE, intent(out) :: outvalues(lanblocknum)  , outvectors(lansize,lanblocknum)
  real*8 :: stopsum,  error(lanblocknum), temprealvectors(lansize,lanblocknum)
  DATATYPE ::   initvectors(lansize,lanblocknum), values(order*lanblocknum), &
       dvalues_dreal(order*lanblocknum), &
       dvalues_dimag(order*lanblocknum), &
       dinvvalues_dreal(order*lanblocknum), &
       dinvvalues_dimag(order*lanblocknum), &
       laneigvects(lanblocknum*order,order*lanblocknum),&
       lefteigvects(lanblocknum*order,order*lanblocknum),&
       laneigmultvects_a(lanblocknum*order,order*lanblocknum),&
       laneigmultvects_anti(lanblocknum*order,order*lanblocknum),&
       laneigmultvects_herm(lanblocknum*order,order*lanblocknum),&
       bdots(lanblocknum*order), adots(lanblocknum*order), covldots(lanblocknum*order),&
       antihamdots(lanblocknum*order),hermhamdots(lanblocknum*order),&
       expects(lanblocknum*order),epredict(order*lanblocknum), &
       lanhamtot(lanblocknum*order,lanblocknum*order),&
       lanhermham(lanblocknum*order,lanblocknum*order),&
       lanantiham(lanblocknum*order,lanblocknum*order),&
       lanvects(lansize,lanblocknum,order), &
       tempvectors(lansize,lanblocknum), &
       lastvalues(lanblocknum), csums(1,1), &
       AMAT(lanblocknum*order,lanblocknum*order),BMAT(lanblocknum*order,lanblocknum*order)
  real*8 :: vardots(lanblocknum*order),quaddots(lanblocknum*order),ovldots(lanblocknum*order)
  real*8 ::        dreal(order*lanblocknum),       dimag(order*lanblocknum),&
       distance(order*lanblocknum)

  real*8 :: harmvalues(order*lanblocknum)
  integer ::  flag,j,i, thisdim,restarts,sortorder(order*lanblocknum),printflag
  real*8 :: rhs(2), mmat(2,2),mydist,olddist
  integer :: ipivtwo(2),info
  real*8 :: ereal, eimag

  integer :: iiter,iilow,iihigh,nmult
  DATATYPE :: templaneigvects(lanblocknum*order,lanblocknum*order),templefteigvects(lanblocknum*order,lanblocknum*order),&
       tempvalues(lanblocknum*order)


  external :: inprepmultsub,ingenmultsub,inh0multsub,indotsub

  if (abs(etarget).eq.0d0) then
     OFLWR "ETARGET CANNOT BE ZERO."; CFLST
  endif

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
  lanhamtot=0;  lanantiham=0; lanhermham=0
  laneigvects=0; lanvects=0d0;

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

  olddist=1d10
  flag=0

  do while (flag.eq.0)

     call fast_harmritz_setup_3( &
          lanblocknum, lansize,order,maxiter,initvectors,&
          inprintflag,&
          ingenmultsub,&
          inh0multsub,&
          h0multsub_herm_shell, & 
          h0multsub_antiherm_shell, & 
          indotsub,&
          lanhermham,lanantiham,&
          harmvalues,lanvects)

     etargettest=etarget


!! ITERATION ZERO - do HR at etarget.  Then refine, based on epredict.


     do iiter=0,thisdim

        print *, "IITER",iiter

        if (iiter.eq.0) then
           iilow=1
           iihigh=thisdim
           nmult=thisdim
        else
           iilow=iiter
           iihigh=iiter
           nmult=1
           etargettest=epredict(iiter)
        endif

        ereal=real(etargettest,8)
        eimag=imag(etargettest)

        BMAT(:,:)=lanhermham(:,:) - (0d0,1d0) * lanantiham(:,:)
        do i=1,thisdim
           BMAT(i,i)=BMAT(i,i)  - ereal + (0d0,1d0)*eimag
        enddo
        
        AMAT(:,:)=((-2)*ereal*lanhermham(:,:) + (-2)*eimag*lanantiham)
        
        do i=1,thisdim
           AMAT(i,i)=AMAT(i,i)+harmvalues(i) + (ereal**2+eimag**2)
        enddo
        

#if REALFLAG == 0
        call get_eigen_c_gen(BMAT,AMAT,thisdim,order*lanblocknum,templaneigvects,templefteigvects,tempvalues)
#else
        call get_eigen_gen(BMAT,AMAT,thisdim,order*lanblocknum,laneigvects,values)
#endif
        if (iiter.eq.0) then
           laneigvects(:,:)=templaneigvects(:,:)
           lefteigvects(:,:)=templefteigvects(:,:)
           values(:)=tempvalues(:)
        else

           call mycsort_getorder(1/tempvalues,thisdim,(0d0,0d0),sortorder)

           laneigvects(:,iiter)=templaneigvects(:,sortorder(1))
           lefteigvects(:,iiter)=templefteigvects(:,sortorder(1))
           values(iiter)=tempvalues(sortorder(1))

        endif

#if REALFLAG == 0
        call ZGEMM('N','N',thisdim,nmult,thisdim,(1d0,0d0),AMAT(:,:),order*lanblocknum,&
             laneigvects(:,iilow),order*lanblocknum,(0d0,0d0),laneigmultvects_a(:,iilow),order*lanblocknum)
        call ZGEMM('N','N',thisdim,nmult,thisdim,(1d0,0d0),lanhermham(:,:),order*lanblocknum,&
             laneigvects(:,iilow),order*lanblocknum,(0d0,0d0),laneigmultvects_herm(:,iilow),order*lanblocknum)
        call ZGEMM('N','N',thisdim,nmult,thisdim,(1d0,0d0),lanantiham(:,:),order*lanblocknum,&
             laneigvects(:,iilow),order*lanblocknum,(0d0,0d0),laneigmultvects_anti(:,iilow),order*lanblocknum)
#else
        call DGEMM('N','N',thisdim,nmult,thisdim,1d0,AMAT(:,:),order*lanblocknum,&
             laneigvects(:,iilow),order*lanblocknum,0d0,laneigmultvects_a(:,iilow),order*lanblocknum)
        call DGEMM('N','N',thisdim,nmult,thisdim,1d0,lanhermham(:,:),order*lanblocknum,&
             laneigvects(:,iilow),order*lanblocknum,0d0,laneigmultvects_herm(:,iilow),order*lanblocknum)
        call DGEMM('N','N',thisdim,nmult,thisdim,1d0,lanantiham(:,:),order*lanblocknum,&
             laneigvects(:,iilow),order*lanblocknum,0d0,laneigmultvects_anti(:,iilow),order*lanblocknum)
#endif

        do i=iilow,iihigh


           quaddots(i)=DOT_PRODUCT(&   !! OK CONV
                laneigvects(:,i),laneigvects(:,i)*harmvalues(:)) 

           ovldots(i)=DOT_PRODUCT(&   !! OK CONV
                laneigvects(:,i),laneigvects(:,i))
           
           
           hermhamdots(i)=DOT_PRODUCT(laneigvects(:,i),laneigmultvects_herm(:,i))
           antihamdots(i)=DOT_PRODUCT(laneigvects(:,i),laneigmultvects_anti(:,i))
           
           expects(i)=(hermhamdots(i)+(0d0,1d0)*antihamdots(i))/ovldots(i)
           
           vardots(i)=abs(quaddots(i)/ovldots(i)-abs(expects(i))**2)

!!!!           lefteigvects(:,:)=laneigvects(:,:)
           
           covldots(i)=DOT_PRODUCT(&
                lefteigvects(:,i),laneigvects(:,i))
           
           hermhamdots(i)=DOT_PRODUCT(lefteigvects(:,i),laneigmultvects_herm(:,i))
           antihamdots(i)=DOT_PRODUCT(lefteigvects(:,i),laneigmultvects_anti(:,i))
           
           bdots(i)= hermhamdots(i) - (0d0,1d0)* antihamdots(i) -  (ereal - (0d0,1d0)*eimag)*covldots(i)
           
           adots(i)=DOT_PRODUCT(lefteigvects(:,i),laneigmultvects_a(:,i))

!! H = H_H + I * H_A
!!
!! A = H_H^2 + H_A^2 - 2Er H_R + 2 Ei H_A + Er^2 + Ei^2

!! A is overlap

           
           dinvvalues_dreal(i)=( (-1)*covldots(i)  - values(i)* ((-2)*hermhamdots(i)+2*ereal*covldots(i)) ) / adots(i)
           
           dinvvalues_dimag(i)=( (0d0,1d0)*covldots(i) - values(i)*(2*antihamdots(i)+2*eimag*covldots(i)) ) / adots(i)
           dvalues_dreal(i)=(-1)/values(i)**2 * dinvvalues_dreal(i)
           dvalues_dimag(i)=(-1)/values(i)**2 * dinvvalues_dimag(i)
           
           values(i)=1d0/values(i)

!! looking for value 0

!! solve 0 = values(i) +  dvalues_dreal(i) * dreal + dvalues_dimag(i) * dimag
!!   for real dreal and dimag
!!

           rhs(:)=(/ -real(values(i),8), -imag(values(i))/)
           
           mmat(:,:)=RESHAPE( (/ &
                real(dvalues_dreal(i),8), imag(dvalues_dreal(i)), &
                real(dvalues_dimag(i),8), imag(dvalues_dimag(i)) /), (/2,2/))
           call DGESV(2,1,mmat,2,ipivtwo,rhs,2,info)
           if (info.ne.0) then
              OFLWR "OOGABLAH DGESV INFO ", info; stop
           endif
           dreal(i)=rhs(1); dimag(i)=rhs(2)
           
           values(i)=values(i)+etargettest
           
           epredict(i)=etargettest+dreal(i)+(0d0,1d0)*dimag(i)

        enddo

        print *, "OKAY DONE LOOP IITER=",iiter



        if (iiter.eq.0.or.iiter.eq.thisdim) then

           distance(:)=abs(etarget - epredict(:))

           call myrsort_getorder(vardots,thisdim,0d0,sortorder)
           
           do i=1,lanblocknum
              mydist=distance(sortorder(i))
              if (mydist.lt.olddist.and.mydist.lt.improvethresh  .and.  1==0) then
                 etargettest=epredict(sortorder(i))
                 OFLWR "SETTING ETARGET ", etargettest; CFL
                 olddist=mydist
              endif
           enddo
           
           
           call myrsort_vals(distance,thisdim,sortorder)
           call myrsort_vals(vardots,thisdim,sortorder)
           
#if REALFLAG == 0 
           call mycsort_vecs(laneigvects,thisdim,order*lanblocknum,sortorder)
           call mycsort_vals(values,thisdim,sortorder)
           call mycsort_vals(expects,thisdim,sortorder)
           call mycsort_vals(epredict,thisdim,sortorder)
           call ZGEMM('N','N',lansize,lanblocknum,thisdim,(1d0,0d0),lanvects(:,:,:),lansize,laneigvects,order*lanblocknum,&
                (0d0,0d0),outvectors,lansize)
#else
           call myrsort_vecs(laneigvects,thisdim,order*lanblocknum,sortorder)
           call myrsort_vals(values,thisdim,sortorder)
           call myrsort_vals(expects,thisdim,sortorder)
           call myrsort_vals(epredict,thisdim,sortorder)
           call DGEMM('N','N',lansize,lanblocknum,thisdim,1d0,lanvects(:,:,:),lansize,laneigvects,order*lanblocknum,&
          0d0,outvectors,lansize)
#endif


           do  j=1, lanblocknum
              call indotsub(outvectors(:,j:j),outvectors(:,j:j),lansize,1,1,csums)
              outvectors(:,j)=outvectors(:,j)/sqrt(csums(1,1))
           enddo
           
           call inh0multsub(outvectors(:,:),tempvectors(:,:),lanblocknum,lansize)
           
           
           do  j=1, lanblocknum
              call indotsub(outvectors(:,j:j),tempvectors(:,j:j),lansize,1,1,outvalues(j))
              tempvectors(:,j) = tempvectors(:,j) - outvalues(j) * outvectors(:,j)            
              call indotsub(tempvectors(:,j:j),tempvectors(:,j:j),lansize,1,1,csums)
              error(j)=sqrt(abs(csums(1,1)))
           enddo
     
           if (printflag.ne.0) then
              
              OFL; write(mpifileptr,'(A15,1000F15.10)')        "EIGVALS NOW   ", values(1:thisdim); CFL
              OFL; write(mpifileptr,'(A15,1000F15.10)')        "EPREDICT      ", epredict(1:thisdim); CFL
              OFL; write(mpifileptr,'(A15,1000F15.10)')        "EXPECT0       ", expects(1:thisdim); CFL
              OFL; write(mpifileptr,'(A15,F15.10,1000F30.10)') "DEV           ", sqrt(vardots(1:thisdim)); CFL

!        OFL; write(mpifileptr,'(A25,1000F19.12)')        "EIGVALS NOW           ", values(1:lanblocknum); CFL
!        OFL; write(mpifileptr,'(A25,1000F19.12)')        "EPREDICT              ", epredict(1:lanblocknum); CFL
!        if (lanblocknum.gt.1) then
!           OFL; write(mpifileptr,'(A25,1000F19.12)')        "EXPECT0               ", expects(1:lanblocknum); CFL
!           OFL; write(mpifileptr,'(A25,F19.12,1000F38.12)') "DEV                   ", sqrt(vardots(1:lanblocknum)); CFL
!        endif

!        OFL; write(mpifileptr,'(A25,1000F19.12)')        "EXPECTATION VALS  ", outvalues(:); CFL
!        OFL; write(mpifileptr,'(A25,F19.12,1000F38.12)') "ERRORS            ", error(:); CFL
           endif

        endif  !! IITER EQ 0 OR THISDIM

     enddo !! IITER

     etargettest=etarget

     call myrsort_getorder(error(1:lanblocknum),lanblocknum,0d0,sortorder)
     
     stopsum=error(sortorder(numvects))
     
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
        OFLWR "  Not converged. restarting.", etargettest;
        WRFL  "     ", stopsum,lanthresh,restarts; CFL
        restarts=restarts+1
     endif
     
     initvectors(:,:)=outvectors(:,:)
     
  enddo


contains

#include "shell.h.f90"
  
end subroutine betterharmonic






end module hr_coreblocklanmod




subroutine harmonic_ritz( &
     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
     inprintflag,guessflag,lanthresh,&
     inh0multsub,indotsub,   inetarget,inrange,in_improverange,in_improvethresh, &
     in_invmult_order, in_invexpmult_order, in_expmult_order)

  use hr_coreblocklanmod
  use hr_etargetmod
  use hr_multmod
  implicit none 

  integer,intent(in) :: lansize,maxiter,lanblocknum,inprintflag,order,numvects,guessflag, &
       in_invmult_order, in_invexpmult_order, in_expmult_order
  DATATYPE, intent(in) :: inetarget
  real*8, intent(in) :: lanthresh,inrange,in_improvethresh, in_improverange
  DATATYPE, intent(out) :: outvalues(lanblocknum)  , outvectors(lansize,lanblocknum)
  external :: inh0multsub,indotsub

  etarget=inetarget
  etargettest=etarget
  range=inrange
  improverange=in_improverange
  improvethresh=in_improvethresh
  invmult_order=in_invmult_order
  invexpmult_order=in_invexpmult_order
  expmult_order=in_expmult_order

!  print *, "xCHECK:", &
!     lanblocknum, numvects, lansize,order, maxiter
!  print *, etarget,range,invmult_order,invexpmult_order,expmult_order; stop
  

  call betterharmonic( &
     lanblocknum, numvects, lansize,order, maxiter, outvectors, outvalues,&
     inprintflag,guessflag,lanthresh,&
!     expmultsub,invmultsub3,inh0multsub,indotsub)
     expmultsub,invmultsub,inh0multsub,indotsub)
!     expmultsub,invexpmultsub,inh0multsub,indotsub)
!     expmultsub,inh0multsub,inh0multsub,indotsub)



end subroutine harmonic_ritz
















