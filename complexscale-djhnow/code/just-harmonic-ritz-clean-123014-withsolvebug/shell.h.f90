
!!include me

subroutine h0multsub_herm_shell(in,out,howmany,checksize) 
  implicit none
  integer, intent(in) :: howmany,checksize
  complex*16 :: in(checksize,howmany),out(checksize,howmany), realtemp(checksize,howmany),&
       realtemp2(checksize,howmany)

!! hermitian part is real part.
  realtemp(:,:)=real(in(:,:),8)  !! ok conv

  call inh0multsub(realtemp,realtemp2,howmany,checksize)

  out(:,:)=real(realtemp2(:,:),8)  ! ok conv
  
  realtemp(:,:)=real(in(:,:)/(0d0,1d0),8)  !! ok conv
  
  call inh0multsub(realtemp,realtemp2,howmany,checksize)
  
  out(:,:)=out(:,:)+(0d0,1d0)*real(realtemp2(:,:),8)

end subroutine h0multsub_herm_shell



subroutine h0multsub_antiherm_shell(in,out,howmany,checksize)
  implicit none
  integer, intent(in) :: howmany,checksize
  complex*16 :: in(checksize,howmany),out(checksize,howmany), realtemp(checksize,howmany),&
       realtemp2(checksize,howmany)

!! antihermitian part is imag part.
  realtemp(:,:)=real(in(:,:),8)  !! ok conv

  call inh0multsub(realtemp,realtemp2,howmany,checksize)


  out(:,:)=real(realtemp2(:,:)/(0d0,1d0),8)  ! ok conv

  
  realtemp(:,:)=real(in(:,:)/(0d0,1d0),8)  !! ok conv
  
  call inh0multsub(realtemp,realtemp2,howmany,checksize)


  out(:,:)=out(:,:)+(0d0,1d0)*real(realtemp2(:,:)/(0d0,1d0),8)


end subroutine h0multsub_antiherm_shell
