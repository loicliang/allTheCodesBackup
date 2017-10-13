 module Fatigue_Matake
 
 contains
 
    subroutine Matake_Macro(sigma,alpha,crit)
    use Tensor
    use mod_icosahedron
    use sphere_fibonacci_grid
    implicit none

    type(Tensor2s), intent(in) :: sigma(:)
    real(8), intent(in) :: alpha
    real(8), intent(out) :: crit
    type(Tensor1), allocatable :: nor(:)

    allocate(nor(2000))
    call sphere_fibonacci_grid_points(size(nor),nor)
    call Matake_Homo(sigma,nor,alpha,crit)
    end subroutine Matake_Macro

    subroutine Matake_Homo(sigma,nor,alpha,crit)
    use Tensor
    use MinimumCircumscribedCircle
    implicit none
    type(Tensor1), intent(in) :: nor(:)
    type(Tensor2s), intent(in) :: sigma(:)
    real(8), intent(in) :: alpha
    real(8), intent(out) :: crit

    integer :: i,j,k
    type(tensor1), allocatable :: tau(:)
    type(tensor1) :: taum,sigNor

    real(8) :: taua,maxsignor,maxtaua
    integer :: ninc,nnor,slipplane
   
    ninc=size(sigma)
    nnor=size(nor)
    allocate(tau(ninc))
   
    
    maxtaua=0.0d0
    do i=1,nnor
        do j=1,ninc
            sigNor=(nor(i).dot.sigma(j).dot.nor(i))*nor(i)   
            tau(j)=(sigma(j).dot.nor(i))-sigNor
        end do
        call mcc_randomised_algorithm2(tau,taum,taua)
        if(taua>maxtaua) then
            maxtaua=taua;slipplane=i;
        endif
    end do
    
    maxsignor=0.0d0
    do j=1,ninc
        sigNor=(nor(slipplane).dot.sigma(j).dot.nor(slipplane))*nor(slipplane)
        maxsignor=max(maxsignor,norm(sigNor))
    end do
    
    crit=maxsignor+maxtaua
    
  end subroutine Matake_Homo
  
    subroutine F2C_Matake_Macro(ninc,sig,alpha,crit)  BIND(C,NAME="F2C_Matake_Macro")
    use, intrinsic :: ISO_C_BINDING
    use Tensor
    implicit none
    integer(kind=C_INT), intent(in), value :: ninc
    real(kind=C_DOUBLE), dimension(6,int(ninc)), intent(in) :: sig
    real(kind=C_DOUBLE),  intent(in), value :: alpha
    real(kind=C_DOUBLE),  intent(out) :: crit

    integer :: i
    type(Tensor2s), allocatable :: sigma(:)
    real(8) :: alp,cr

    allocate(sigma(int(ninc)))

    do i=1,int(ninc)
        sigma(i)%val(:)=sig(:,i)
    end do

    alp=alpha
    call Matake_Macro(sigma,alp,cr)
    crit=cr

    end subroutine F2C_Matake_Macro
  
    subroutine F2C_Matake_CRYSTAL(ninc,nsys,sig,n,l,alpha,crit)  BIND(C,NAME="F2C_Matake_CRYSTAL")
    use, intrinsic :: ISO_C_BINDING
    use Tensor
    implicit none
    integer(kind=C_INT), intent(in), value :: ninc
    integer(kind=C_INT), intent(in), value :: nsys
    real(kind=C_DOUBLE), dimension(6,int(ninc)), intent(in) :: sig
    real(kind=C_DOUBLE), dimension(3,int(nsys)),  intent(in) :: n
    real(kind=C_DOUBLE), dimension(3,int(nsys)),  intent(in) :: l
    real(kind=C_DOUBLE),  intent(in), value :: alpha
    real(kind=C_DOUBLE),  intent(out) :: crit

    integer :: i,j,k
    type(Tensor1), allocatable :: nor(:),dir(:)
    type(Tensor2s), allocatable :: sigma(:)
    real(8) :: alp,cr

    allocate(nor(int(nsys)))
    allocate(dir(int(nsys)))
    allocate(sigma(int(ninc)))
    do i=1,int(nsys)
        nor(i)%val(:)=n(:,i)
        dir(i)%val(:)=l(:,i)
    end do

    do i=1,int(ninc)
        sigma(i)%val(:)=sig(:,i)
    end do

    alp=alpha
    call Matake_CRYSTAL(sigma,nor,dir,alp,cr)
    crit=cr
   
   end subroutine F2C_Matake_CRYSTAL

    subroutine Matake_CRYSTAL(sigma,nor,dir,alpha,crit)
    use Tensor
    implicit none

    type(Tensor1), intent(in) :: nor(:)
    type(Tensor1), intent(in) :: dir(:)
    type(Tensor2s), intent(in) :: sigma(:)
    real(8), intent(in) :: alpha
    real(8), intent(out) :: crit

    integer :: i,j,k
    type(tensor2s), allocatable :: schmid(:)
    real(8), allocatable :: taum(:),taumin(:),taumax(:)
    real(8) :: tau,shydro,val

    integer:: ninc,nsys

!     ninc=size(sigma)
!     nsys=size(dir)
! 
!     allocate(schmid(nsys))
!     do i=1,nsys
!         schmid(i)=sym(dir(i)*nor(i))
!     end do
! 
! 
!     allocate(taum(nsys))
!     allocate(taumin(nsys))
!     allocate(taumax(nsys))
! 
!     do i=1,ninc
!         do j=1,nsys
!         tau=schmid(j).ddot.sigma(i)
!         if (i==1) then
!             taumin(j)=tau
!             taumax(j)=tau
!         else
!             taumin(j)=min(taumin(j),tau)
!             taumax(j)=max(taumax(j),tau)
!         end if
!         end do
!     end do
!     taum(:)=0.5d0*(taumax(:)+taumin(:))
!     !print*, taum
!     crit=0.0d0
!     do i=1,ninc
!         shydro=sigma(i)%I1()/3.0d0
!         do j=1,nsys
!             tau=schmid(j).ddot.sigma(i)
!      !if (i==25) print*, tau
!             val=abs(tau-taum(j))+alpha*shydro
!             if ((i==1).and.(j==1)) then
!                 crit=val
!             else
!                 crit=max(crit,val)
!             end if
!         end do
!     end do
    
    end subroutine Matake_CRYSTAL


  
end module Fatigue_Matake  
