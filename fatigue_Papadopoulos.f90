
module Fatigue_Papadopoulos

 contains

  subroutine F2C_Papadopoulos_POLYCRYSTAL(ninc,ngrains,nsys,sig,n,l,frac,alpha,crit,gcrit)  &
             BIND(C,NAME="F2C_Papadopoulos_POLYCRYSTAL")
   use, intrinsic :: ISO_C_BINDING
   use Tensor
   implicit none
   integer(kind=C_INT), intent(in), value :: ninc
   integer(kind=C_INT), intent(in), value :: ngrains
   integer(kind=C_INT), intent(in), value :: nsys
   real(kind=C_DOUBLE), dimension(6,int(ninc),int(ngrains)), intent(in) :: sig
   real(kind=C_DOUBLE), dimension(3,int(nsys),int(ngrains)),  intent(in) :: n
   real(kind=C_DOUBLE), dimension(3,int(nsys),int(ngrains)),  intent(in) :: l
   real(kind=C_DOUBLE), dimension(int(ngrains)),  intent(in) :: frac
   real(kind=C_DOUBLE),  intent(in), value :: alpha
   real(kind=C_DOUBLE),  intent(out) :: crit
   real(kind=C_DOUBLE), dimension(int(ngrains)), intent(out) :: gcrit

   integer :: i,g,k
   type(Tensor1), allocatable :: nor(:,:),dir(:,:)
   type(Tensor2s), allocatable :: sigma(:,:)
   real(8) :: alp,cr
   real(8), allocatable :: gcr(:),fvol(:)

   allocate(nor(int(nsys),int(ngrains)))
   allocate(dir(int(nsys),int(ngrains)))
   allocate(sigma(int(ninc),int(ngrains)))
   allocate(gcr(int(ngrains)))
   allocate(fvol(int(ngrains)))

   do g=1,int(ngrains)
    do i=1,int(nsys)
     nor(i,g)%val(:)=n(:,i,g)
     dir(i,g)%val(:)=l(:,i,g)
    end do
    do i=1,int(ninc)
     sigma(i,g)%val(:)=sig(:,i,g)
    end do
   end do

   alp=alpha
   fvol(:)=frac(:)
   call Papadopoulos_POLYCRYSTAL(sigma,nor,dir,fvol,alp,cr,gcr)
   crit=cr
   gcrit(:)=gcr(:)

  end subroutine F2C_Papadopoulos_POLYCRYSTAL

  subroutine F2C_Papadopoulos_Macro(ninc,sig,alpha,crit)  &
             BIND(C,NAME="F2C_Papadopoulos_Macro")
   use, intrinsic :: ISO_C_BINDING
   use Tensor
   implicit none
   integer(kind=C_INT), intent(in), value :: ninc
   real(kind=C_DOUBLE), dimension(6,int(ninc)), intent(in) :: sig
   real(kind=C_DOUBLE),  intent(in), value :: alpha
   real(kind=C_DOUBLE),  intent(out) :: crit

   integer :: i,g,k
   type(Tensor1), allocatable :: nor(:,:),dir(:,:)
   type(Tensor2s), allocatable :: sigma(:)
   real(8) :: alp,cr

   allocate(sigma(int(ninc)))


    do i=1,int(ninc)
     sigma(i)%val(:)=sig(:,i)
    end do


   alp=alpha
   call Papadopoulos_MACRO(sigma,alp,cr)
   crit=cr

  end subroutine F2C_Papadopoulos_Macro

  
  
  
  subroutine Papadopoulos_POLYCRYSTAL(sigma,nor,dir,fvol,alpha,crit,gcrit)
   use Tensor
   implicit none
   type(Tensor1), intent(in) :: nor(:,:)
   type(Tensor1), intent(in) :: dir(:,:)
   type(Tensor2s), intent(in) :: sigma(:,:)
   real(8), intent(in) :: fvol(:)
   real(8), intent(in) :: alpha
   real(8), intent(out) :: crit
   real(8), intent(out), optional :: gcrit(:)

   integer :: i,j,k,g,nsys,ngrains,ninc
   type(tensor2s) :: schmid
   real(8) :: taumin,taumax,val,tau,taua
   real(8), allocatable :: var(:)
   type(Tensor2s), allocatable :: sigma2(:,:)


   nsys=size(nor,dim=1)
   ngrains=size(sigma,dim=2)
   ninc=size(sigma,dim=1)
   allocate(var(ngrains))






   sigma2=transpose(sigma)
   do i=1,ninc
    do g=1,ngrains
     var(g)=0.0d0
     do j=1,nsys
      !var(g)=var(g)+nor(j,g).dot.sigma(i,g).dot.nor(j,g)
      var(g)=var(g)+((nor(j,g).dot.sigma2(g,i)).dot.nor(j,g))
     end do
     if (present(gcrit)) then
      if (i==1) then
       gcrit(g)=var(g)/real(nsys,8)
      else
       gcrit(g)=max(var(g)/real(nsys,8),gcrit(g))
      end if
     end if
     var(g)=fvol(g)*var(g)
    end do
    val=sum(var)/real(nsys,8)
    if (i==1) then
     crit=val
    else
     crit=max(crit,val)
    end if
   end do





   do g=1,ngrains
    var(g)=0.0d0
    do i=1,nsys
     schmid=sym(dir(i,g)*nor(i,g))
     do j=1,ninc
      tau=schmid.ddot.sigma(j,g)
      if (j==1) then 
       taumax=tau
       taumin=tau
      else
       taumax=max(tau,taumax)
       taumin=min(tau,taumin)
      end if
     end do
     var(g)=var(g)+(0.5d0*(taumax-taumin))**2
    end do
    if (present(gcrit)) then
     gcrit(g)=alpha*gcrit(g)+sqrt(5.0d0*var(g)/real(nsys,8))
    end if
    var(g)=fvol(g)*var(g)
   end do
   crit=alpha*crit+sqrt(5.0d0*sum(var)/real(nsys,8))




  end subroutine Papadopoulos_POLYCRYSTAL



  subroutine Papadopoulos_MACRO(sigma,alpha,crit)
   use Tensor
   use sphere_lebedev_rule
   implicit none
   type(Tensor2s), intent(in) :: sigma(:)
   real(8), intent(in) :: alpha
   real(8), intent(out) :: crit

   integer :: n,m,order,available,k,ndir,i
   real(8), dimension(:), allocatable :: x,y,z,w,theta,phi,NN
   real ( kind = 8 ),  parameter :: pi = 3.14159265358979323846D+00
   real(8) :: xsi,integrale1,integrale2,TT
   type(Tensor1) :: nor,dir
   type(tensor2s) :: schmid
   real(8) :: tau,taumin,taumax,taua,mat(3,3),toto
   
   allocate(NN(size(sigma)))

   !toto=0.0d0

   ndir=100
   n=20!65
   available = available_table ( n )
   if ( available == 1 ) then
    order = order_table ( n )
    !print*, order
    allocate(x(order))
    allocate(y(order))
    allocate(z(order))
    allocate(w(order))
    allocate(theta(order))
    allocate(phi(order))
    call ld_by_order ( order, x, y, z, w )
    do m = 1, order
     call xyz_to_tp ( x(m), y(m), z(m), theta(m), phi(m) )
    end do
    theta(:)=theta(:)*pi/180.0d0
    phi(:)=phi(:)*pi/180.0d0
    deallocate(x)
    deallocate(y)
    deallocate(z)
    TT=0.0d0
    NN(:)=0.0d0
    do m = 1, order
     mat(1,1)=sin(phi(m))*cos(theta(m))
     mat(2,1)=sin(phi(m))*sin(theta(m))
     mat(3,1)=cos(phi(m))
     mat(1,2)=cos(phi(m))*cos(theta(m))
     mat(2,2)=cos(phi(m))*sin(theta(m))
     mat(3,2)=-sin(phi(m))
     mat(1,3)=-sin(theta(m))
     mat(2,3)=cos(theta(m))
     mat(3,3)=0.0d0
     nor%val(1:3)=mat(1:3,1)
     integrale1=0.0d0
     do k=1,ndir
      xsi=2.0D+00*pi*real(k-1,kind=8)/real(ndir,kind=8)
      dir%val(1)=0.0d0
      dir%val(2)=cos(xsi)
      dir%val(3)=sin(xsi)
      dir%val=matmul(mat,dir%val)
      schmid=sym(dir*nor)
      do i=1,size(sigma)
       tau=schmid.ddot.sigma(i)
       if (i==1) then
        taumin=tau
        taumax=tau
       else
        taumin=min(taumin,tau)
        taumax=max(taumax,tau)
       end if
      end do
      taua=0.5d0*(taumax-taumin)
      !toto=max(toto,taua)
      integrale1=integrale1+(1.0D+00/real(ndir,kind=8))*taua*taua
     end do
     TT=TT+w(m)*integrale1
     do i=1,size(sigma)
      NN(i)=NN(i)+w(m)*(nor.dot.sigma(i).dot.nor)
     end do
    end do
    TT=sqrt(5.0d0*TT)

    crit=0.0d0
    do i=1,size(sigma)
     crit=max(crit,TT+alpha*NN(i))
    end do
   end if
   !print*, 'ee', toto

  end subroutine Papadopoulos_MACRO







end module Fatigue_Papadopoulos

