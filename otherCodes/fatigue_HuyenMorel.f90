
module Fatigue_HuyenMorel


 contains


  subroutine F2C_HuyenMorel_POLYCRYSTAL(ninc,ngrains,nsys,sig,n,l,alpha,gamma,tau0p,m,crit,gcrit)  &
             BIND(C,NAME="F2C_HuyenMorel_POLYCRYSTAL")
   use, intrinsic :: ISO_C_BINDING
   use Tensor
   implicit none
   integer(kind=C_INT), intent(in), value :: ninc
   integer(kind=C_INT), intent(in), value :: ngrains
   integer(kind=C_INT), intent(in), value :: nsys
   real(kind=C_DOUBLE), dimension(6,int(ninc),int(ngrains)), intent(in) :: sig
   real(kind=C_DOUBLE), dimension(3,int(nsys),int(ngrains)),  intent(in) :: n
   real(kind=C_DOUBLE), dimension(3,int(nsys),int(ngrains)),  intent(in) :: l
   real(kind=C_DOUBLE),  intent(in), value :: alpha,gamma,tau0p
   integer(kind=C_INT),  intent(in), value :: m
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

   do g=1,int(ngrains)
    do i=1,int(nsys)
     nor(i,g)%val(:)=n(:,i,g)
     dir(i,g)%val(:)=l(:,i,g)
    end do
    do i=1,int(ninc)
     sigma(i,g)%val(:)=sig(:,i,g)
    end do
   end do


   call HuyenMorel_POLYCRYSTAL(sigma,nor,dir,real(alpha,8),real(gamma,8),real(tau0p,8),int(m),cr,gcr)

   crit=cr
   gcrit(:)=gcr(:)

  end subroutine F2C_HuyenMorel_POLYCRYSTAL







  subroutine F2C_HuyenMorel_MACRO(ninc,nelt,sig,volume,alpha,gamma,tau0p,m,crit,ecrit)  &
             BIND(C,NAME="F2C_HuyenMorel_MACRO")
   use, intrinsic :: ISO_C_BINDING
   use Tensor
   implicit none
   integer(kind=C_INT), intent(in), value :: ninc
   integer(kind=C_INT), intent(in), value :: nelt
   real(kind=C_DOUBLE), dimension(6,int(ninc),int(nelt)), intent(in) :: sig
   real(kind=C_DOUBLE), dimension(int(nelt)), intent(in) :: volume
   real(kind=C_DOUBLE),  intent(in), value :: alpha,gamma,tau0p
   integer(kind=C_INT), intent(in), value :: m
   real(kind=C_DOUBLE), intent(out) :: crit
   real(kind=C_DOUBLE), dimension(int(nelt)), intent(out) :: ecrit

   type(Tensor2s), allocatable :: sigma(:,:)
   real(8) :: PF
   real(8), allocatable :: PFe(:)
   integer :: i,j


   allocate(sigma(int(ninc),int(nelt)))
   allocate(PFe(int(nelt)))

   do i=1,nelt
    do j=1,ninc
     sigma(j,i)%val(:)=sig(:,j,i)
    end do
   end do
   call HuyenMorel_MACRO(sigma,real(volume,8),real(alpha,8),real(gamma,8),real(tau0p,8),int(m),PF,PFe)
   crit=PF
   ecrit(:)=PFe(:)


 end subroutine F2C_HuyenMorel_MACRO








  subroutine HuyenMorel_POLYCRYSTAL(sigma,nor,dir,alpha,gamma,tau0p,m,crit,gcrit)
   use Tensor
   implicit none
   type(Tensor1), intent(in) :: nor(:,:)
   type(Tensor1), intent(in) :: dir(:,:)
   type(Tensor2s), intent(in) :: sigma(:,:)
   real(8), intent(in) :: alpha,gamma,tau0p
   integer, intent(in) :: m
   real(8), intent(out) :: crit
   real(8), intent(out), optional :: gcrit(:)

   integer :: i,j,k,g,nsys,ngrains,ninc
   type(tensor2s), allocatable :: schmid(:)
   real(8) :: var,tau,sn
   real(8), allocatable, dimension(:) :: taumin,taumax,snmin,snmax,taua,snm,sna,PFG


   nsys=size(nor,dim=1)
   ngrains=size(sigma,dim=2)
   ninc=size(sigma,dim=1)
   allocate(PFG(ngrains))





   allocate(schmid(nsys))
   allocate(taua(nsys))
   allocate(taumin(nsys))
   allocate(taumax(nsys))
   allocate(sna(nsys))
   allocate(snm(nsys))
   allocate(snmin(nsys))
   allocate(snmax(nsys))

   do g=1,ngrains
    do i=1,nsys
     schmid(i)=sym(dir(i,g)*nor(i,g))
    end do

    do i=1,ninc
     do j=1,nsys
      tau=schmid(j).ddot.sigma(i,g)
      sn=nor(j,g).dot.sigma(i,g).dot.nor(j,g)
      if (i==1) then
       taumin(j)=tau
       taumax(j)=tau
       snmin(j)=sn
       snmax(j)=sn
      else
       taumin(j)=min(taumin(j),tau)
       taumax(j)=max(taumax(j),tau)
       snmin(j)=min(snmin(j),sn)
       snmax(j)=max(snmin(j),sn)
      end if
     end do
    end do
    taua(:)=0.5d0*(taumax(:)-taumin(:))
    snm(:)=0.5d0*(taumax(:)+taumin(:))
    sna(:)=0.5d0*(snmax(:)-snmin(:))

    do i=1,nsys
     if (abs(taua(i))>1.0d-4)then
      var=tau0p*(1.0d0-gamma*snm(i))/(1.0d0+alpha*(sna(i)/taua(i)))
     else
      var=0.0d0
     end if

     !print*, var
     var=taua(i)/var
     !print*, var,exp(-var**m),taua(i)
     !print*


  !print*, 'T0', tau0p*(1.0d0-gamma*snm(i))/(1.0d0+alpha*(sna(i)/taua(i)))
  !print*, 'Ta', taua(i)
     !end if
   
     !var=taua(i)+alpha*sna(i)
     !if (abs(var)<1.0d-1) then
     ! var=0.0d0
     !else
     ! var=tau0p*(1.0d0-gamma*snm(i))/var
     !end if
     !print*, var,var**m
     if (var>0.0d0) then
      var=1.0d0-exp(-var**m)
     else
      var=1.0d0
     end if
     !print*, var
     !print*, '--'
     if (i==1) then
      PFG(g)=var
     else
      PFG(g)=max(var,PFG(g))
     end if
    end do
    !print*, PFG(g)
   end do
   if (present(gcrit)) gcrit(:)=PFG(:)
   PFG(:)=1.0d0-PFG(:)
   crit=1.0d0-product(PFG)
   !print*, 'aa', crit



  end subroutine HuyenMorel_POLYCRYSTAL


  subroutine HuyenMorel_MACRO(sigma,volume,alpha,gamma,tau0p,m,PF,PFe)
   use Tensor
   use sphere_lebedev_rule
   implicit none

   type(Tensor2s), intent(in) :: sigma(:,:)
   real(8), intent(in) :: volume(:)
   real(8), intent(in) :: alpha,gamma,tau0p
   integer, intent(in) :: m
   real(8), intent(out) :: PF
   real(8), intent(out), optional :: PFe(:)

   real ( kind = 8 ),  parameter :: pi = 3.14159265358979323846D+00

   type(Tensor1) :: nor, dir
   type(Tensor2s) :: schmid
   real(8) :: sn,snmin,snmax,sna,snm,tau,taumin,taumax,taua,xsi
   real(8) :: mat(3,3),integrale1,var
   integer :: i,k,e,n,ndir,order,available
   real(8), allocatable :: x(:),y(:),z(:),w(:),theta(:),phi(:)









   ndir=100
   n=20!65
   available = available_table ( n )
   write(*,*) "Total Volume",sum(volume(:))
   if ( available == 1 ) then
   write(*,*) "HuyenMorel Analysis Begins"
    order = order_table ( n )
!write(*,*) "order value",order
    allocate(x(order))
    allocate(y(order))
    allocate(z(order))
    allocate(w(order))
    allocate(theta(order))
    allocate(phi(order))
    call ld_by_order ( order, x, y, z, w )
    do n = 1, order
     call xyz_to_tp ( x(n), y(n), z(n), theta(n), phi(n) )
    end do
    theta(:)=theta(:)*pi/180.0d0
    phi(:)=phi(:)*pi/180.0d0
    deallocate(x)
    deallocate(y)
    deallocate(z)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do n=1,order
!write(*,*) w(n)
!end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

    PF=0.0d0
    do e=1,size(sigma,dim=2)
     do n = 1, order
      mat(1,1)=sin(phi(n))*cos(theta(n))
      mat(2,1)=sin(phi(n))*sin(theta(n))
      mat(3,1)=cos(phi(n))
      mat(1,2)=cos(phi(n))*cos(theta(n))
      mat(2,2)=cos(phi(n))*sin(theta(n))
      mat(3,2)=-sin(phi(n))
      mat(1,3)=-sin(theta(n))
      mat(2,3)=cos(theta(n))
      mat(3,3)=0.0d0
      nor%val(1:3)=mat(1:3,1)
      integrale1=0.0d0

      do i=1,size(sigma,dim=1)
       sn=nor.dot.sigma(i,e).dot.nor
       if (i==1) then
        snmin=sn
        snmax=sn
       else
        snmin=min(snmin,sn)
        snmax=max(snmax,sn)
       end if
      end do
      sna=0.5d0*(snmax-snmin)
      snm=0.5d0*(snmax+snmin)


      do k=1,ndir
       xsi=2.0D+00*pi*real(k-1,kind=8)/real(ndir,kind=8)
       dir%val(1)=0.0d0
       dir%val(2)=cos(xsi)
       dir%val(3)=sin(xsi)
       dir%val=matmul(mat,dir%val)
       schmid=sym(dir*nor)


       do i=1,size(sigma,dim=1)
        tau=schmid.ddot.sigma(i,e)
        if (i==1) then
         taumin=tau
         taumax=tau
        else
         taumin=min(taumin,tau)
         taumax=max(taumax,tau)
        end if
       end do
       taua=0.5d0*(taumax-taumin)

       if (abs(taua)>1.0d-5)then
        var=tau0p*(1.0d0-gamma*snm)/(1.0d0+alpha*(sna/taua))
        var=taua/var
       else
        var=alpha*sna/(tau0p*(1.0d0-gamma*snm))
       end if

      
       if (var>0.0d0) then
        var=var**m
       else
        var=0.0d0
       end if



       !toto=max(toto,taua)
       integrale1=integrale1+(1.0D+00/real(ndir,kind=8))*var
      end do
      PF=PF+w(n)*integrale1*volume(e)
      if (present(PFe)) PFe(e)=PFe(e)+w(n)*integrale1
     end do
    end do
!write(*,*) PF    
    PF=1.0d0-exp(-PF/sum(volume(:)))
   end if

   if (present(PFe)) PFe(:)=1.0d0-exp(-PFe(:))

   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (present(PFe)) then
!write(*,*) "presented"
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
   write(*,*) "HuyenMorel Analysis Ends"






  end subroutine HuyenMorel_MACRO







end module Fatigue_HuyenMorel

