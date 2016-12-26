
!---------- 1-D Time Dependent Diffusion equation solver ---------------
!-------------------- Last Update 2015/10/20 ---------------------------

subroutine ElectronDiffusion(nele,qion,ugrd,res_nele)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i
   double precision,dimension(1:NX),intent(inout) :: nele               ![m-3]    Electron density
   double precision,dimension(1:NX),intent(in)    :: qion               ![m-3s-1] Ion production rate
   double precision                ,intent(inout) :: res_nele           ![-]      Normalized difference (residual)
   double precision                ,intent(in)    :: ugrd               ![ms-1]   Uniform grid speed
   double precision,dimension(-1:NX+2)            :: var                !Variable including BC
   double precision,dimension(-1:NX+2)            :: adv                !Coefficient for advection term including BC
   double precision,dimension(0:NX+1)             :: dif                !Coefficient for diffusion term including BC
   double precision,dimension(1:NX)               :: sce                !Source term
   double precision,dimension(1:NX)               :: deltaele           !Delta of electron number density
   double precision                               :: dt                 !Time step interval

   !Coefficient setting
   do i = 1,NX
      adv(i) =-ugrd
      var(i) = nele(i)
      dif(i) = DAMB0
      sce(i) = qion(i)
   enddo

   !Boundary condition
   var(0)    = 2.0d0*var(1)   -var(2)
   var(NX+1) = 2.0d0*NELE0    -var(NX)!*dexp(-dble(it)**2.0d0/(0.1d0*ITMAX)**2.0d0)
   var(-1)   = 2.0d0*var(0)   -var(1)
   var(NX+2) = 2.0d0*var(NX+1)-var(NX)
   dif(0)    = dif(1)
   dif(NX+1) = dif(NX)
   adv(0)    = adv(1)
   adv(NX+1) = adv(NX)
   adv(-1)   = 2.0d0*adv(0)-adv(1)
   adv(NX+2) = 2.0d0*adv(NX+1)-adv(NX)
   dt = DTNOD

   !Call Advection-Diffusion Solver
   call AdvectionDiffusion1D(NX,dt,DXL,var,adv,dif,sce,deltaele)

   !Update variable
   do i = 1,NX
      var(i) = var(i)+deltaele(i)
   enddo

   !Update electron number density (which must be greater than 0)
   do i = 1,NX
      nele(i) = dmax1(var(i),0.0d0)
   enddo

   !Calculate normalized difference
   if(mod(it,ISMP).eq.0) then
      res_nele = 0.0d0
      do i = 1,NX
         res_nele = res_nele+(dabs(deltaele(i))/(dabs(var(i))+1.0d15))**2.0d0
      enddo
      res_nele = dsqrt(res_nele/dble(NX))
   endif

   return
endsubroutine







