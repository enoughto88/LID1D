
!---------- 1-D Time Dependent Diffusion equation solver ---------------
!-------------------- Last Update 2015/10/20 ---------------------------

subroutine EnergyAdvectionDiffusion(nele,Tele,Ilas,qsdm,ugrd,res_Tele)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i
   double precision,dimension(1:NX),intent(in)    :: nele               ![m-3]    Electron density
   double precision,dimension(1:NX),intent(inout) :: Tele               ![eV]     Electron temperature
   double precision,dimension(1:NX),intent(in)    :: Ilas               ![Wm-2]   Laser intensity
   double precision,dimension(1:NX),intent(in)    :: qsdm               ![m-3s-1] Ion production rate by streamer discharge model
   double precision                ,intent(inout) :: res_Tele           ![-]      Normalized difference (residual)
   double precision,intent(in)                    :: ugrd               ![ms-1]   Uniform grid speed
   double precision,dimension(-1:NX+2)            :: var                !Variable including BC
   double precision,dimension(-1:NX+2)            :: adv                !Coefficient for advection term including BC
   double precision,dimension(0:NX+1)             :: dif                !Coefficient for diffusion term including BC
   double precision,dimension(1:NX)               :: sce                !Source term
   double precision,dimension(1:NX)               :: deltaene           !Delta of electron internal energy
   double precision                               :: dt                 !Time step interval

   !Coefficient setting
   do i = 1,NX
      if     (i.eq.1) then
         adv(i) =-ugrd+5.0d0/3.0d0*(DELE0-DAMB0)*(dlog(nele(i+1))-dlog(nele(i  )))/DXL
      else if(i.ge.2 .and. i.le.NX-1) then
         adv(i) =-ugrd+5.0d0/3.0d0*(DELE0-DAMB0)*(dlog(nele(i+1))-dlog(nele(i-1)))/2.0d0/DXL
      else if(i.eq.NX) then
         adv(i) =-ugrd+5.0d0/3.0d0*(DELE0-DAMB0)*(dlog(nele(i  ))-dlog(nele(i-1)))/DXL
      endif
      var(i) = 3.0d0/2.0d0*ECH*nele(i)*Tele(i)
      dif(i) = 5.0d0/3.0d0*DELE0
      sce(i) = nele(i)*KABS0*Ilas(i)-ECH*EION*qsdm(i)!-2.0d0*nele(i)*nele(i)*ALPHA0!-ECH*EION*qion(i)
      !sce(i) = nele(i)*KABS0*Ilas(i)-ECH*EION*qion(i)!-2.0d0*nele(i)*nele(i)*ALPHA0!-ECH*EION*qion(i)
   enddo

   !Boundary condition
   var(0)    = var(1)
   var(-1)   = 2.0d0*var(0)-var(1)
   var(NX+1) = 2.0d0*3.0d0/2.0d0*ECH*nele(NX)*TELE0-var(NX)
   var(NX+2) = 2.0d0*var(NX+1)-var(NX)
   dif(0)    = dif(1)
   dif(NX+1) = dif(NX)
   adv(0)    = adv(1)
   adv(NX+1) = adv(NX)
   adv(-1)   = 2.0d0*adv(0)-adv(1)
   adv(NX+2) = 2.0d0*adv(NX+1)-adv(NX)
   dt = DTNOD

   !Call Advection-Diffusion Solver
   call AdvectionDiffusion1D(NX,dt,DXL,var,adv,dif,sce,deltaene)

   !Update variable
   do i = 1,NX
      var(i) = var(i)+deltaene(i)
   enddo

   !Update electron temperature (which must be greater than 0)
   do i = 1,NX
      Tele(i) = dmax1(var(i)/(3.0d0/2.0d0*ECH*nele(i)),0.0d0)
   enddo

   !Calculate normalized difference
   if(mod(it,ISMP).eq.0) then
      res_Tele = 0.0d0
      do i = 1,NX
         res_Tele = res_Tele+(dabs(deltaene(i))/(dabs(var(i))+ECH*1.0d12))**2.0d0
      enddo
      res_Tele = dsqrt(res_Tele/dble(NX))
   endif

   return
endsubroutine






