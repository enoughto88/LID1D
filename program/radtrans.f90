!***********************************************************************
!*****             Radiative Transport Equation Solver             *****
!*****         Input: nele,nneu   Output : Ilas,IuvR,IuvL          *****
!***********************************************************************

subroutine RadiativeTransport(nele,nneu,Ilas,IuvR,IuvL)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,k
   double precision,dimension(1:NX),intent(in)    :: nele               !Electron number density
   double precision,dimension(1:NX),intent(in)    :: nneu               !Neutral number density
   double precision,dimension(1:NX),intent(out)   :: Ilas               !Laser intensity
   double precision,dimension(1:NX,8),intent(out) :: IuvR               !Right-running UV light intensity
   double precision,dimension(1:NX,8),intent(out) :: IuvL               !Left-running UV light intensity
   double precision,dimension(1:NX)               :: Ieps
   double precision,dimension(1:NX)               :: aa_Ilas,aa_IuvR,aa_IuvL
   double precision,dimension(1:NX)               :: bb_Ilas,bb_IuvR,bb_IuvL
   double precision              :: ap3,ap2,ap1,ap0,bp3,bp2,bp1,bp0,cp0,eps
   double precision,dimension(8) :: kphoto,kexp

   ap3 = 4.144293d-04
   ap2 =-4.231997d-02
   ap1 = 1.223061d+00
   ap0 =-7.112518d+00
   bp3 =-5.528983d-04
   bp2 = 8.594670d-02
   bp1 =-4.381191d+00
   bp0 = 7.360449d+01
   cp0 = 1.101049d-01

   do k = 1,8
      if(k.eq.1) then
         eps = EION
      else
         eps = 20.0d0+dble(k-2)*5.0d0
      endif
      if    (eps.gt.15.0d0 .and. eps.le.37.0d0) then
         kphoto(k) = 1.0d-21*(ap3*eps**3.0d0+ap2*eps**2.0d0+ap1*eps+ap0)
      elseif(eps.gt.37.0d0 .and. eps.le.45.0d0) then
         kphoto(k) = 1.0d-21*(bp3*eps**3.0d0+bp2*eps**2.0d0+bp1*eps+bp0)
      elseif(eps.gt.45.0d0) then
         kphoto(k) = 1.0d-21*(cp0)
      endif
      !if(k.eq.1) then
      !   kexp(k) = dexp(-EION/2.0d0)-dexp(-20.0d0/2.0d0)
      !elseif(k.ge.2 .and. k.le.7) then
      !   kexp(k) = dexp(-eps/2.0d0)-dexp(-(eps+5.0d0)/2.0d0)
      !elseif(k.eq.8) then
      !   kexp(k) = dexp(-50.0d0/2.0d0)
      !endif
      if(k.eq.1) then
         kexp(k) = dexp(-EION/TeUV)-dexp(-20.0d0/TeUV)
      elseif(k.ge.2 .and. k.le.7) then
         kexp(k) = dexp(-eps/TeUV)-dexp(-(eps+5.0d0)/TeUV)
      elseif(k.eq.8) then
         kexp(k) = dexp(-50.0d0/TeUV)
      endif
   enddo

   !Ilas is chosen for variable
   do i = 1,NX
      !aa_Ilas(i) =-nele(i)*KABS0
      aa_Ilas(i) =-nele(i)*(KABS0+nele(i)*KABS1)
      bb_Ilas(i) = 0.0d0
   enddo
   call RTESolver(NX,DXL,'L',aa_Ilas,bb_Ilas,ILAS0,Ilas)

   do k = 1,8
      do i = 1,NX
         aa_IuvR(i) =-nneu(i)*kphoto(k)
         bb_IuvR(i) = 3.49d-37*nele(i)*nele(i)*dsqrt(TeUV)*kexp(k)
      enddo
      call RTESolver(NX,DXL,'R',aa_IuvR,bb_IuvR,0.0d0,Ieps)
      do i = 1,NX
         IuvR(i,k) = Ieps(i)
      enddo

      do i = 1,NX
         aa_IuvL(i) =-nneu(i)*kphoto(k)
         bb_IuvL(i) = 3.49d-37*nele(i)*nele(i)*dsqrt(TeUV)*kexp(k)
      enddo
      call RTESolver(NX,DXL,'L',aa_IuvL,bb_IuvL,0.0d0,Ieps)
      do i = 1,NX
         IuvL(i,k) = Ieps(i)
      enddo
   enddo

   return
endsubroutine


!************************ RTE Solver ***********************************

subroutine RTESolver(NX,dx,direction,aa,bb,Ilas0,Ilas)

   implicit none
   integer :: i,k
   integer                            ,intent(in)   :: NX               !Number of cells
   double precision                   ,intent(in)   :: dx               !Grid sizing
   character(len=1)                   ,intent(in)   :: direction        !Beam direction; R:right, L:left
   double precision,dimension(1:NX)   ,intent(in)   :: aa               !Source term coefficient
   double precision,dimension(1:NX)   ,intent(in)   :: bb               !Source term coefficient
   double precision                   ,intent(in)   :: Ilas0            !Light intensity at boundary
   double precision,dimension(1:NX)   ,intent(out)  :: Ilas             !Light intensity
   double precision,dimension(-1:NX+2)              :: var              !Variable including ghost cells
   double precision,dimension(1:NX,3)               :: coef             !Coefficients for direct method


   if(direction.eq.'R') then
      !Boundary condition
      var(0) = Ilas0
      if(dabs(1.0d0-dx*aa(1)).lt.1.0d-6) then
         write(*,*) 'Warning in RTE...',(1.0d0-dx*aa(1))
         var(-1) = var(0)
      else
         var(-1) = (1.0d0-2.0d0*dx*aa(1))/(1.0d0-dx*aa(1))*var(0)-dx*bb(1)/(1.0d0-dx*aa(1))
      endif

      !Coefficient calculation
      do i = 1,NX
         if(dabs(3.0d0-2.0d0*dx*aa(i)).lt.1.0d-6) then
            write(*,*) 'Warning in RTE...',(3.0d0-2.0d0*dx*aa(i))
            coef(i,1) = 1.0d0/(1.0d0-dx*aa(i))
            coef(i,2) = 0.0d0
            coef(i,3) = dx*bb(i)/(1.0d0-dx*aa(i))
         else
            coef(i,1) = 4.0d0/(3.0d0-2.0d0*dx*aa(i))
            coef(i,2) =-1.0d0/(3.0d0-2.0d0*dx*aa(i))
            coef(i,3) = 2.0d0*dx*bb(i)/(3.0d0-2.0d0*dx*aa(i))
         endif
      enddo

      !Direct method for radiative transfer equation
      do i = 1,NX
         var(i) = coef(i,1)*var(i-1)+coef(i,2)*var(i-2)+coef(i,3)
      enddo
      do i = 1,NX
         Ilas(i) = dmax1(var(i),0.0d0)
      enddo
   else if(direction.eq.'L') then
      !Boundary condition
      var(NX+1) = Ilas0
      if(dabs(1.0d0-dx*aa(NX)).lt.1.0d-6) then
         write(*,*) 'Warning in RTE...',(1.0d0-dx*aa(NX))
         var(NX+2) = var(NX+1)
      else
         var(NX+2) = (1.0d0-2.0d0*dx*aa(NX))/(1.0d0-dx*aa(NX))*var(NX+1)-dx*bb(NX)/(1.0d0-dx*aa(NX))
      endif

      !Coefficient calculation
      do i = 1,NX
         if(dabs(3.0d0-2.0d0*dx*aa(i)).lt.1.0d-6) then
            write(*,*) 'Warning in RTE...',(3.0d0-2.0d0*dx*aa(i))
            coef(i,1) = 1.0d0/(1.0d0-dx*aa(i))
            coef(i,2) = 0.0d0
            coef(i,3) = dx*bb(i)/(1.0d0-dx*aa(i))
         else
            coef(i,1) = 4.0d0/(3.0d0-2.0d0*dx*aa(i))
            coef(i,2) =-1.0d0/(3.0d0-2.0d0*dx*aa(i))
            coef(i,3) = 2.0d0*dx*bb(i)/(3.0d0-2.0d0*dx*aa(i))
         endif
      enddo

      !Direct method for radiative transfer equation
      do k = 1,NX
         i = NX+1-k
         var(i) = coef(i,1)*var(i+1)+coef(i,2)*var(i+2)+coef(i,3)
      enddo
      do i = 1,NX
         Ilas(i) = dmax1(var(i),0.0d0)
      enddo
   else
      write(*,*) 'Error in RTE direction...',direction
      stop
   endif

   return
endsubroutine



