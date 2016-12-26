!***********************************************************************
!*****         1-D Laser-induced Ionization-wave Simulator         *****
!*****         Developed by Rei Kawashima (Univ. of Tokyo)         *****
!*****                 Last Update 2016/06/10                      *****
!***********************************************************************

program main
   use parameters_mod
   use global_mod
   implicit none
   integer :: time_ini,time_end,i
   double precision,dimension(1:NX)     :: nele                         ![m-3]    Electron number density (Total)
   double precision,dimension(1:NX)     :: npho                         ![m-3]    Electron number density (Photoionization)
   double precision,dimension(1:NX)     :: ncol                         ![m-3]    Electron number density (Collisional ionization)
   double precision,dimension(1:NX)     :: uele                         ![ms-1]   Electron velocity
   double precision,dimension(1:NX)     :: Tele                         ![eV]     Electron temperature
   double precision,dimension(1:NX)     :: Ilas                         ![Wm-2]   Laser intensity
   double precision,dimension(1:NX,8)   :: IuvR                         ![Wm-2]   Right-running UV light intensity
   double precision,dimension(1:NX,8)   :: IuvL                         ![Wm-2]   Left-running UV light intensity
   double precision,dimension(1:NX)     :: nneu                         ![m-3]    Neutral number density
   double precision,dimension(1:NX)     :: qcol                         ![m-3s-1] Ion production rate by collisional ionization
   double precision,dimension(1:NX)     :: qpho                         ![m-3s-1] Ion production rate by photoionization
   double precision,dimension(1:NX)     :: qsdm                         ![m-3s-1] Ion production rate by streamer discharge model
   double precision                     :: ugrd                         ![ms-1]   Uniform grid speed (corresponding to ionization wave velocity)
   double precision                     :: xwav                         ![m]      Wave position
   double precision                     :: res_ncol,res_npho            ![-]      Normalized difference of nele
   double precision                     :: res_Tele                     ![-]      Normalized difference of Tele
   double precision                     :: res_nneu                     ![-]      Normalized difference of nneu
   double precision                     :: res_ugrd                     ![-]      Normalized difference of ugrd

   write(*,*) 'Program start...'

   call system_clock(count = time_ini)

   !Initial setups
   call InitialCondition(nele,uele,Tele,nneu,ugrd)

   do i = 1,NX
      ncol(i) = nele(i)
      npho(i) = 0.0d0
   enddo

   MARCH: do it = 1,ITMAX
      if(mod(it,ISMP).eq.0) nstp = nstp+1

      !Radiative transport equation solver
      call RadiativeTransport(nele,nneu,Ilas,IuvR,IuvL)

      !Ion production rate calc.
      call IonProduction(nele,Tele,Ilas,IuvR,IuvL,nneu,qcol,qpho,qsdm)

      !Electron advection-diffusion equation solver
      call ElectronDiffusion(ncol,qsdm,ugrd,res_ncol)
      call ElectronDiffusion(npho,qpho,ugrd,res_npho)
      do i = 1,NX
         nele(i) = ncol(i)+npho(i)
      enddo

      !Energy advection-diffusion equation solver
      call EnergyAdvectionDiffusion(nele,Tele,Ilas,qsdm,ugrd,res_Tele)

      !Neutral steady calculation
      call NeutralSteady(ugrd,nneu,qsdm,res_nneu)

      !Grid speed calculation
      call GridSpeed(nele,ugrd,xwav,res_ugrd)

      !Output distribution
      if(mod(it,ISMP).eq.0) then
         call Output(nele,ncol,npho,Tele,Ilas,IuvR,IuvL,nneu,qcol,qpho,qsdm,ugrd,xwav,res_ncol,res_npho,res_Tele,res_nneu,res_ugrd)

         !Information for monitoring
         write(*,*) '***********************************************'
         write(*,'(A10,I10,A)')   '          ',nstp,'-th Step'
         write(*,'(A20,E12.3,A)') '         Time Step =',DTNOD          ,' [s]'
         write(*,'(A20,E12.3,A)') '        Grid Speed =',ugrd           ,' [m s-1]'
         write(*,'(A20,E12.3,A)') '   Norm. Wave Pos. =',xwav/XL        ,' [-]'
         write(*,'(A20,E12.3,A)') '     ncol Residual =',res_ncol       ,' [-]'
         write(*,'(A20,E12.3,A)') '     npho Residual =',res_npho       ,' [-]'
         write(*,'(A20,E12.3,A)') '     Tele Residual =',res_Tele       ,' [-]'
         write(*,'(A20,E12.3,A)') '     nneu Residual =',res_nneu       ,' [-]'
         write(*,'(A20,E12.3,A)') '     ugrd Residual =',res_ugrd       ,' [-]'

         if(res_ncol.lt.UPS .and. res_Tele.lt.UPS .and. res_nneu.lt.UPS .and. res_ugrd.lt.UPS) exit MARCH
      endif
   enddo MARCH

   call system_clock(count = time_end)

   write(*,'(A44)')         ' The calculation was successfully converged.'
   write(*,'(A21,I11,A)')   ' Calculation time was',time_end-time_ini,' [cpus]'
   write(*,*) 'Program end...'

   stop
end program main


!***********************************************************************
!*****                 Initial condition setting                   *****
!***********************************************************************

subroutine InitialCondition(nele,uele,Tele,nneu,ugrd)

   use parameters_mod
   use global_mod
   implicit none
   character(len=30) :: dname
   integer :: i
   double precision,dimension(1:NX),intent(out)     :: nele             ![m-3]    Electron density
   double precision,dimension(1:NX),intent(out)     :: uele             ![ms-1]   Electron velocity
   double precision,dimension(1:NX),intent(out)     :: Tele             ![eV]     Electron temperature
   double precision,dimension(1:NX),intent(out)     :: nneu             ![m-3]    Neutral number density
   double precision                ,intent(out)     :: ugrd             ![ms-1]   Uniform grid speed

   !Initial condition for variables
   do i = 1,NX
      if(i.le.NX/2) then
         nele(i) = 3.0d24
      else
         nele(i) = 3.0d24*(1.0d0-2.0d0*(DXL*dble(i)-0.5d0*XL)/XL)
      endif
      uele(i) = 0.0d0
      Tele(i) = 2.0d0
      nneu(i) = NNEU0
   enddo

   ugrd = UGRD0

   dname = '/output/datafile/'
   open(unit=12,file=trim(TOPDIR)//trim(dname)//'input.dat',&
      form='formatted',status='unknown')
      do i = 1,NX
         write(12,'(1I5,4E15.5)') &
            i,nele(i),uele(i),Tele(i),nneu(i)
      enddo
   close(12)
   write(*,*) 'Wrote input.dat...'

   return
end subroutine


!***********************************************************************
!*****                    Output file generator                    *****
!***********************************************************************

subroutine Output(nele,ncol,npho,Tele,Ilas,IuvR,IuvL,nneu,qcol,qpho,qsdm,ugrd,xwav,res_ncol,res_npho,res_Tele,res_nneu,res_ugrd)

   use parameters_mod
   use global_mod
   implicit none
   character(len=6)  :: cname
   character(len=30) :: dname
   integer :: i
   double precision,dimension(1:NX)  ,intent(in) :: nele                !Electron density
   double precision,dimension(1:NX)  ,intent(in) :: ncol                !Electron density
   double precision,dimension(1:NX)  ,intent(in) :: npho                !Electron density
   double precision,dimension(1:NX)  ,intent(in) :: Tele                !Electron temperature
   double precision,dimension(1:NX)  ,intent(in) :: Ilas                !Laser intensity
   double precision,dimension(1:NX,8),intent(in) :: IuvR                !Right-running UV light intensity
   double precision,dimension(1:NX,8),intent(in) :: IuvL                !Left-running UV light intensity
   double precision,dimension(1:NX)  ,intent(in) :: nneu                !Neutral number density
   double precision,dimension(1:NX)  ,intent(in) :: qcol                ![m-3s-1] Ion production rate by collisional ionization
   double precision,dimension(1:NX)  ,intent(in) :: qpho                ![m-3s-1] Ion production rate by photoionization
   double precision,dimension(1:NX)  ,intent(in) :: qsdm                ![m-3s-1] Ion production rate by streamer discharge model
   double precision,intent(in)                   :: ugrd                !Uniform grid speed
   double precision,intent(in)                   :: xwav                !Wave position
   double precision,intent(in)                   :: res_ncol            !Normalized difference of nele
   double precision,intent(in)                   :: res_npho            !Normalized difference of nele
   double precision,intent(in)                   :: res_Tele            !Normalized difference of Tele
   double precision,intent(in)                   :: res_nneu            !Normalized difference of nneu
   double precision,intent(in)                   :: res_ugrd            !Normalized difference of ugrd

   write(cname,'(I6.6)') nstp
   dname = '/output/distribution/'
   open(unit=20,file=trim(TOPDIR)//trim(dname)//'distribution.'//&
      cname//'.dat',form='formatted',status='unknown')
      do i = 1,NX
         write(20,'(1I8,6E15.5)') &
            i,nele(i),ncol(i),npho(i),Tele(i),Ilas(i),nneu(i)
      enddo
   close(20)
   open(unit=21,file=trim(TOPDIR)//trim(dname)//'ionization.'//&
      cname//'.dat',form='formatted',status='unknown')
      do i = 1,NX
         write(21,'(1I8,3E15.5)') &
            i,qcol(i),qpho(i),qsdm(i)
      enddo
   close(21)
   open(unit=22,file=trim(TOPDIR)//trim(dname)//'uvlight.'//&
      cname//'.dat',form='formatted',status='unknown')
      do i = 1,NX
         write(22,'(1I8,16E12.3)') &
            i,IuvR(i,1),IuvR(i,2),IuvR(i,3),IuvR(i,4),IuvR(i,5),IuvR(i,6),IuvR(i,7),IuvR(i,8),&
              IuvL(i,1),IuvL(i,2),IuvL(i,3),IuvL(i,4),IuvL(i,5),IuvL(i,6),IuvL(i,7),IuvL(i,8)
      enddo
   close(22)

   dname = '/output/datafile/'
   open(unit=31,file=trim(TOPDIR)//trim(dname)//'history.dat', &
      form='formatted',status='unknown',position='append')
      write(31,'(1I7,7E15.5)') nstp,ugrd,xwav,res_ncol,res_npho,res_Tele,res_nneu,res_ugrd
   close(31)

   return
endsubroutine




!***********************************************************************
!*****                 Neutral Number Density Solver               *****
!*****                    Last Update 2016/01/15                   *****
!***********************************************************************

subroutine NeutralSteady(ugrd,nneu,qion,res_nneu)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i
   double precision,intent(in)                     :: ugrd              !Uniform grid speed
   double precision,dimension(1:NX),intent(inout)  :: nneu              !Neutral number density
   double precision,dimension(1:NX),intent(in)     :: qion              !Ion production rate
   double precision                ,intent(inout)  :: res_nneu          !Normalized difference of nneu
   double precision,dimension(-1:NX+2)             :: var               !Variable
   double precision,dimension(-1:NX+2)             :: adv               !Coefficient for advection term including BC
   !double precision,dimension(1:NX)                :: rhs               !Right-hand side of equations
   double precision,dimension(1:NX)                :: deltaneu          !delta
   double precision,dimension(1:NX)                :: S                 !Source term
   !double precision,dimension(0:NX)                :: Fadv              !Advection flux
   double precision                                :: dt

   !Variable calculation
   do i = 1,NX
      var(i) = nneu(i)
      adv(i) =-ugrd
   enddo

   !Boundary conditions
   var(0)    = 2.0d0*var(1)   -var(2)
   var(-1)   = 2.0d0*var(0)   -var(1)
   var(NX+1) = 2.0d0*NNEU0    -var(NX)
   var(NX+2) = 2.0d0*var(NX+1)-var(NX)
   adv(0)    =-ugrd
   adv(-1)   =-ugrd
   adv(NX+1) =-ugrd
   adv(NX+2) =-ugrd
   !Source term calculation
   do i = 1,NX
      S(i) =-qion(i)
   enddo
   dt = DTNOD

   call Advection1D(NX,dt,DXL,var,adv,S,deltaneu)

   do i = 1,NX
      var(i) = var(i)+deltaneu(i)
   enddo
   !Variable calculation
   do i = 1,NX
      nneu(i) = var(i)
      if(nneu(i).lt.0.0d0) then
         !write(*,*) 'Warning...No neutral particle...',it,i,nneu(i)
         nneu(i) = 0.0d0
      endif
      !if(mod(it,ISMP).eq.0 .and. mod(i,10).eq.0) write(*,*) it,i,nneu(i)
   enddo

   !Calculation of normalized difference
   if(mod(it,ISMP).eq.0) then
      res_nneu = 0.0d0
      do i = 1,NX
         res_nneu = res_nneu+(dabs(deltaneu(i))/(dabs(var(i))+1.0d12))**2.0d0
      enddo
      res_nneu = dsqrt(res_nneu/dble(NX))
   endif


   return
endsubroutine




!***********************************************************************
!*****               Ion production rate calculation               *****
!***********************************************************************

subroutine IonProduction(nele,Tele,Ilas,IuvR,IuvL,nneu,qcol,qpho,qsdm)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,k
   double precision,dimension(1:NX),intent(in)     :: nele              ![m-3]    Electron number density
   double precision,dimension(1:NX),intent(in)     :: Tele              ![eV]     Electron temperature
   double precision,dimension(1:NX),intent(in)     :: Ilas              ![Wm-2]   Laser intensity
   double precision,dimension(1:NX,8),intent(in)   :: IuvR              ![Wm-2]   Right-running UV light intensity
   double precision,dimension(1:NX,8),intent(in)   :: IuvL              ![Wm-2]   Left-running UV light intensity
   double precision,dimension(1:NX),intent(in)     :: nneu              ![m-3]    Neutral number density
   double precision,dimension(1:NX),intent(out)    :: qcol              ![m-3s-1] Ion production rate by collisional ionization
   double precision,dimension(1:NX),intent(out)    :: qpho              ![m-3s-1] Ion production rate by photoionization
   double precision,dimension(1:NX),intent(out)    :: qsdm              ![m-3s-1] Ion production rate by streamer discharge model
   double precision  :: a6,a5,a4,a3,a2,a1,a0,te,rate
   double precision  :: ap3,ap2,ap1,ap0,bp3,bp2,bp1,bp0,cp0,eps
   double precision,dimension(8) :: kphoto

   a6  = 1.7636941284d-05
   a5  =-4.4105941303d-04
   a4  = 3.7169449647d-03
   a3  =-9.7065537389d-03
   a2  = 8.9385677620d-03
   a1  =-2.0345182611d-03
   a0  =-1.0350512923d-04

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
   enddo

   do i = 1,NX
      te = Tele(i)
      if(te.lt.2.0d0) then
         rate = 0.0d0
      else
         rate = 1.0d-14*(a6*te**6.0d0+a5*te**5.0d0+a4*te**4.0d0+a3*te**3.0d0+a2*te**2.0d0+a1*te+a0)
      endif
      qcol(i) = nele(i)*nneu(i)*rate                                    !Collision
      !qsdm(i) = ERATIO*1.0d0/(ECH*EION)*(nele(i)*KABS0*Ilas(i))         !Streamer discharge model
      qsdm(i) = ERATIO*1.0d0/(ECH*EION)*(nele(i)*(KABS0+nele(i)*KABS1)*Ilas(i)) !Streamer discharge model (E-N and E-I IB absorption)
   enddo
   !Photoionization
   do i = 1,NX
      qpho(i) = 0.0d0
      do k = 1,8
         qpho(i) = qpho(i)+1.0d0/(ECH*EION)*(nneu(i)*kphoto(k)*(IuvR(i,k)+IuvL(i,k)))
      enddo
   enddo

   return
endsubroutine


!***********************************************************************
!*****                    Grid speed calculation                   *****
!***********************************************************************

subroutine GridSpeed(nele,ugrd,xwav,res_ugrd)

   use parameters_mod
   use global_mod
   implicit none
   integer :: i,k
   double precision,dimension(1:NX),intent(in)     :: nele              !Electron density
   double precision,intent(inout)                  :: ugrd              !Uniform grid speed
   double precision,intent(inout)                  :: xwav              !Wave position
   double precision,intent(inout)                  :: res_ugrd          !Normalized difference of ugrd
   double precision :: delta,xpreset

   xwavp = xwav
   do k = 1,NX
      i = NX+1-k
      if(nele(i).gt.1.0d20) then
         if(i.eq.NX) then
            xwav = DXL*(dble(i)-0.5d0)
         else
            xwav = DXL*(dble(i)-0.5d0)+DXL*(nele(i)-1.0d20)/(nele(i)-nele(i+1))
         endif
         exit
      endif
   enddo
   if(k.eq.NX+1) xwav = 0.0d0
   xpreset = XL*0.5d0
   delta   = DTNOD*ALPHA0*(xwav-xpreset)+ALPHA1*(xwav-xwavp)
   ugrd    = ugrd+delta

   !Calculation of normalized difference
   if(mod(it,ISMP).eq.0) then
      res_ugrd = dabs(delta)/(dabs(ugrd)+1.0d0)
   endif

   return
endsubroutine




