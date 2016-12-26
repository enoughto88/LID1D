
!----------------- 1-D Advection Equation Solver -----------------------
!      Input : NX,dt,dx,u(-1:NX+2),C(-1:NX+2),S(1,NX)
!      Output: delta(1:NX)
!-----------------------------------------------------------------------

subroutine Advection1D(NX,dt,dx,u,C,S,delta)

   implicit none
   integer :: i
   integer                            ,intent(in)    :: NX              !Number of cells
   double precision                   ,intent(in)    :: dt              !Time interval
   double precision                   ,intent(in)    :: dx              !Space interval
   double precision,dimension(-1:NX+2),intent(in)    :: u               !Variable including BC
   double precision,dimension(-1:NX+2),intent(in)    :: C               !Advection coefficient including BC
   double precision,dimension(1:NX)   ,intent(in)    :: S               !Source term
   double precision,dimension(1:NX)   ,intent(out)   :: delta           !Delta of variable (n+1 - n)
   double precision,dimension(1:NX)                  :: rhs             !Right-hand side in explicit form
   double precision,dimension(1:NX,5)                :: ccoef
   double precision,dimension(1:NX,5)                :: coef,LU

   !Coefficient calculation of advection term
   call AdvectionCoefficient1D(NX,dt,dx,C,u,ccoef)

   !Element calculation in 5-band matrix
   do i = 1,NX
      coef(i,1) = ccoef(i,1)
      coef(i,2) = ccoef(i,2)
      coef(i,3) = ccoef(i,3)+1.0d0
      coef(i,4) = ccoef(i,4)
      coef(i,5) = ccoef(i,5)
   enddo

   !Calculate right-hand side
   do i = 1,NX
      rhs(i)    =-(coef(i,1)*u(i-2)+coef(i,2)*u(i-1)+(coef(i,3)-1.0d0)*u(i)&
                  +coef(i,4)*u(i+1)+coef(i,5)*u(i+2))+dt*S(i)
   enddo

   !Matrix inversion calculation using LU factorization for 5-band matrix
   ! 2*W+1 = 5  =>  W = 2
   call LUfactorization(NX,2,coef,LU)
   call ForBackSubstitution(NX,2,LU,rhs,delta)

   return
endsubroutine



!-------------------- 1D Diffusion Equation Solver ---------------------
!          Input : NX,dt,dx,u(0:NX+1),D(0:NX+1),S(1,NX)
!          Output: delta(1:NX)
!-----------------------------------------------------------------------

subroutine Diffusion1D(NX,dt,dx,u,D,S,delta)

   implicit none
   integer :: i
   integer,intent(in)                                :: NX              !Number of cells
   double precision,intent(in)                       :: dt              !Time step
   double precision,intent(in)                       :: dx              !Space interval
   double precision,dimension(0:NX+1),intent(in)     :: u               !Variable including BC
   double precision,dimension(0:NX+1),intent(in)     :: D               !Coefficient for diffusion term including BC
   double precision,dimension(1:NX),intent(in)       :: S               !Source term coefficient
   double precision,dimension(1:NX),intent(out)      :: delta           !Delta of variable in implicit form
   double precision,dimension(1:NX)                  :: rhs             !Right-hand side in explicit form
   double precision,dimension(1:NX,3)                :: dcoef           !Diffusion coefficients for LHS (1:i-1,2:i,3:i+1)
   double precision,dimension(1:NX,3)                :: lcoef           !LHS coefficients (1:i-1,2:i,3:i+1)


   !Diffusion term (dt*dE^{dif}/dx) calculation for RHS
   call DiffusionCoefficient1D(NX,dt,dx,D,dcoef)

   !LHS coefficient calculation
   do i = 1,NX
      lcoef(i,1) = dcoef(i,1)
      lcoef(i,2) = 1.0d0+dcoef(i,2)
      lcoef(i,3) = dcoef(i,3)
   enddo

   !RHS calculation
   do i = 1,NX
      rhs(i) =-(dcoef(i,1)*u(i-1)+dcoef(i,2)*u(i  )+dcoef(i,3)*u(i+1))+dt*S(i)
   enddo

   !Matrix inversion calculation using Thomas' algorithm
   call Thomas(NX,lcoef,rhs,delta)

   return
endsubroutine


!---------- 1D Diffusion Equation Solver with Internal Loop ------------
!          Input : NX,dt,dx,u(0:NX+1),D(0:NX+1),S(1,NX),LPMAX
!          Output: delta(1:NX)
!-----------------------------------------------------------------------

subroutine DiffusionLpBD1D(NX,dt,dx,u,D,S,delta,LPMAX)

   implicit none
   integer :: i,k
   integer,intent(in)                                :: NX              !Number of cells
   integer,intent(in)                                :: LPMAX           !Maximum iteration number for internal loop
   double precision                                  :: dt              !Time step
   double precision,intent(in)                       :: dx              !Space interval
   double precision,dimension(0:NX+1),intent(in)     :: u               !Variable including BC
   double precision,dimension(0:NX+1),intent(in)     :: D               !Coefficient for diffusion term including BC
   double precision,dimension(1:NX),intent(in)       :: S               !Source term coefficient
   double precision,dimension(0:NX+1)                :: uk              !Variable during internal loop
   double precision,dimension(1:NX)                  :: rhs             !Right-hand side in explicit form
   double precision,dimension(1:NX),intent(inout)    :: delta           !Delta of variable in implicit form
   double precision,dimension(1:NX,3)                :: dcoef           !Diffusion coefficients for LHS (1:i-1,2:i,3:i+1)
   double precision,dimension(1:NX,3)                :: lcoef           !LHS coefficients (1:i-1,2:i,3:i+1)


   dt = 2.0d0/3.0d0*dt

   !Diffusion term (dt*dE^{dif}/dx) calculation for RHS
   call DiffusionCoefficient1D(NX,dt,dx,D,dcoef)

   !LHS coefficient calculation
   do i = 1,NX
      lcoef(i,1) = dcoef(i,1)
      lcoef(i,2) = 1.0d0+dcoef(i,2)
      lcoef(i,3) = dcoef(i,3)
   enddo

   !Variable during internal loop
   do i = 0,NX+1
      uk(i) = u(i)
   enddo
   do k = 1,LPMAX
      !RHS calculation
      do i = 1,NX
         rhs(i) =-(lcoef(i,1)*uk(i-1)+lcoef(i,2)*uk(i  )+lcoef(i,3)*uk(i+1))&
                 +dt*S(i)+u(i)+1.0d0/3.0d0*delta(i)
      enddo

      !Matrix inversion calculation using Thomas' algorithm
      call Thomas(NX,lcoef,rhs,delta)

      !Update
      do i = 1,NX
         uk(i) = uk(i)+delta(i)
      enddo
   enddo
   do i = 1,NX
      delta(i) = uk(i)-u(i)
   enddo

   return
endsubroutine


!------------ 1-D Advection-Diffusion Equation Solver ------------------
!      Input : NX,dt,dx,u(-1:NX+2),C(-1:NX+2),D(0:NX+1),S(1,NX)
!      Output: delta(1:NX)
!-----------------------------------------------------------------------

subroutine AdvectionDiffusion1D(NX,dt,dx,u,C,D,S,delta)

   implicit none
   integer :: i
   integer                            ,intent(in)    :: NX              !Number of cells
   double precision                   ,intent(in)    :: dt              !Time interval
   double precision                   ,intent(in)    :: dx              !Space interval
   double precision,dimension(-1:NX+2),intent(in)    :: u               !Variable including BC
   double precision,dimension(-1:NX+2),intent(in)    :: C               !Advection coefficient including BC
   double precision,dimension(0:NX+1) ,intent(in)    :: D               !Diffusion coefficient including BC
   double precision,dimension(1:NX)   ,intent(in)    :: S               !Source term
   double precision,dimension(1:NX)   ,intent(out)   :: delta           !Delta of variable (n+1 - n)
   double precision,dimension(1:NX)                  :: rhs             !Right-hand side in explicit form
   double precision,dimension(1:NX,5)                :: ccoef
   double precision,dimension(1:NX,3)                :: dcoef
   double precision,dimension(1:NX,5)                :: lcoef,LU

   !Coefficient calculation of advection and diffusion terms
   call AdvectionCoefficient1D(NX,dt,dx,C,u,ccoef)
   call DiffusionCoefficient1D(NX,dt,dx,D,dcoef)

   !Element calculation in 5-band matrix
   do i = 1,NX
      lcoef(i,1) = ccoef(i,1)
      lcoef(i,2) = ccoef(i,2)+dcoef(i,1)
      lcoef(i,3) = ccoef(i,3)+dcoef(i,2)+1.0d0
      lcoef(i,4) = ccoef(i,4)+dcoef(i,3)
      lcoef(i,5) = ccoef(i,5)
   enddo

   !Calculate right-hand side
   do i = 1,NX
      rhs(i)    =-(lcoef(i,1)*u(i-2)+lcoef(i,2)*u(i-1)+(lcoef(i,3)-1.0d0)*u(i)&
                  +lcoef(i,4)*u(i+1)+lcoef(i,5)*u(i+2))+dt*S(i)
   enddo

   !Matrix inversion calculation using LU factorization for 5-band matrix
   ! 2*W+1 = 5  =>  W = 2
   call LUfactorization(NX,2,lcoef,LU)
   call ForBackSubstitution(NX,2,LU,rhs,delta)

   return
endsubroutine




!------------ 1-D Advection-Diffusion Equation Solver ------------------
!-------------------- Last Update 2015/10/27 ---------------------------
!      Input : NX,dt,dx,u(-1:NX+2),C(-1:NX+2),D(0:NX+1),S(1,NX)
!      Output: delta(1:NX)
!-----------------------------------------------------------------------

subroutine AdvectionDiffusionLpBD1D(NX,dt,dx,u,C,D,S,delta,kMAX)

   implicit none
   integer :: i,k
   integer,intent(in)                                :: NX              !Number of cells
   integer,intent(in)                                :: kMAX            !Maximum iteration number for internal loop
   double precision                                  :: dt              !Time interval
   double precision,intent(in)                       :: dx              !Space interval
   double precision,dimension(-1:NX+2),intent(in)    :: u               !Variable including BC
   double precision,dimension(-1:NX+2),intent(in)    :: C               !Advection coefficient including BC
   double precision,dimension(0:NX+1),intent(in)     :: D               !Diffusion coefficient including BC
   double precision,dimension(1:NX),intent(in)       :: S               !Source term
   double precision,dimension(1:NX),intent(inout)    :: delta           !Delta of variable in implicit form
   double precision,dimension(1:NX)                  :: rhs             !Right-hand side in explicit form
   double precision,dimension(-1:NX+2)               :: uk              !Variable during internal loop
   double precision,dimension(1:NX)                  :: deltak          !Delta during internal loop
   double precision,dimension(1:NX,5)                :: ccoef
   double precision,dimension(1:NX,3)                :: dcoef
   double precision,dimension(1:NX,5)                :: coef,LU

   dt = 2.0d0/3.0d0*dt

   !Coefficient calculation of advection and diffusion terms
   call AdvectionCoefficient1D(NX,dt,dx,C,u,ccoef)
   call DiffusionCoefficient1D(NX,dt,dx,D,dcoef)

   !Element calculation in 5-band matrix
   do i = 1,NX
      coef(i,1) = ccoef(i,1)
      coef(i,2) = ccoef(i,2)+dcoef(i,1)
      coef(i,3) = ccoef(i,3)+dcoef(i,2)+1.0d0
      coef(i,4) = ccoef(i,4)+dcoef(i,3)
      coef(i,5) = ccoef(i,5)
   enddo
   !LU factorization of 5-band matrix
   call LUfactorization(NX,2,coef,LU)

   !Variable during internal loop
   do i =-1,NX+2
      uk(i) = u(i)
   enddo

   !Internal loop
   do k = 1,kMAX
      do i = 1,NX
         rhs(i) =-(coef(i,1)*uk(i-2)+coef(i,2)*uk(i-1)+coef(i,3)*uk(i)&
                  +coef(i,4)*uk(i+1)+coef(i,5)*uk(i+2))+dt*S(i)+u(i)+1.0d0/3.0d0*delta(i)
      enddo
      !Forward-Backward substitution using LU
      call ForBackSubstitution(NX,2,LU,rhs,deltak)
      !Update
      do i = 1,NX
         uk(i) = uk(i)+deltak(i)
      enddo
   enddo

   !Finally calculate delta
   do i = 1,NX
      delta(i) = uk(i)-u(i)
   enddo

   return
endsubroutine


!--------------- Coefficient for 1D advection term ---------------------
!      Input:  dt,dx,c(-1:NX+2), u(-1:NX+2)
!      Output: a(1:NX,5)
!      dt(d(Cu)/dx) = a(i,1)u(i-2)+a(i,2)u(i-1)+a(i,3)u(i)
!                    +a(i,4)u(i+1)+a(i,5)u(i+2)
!-----------------------------------------------------------------------

subroutine AdvectionCoefficient1D(NX,dt,dx,c,u,a)

   implicit none
   integer :: i
   integer,intent(in)                               :: NX
   double precision,intent(in)                      :: dt,dx
   double precision,dimension(-1:NX+2),intent(in)   :: c
   double precision,dimension(-1:NX+2),intent(in)   :: u
   double precision,dimension(1:NX,5),intent(out)   :: a
   double precision,dimension(0:NX+1)               :: sc1,sc2,su1,su2
   double precision,dimension(1:NX+1)               :: cL,cR,machL,machR
   double precision,dimension(1:NX+1)               :: cLp,cRm
   double precision :: db,df,CNORM

   !Determine switch for c
   do i = 0,NX+1
      db = c(i)-c(i-1); df = c(i+1)-c(i)
      if(db*df.lt.0.0d0) then
         sc1(i) = 0.0d0; sc2(i) = 0.0d0
      else if(dabs(db).le.dabs(df)) then
         sc1(i) = 1.0d0; sc2(i) = 0.0d0
      else if(dabs(db).gt.dabs(df)) then
         sc1(i) = 0.0d0; sc2(i) = 1.0d0
      endif
   enddo
   !Determine switch for u
   do i = 0,NX+1
      db = u(i)-u(i-1); df = u(i+1)-u(i)
      if(db*df.lt.0.0d0) then
         su1(i) = 0.0d0; su2(i) = 0.0d0
      else if(dabs(db).le.dabs(df)) then
         su1(i) = 1.0d0; su2(i) = 0.0d0
      else if(dabs(db).gt.dabs(df)) then
         su1(i) = 0.0d0; su2(i) = 1.0d0
      endif
   enddo
   !Calculate interpolated c
   do i = 1,NX+1
      cL(i) = c(i-1)+0.5d0*(sc1(i-1)*(c(i-1)-c(i-2))+sc2(i-1)*(c(i  )-c(i-1)))
      cR(i) = c(i  )-0.5d0*(sc1(i  )*(c(i  )-c(i-1))+sc2(i  )*(c(i+1)-c(i  )))
   enddo
   CNORM = 1.0d0
   !Calculated upwinded c
   do i = 1,NX+1
      machL(i) = cL(i)/CNORM
      machR(i) = cR(i)/CNORM
      !if(dabs(machL(i)).ge.1.0d0) then
      !   cLp(i) = CNORM*0.5d0*(machL(i)+dabs(machL(i)))
      !else
      !   cLp(i) = CNORM*0.25d0*(machL(i)+1.0d0)**2.0d0
      !endif
      !if(dabs(machR(i)).ge.1.0d0) then
      !   cRm(i) = CNORM*0.5d0*(machR(i)-dabs(machR(i)))
      !else
      !   cRm(i) =-CNORM*0.25d0*(machR(i)-1.0d0)**2.0d0
      !endif
      cLp(i) = 0.5d0*(cL(i)+dabs(cL(i)))
      cRm(i) = 0.5d0*(cR(i)-dabs(cR(i)))
   enddo
   !Calculate coefficients a
   do i = 1,NX
      a(i,1) = dt/2.0d0/dx*( cLp(i  )*su1(i-1))
      a(i,2) = dt/2.0d0/dx*(-cLp(i+1)*su1(i  )&
                            -cRm(i  )*su1(i  )&
                            -cLp(i  )*(2.0d0+su1(i-1)-su2(i-1)))
      a(i,3) = dt/2.0d0/dx*( cLp(i+1)*(2.0d0+su1(i  )-su2(i  ))&
                            +cRm(i+1)*su1(i+1)&
                            -cLp(i  )*su2(i-1)&
                            -cRm(i  )*(2.0d0-su1(i  )+su2(i  )))
      a(i,4) = dt/2.0d0/dx*( cLp(i+1)*su2(i  )&
                            +cRm(i+1)*(2.0d0-su1(i+1)+su2(i+1))&
                            +cRm(i  )*su2(i  ))
      a(i,5) = dt/2.0d0/dx*(-cRm(i+1)*su2(i+1))
   enddo

   return
endsubroutine



!--------------- Coefficient for 1D diffusion term ---------------------
!      Input : NX,dt,dx,D(0:NX+1),u(0:NX+1)
!      Output: a(1:NX,3)
!      dt(d/dx(-Dd/dx(u))) = a(i,1)u(i-1)+a(i,2)u(i)+a(i,3)u(i+1)
!-----------------------------------------------------------------------

subroutine DiffusionCoefficient1D(NX,dt,dx,D,a)

   implicit none
   integer :: i
   integer,intent(in)                               :: NX
   double precision,intent(in)                      :: dt,dx
   double precision,dimension(0:NX+1),intent(in)    :: D
   double precision,dimension(1:NX,3),intent(out)   :: a

   !Calculate coefficients a
   do i = 1,NX
      a(i,1) =-dt/2.0d0/dx/dx*(D(i  )+D(i-1))
      a(i,2) = dt/2.0d0/dx/dx*(D(i+1)+2.0d0*D(i  )+D(i-1))
      a(i,3) =-dt/2.0d0/dx/dx*(D(i+1)+D(i  ))
   enddo

   return
endsubroutine




!------------- Thomas algorithm for tridiagonal matrix -----------------
! Input  NX: Size of matrix, coef: tri-diagonal elements in matrix
!        dd: Right-hand side
! Output xx: Solution
!-----------------------------------------------------------------------

subroutine Thomas(NX,coef,dd,xx)

   implicit none
   integer,intent(in)                            :: NX
   double precision,dimension(1:NX,3),intent(in) :: coef
   double precision,dimension(1:NX),intent(in)   :: dd
   double precision,dimension(1:NX),intent(out)  :: xx
   double precision,dimension(1:NX) :: aa,bb,cc
   double precision,dimension(0:NX) :: ee,ff
   integer :: i,k,l


   do i = 1,NX
      aa(i) = coef(i,1)
      bb(i) = coef(i,2)
      cc(i) = coef(i,3)
   enddo
   ee(0) = 0.0d0
   ff(0) = 0.0d0
   do k = 1,NX
      ee(k) =-cc(k)/(bb(k)+aa(k)*ee(k-1))
      ff(k) = (dd(k)-aa(k)*ff(k-1))/(bb(k)+aa(k)*ee(k-1))
   enddo
   xx(NX) = ff(NX)
   do k = 1,NX-1
      l = NX-k
      xx(l) = ee(l)*xx(l+1)+ff(l)
   enddo


   return
endsubroutine


!---------- LU factorization in Matrix Inversion -----------------------
! Input  N: Size of matrix, 2*W+1: Width of band matrix, A: Band matrix
! Output LU: LU-factorized matrix
!-----------------------------------------------------------------------

subroutine LUfactorization(N,W,A,LU)

   implicit none
   integer,intent(in) :: N,W
   double precision,dimension(N,2*W+1),intent(in)  :: A
   double precision,dimension(N,2*W+1),intent(out) :: LU
   integer :: i,j,k
   double precision  :: dtmp,dakj

   do j = 1,2*W+1
      do i = 1,N
         LU(i,j) = A(i,j)
      enddo
   enddo
   !LU factorization based on Gauss method
   do k = 1,N
      if(LU(k,W+1).eq.0.0d0) write(*,*) 'Error in LU factorization...'
      dtmp = 1.0d0/LU(k,W+1)
      if(k+W.le.N) then
         do i = k+1,k+W
            LU(i,W+1+k-i) = LU(i,W+1+k-i)*dtmp
         enddo
      else
         do i = k+1,N
            LU(i,W+1+k-i) = LU(i,W+1+k-i)*dtmp
         enddo
      endif
      if(k+W.le.N) then
         do j = k+1,k+W
            dakj = LU(k,W+1+j-k)
            do i = k+1,k+W
               LU(i,W+1+j-i) = LU(i,W+1+j-i)-LU(i,W+1+k-i)*dakj
            enddo
         enddo
      else
         do j = k+1,N
            dakj = LU(k,W+1+j-k)
            do i = k+1,N
               LU(i,W+1+j-i) = LU(i,W+1+j-i)-LU(i,W+1+k-i)*dakj
            enddo
         enddo
      endif
   enddo

   return
endsubroutine


!------- Forward-Backward Substitution in Matrix Inversion -------------
! Input  N: Size of matrix, 2*W+1: Width of band matrix,
!        LU: LU-factorized matrix, b: Right-hand side
! Output x: Solution of linear system
!-----------------------------------------------------------------------

subroutine ForBackSubstitution(N,W,LU,b,x)

   implicit none
   integer,intent(in) :: N,W
   double precision,dimension(N,2*W+1),intent(in)  :: LU
   double precision,dimension(N),intent(in)    :: b
   double precision,dimension(N),intent(out)   :: x
   double precision,dimension(N)               :: y
   integer :: i,j,k

   ! Forward substitution
   y = 0.0d0
   do i = 1,N
      if(i.eq.1) then
         y(i) = b(i)
      else if(i.le.W) then
         do j = 1,i-1
            y(i) = y(i)+LU(i,W+1+j-i)*y(j)
         enddo
         y(i) = b(i)-y(i)
      else
         do j = i-W,i-1
            y(i) = y(i)+LU(i,W+1+j-i)*y(j)
         enddo
         y(i) = b(i)-y(i)
      endif
   enddo
   ! Backward substitution
   x = 0.0d0
   do i = 1,N
      k = N-i+1
      if(LU(k,W+1).eq.0.0d0) then
         write(*,*) 'Error in ForBackSubs...'
         stop
      endif
      if(i.eq.1) then
         x(k) = y(k)/LU(k,W+1)
      else if(i .le. W) then
         do j = k+1,N
            x(k) = x(k)+LU(k,W+1+j-k)*x(j)
         enddo
         x(k) = (y(k)-x(k))/LU(k,W+1)
      else
         do j = k+1,k+W
            x(k) = x(k)+LU(k,W+1+j-k)*x(j)
         enddo
         x(k) = (y(k)-x(k))/LU(k,W+1)
      endif
   enddo

   return
endsubroutine



