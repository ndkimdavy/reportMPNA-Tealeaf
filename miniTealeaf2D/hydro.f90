module hydro
   use config
   implicit none
   private
   public :: init, kernel

contains

   subroutine init()
      integer :: i, j
      real(8) :: v0

      rho(:,:)  = 1.0d0
      Etot(:,:) = 1.0d0

      v0 = 0.5d0
      do j = 1, ny
         do i = 1, nx
            momentumX(i,j) = rho(i,j) * v0
            momentumY(i,j) = 0.0d0
         end do
      end do
   end subroutine init


   subroutine kernel(dt)
      real(8), intent(in) :: dt
      integer :: i, j
      real(8) :: invdx, invdy
      real(8) :: U(4), Un(4), Fm(4), Fp(4), Gm(4), Gp(4)
      real(8), allocatable :: rhoN(:,:), mxN(:,:), myN(:,:), EN(:,:)

      invdx = 1.0d0/dx
      invdy = 1.0d0/dy

      allocate(rhoN(nx,ny), mxN(nx,ny), myN(nx,ny), EN(nx,ny))

      rhoN(:,:) = rho(:,:)
      mxN(:,:)  = momentumX(:,:)
      myN(:,:)  = momentumY(:,:)
      EN(:,:)   = Etot(:,:)

      do j = 2, ny-1
         do i = 2, nx-1

            U  = packU(i,j)

            Fm = rusanov_flux(1, packU(i-1,j), packU(i,  j))
            Fp = rusanov_flux(1, packU(i,  j), packU(i+1,j))
            Gm = rusanov_flux(2, packU(i,j-1), packU(i,  j))
            Gp = rusanov_flux(2, packU(i,  j), packU(i,j+1))

            Un = U - dt*invdx*(Fp - Fm) - dt*invdy*(Gp - Gm)

            rhoN(i,j) = Un(1)
            mxN(i,j)  = Un(2)
            myN(i,j)  = Un(3)
            EN(i,j)   = Un(4)
         end do
      end do

      rho(:,:)       = rhoN(:,:)
      momentumX(:,:) = mxN(:,:)
      momentumY(:,:) = myN(:,:)
      Etot(:,:)      = EN(:,:)

      deallocate(rhoN, mxN, myN, EN)

   contains

      function packU(i,j) result(U)
         integer, intent(in) :: i, j
         real(8) :: U(4)
         U(1) = rho(i,j)
         U(2) = momentumX(i,j)
         U(3) = momentumY(i,j)
         U(4) = Etot(i,j)
      end function packU

      pure subroutine prim(U, vx, vy, p, cs)
         real(8), intent(in)  :: U(4)
         real(8), intent(out) :: vx, vy, p, cs
         real(8) :: r, mx, my, E
         r = U(1); mx = U(2); my = U(3); E = U(4)
         vx = mx/r
         vy = my/r
         p  = (gamma-1.0d0)*(E - 0.5d0*(mx*mx + my*my)/r)
         cs = sqrt(gamma*p/r)
      end subroutine prim

      pure function flux(dir, U) result(F)
         integer, intent(in) :: dir   ! 1:x, 2:y
         real(8), intent(in) :: U(4)
         real(8) :: F(4), vx, vy, p, cs
         call prim(U, vx, vy, p, cs)
         if (dir == 1) then
            F(1)=U(2);        F(2)=U(2)*vx + p;  F(3)=U(3)*vx;       F(4)=(U(4)+p)*vx
         else
            F(1)=U(3);        F(2)=U(2)*vy;      F(3)=U(3)*vy + p;   F(4)=(U(4)+p)*vy
         end if
      end function flux

      pure function rusanov_flux(dir, UL, UR) result(Fh)
         integer, intent(in) :: dir
         real(8), intent(in) :: UL(4), UR(4)
         real(8) :: Fh(4)
         real(8) :: vxL, vyL, pL, csL, vxR, vyR, pR, csR, smax
         call prim(UL, vxL, vyL, pL, csL)
         call prim(UR, vxR, vyR, pR, csR)
         if (dir == 1) then
            smax = max(abs(vxL)+csL, abs(vxR)+csR)
         else
            smax = max(abs(vyL)+csL, abs(vyR)+csR)
         end if
         Fh = 0.5d0*(flux(dir, UL) + flux(dir, UR)) - 0.5d0*smax*(UR - UL)
      end function rusanov_flux

   end subroutine kernel

end module hydro
