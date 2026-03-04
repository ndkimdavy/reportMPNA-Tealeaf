module radiation
   use config
   implicit none
   private
   public :: init, kernel

contains

   subroutine init()
      integer :: i, j
      real(8) :: x, y, x0, y0, sig

      x0  = 0.5d0*nx*dx
      y0  = 0.5d0*ny*dy
      sig = 0.10d0*nx*dx

      do j = 1, ny
         y = (j-0.5d0)*dy
         do i = 1, nx
            x = (i-0.5d0)*dx
            Er(i,j) = exp(-((x-x0)**2 + (y-y0)**2)/(sig*sig))
         end do
      end do
   end subroutine init


   subroutine kernel(dt)
      real(8), intent(in) :: dt
      integer :: i, j
      real(8) :: dErx, dEry
      real(8) :: r, mx, my, E, kin, u, T
      real(8) :: sigma, Er_new, dEr

      real(8), allocatable :: b(:,:), x(:,:)

      !-------------------------------------------------
      ! (1) Radiation pressure force: div P = (1/3) grad(Er)
      ! momentum -= dt * (1/3) * grad(Er)
      !-------------------------------------------------
      do j = 2, ny-1
         do i = 2, nx-1
            dErx = (Er(i+1,j) - Er(i-1,j)) / (2.0d0*dx)
            dEry = (Er(i,j+1) - Er(i,j-1)) / (2.0d0*dy)
            momentumX(i,j) = momentumX(i,j) - dt*(1.0d0/3.0d0)*dErx
            momentumY(i,j) = momentumY(i,j) - dt*(1.0d0/3.0d0)*dEry
         end do
      end do

      !-------------------------------------------------
      ! (2) Local source: S = c*rho*Kop*(Ag*T^4 - Er)
      ! Semi-implicit in Er (T frozen from current Etot):
      ! Er^{n+1} = (Er + dt*sigma*Ag*T^4) / (1 + dt*sigma)
      ! with sigma = c*rho*Kop
      ! Energy conservation locally: Etot -= (Er^{n+1}-Er)
      !-------------------------------------------------
      do j = 2, ny-1
         do i = 2, nx-1
            r  = rho(i,j)
            mx = momentumX(i,j)
            my = momentumY(i,j)
            E  = Etot(i,j)

            kin = 0.5d0*(mx*mx + my*my)/r
            u   = (E - kin)/r
            if (u < 0.0d0) u = 0.0d0
            T   = u / cv

            sigma  = c * r * Kop
            Er_new = (Er(i,j) + dt*sigma*Ag*(T**4)) / (1.0d0 + dt*sigma)

            dEr       = Er_new - Er(i,j)
            Er(i,j)   = Er_new
            Etot(i,j) = Etot(i,j) - dEr
         end do
      end do

      !-------------------------------------------------
      ! (3) Implicit diffusion: (I - dt L) Er^{n+1} = Er*
      ! L(x) = div( D grad x ), D = c/(3 rho Kop)
      ! Solve by CG on interior cells.
      !-------------------------------------------------
      allocate(b(nx,ny), x(nx,ny))
      b(:,:) = Er(:,:)
      x(:,:) = Er(:,:)

      call solverCG(dt, b, x)

      Er(:,:) = x(:,:)

      deallocate(b, x)
   end subroutine kernel


   !==========================================================
   ! CG solver for: A(x) = b, with A(x)=x - dt*L(x)
   ! L(x)=div(D grad x), D=c/(3 rho Kop)
   !==========================================================
   subroutine solverCG(dt, b, x)
      real(8), intent(in)    :: dt
      real(8), intent(in)    :: b(nx,ny)
      real(8), intent(inout) :: x(nx,ny)

      integer :: k, i, j
      real(8) :: alpha, beta
      real(8) :: r2, r2new, pAp
      real(8) :: bn, res

      real(8), allocatable :: r(:,:), p(:,:), q(:,:)

      allocate(r(nx,ny), p(nx,ny), q(nx,ny))

      call A(dt, x, q)
      r(:,:) = b(:,:) - q(:,:)
      p(:,:) = r(:,:)

      r2 = dot(r, r)
      bn = sqrt(max(dot(b,b), 1.0d-300))

      do k = 1, ITmaxCG

         call A(dt, p, q)
         pAp = dot(p, q)
         if (abs(pAp) < 1.0d-300) exit

         alpha = r2 / pAp

         do j = 2, ny-1
            do i = 2, nx-1
               x(i,j) = x(i,j) + alpha*p(i,j)
               r(i,j) = r(i,j) - alpha*q(i,j)
            end do
         end do

         r2new = dot(r, r)
         res   = sqrt(r2new) / bn
         if (res < CGtol) exit

         beta = r2new / r2
         r2   = r2new

         do j = 2, ny-1
            do i = 2, nx-1
               p(i,j) = r(i,j) + beta*p(i,j)
            end do
         end do

      end do

      deallocate(r, p, q)

   contains

      subroutine A(dt, x, Ax)
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: x(nx,ny)
         real(8), intent(out) :: Ax(nx,ny)
         real(8) :: Lx(nx,ny)

         call L(x, Lx)
         Ax(:,:) = x(:,:) - dt*Lx(:,:)
      end subroutine A


      subroutine L(x, Lx)
         real(8), intent(in)  :: x(nx,ny)
         real(8), intent(out) :: Lx(nx,ny)

         integer :: i, j
         real(8) :: D, De, Dw, Dn, Ds

         Lx(:,:) = 0.0d0

         do j = 2, ny-1
            do i = 2, nx-1

               D  = c / (3.0d0 * rho(i,j) * Kop)

               De = 0.5d0*(D + c/(3.0d0*rho(i+1,j)*Kop))
               Dw = 0.5d0*(D + c/(3.0d0*rho(i-1,j)*Kop))
               Dn = 0.5d0*(D + c/(3.0d0*rho(i,j+1)*Kop))
               Ds = 0.5d0*(D + c/(3.0d0*rho(i,j-1)*Kop))

               Lx(i,j) = ( De*(x(i+1,j) - x(i,j)) - Dw*(x(i,j) - x(i-1,j)) )/(dx*dx) &
                  + ( Dn*(x(i,j+1) - x(i,j)) - Ds*(x(i,j) - x(i,j-1)) )/(dy*dy)
            end do
         end do
      end subroutine L


      pure function dot(a, b) result(s)
         real(8), intent(in) :: a(nx,ny), b(nx,ny)
         real(8) :: s
         integer :: i, j
         s = 0.0d0
         do j = 2, ny-1
            do i = 2, nx-1
               s = s + a(i,j)*b(i,j)
            end do
         end do
      end function dot

   end subroutine solverCG

end module radiation
