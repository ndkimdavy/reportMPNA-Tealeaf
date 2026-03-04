program miniHades2D
   use config
   use hydro,     only: hydro_init => init, hydro_step => kernel
   use radiation, only: rad_init   => init, rad_step   => kernel
   use io,        only: dump
   implicit none

   integer :: i, j
   real(8) :: vx, vy, p, cs, r, mx, my, E, smax, kin

   allocate(rho(nx,ny), momentumX(nx,ny), momentumY(nx,ny), Etot(nx,ny), Er(nx,ny))

   call hydro_init()
   call rad_init()

   write(*,*)
   write(*,*) "Simulation running..."
   write(*,*) "Run ./plot (same directory) for live visualization."
   write(*,*)

   t  = 0.0d0
   it = 0

   call dump(Er)

   do while (t < tf .and. it < ITmax)

      smax = 0.0d0

      do j = 2, ny-1
         do i = 2, nx-1

            r = rho(i,j)
            if (r <= 0.d0) cycle

            mx = momentumX(i,j)
            my = momentumY(i,j)
            E  = Etot(i,j)

            vx  = mx / r
            vy  = my / r
            kin = 0.5d0 * (mx*mx + my*my) / r

            p = (gamma - 1.d0) * (E - kin)
            if (p < 1.d-14) p = 1.d-14

            cs = sqrt(gamma * p / r)

            smax = max(smax, abs(vx) + cs, abs(vy) + cs)

         end do
      end do

      if (smax > 0.d0) then
         dt = CFL * min(dx,dy) / smax
      else
         dt = 0.d0
      end if

      if (t + dt > tf) dt = tf - t
      if (dt <= 0.d0) exit

      call hydro_step(dt)
      call rad_step(dt)

      t  = t + dt
      it = it + 1

      if (mod(it, dumpEvery) == 0) call dump(Er)

   end do

   call dump(Er)

   write(*,*) "Simulation finished."

   deallocate(rho, momentumX, momentumY, Etot, Er)

end program miniHades2D
