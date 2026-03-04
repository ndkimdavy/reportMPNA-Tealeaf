module io
   use config
   implicit none
   private
   public :: dump

contains

   subroutine dump(field)
      real(8), intent(in) :: field(nx,ny)
      integer :: i, j, u
      character(len=128) :: fname, tmp
      real(8) :: x, y

      write(fname,'(A,".dat")') trim(filename)
      write(tmp ,'(A,".tmp")') trim(filename)

      open(newunit=u, file=tmp, status="replace", action="write")
      do j = 1, ny
         y = (j-1)*dy
         do i = 1, nx
            x = (i-1)*dx
            write(u,'(3(ES16.8,1X))') x, y, field(i,j)
         end do
         write(u,*)
      end do
      flush(u)
      close(u)

      call execute_command_line("mv -f " // trim(tmp) // " " // trim(fname), wait=.true.)
   end subroutine dump

end module io
