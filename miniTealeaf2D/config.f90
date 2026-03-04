module config
   implicit none

   !--------------------------------
   ! Domain size (2D)
   !--------------------------------
   integer :: nx = 200
   integer :: ny = 200
   real(8) :: dx = 1.0d0
   real(8) :: dy = 1.0d0

   !--------------------------------
   ! Hydrodynamics parameters
   !--------------------------------
   real(8) :: gamma = 1.4d0
   real(8) :: CFL   = 0.4d0
   real(8) :: cv    = 1.0d0
   real(8) :: V0    = 1.5d0     ! Initial advection speed

   !--------------------------------
   ! Radiation parameters
   !--------------------------------
   real(8) :: c   = 1.0d0       ! Radiation speed
   real(8) :: Ag  = 1.0d0
   real(8) :: Kop = 1.0d0       ! Moderate opacity

   !--------------------------------
   ! Time control
   !--------------------------------
   real(8) :: tf    = 1000.0d0   ! Final time
   integer :: ITmax = 10000    ! Max iterations

   real(8) :: t  = 0.0d0         ! Time
   real(8) :: dt = 0.0d0         ! Time step
   integer :: it = 0             ! Iteration


   !--------------------------------
   ! Conjugate Gradient parameters
   !--------------------------------
   real(8) :: CGtol   = 1.0d-10
   integer :: ITmaxCG = 1000

   !--------------------------------
   ! Allocated fields
   !--------------------------------
   real(8), allocatable :: rho(:,:)
   real(8), allocatable :: momentumX(:,:)
   real(8), allocatable :: momentumY(:,:)
   real(8), allocatable :: Etot(:,:)
   real(8), allocatable :: Er(:,:)


   !--------------------------------
   ! Output configuration
   !--------------------------------
   character(len=64) :: filename = "out"
   integer :: dumpEvery = 1

end module config
