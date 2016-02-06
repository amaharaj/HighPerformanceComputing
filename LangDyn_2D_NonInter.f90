module globals
use omp_lib
! Global variables
implicit none
integer :: n=500                                ! number of particles
double precision :: L=1                         ! domain size with periodic boundary conditions
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
double precision :: dt,kT,g,m                   ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy
!integer, allocatable :: seed(:)
!allocate(seed(1))
contains
subroutine set_parameters

! Set time step and physical parameters
dt=0.01d0
kT=1.0d0
g=1.0d0
m=1.0d0

! Set auxiliary parameters
pref1=g
pref2=sqrt(24*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles
integer :: i
double precision :: ran1,ran2,gr1,gr2

! Give particles initial position and velocity
do i=1,n
   x(i)=0d0
   y(i)=0d0
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2)
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx(i)=gr1
   vy(i)=gr2
end do

end subroutine initialize_particles
end module Langevin

program main
use globals
use Langevin
implicit none
integer :: i,nt,myt,threadnum
double precision :: t,t_max,ran1,ran2,m1,m2,begin,end,wtime

! Open files
open(11,file='trajectories')
open(12,file='means')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n))

t=0d0
t_max=1d+00

call set_parameters
call initialize_particles


! OMP implementation of CPU clock
begin=omp_get_wtime()


do while(t.lt.t_max)

!$omp parallel private(i,ran1,ran2)

!$omp do
   do i=1,n
      vhx(i)=vx(i)+0.5*ax(i)*dt
      vhy(i)=vy(i)+0.5*ay(i)*dt
      x(i)=x(i)+vhx(i)*dt
      y(i)=y(i)+vhy(i)*dt
      ax(i)=0d0                   ! Add forces here if any
      ay(i)=0d0                   ! Add forces here if any
      call random_number(ran1)
      ran1=ran1-0.5d0
      call random_number(ran2)
      ran2=ran2-0.5d0
      ax(i)=ax(i)-pref1*vhx(i)+pref2*ran1
      ay(i)=ay(i)-pref1*vhy(i)+pref2*ran2
         
      vx(i)=vhx(i)+0.5*ax(i)*dt
      vy(i)=vhy(i)+0.5*ay(i)*dt



      if(x(i).lt.-L/2d0) x(i)=x(i)+L
      if(x(i).gt.L/2d0) x(i)=x(i)-L
      if(y(i).lt.-L/2d0) y(i)=y(i)+L
      if(y(i).gt.L/2d0) y(i)=y(i)-L
 

   end do
!$omp end do
!$omp end parallel

   ! Move mean calculations outside of the loop over particles
   m1=sum(x)/dfloat(n)
   m2=sum(x**2)/dfloat(n) 
   ! Need to write in serial, otherwise can not follow path of particle
   do i=1,n
      write(11,'(2(1x,e10.3))') x(i),y(i)
   end do
   ! Only need to calculate mean once per time step
   write(12,'(2(1x,e10.3))') m1,m2

   write(11,*) ' '
   t=t+dt
end do

end=omp_get_wtime()
wtime=end-begin

print *,n,wtime

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay)
! Close files
close(11)
close(12)

end program main
