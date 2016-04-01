module globals
use omp_lib
!use mt95
! Global 1(I4),2(1x,e10.31(I4),2(1x,e10.3)bles
implicit none
integer :: n,t_max=100                                ! number of particles
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
double precision, allocatable, dimension(:) :: rnds
contains
subroutine set_parameters

! Set time step and physical parameters
dt=1.0d0
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
double precision :: t,ran1,ran2,m1,m2,begin,end,wtime
integer, allocatable :: seed(:)
integer size
character(len=32) :: arg

call getarg(1,arg)
read (arg(1:20),'(i20)') n

! Open files
open(11,file='trajectories')
open(12,file='means')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n))
allocate(rnds(2*n))

t=0d0
! t_max=1d+00

call set_parameters
call initialize_particles

! OMP implementation of CPU clock
begin=omp_get_wtime()

do while(t.lt.t_max)
    
  do i=1,n
     call random_number(rnds(i))
     call random_number(rnds(2*i+1))
  end do
  print *,rnds(1),t

!$omp parallel private(i,threadnum,ran1,ran2)

! Check if there are multiple threads being used
! !$ threadnum=omp_get_num_threads()
! !$ print *,'Thread num:',threadnum


!$omp do
   do i=1,n
      vhx(i)=vx(i)+0.5*ax(i)*dt
      vhy(i)=vy(i)+0.5*ay(i)*dt
      x(i)=x(i)+vhx(i)*dt
      y(i)=y(i)+vhy(i)*dt
      ax(i)=0d0                   ! Add forces here if any
      ay(i)=0d0                   ! Add forces here if any
      
      ran1=rnds(i) - 0.5d0
      ran2=rnds(2*i+1) - 0.5d0
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
!$omp barrier
!$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!
!!DO THE MEAN IN PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!

m1=0d0
m2=0d0

!$omp parallel private(i)
!$omp do reduction (+ : m1,m2) 
do i=1,n
   m1 = m1 + x(i)
   m2 = m2 + (x(i))**2 
end do
!$omp end do
!$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!
!!DO THE MEAN IN SERIAL
!!!!!!!!!!!!!!!!!!!!!!!!!
!   m1=sum(x)/dfloat(n)
!   m2=sum(x**2)/dfloat(n) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!WRITE TO FILE !!!!!!!!!!!
   do i=1,n
      write(11,'(1(I4),2(1x,e10.3) )') i,x(i),y(i)
   end do
   write(12,'(2(1x,e10.3))') m1,m2
!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(11,*) ' '
   t=t+dt
   print *,'t_max: ',t_max,'     ',' n: ', n
end do

end=omp_get_wtime()
wtime=end-begin

!$ threadnum=omp_get_num_threads()
!$ print *,'Thread num:',threadnum

! print *,n,wtime,threadnum

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,rnds)
! Close files
close(11)
close(12)

end program main
