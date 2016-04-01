module globals
! Latest version: Feb 29 2016; starting point of assignment 2.
! Global variables
use omp_lib
implicit none
integer :: n, d                                ! number of particles
double precision :: L=1                         ! domain size with periodic boundary conditions
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
double precision :: dt,kT,g,m,sigma,eps,rc      ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy
integer :: neighbours(8)
!integer, allocatable, dimension(:) :: neighbours
contains
subroutine set_parameters

! Set time step and physical parameters
dt=1d-5
kT=1.0d0
g=1.0d0
m=1.0d0

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

! Set Lennard-Jones parameters
sigma=1d-4
eps=1d0
rc=sigma*2d0**(1d0/6d0)

end subroutine set_parameters
subroutine initialize_particles
integer :: i
double precision :: ran1,ran2,gr1,gr2

! Give particles initial position and velocity.
! In this version, the particles a randomly placed in a square of length L/20,
! you may want to change this to take into account the number and size of the particles. 
do i=1,n
   call random_number(ran1)
   call random_number(ran2)
   x(i)=(ran1-0.5d0)*L/10d0
   y(i)=(ran2-0.5d0)*L/10d0
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2d0*log(ran1))*cos(2d0*pi*ran2)
   gr2=sqrt(kT/(m))*sqrt(-2d0*log(ran1))*sin(2d0*pi*ran2)
   vx(i)=gr1
   vy(i)=gr2
   vhx(i)=0d0
   vhy(i)=0d0
end do

end subroutine initialize_particles
end module Langevin

include 'qsort_deb.f95'

module geometry
use globals
use qsort_c_module
implicit none
integer :: ndim=2,nn=9   ! Number of spatial dimensions, sub domains per spatial dimension, number of neighbours for each cell (including the cell itself).
contains
subroutine periodic_bc(x,y)
! Implement periodic boundary conditions.
! Note, that if aparticle jumps farther than distance L in one time step, it is excluded from further computation.
double precision :: x,y
if(x.lt.-L/2d0) x=x+L
if(x.gt.L/2d0) x=x-L
if(y.lt.-L/2d0) y=y+L
if(y.gt.L/2d0) y=y-L
return
end subroutine periodic_bc

subroutine dist(x1,y1,x2,y2,md,xg,yg)
integer :: i,in
double precision :: x1,x2,y1,y2,md,xg,yg,dst(nn),X(nn,ndim) ! Distance renamed to "dst" and given dimension nn.

X(1,:)=(/ x2, y2 /)
X(2,:)=(/ x2+L, y2 /)
X(3,:)=(/ x2-L, y2 /)
X(4,:)=(/ x2, y2+L /)
X(5,:)=(/ x2, y2-L /)
X(6,:)=(/ x2+L, y2+L /)
X(7,:)=(/ x2-L, y2+L /)
X(8,:)=(/ x2+L, y2-L /)
X(9,:)=(/ x2-L, y2-L /)
do i=1,9
   dst(i)=norm2((/ x1, y1/)-X(i,:))
end do
in=minloc(dst,1)
md=dst(in)
xg=X(in,1)
yg=X(in,2)

return
end subroutine dist

subroutine divide(x,y,in)
integer :: in,i1,i2
double precision :: x,y

if((abs(x).le.L/2).and.(abs(y).le.L/2)) then  ! Particle in domain, assign domain index.
   i1=floor((x+L/2d0)/(L/float(d)))
   i2=floor((y+L/2d0)/(L/float(d)))
   in=i1+i2*d
else
   in=d*d                                     ! Particle outside bounds, put in dumpster.
end if

return
end subroutine divide

subroutine sort(x,y,vhx,vhy,in,lin,nt)
! Input: unsorted list of positions and domain indices, nr. of threads.
! Output: sorted lists and list of loop indeces so that thread k processes lin(k,1) to lin(k,2).
integer :: i,nt
integer :: in(n)
integer :: lin(0:nt,2)
double precision :: x(n),y(n),vhx(n),vhy(n)

call QsortC(in,x,y,vhx,vhy) ! quicksort lists (library routine adjusted to change the order of in, x, y, vhx and vhy simultaneously)

! Now go over the sorted list and find the first and last index for each block of domain index.
lin=0
lin(in(1),1)=1

do i=2,n
   if(in(i).gt.in(i-1)) then
      lin(in(i-1),2)=i-1
      lin(in(i),1)=i
   end if
end do
lin(in(n),2)=n
! Note: the highest sector/processor number is d*d-1, index d*d is used for "escaped" particles.

return
end subroutine sort

subroutine nnCheck(id,neighbours)
integer :: i,n , id, val
integer :: tr, tm, tl, mr, mm, ml, br, bm, bl
integer, dimension(0:(d**2-1)) :: a, shift_left, shift_right, test
integer :: neighbours(8)


! Define main array and shifted arrays
n = (d**2 - 1)
do i=0,n
   a(i) = i
end do 
shift_left = cshift(a, shift=-1, dim=1)
shift_right = cshift(a, shift=1, dim=1)

test = shift_left
val = d-3

mm = id 

! Check right edge, else move one space
   if (modulo(id,d).eq.(d-1)) then
      tr = shift_left(modulo((id+d-1),size(a))) - val
      br = shift_left(modulo((id-d-1),size(a))) - val
      mr = shift_left(modulo((id-1),size(a))) - val  
   else  
      mr = a(modulo((id+1),size(a))) 
      tr = a(modulo((id+d+1),size(a)))
      br = a(modulo((id-d+1),size(a)))
   end if

    ! Check left edge, else move one space
    if (modulo(id,d).eq.0) then 
        tl = shift_right(modulo((id+d+1), size(a))) + val
        bl = shift_right(modulo((id-d+1),size(a))) + val
        ml = shift_right(modulo((id+1),size(a))) + val
    else 
        ml = a(modulo((id-1),size(a)))
        tl = a(modulo((id+d-1),size(a)))
        bl = a(modulo((id-d-1),size(a)))
    end if 
    
    ! Check top edge
    if (dble(id)/dble(d).ge.2) then 
        tm = shift_right(modulo((id+d-1),size(a)))
    else 
        tm = a((id+d)) 
    end if
        
    ! Check bottom edge
    if (dble(id)/dble(d).le.1) then 
        bm = shift_left(modulo((id-(d-1)),size(a)))
    else 
        bm = a((id-d)) 
    end if 

neighbours = (/ tl, tm, tr, ml, mr, bl, bm, br /)

!print *, neighbours(1)
 
!print *, neighbours
!print *, "************"
!print *, "  "
!print *, 'tl:', tl, 'tm:', tm, 'tr:', tr
!print *, 'ml:', ml, 'mm:', mm, 'mr:', mr
!print *, 'bl:', bl, 'bm:', bm, 'br:', br
!print *, "  "

return

end subroutine


end module geometry

program main
use globals
use Langevin
use geometry
implicit none
integer :: i,j,p,neigh,nt,mytn,step,nstep,left
double precision :: t,t_max,m1,m2,msd,begin,end,wtime,r,F,xg,yg
double precision, allocatable, dimension(:) :: xc,yc
double precision, allocatable, dimension(:,:) :: ran
integer, allocatable, dimension(:) :: in,inc
integer, allocatable, dimension(:,:) :: lin

character(len=32) :: arg
call getarg(1,arg)
read (arg(1:20),'(i20)') n
call getarg(2,arg)
read (arg(1:20),'(i20)') d


! Open files
open(11,file='trajectories')
open(12,file='means')
open(13,file='walltimes')

! Allocate arrays
nt=1 ! default in absence of -fopenmp flag
!$ nt=d*d
!$ call omp_set_num_threads(nt)
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),in(n),lin(0:nt,2),ran(n,2),xc(n),yc(n),inc(n))

left=0
t=0d0
t_max=1d-01

call set_parameters
call initialize_particles

nstep=1000
!nstep=int(t_max/dt)
!nstep=100
print *,'Taking ',nstep,' time steps of size ',dt
print *,'Number of sections per dimension=',d

do i=1,n
   call divide(x(i),y(i),in(i))
end do

call sort(x,y,vhx,vhy,in,lin,nt)
xc=x                                                                  ! Copy positions for writing to disk.
yc=y
inc=in

! CPU clock
call cpu_time(begin)
!$ begin=omp_get_wtime()

write(11,'(a7,i6)') '# Step ',0
do i=1,n
   write(11,'(2(1x,e14.7),2i6)') x(i),y(i),i,in(i)
end do
write(11,*) ' '

!$omp parallel private(i,j,p,mytn,step,r,xg,yg,F,neighbours)
!$ mytn=omp_get_thread_num()

do step=1,nstep
! Update particle positions and domain indices.
   if(lin(mytn,2).gt.lin(mytn,1)) then ! If there are particles in my domain..
      do i=lin(mytn,1),lin(mytn,2)
         vhx(i)=vx(i)+0.5*ax(i)*dt
         vhy(i)=vy(i)+0.5*ay(i)*dt
         x(i)=x(i)+vhx(i)*dt
         y(i)=y(i)+vhy(i)*dt
         call periodic_bc(x(i),y(i))
         call divide(x(i),y(i),in(i))
      end do
   end if

!$omp single
   t=t+dt
   !print *, t
!$omp end single
! Barrier implied by "end single", threads must be synchronized prior to sorting.

! Sort the particle lists, prefetch random numbers and write to disk.
! The attributes for "left" ensure that the process wiriting to disk uses the value for the previous
! time step, while the process sorting computes the new value.
!$omp sections firstprivate(left) lastprivate(left)
!$omp section
   call random_number(ran)
   ran=ran-0.5d0
!$omp section 
!   if(mod(step,2).eq.0) then                                ! Magic number alert: data written every other step, you may want to change this.
      write(11,'(a5,i6)') 'Step ',step
      do i=1,n
         write(11,'(2(1x,e14.7),2i6)') xc(i),yc(i),i,inc(i) ! Note that i and inc(i) are here for debuggin only.
      end do
      write(11,*) ' '
 !  else
      m1=sum(xc(1:n-left))/dfloat(n-left)
      m2=sum(xc(1:n-left)**2)/dfloat(n-left)
      msd=m2-m1**2
      write(12,'(3(1x,e10.3))') t,m1,m2,msd
 !  end if
!$omp section
   call sort(x,y,vhx,vhy,in,lin,nt)
   if(lin(nt,2).gt.lin(nt,1)) left=lin(nt,2)-lin(nt,1)+1                 ! Number of particles that escaped.
!$omp end sections 
   xc=x                                                                  ! Copy positions for writing to disk.
   yc=y
   inc=in

! A barrier is implied by the "end sections" directive, last value of "left" has now been shared by the sorting process.

   if(lin(mytn,2).gt.lin(mytn,1)) then             ! If there are particles in my domain..
      do i=lin(mytn,1),lin(mytn,2)                 ! Add external forces here if any
         ax(i)=0d0
         ay(i)=0d0                  
      end do
   end if

   if(lin(mytn,2).gt.lin(mytn,1)) then             ! If there are particles in my domain..
      do i=lin(mytn,1),lin(mytn,2)                 ! check all particle pairs inside my domain
         do j=i+1,lin(mytn,2)
            call dist(x(i),y(i),x(j),y(j),r,xg,yg)
            if(r.lt.rc) then
               F=-4d0*eps*( -12d0*sigma**12/r**13 + 6D0* sigma**6/r**7 );
               ax(i)=ax(i)+F*(x(i)-xg)/(r*m)
               ay(i)=ay(i)+F*(y(i)-yg)/(r*m)
               ax(j)=ax(j)+F*(xg-x(i))/(r*m)
               ay(j)=ay(j)+F*(yg-y(i))/(r*m)
             !  print *, ax(i)
!               print *,'pair ',i,j,F*dt**2       ! Enable to detect large interactions (F*dt^2>L)

            end if
         end do
          
       !  do neigh=0,(d**2-1)
       !     call nnCheck(neigh ,neighbours)
       !     print *, "Neighbours: ", neighbours
       !     print *, "mytn: ", neigh
       !  end do 


         call nnCheck(mytn, neighbours)


         do p=1,(nn-1) 
            do j=lin(neighbours(int(p)), 1),lin(neighbours(int(p)) ,2)
                call dist(x(i),y(i),x(j),y(j),r,xg,yg)
     !           print *, "r: ", r
                if(r.lt.rc) then
                   F=-4d0*eps*( -12d0*sigma**12/r**13 + 6D0* sigma**6/r**7 );
                   ax(i)=ax(i)+F*(x(i)-xg)/(r*m)
                   ay(i)=ay(i)+F*(y(i)-yg)/(r*m)
                !   print *, "Neighbours: ", neighbours
                !   print *, "mytn: ", mytn
                !   print *, "r: ", r
                end if
            end do     
         end do

         ax(i)=ax(i)-pref1*vhx(i)+pref2*ran(i,1)
         ay(i)=ay(i)-pref1*vhy(i)+pref2*ran(i,2)
         vx(i)=vhx(i)+0.5*ax(i)*dt
         vy(i)=vhy(i)+0.5*ay(i)*dt
   !      print *, x(i), "t: ", t/dt
 
      end do
   end if

end do
!$omp end parallel
call cpu_time(end)
!$ end=omp_get_wtime()
wtime=end-begin

print *,'Complete... lost ',left,' particles.'
print *,'Initial number of particles was ',n,' wall time was ',wtime
write(13,'(3(1x,e10.3))') wtime

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,vhx,vhy,in,lin,ran,xc,yc,inc)!,neighbours)
! Close files
close(11)
close(12)
close(13)
end program main
