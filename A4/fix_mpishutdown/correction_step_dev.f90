subroutine correction_step(z,z0,T,allocate_err)
! Sinlge correction step for the BVP, uses finite differencing
use globals
implicit none
integer :: k,l,piv(N+1),lapack_flag,allocate_err
double precision :: z(N+1),z0(N+1),T(N+1),eps,re(N+1),r(N+1),y(N+1)
double precision , allocatable, dimension (:,:) :: J 


!apple = 70000000
allocate(J(N+1, N+1), stat=allocate_err)
if (allocate_err /= 0 ) then
print *, "Allocation failed"
else 
print *, "Allocation works!"
end if  

eps=1d-4
call residual(z,z0,T,r)

! Build Jacobian matrix
do k=1,N+1
   y=z
   y(k)=y(k)+eps
   call residual(y,z0,T,re)
   J(:,k)=(re-r)/eps
end do

! Solve for the Newton-Raphson update step
r=-r
call DGESV(N+1,1,J,N+1,piv,r,N+1,lapack_flag)
if(lapack_flag.ne.0) print *,'LAPACK warning ',lapack_flag

! Apply update step
z=z+r

deallocate(J)

return
end subroutine correction_step
