subroutine residual(z,z0,T,r)
! Compute the residual of the BVP.
! -c w' + w w' + w'' + lambda w'''' - A sin(w)=0 with periodic BC
! In: vector of unknowns z=(Re(a_0),Re(a_1),Im(a_1),...,Re(a_{m/2-1}),Im(a_{m/2-1}),Re(a_m/2),c,\lambda)
! In: initial guess for sequence of Newton-Raphson iterations, z0; approximate tangent T to the continuation curve z(s).
! Out: residual vector r of length N+1=m+2.
integer :: fft_status,k
double precision :: z0(N+1),T(N+1),z(N+1),r(N+1),w(m),wx(m),am(m)
double precision :: lambda,c,c0
c0=z0(N)
c=z(N)
lambda=z(N+1)

! Evaluate linear part
r(1:m)=-c*sd(z(1:m))            ! -c w'
r(1:m)=r(1:m)+sl(z(1:m),lambda) ! +w''+lambda w''''

! Evaluate advection term
fft_status = DftiComputeBackward(hand, z(1:m), w)
fft_status = DftiComputeBackward(hand, sd(z(1:m)), wx)
wx=w*wx
fft_status = DftiComputeForward(hand, wx, am)
r(1:m)=r(1:m)+am

! Evaluate sin(w)
wx=-A*sin(w)
fft_status = DftiComputeForward(hand, wx, am)
r(1:m)=r(1:m)+am

! First phase condition
r(N)=dot_product(z(1:m)-z0(1:m),sd(z0(1:m)))+(c-c0)

! Second phase condition
r(N+1)=dot_product((z-z0),T)

return
end subroutine residual

function sd(a)
! Spectral derivative
integer :: k
double precision :: sd(m),a(m)

sd(1)=0d0
do k=1,m/2-1
   sd(2*k)=-a(2*k+1)*real(k,8)
   sd(2*k+1)=a(2*k)*real(k,8)
end do
sd(m)=0d0

return
end function sd

function sl(a,lambda)
! Symmetric linear part
integer k
double precision :: sl(m),a(m),q,lambda

sl(1)=0d0
do k=1,m/2-1
   q=lambda*real(k,8)**4 - real(k,8)**2
   sl(2*k)=q*a(2*k)
   sl(2*k+1)=q*a(2*k+1)
end do
q=lambda*real(m/2,8)**4 - real(m/2,8)**2
sl(m)=q*a(m)

return
end function sl
