subroutine master(numprocs)
! Go over the pseudo-arclength continuation loop, distributing concurrent correction steps over processes.
use globals
implicit none
logical :: init_fail,fail,working,hang,allocation_error
integer :: i,i_max,k,ierr,numprocs,p,nnwt,nnwt_thresh,nkill,accepted,AllDiverge,AllConverge,AllWillConverge,OneSolution,Unsure,ConvergeWillConverge,CodeHangs
integer, allocatable, dimension(:) :: list,conv,willconv
double precision :: z(N+1),T(N+1),th(3),z_prev(N+1),xs(m),alpha,maxdelconv,maxdelwillconv,maxlocation(2),prodconv(numprocs-1),prodwillconv(numprocs-1),hmin
double precision, allocatable, dimension(:) :: del,r,rnew,rtemp,rold!,prodconv,prodwillconv
double precision, allocatable, dimension(:,:) :: pred,cur

allocate(del(numprocs-1),r(numprocs-1),rnew(numprocs-1),rtemp(numprocs-1),rold(numprocs-1),list(numprocs-1),conv(numprocs-1),willconv(numprocs-1))
allocate(pred(N+1,numprocs-1),cur(N+1,numprocs-1))

print *,'Master process starting...'
print *,'Nr. of slaves=',numprocs-1
! Initializations
init_fail=.false.
! Open file for output
open(110,file='continuation_curve')
! Read initial data from disk
open(11,file='z0_12',err=2)
open(111,file='Step_Stats')
read(11,*,err=2) z
close(11)
goto 1
2 print *,'Trouble reading intial data, exiting.'
init_fail=.true.
1 continue


! Set initial tangent vector: vary only wave speed c
T=0d0
T(N)=1d0
conv=0d0
willconv=0d0
prodconv=0d0
prodwillconv=0d0
! Set up arrays for residual comparisons
rold=0d0
rnew=0d0
rtemp=0d0
! Set the initial distribution of step sizes (experiment with this)
del(1)=5d0
alpha=1.9d0
! Initialize Counters
AllDiverge=0
AllConverge=0
AllWillConverge=0
OneSolution=0
Unsure=0
ConvergeWillConverge=0
CodeHangs=0
! Set threshold for number of newton iterations
nnwt_thresh=10

! Initialize Delta
if(numprocs.gt.2) then
   do i=2,numprocs-1
      del(i)=alpha*del(i-1)
   end do
end if
i_max=200

! Set the thresholds that determine whether an approximate solution is diverging,
! will likely converge after one more step or is converged
th(1)=500
th(2)=1d-4
th(3)=1d-9

! Everybody verifies that no initialization errors occured
call mpi_allreduce(init_fail,fail,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
! End initializations

if(.not.(fail)) then
   i=0
! Note: naming the loops (here inner and outter) allows us to cycle and exit either one)
outter:   do while(i.le.i_max)

      ! Predict using multiple step sizes
10    do k=1,numprocs-1
         pred(:,k)=z+del(k)*T
      end do
      cur=pred ! Set current approximations to initial guesses

!     Set up the list of processors to send data to, initially all of them
      do k=1,numprocs-1
         list(k)=k
      end do
      z_prev=z
      working=.true.
      nnwt=0
      ! Start correction step

inner: do while(working)
        ! set the previous residuals to compare against
        conv=0d0
        willconv=0d0       
        if ((nnwt.ge.0).AND.(i.gt.1)) then 
           do k=1,numprocs-1
              rtemp(k)=rold(k)
              r(k)=rold(k)
           end do 
        end if
         ! Send the predictions to the slaves
         do k=1,numprocs-1
            if(list(k).lt.0) cycle
            p=list(k)
            call MPI_SEND(pred(:,p),N+1,MPI_DOUBLE_PRECISION,p, DEFAULT_TAG, MPI_COMM_WORLD, ierr)
            call MPI_SEND(cur(:,p),N+1,MPI_DOUBLE_PRECISION,p, DEFAULT_TAG, MPI_COMM_WORLD, ierr)
            call MPI_SEND(T,N+1,MPI_DOUBLE_PRECISION,p,DEFAULT_TAG,MPI_COMM_WORLD,ierr)
         end do
         ! Check for allocation errors
         do k=1,numprocs-1
            if(list(k).lt.0) cycle
            p=list(k)
            call MPI_RECV(allocation_error,1,MPI_LOGICAL,p,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            if (allocation_error.eq.(.true.)) then
               exit outter
            end if
         end do
         ! Receive the results
         do k=1,numprocs-1
            if(list(k).lt.0) cycle
            p=list(k)
            call MPI_RECV(cur(:,p),N+1,MPI_DOUBLE_PRECISION,p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            call MPI_RECV(r(p),1,MPI_DOUBLE_PRECISION, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         end do
         nnwt=nnwt+1
         do k=1,numprocs-1
            rnew(k)=r(k)
         end do 
      ! Make a decision (experiment and expand this, now written for two slaves only)
      !  First eliminate diverging approximate solutions
         nkill=0
         rold=0d0
         do k=1,numprocs-1
            rold(k) = r(k)
            if ((nnwt.gt.0).AND.(i.gt.1).AND.(hang.eq..True.)) then 
               if (rnew(k).gt.rtemp(k)) then
                  ! kill if residual increases between steps
                  r(k) = th(1) + 100
                  nkill=nkill+1
               end if 
            end if 
            hang=.False.
            if(r(k).gt.th(1)) then
               list(k)=-1
               nkill=nkill+1
            else if (r(k).lt.th(3)) then
               conv(k) = 1
            else if ((r(k).lt.th(2)).AND.(r(k).gt.th(3))) then
               willconv(k) = 1
            end if
         end do
      ! If all approximate solutions diverged, retry the step with smaller step sizes
         if(nkill.eq.numprocs-1) then
            AllDiverge = AllDiverge + 1
            del=del/alpha**(numprocs-1) 
            cycle outter ! Go to the next iteration of the outter loop without adding to point counter i
         end if
      ! If only one approximate solution is left then ...
         if(numprocs-1-nkill.eq.1) then
            OneSolution = OneSolution + 1
            k=minloc(r,1)
      ! ... if it has not yet converged, take another Newton-Raphson step,
            if(r(k).gt.th(3).AND.nnwt.lt.nnwt_thresh) then
               cycle inner
            else
      ! ... if it has converged, accept it as the next point on the continuation curve.
               accepted=k
               working=.false.
               cycle inner
            end if
         end if

      !  ... if one has a residual we are not sure about (i.e. greater than th(2)), take another Newton-Raphson step,
            
            if(maxval(r).gt.th(2).AND.(maxval(r).lt.th(1))) then
               hang=.True.
               if (hang.eq..True.) then 
                  CodeHangs = CodeHangs + 1
               end if 
               cycle inner
            end if

      ! If neither approximate solution has diverged, then ...
      !   if(nkill.eq.0) then
      ! ... if ALL have converged, select the one with the largest step size,
            if(maxval(r).lt.th(3)) then
               AllConverge = AllConverge + 1
               accepted=maxloc(del,1)
               del(1)=del(accepted)               ! Increase all step sizes
               do k=2,numprocs-1
                  del(k)=alpha*del(k-1)
               end do  
               working=.false.
               cycle inner
            end if
     !  ... Comments
            prodconv=0d0
            prodwillconv=0d0
            prodconv = conv*del
            prodwillconv = willconv*del
            maxdelconv = maxval(prodconv,1)
            maxdelwillconv = maxval(prodwillconv,1)
            maxlocation(1) = maxloc(prodconv,1)
            maxlocation(2) = maxloc(prodwillconv,1)
            if (sum(willconv).eq.(numprocs-1)) then
               AllWillConverge = AllWillConverge + 1
               accepted=maxloc(prodwillconv,1)
               del(1)=del(accepted)               ! Increase all step sizes
               do k=2,numprocs-1
                  del(k)=alpha*del(k-1)
               end do 
               working=.false.
               cycle inner
            end if
            if ((sum(conv)).AND.(sum(willconv)).gt.(0)) then 
      !  ... if one has converged and one will likely converge at the next step, select the one the gives the greatest
               ConvergeWillConverge = ConvergeWillConverge + 1      
      !      increase in arclength.
                  if(maxdelconv/real(nnwt,8).gt.maxdelwillconv/real(nnwt+1,8)) then
                     accepted=maxloc(prodconv,1)
                     del(1)=del(accepted)               ! Increase all step sizes
                     do k=2,numprocs-1
                        del(k)=alpha*del(k-1)
                     end do 
                     working=.false.
                     cycle inner
                  else if (maxdelconv/real(nnwt,8).lt.maxdelwillconv/real(nnwt+1,8)) then
                     accepted=maxloc(prodwillconv,1)
                     del(1)=del(accepted)               ! Increase all step sizes
                     do k=2,numprocs-1
                        del(k)=alpha*del(k-1)
                     end do
                     working=.false.
                     cycle inner
                  else
                     accepted=maxloc(maxlocation,1)
                     del(1)=del(accepted)               ! Increase all step sizes
                     do k=2,numprocs-1
                        del(k)=alpha*del(k-1)
                     end do
                     working=.false.
                     cycle inner
                  end if
            else 
               working=.false.
               cycle inner
            end if 
      end do inner

      ! Accept step
      z=cur(:,accepted)
      print *,'Step accepted, lambda=',z(N+1)
!     Write lambda, c, solution norm and step size
      write(110,'(i5,4(1x,e17.10))') i,z(N+1),z(N),sqrt(sum(z(1:m)**2)),del(accepted)
      write(111,'(8(2I4))') i,AllDiverge,AllConverge,AllWillConverge,OneSolution,Unsure,ConvergeWillConverge,Hang
      ! Recompute T
      T=(z-z_prev)/sqrt(sum((z-z_prev)**2))
      ! Add to point counter
      i=i+1
   end do outter
   call stop_slaves(numprocs)
end if

close(110)
close(111)
deallocate(del,pred,r,cur,list,conv,willconv,rtemp,rold,rnew)
end subroutine master

subroutine stop_slaves(numprocs)
use globals
implicit none
integer :: k,numprocs,ierr

do k=1,numprocs-1
   call MPI_SEND(void,1,MPI_DOUBLE_PRECISION,k, EXIT_TAG, MPI_COMM_WORLD, ierr)
end do

return
end subroutine stop_slaves
