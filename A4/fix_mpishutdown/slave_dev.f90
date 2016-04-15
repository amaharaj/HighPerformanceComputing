subroutine slave
use globals
use fft
use cont_problem
implicit none
logical :: init_fail,fail
integer :: received_tag,ierr,allocate_err
integer, dimension(mpi_status_size) :: status
double precision :: z0(N+1),z(N+1),T(N+1),r(N+1),res



! Initialize the FFT
init_fail=.false.
!666 if (allocate_err.ne.0) then
!       print *, "ALLOCATION ERRORRRR"
!       init_fail=.true.
!       call MPI_SEND(init_fail,1,MPI_INTEGER,0,DEFAULT_TAG,MPI_COMM_WORLD,ierr)
!    end if 
call setup_fft(ierr)
if(ierr.ne.0) then
   print *,'MKL FFT initialization trouble, exiting.'
   init_fail=.true.
end if

!call MPI_SEND(init_fail,1,MPI_INTEGER,0,DEFAULT_TAG,MPI_COMM_WORLD,ierr)
call mpi_allreduce(init_fail,fail,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
! End initialization

if(.not.(fail)) then
   print *,'Slave process starting...'
   do while(.true.)
      ! Receive initial guess for my step size
      call MPI_RECV(z0,N+1,MPI_DOUBLE_PRECISION,ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      received_tag = status (MPI_TAG)
      if(received_tag.eq.EXIT_TAG)  exit ! Exit when root tells me to
      ! Receive current approximation for my step size
      call MPI_RECV(z,N+1,MPI_DOUBLE_PRECISION,ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      call MPI_RECV(T,N+1,MPI_DOUBLE_PRECISION,ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      call correction_step(z,z0,T,allocate_err)
      if (allocate_err.eq.(0)) then 
        print *, "!!!!" 
         call MPI_SEND(init_fail,1,MPI_LOGICAL,0,DEFAULT_TAG,MPI_COMM_WORLD,ierr)
      else
         goto 666  
        ! call MPI_SEND(allocate_err,1,MPI_INTEGER,0,DEFAULT_TAG,MPI_COMM_WORLD,ierr)
      end if 
      call residual(z,z0,T,r)
      call MPI_SEND(z,N+1,MPI_DOUBLE_PRECISION,0, DEFAULT_TAG, MPI_COMM_WORLD, ierr)
      res=sqrt(sum(r**2))                ! 2-norm of the residual vector
      call MPI_SEND(res,1,MPI_DOUBLE_PRECISION,0, DEFAULT_TAG, MPI_COMM_WORLD, ierr)
   end do
end if

666 if (allocate_err.ne.0) then
       print *, "ALLOCATION ERRORRRR"
       init_fail=.true.
       call MPI_SEND(init_fail,1,MPI_LOGICAL,0,DEFAULT_TAG,MPI_COMM_WORLD,ierr)
    end if


call teardown_fft
return
end subroutine slave