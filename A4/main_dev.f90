module globals
implicit none
include 'mpif.h'
integer, parameter :: m=2**12,N=m+1        ! Numbers of grid points and equations (there are N+1 unknowns)
integer :: DEFAULT_TAG=0,EXIT_TAG=1,ROOT=0 ! Default tag and exit tag for mpi_send, root rank
double precision :: VOID=0d0               ! Empty message for sending exit tag
end module globals

! Module for MKL FFT implementation, used by the slaves for the evaluation of the BVP
module fft
use MKL_DFTI
use globals
implicit none
integer :: fft_status = 0, ignored_status
type(DFTI_DESCRIPTOR), POINTER :: hand
contains
include 'setup_teardown_fft.f90'
end module fft

! Module containing the constants and subroutined for evaluation of the BVP
module cont_problem
use globals
use fft
implicit none
double precision :: A=8.09d0                        ! Parameter in the BVP we do not change
double precision,parameter :: PI=3.14159265358979d0 
complex(kind=8), parameter :: AI=(0.d0,1.d0)
contains
include 'residual_dev.f90'
include 'correction_step_dev.f90'
end module cont_problem

! Main program, only does MPI initalization/finalization and separates master and slaves
program main
use globals
implicit none
integer :: ierr,numprocs,proc_num

call mpi_initialize(numprocs,proc_num,ierr)

if(ierr.eq.0) then 
   if(proc_num.eq.ROOT) then
      call master(numprocs)
   else
      call slave
   end if
else
   print *,'MPI initialization failed, exiting.'
end if

call mpi_shut_down
end program main

subroutine mpi_initialize(numprocs,proc_num,ierr)
use globals
implicit none
integer ierr,numprocs,proc_num
! MPI initialization - can be tucked away in a subroutine.
    call mpi_init(ierr)
    if(ierr.ne.0) then
       print *,'mpi_init failed.'
       goto 1
    end if
    call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
     if(ierr.ne.0) then
       print *,'mpi_comm_size failed, ierr=',ierr
       goto 1
    end if
    call mpi_comm_rank(MPI_COMM_WORLD, proc_num, ierr)
    if(ierr.ne.0) then
       print *,'mpi_comm_rank failed, ierr=',ierr
       goto 1
    end if 
! Check if we have at least one slave
    if(numprocs.eq.1) then
       print *,'We need more than 1 process!'
       ierr=1
    end if
1 continue

  end subroutine mpi_initialize

subroutine mpi_shut_down
implicit none
logical :: init
integer :: ierr

call mpi_initialized(init,ierr)
if(init) call mpi_finalize(ierr)

return
end subroutine mpi_shut_down

! For easy compiling, include the files with the master and alve code
include 'master_dev.f90'
include 'slave_dev.f90'
