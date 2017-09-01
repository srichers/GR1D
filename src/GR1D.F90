!-*-f90-*-
program GR1D	
    	  
  use GR1D_module
  implicit none

  include 'mpif.h'
  integer ierr

  call MPI_INIT ( ierr )
  call MPI_COMM_RANK (MPI_COMM_WORLD, myID, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, Nprocs, ierr)
  
  !Welcome to GR1D
  if(myID==0) then
  write(*,*) "#################################################"
  write(*,*) "#################################################"
  write(*,*) "########### GR1D SPHERICAL HYDRO v2 #############"
  write(*,*) "######### Now with Neutrino Transport ###########"
  write(*,*) "################# Nov ??, 2014 ##################"
  write(*,*) "#################################################"
  endif

  ! Call problem setup and allocate/initialize variables 
  call start
  if(myID==0) then
  write(*,*) "Done with initial data :-)"

  write(*,*) "Begin time integration loop:"
  endif
  IntegrationLoop: do 

     call SetTimeStep

     call handle_output
     
!!   Integrate
     call Step(dt)

     call postStep_analysis
     
  enddo IntegrationLoop
      
  if(myID==0) then
  write(*,*) "Shutting down!"
  write(*,*) " "
  endif

  call MPI_FINALIZE ( ierr )
  
end program GR1D
