module acquisition

  ! External modules
  use precision
  use grid_patch
  use parallel
  use state

  implicit none

  real(WP),allocatable,dimension(:,:) :: acqData
  real(WP),allocatable,dimension(:) :: acqTime
  type(t_Patch), pointer :: acqPatch
  integer :: acqSubarrayType, acqLocalSize(3), acqGlobalSize(3),acqOffset(3),               &
  nacqPoints, nacqPointsGlobal,acqComm, acqRoot
  integer::currentStep
  integer::nSteps
  integer::startIter
  integer::stopIter
  integer::skipIter

end module acquisition


! ==========================  !
! Setup the acquisition patch !
! ==========================  !
subroutine acquisition_setup(debug)

  ! Internal modules
  use acquisition

  ! External modules
  use parser

  implicit none
 
  integer :: i  
  integer :: starts(3), extent(6)
  real(WP), allocatable :: buffer_WP(:,:)
  integer :: ierror
  logical,optional::debug

   ! Connect the acquisition patch if there is one
  do i = 1, nPatches
     if (patches(i)%patchType .eq. ACQUIRE_SUBREGION) then
        allocate(acqPatch)
        acqPatch => patches(i)
        exit
      end if
  end do

  if (.not. associated(acqPatch)) then
    print *,"No ACQUIRE_SUBREGION patch exists: Exiting"
    return
  end if

  ! Set the visualization data environment
  if (associated(acqPatch)) then
     acqLocalSize = acqPatch%localSize
     acqGlobalSize = acqPatch%globalSize
     acqOffset = acqPatch%offset
     acqComm = acqPatch%comm
     acqRoot =  acqPatch%masterRank
  end if

  nAcqPoints = product(acqLocalSize)
  nAcqPointsGlobal = product(acqGlobalSize)

  call parser_read('acquire subregion start', startIter, 0)
  call parser_read('acquire subregion stop', stopIter, 0)
  call parser_read('acquire subregion skip', skipIter, 0)
  nSteps=abs(startIter-stopIter)/skipIter+1

  if (.not. present(debug)) debug=.false.
  if (debug) then
     print *,"acquire subregion start: ",startIter
     print *,"acquire subregion start: ",stopIter
     print *,"acquire subregion start: ",skipIter
     print *,"total steps: ",nSteps
  end if


  if ( nAcqPoints .gt. 0) then
     if (associated(acqPatch)) then
        extent = (/ acqPatch%iMin, acqPatch%iMax, acqPatch%jMin, acqPatch%jMax,              &
             acqPatch%kMin, acqPatch%kMax /)
        starts = acqPatch%offset - extent(1::2) + 1
        call MPI_TYPE_CREATE_SUBARRAY(3, acqGlobalSize, acqLocalSize, starts,                &
             MPI_ORDER_FORTRAN, MPI_REAL_WP, acqSubarrayType, ierror)
        call MPI_TYPE_COMMIT(acqSubarrayType, ierror)
     else
        call MPI_TYPE_CREATE_SUBARRAY(3, acqGlobalSize, acqLocalSize, acqOffset,             &
             MPI_ORDER_FORTRAN, MPI_REAL_WP, acqSubarrayType, ierror)
        call MPI_TYPE_COMMIT(acqSubarrayType, ierror)
     end if

     if (.not. allocated(acqData)) allocate(acqData(nAcqPoints,nSteps))
     acqData=0._WP
     if (.not. allocated(acqTime)) allocate(acqTime(nSteps))
     acqTime=0._WP

     if (debug) then
     print *,"acquire subregion nAcqPoints: ",nAcqPoints
  end if

     currentStep=1

  end if

  return
end subroutine acquisition_setup


! ========================== !
! Acquire the data on the patch
! ========================== !
subroutine acquire_data

use state
use equation_of_state
use acquisition
use time_info 

implicit none

real(WP),dimension(nGridPoints,1)::tmp

  if (iRank .eq. acqRoot) then
    acqTime(currentStep)=time
  end if

if (.not. associated(acqPatch)) then
    return
end if

if (timestep .ge. startIter .and. timestep .le. stopIter) then
 if (mod(timestep-startIter,skipIter) .eq. 0) then

! ------------------------------------------------------------------------------------------
! LUCAS: in this version of jCODE  variables 'forwardMode' and 'baseState' do not exist
! only consider NONLINEAR case, with conserved variables as is
! select case (trim(forwardMode))
!       case ('LINEARIZED')
!        call compute_dependent_variables(conservedVariables+baseState,pressure=tmp(:,1))
!       case ('NONLINEAR')
!        call compute_dependent_variables(conservedVariables,pressure=tmp(:,1))    
!  end select

  ! LUCAS only consider NONLINEAR case, with conserved variables as is
  call compute_dependent_variables(conservedVariables,pressure=tmp(:,1))    
! ------------------------------------------------------------------------------------------

  if (nAcqPoints .eq. 0) return

  if (associated(acqPatch)) then
   call patch_collect(acqPatch, tmp(:,1), acqData(:,currentStep))
  end if
  currentStep=currentStep+1 !incrementcurrentStep=currentStep+1 !increment
 end if
end if


end subroutine

! ============================ !
! Cleanup the acquisition patches !
! ============================ !
subroutine acquisition_cleanup

  ! Internal modules
  use acquisition
  use time_info

  implicit none

   logical::fileexists
  integer(kind = MPI_OFFSET_KIND) :: offset
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK, nt_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, t_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK
  integer :: ierror, iunit
  character(len = str_medium) :: filename
  integer::t
  logical::fileIsTHere

   integer, dimension(4) :: dims

   if (nAcqPoints .eq. 0) return

      ! Prefix the filename to prevent potential file locking
    filename = trim(mpiiofs) // trim('data.subregion')

    ! Open the file to write
    inquire(file = filename, exist = fileIsThere)
    if (fileIsThere .and. irank.eq.iroot) call MPI_FILE_DELETE(filename, mpiInfo, ierror)
    call MPI_FILE_OPEN(acqComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo,         &
       iunit, ierror)

    ! Write header (root process only)
    if (iRank .eq. acqRoot) then
     dims(1:3) = acqGlobalSize
     dims(4) = nSteps
     call MPI_FILE_WRITE(iunit, dims, 4, MPI_INTEGER, status, ierror)
     call MPI_FILE_WRITE(iunit, acqTime, nSteps, MPI_REAL_WP, status, ierror)
    end if

    ! Resize some integers so MPI can read even the biggest files
    nx_MOK    = int(acqGlobalSize(1), MPI_Offset_kind)
    ny_MOK    = int(acqGlobalSize(2), MPI_Offset_kind)
    nz_MOK    = int(acqGlobalSize(3), MPI_Offset_kind)
    nt_MOK    = int(nSteps, MPI_Offset_kind)
    WP_MOK    = int(WP,               MPI_Offset_kind)

    ! Write the file
    if (nAcqPoints .gt. 0) then
     do t = 1, nSteps
     t_MOK = int(t, MPI_Offset_kind)
     offset = 4*4 + nt_MOK*WP_MOK + nx_MOK*ny_MOK*nz_MOK*WP_MOK*(t_MOK-1)
     call MPI_FILE_SET_VIEW(iunit, offset,MPI_REAL_WP, acqSubarrayType, "native",            &
          mpiInfo, ierror)
     call MPI_FILE_WRITE_ALL(iunit, acqData(:,t), nAcqPoints, MPI_REAL_WP,                  &
          MPI_STATUS_IGNORE, ierror)
     end do
    end if


    ! Close the file.
    call MPI_FILE_CLOSE(iunit, ierror)

  !delete data structures
  if (allocated(acqData)) deallocate(acqData)

  return
end subroutine acquisition_cleanup

