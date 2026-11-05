module acquisition

  ! External modules
  use precision
  use grid_patch
  use parallel
  use state

  implicit none

  real(WP),allocatable,dimension(:,:) :: acqData
  real(WP),allocatable,dimension(:) :: acqTime
  type(t_Patch), pointer :: acqPatch => null()
  integer :: acqSubarrayType, acqLocalSize(3), acqGlobalSize(3),acqOffset(3),               &
  nacqPoints = 0, nacqPointsGlobal = 0, acqComm = MPI_COMM_NULL, acqRoot = -1
  integer :: currentStep = 1
  integer :: nSteps
  integer :: startIter
  integer :: stopIter
  integer :: skipIter

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
  logical, intent(in) :: debug

   ! Connect the acquisition patch if there is one
  do i = 1, nPatches
     if (patches(i)%patchType .eq. ACQUIRE_SUBREGION) then
        ! allocate(acqPatch)
        acqPatch => patches(i)
        exit
      end if
  end do

  if (.not. associated(acqPatch)) then
     if (iRank .eq. iRoot) print *, "No ACQUIRE_SUBREGION patch found: acquisition disabled"
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
  call parser_read('acquire subregion skip', skipIter, 1)
  if (skipIter .lt. 1) call die('acquisition_setup: "acquire subregion skip" must be >= 1')
  if (stopIter .lt. startIter) call die('acquisition_setup: "acquire subregion stop" must be >= "start"')
  nSteps = (stopIter - startIter) / skipIter + 1

  if (debug .and. iRank .eq. iRoot) then
     print *, "acquire subregion start: ", startIter
     print *, "acquire subregion stop:  ", stopIter
     print *, "acquire subregion skip:  ", skipIter
     print *, "total steps: ", nSteps
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

! real(WP),dimension(nGridPoints,1)::tmp
real(WP), allocatable :: tmp(:)

! Ranks with no patch points (or no patch at all) have nothing to do
if (.not. associated(acqPatch)) return
if (nAcqPoints .eq. 0) return

! Sample only inside [start, stop], every skip steps
if (timestep .lt. startIter .or. timestep .gt. stopIter) return
if (mod(timestep - startIter, skipIter) .ne. 0) return

! Never write past the arrays sized in acquisition_setup
if (currentStep .gt. nSteps) return

allocate(tmp(nGridPoints)) ! NOTE : this only returns pressure this way
call compute_dependent_variables(conservedVariables,pressure=tmp)    
call patch_collect(acqPatch, tmp, acqData(:,currentStep))

if (iRank .eq. acqRoot) acqTime(currentStep) = time
deallocate(tmp)

currentStep=currentStep+1

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
    call MPI_BARRIER(acqComm, ierror)
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
    call MPI_TYPE_FREE(acqSubarrayType, ierror)
    
  !delete data structures
  if (allocated(acqData)) deallocate(acqData)
  if (allocated(acqTime)) deallocate(acqTime)

  return
end subroutine acquisition_cleanup

