
program HistToStep
  implicit none
!	use HistToStep
  character (100) :: InputFile, OutputFile 
  integer :: flag,i,j,counter=0,argc
  real(kind(1d0)),dimension(1:4,1:3):: d1 
  logical :: exist=.false.

!  if (COMMAND_ARGUMENT_COUNT().lt.4) then
!    print *, "Give 3 three arguments mid min max out"
!  end if
  argc=command_argument_count()  
  
  do i=0,2
    if( i+1.gt.argc) exit

    call get_command_argument(i+1,value=InputFile,status=flag)
    if(argc.eq. 4) then
      Inquire(FILE=InputFile,Exist=exist)
    end if
       
    if (exist) then
      write(6,*) "Read ",inputfile
      open(unit=10+i, file=InputFile,status='old', action='read')
      counter=counter+1
    end if
  end do
  call get_command_argument(4,value=OutputFile,status=flag)
  open(unit=13, file=OutputFile,status='replace', action='write')
  write(6,*) "Write ", outputfile

! 100 continue
  do j=1,100
  do i=1,counter
    read(10+i-1,*,end=999) d1(:,i)
  end do
  write(13,101) d1(1,1),d1(3,1),d1(3,2),d1(3,3) 
  write(13,101) d1(2,1),d1(3,1),d1(3,2),d1(3,3) 
!  go to 100
  end do
  101 format(e12.5,3(1x,e12.5))

  999 print*, 'End of file'
!  write(13,101) d1(2),d1(3),d2(3),d3(3) 

  close(10)
  close(11)
  close(12)
  close(13)

end program HistToStep
