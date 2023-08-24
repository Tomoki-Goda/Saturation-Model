!module HistToStep
	
!end module HistToStep

program HistToStep
	implicit none
!	use HistToStep
	character (100) :: InputFile, OutputFile 
	integer :: flag,i
  real(kind(1d0)),dimension(1:4):: d1,d2,d3 
  
  if (COMMAND_ARGUMENT_COUNT().lt.4) then
  	print *, "Give 3 three arguments mid min max out"
  end if
  do i=0,2
  call get_command_argument(i+1,value=InputFile,status=flag)
	open(unit=10+i, file=InputFile,status='old', action='read')
!	open(unit=11, file=InputFile,status='old', action='read')
!	open(unit=12, file=InputFile,status='old', action='read')
  end do
	call get_command_argument(4,value=OutputFile,status=flag)
	open(unit=13, file=OutputFile,status='replace', action='write')
  do i=1,100
  read(10,*,end=999) d1
  read(11,*,end=999) d2
  read(12,*,end=999) d3
    write(13,101) d1(1),d1(3),d2(3),d3(3) 
    write(13,101) d1(2),d1(3),d2(3),d3(3) 
  end do
  101 format(e12.5,3(1x,e12.5))

  999 print*, 'End of file'
!  write(13,101) d1(2),d1(3),d2(3),d3(3) 

	close(10)
	close(11)
	close(12)
	close(13)
	
end program HistToStep
	
