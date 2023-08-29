!usage 
!./rescale.out ./gluon-grids/BGKSkt.dat ./gluon-grids/BGKSkt


program HistToStep
	implicit none
!	use HistToStep
	character (100) :: InputFile,OutputFile, OutputFile1, OutputFile2 
	integer :: flag,i
	real(kind(1d0)),dimension(1:4):: d1
	
	call get_command_argument(1,value=InputFile,status=flag)
	open(unit=10, file=InputFile,status='old', action='read')
	call get_command_argument(2,value=OutputFile,status=flag)
	OutputFile1=trim(OutputFile)//"05.dat"  
	OutputFile2=trim(OutputFile)//"2.dat"
	open(unit=11, file=OutputFile1,status='replace', action='write')
	open(unit=12, file=OutputFile2,status='replace', action='write')
 
 101 format(e12.5,3(1x,e12.5))
 100 read(10,*,end=999) d1
	write(11,101) d1(1),d1(2),d1(3)-log(4d0) ,d1(4)
    write(12,101) d1(1),d1(2),d1(3)+log(4d0) ,d1(4)
	go to 100
 999 print*, 'End of file'
 	close(10)
 	close(11)
 	close(12)
end program
 	
