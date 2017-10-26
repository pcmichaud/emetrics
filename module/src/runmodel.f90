program main
	use model
	integer nsim
	character*10 buff
	! gets the first argument of when the function is called to store number of obs
	call get_command_argument(1,buff)
	read(buff,'(I8)') nsim
	write(*,*) 'nsim = ', nsim
	call simulate(nsim) 
end program main