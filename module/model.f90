! shell that tells how to set up a module in fortran

module model 
	! if you do not declare a var and use it, Fortran would assume integer, protect against that by having no implicit
	! code will block at compiling, throw error message
	implicit none

	! declare your global variables (that can be used anywhere in code)
	double precision g_sigma
	! or put them into a header file
	include 'model.h'

	! use allocatable arrays for globals for which sizes will be determined later 
	double precision, allocatable :: g_params(:)

	! use derived types for hybrid constructs
	type person 
		integer id 
		integer age 
		logical married
		double precision consumption
	end type person 
	
	! sets by default that all functions and subroutines are private
	private
	! makes exception for some functions that can be called from outside
	public :: simulate

	contains

	! you can define either subroutines or functions
	! this is the routine that will be called from outside, setting number of simulated individuals
	subroutine simulate(nsim)
		integer nsim
		type (person), allocatable :: pop(:)
		g_nsim = nsim

		call initparams

		allocate(pop(g_nsim*g_nages))
		call dosimulate(pop)

	end subroutine simulate

	! read parameters from file to memory (useful, because no need to recompile if change)
	subroutine initparams
		integer i
		character*80 buffer
		open(1,file='params.info')
			read(1,*) buffer, g_npar
			allocate(g_params(g_npar))
			do i = 1, g_npar, 1
				read(1,*) buffer, g_params(i)
			end do 
			read(1,*) buffer, g_agemin
			read(1,*) buffer, g_agemax
			g_nages = g_agemax - g_agemin + 1
		close(1)
	end subroutine initparams

	! the core simulator, writes simulated data to file
	subroutine dosimulate(pop)
		type (person) pop(g_nsim*g_nages)
		integer i, pp, a
		double precision draw
		pp = 1
		open(1,file='simulatedpop.dat')
		do i = 1, g_nsim, 1
			pop(pp)%id = i
			pop(pp)%age = g_agemin
			call random_number(draw)
			if (draw .lt. 0.5d0) then
				pop(pp)%married = .true.
			else
				pop(pp)%married = .false.
			end if
			pop(pp)%consumption = 10.0d3
			pp = pp + 1
			do a = 2, g_nages, 1
				pop(pp)%id = pop(pp-1)%id
				pop(pp)%age = pop(pp-1)%age + 1
				pop(pp)%married = pop(pp-1)%married
				pop(pp)%consumption = pop(pp-1)%consumption * (1.0d0 + 0.01d0)
				write(1,*) pop(pp)%id, pop(pp)%age, pop(pp)%married, pop(pp)%consumption
				pp = pp + 1
			end do
		end do
		close(1)
		write(*,*) 'number of simulated observations ', pp-1, ' of ', g_nsim * g_nages

	end subroutine dosimulate


end module model
