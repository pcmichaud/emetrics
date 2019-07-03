module chain
	use mcmc
	implicit none

	integer n, k
	double precision, allocatable:: par(:), x(:,:), y(:), xb(:)

	contains

	subroutine initdata
		integer i, j
		double precision q, draw
		! simulate data for a probit model
		n = 10000
		k = 5
		allocate(par(k))
		q = 0.2d0
		
		do i = 1, k, 1
			par(i) = q*dble(i)
			q = - q
		end do

		! generate x
		allocate(x(n,k))
		do i = 1, n, 1
			x(i,1) = 1.0d0
			do j = 2, k, 1
				call random_number(draw)
				x(i,j) = quann(draw)
			end do
		end do

		! generate the probit outcome
		allocate(xb(n), y(n))
		xb = matmul(x,par)
		do i = 1, n, 1
			call random_number(draw)

			if ((xb(i) + quann(draw)) .gt. 0.0d0) then
				y(i) = 1.0d0
			else
				y(i) = 0.0d0
			end if	
		end do
	end subroutine initdata

	subroutine targetf(par, f)
		double precision xb(n), q(n), prob(n), par(k), f
		integer i
		xb = matmul(x,par)
		q = 2.0d0*y - 1.0d0
		do i = 1, n, 1
			prob(i) = probn(q(i)*xb(i))
		end do
		f = sum(dlog(prob))
	end subroutine targetf

	subroutine estimate
		integer nchain, nburn, nkeep
		double precision, allocatable :: initpar(:), chain(:,:)
		call initdata
		allocate(initpar(k))
		initpar(:) = 0.0d0
		nchain = 50000
		nburn = 250
		nkeep = 10000
		allocate(chain(k,nkeep)) 
		call initchain(nburn, nchain, nkeep, k, 0.0d0,0.001d0)
		call dochain(targetf, initpar, chain)
		call dostats(targetf, chain, .true., .true., 2, (/0.025d0,0.975d0/))
	end subroutine estimate

		double precision function probn(value)
			double precision value,  q, bound
			integer ifail, status
			call cdfnor(1,probn,q, value, 0.0d0, 1.0d0, status, bound)
		end function probn

		! wrapper for inverse normal cdf from dcdflib.a
		double precision function quann(prob)
			double precision prob, bound
			integer ifail, status
			ifail = 1
			call cdfnor(2,prob,1.0d0-prob, quann, 0.0d0, 1.0d0, status, bound)
		end function quann



end module chain


program main
use chain

	call estimate

	write(*,*) 'true parameters :'
	do i = 1, k, 1
		write(*,*) i, par(i)
	end do

end program main
