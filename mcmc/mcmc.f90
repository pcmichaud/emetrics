module mcmc
	implicit none
	
	private
		integer nburn, nchain, nkeep, npar
		double precision sigma, eps, kappa
		double precision, allocatable:: initcov(:,:), eye(:,:)
		
	public:: dochain, initchain, dostats
	
	contains

	   subroutine initchain(burn, chain, keep, par, ep, kap)
			integer burn, chain, keep, par, i
			double precision sig, ep, kap
			! assign to globals
			nburn = burn
			nchain = chain
			nkeep = keep
			npar = par
			sigma = ((2.4d0)**2)/dble(npar)
			eps = ep
			kappa = kap
			! initial covariance matrix
			allocate(initcov(npar,npar), eye(npar,npar))
			eye(:,:) = 0.0d0
			do i = 1, npar, 1
				eye(i,i) = 1.0d0
			end do
			initcov = kappa * eye
			write(*,*) '** Initialized parameters'
			write(*,*) '- # chain: ', nchain
			write(*,*) '- # burn-in: ', nburn
			write(*,*) '- # keep: ', nkeep
			write(*,*) '- # parameters: ', npar
			write(*,*) '- s_d (Gelman): ', sigma
			write(*,*) '- epsilon: ', eps
			write(*,*) '- kappa: ', kappa
			write(*,*) ''
	   end subroutine initchain

	   subroutine dochain(targetf, initpar, chain)
	   		double precision initpar(npar), chain(npar,nkeep)
	   		double precision initmean(npar), mean(npar), newmean(npar)
	   		double precision cov(npar,npar), update(npar,npar), newcov(npar,npar)
	   		double precision accept, newpi, pi, prop(npar), state(npar)
	   		double precision unif, x(npar,1), y(npar,1), z(npar,1)
	   		integer i, j, naccept, c
	   		external targetf
	   		write(*,*) '** starting MCMC chain'

	   		! burn in to obtain initial mean and covariance matrix
	   		call doburn(targetf, initmean, initpar)

	   		! initialize chain
	        cov = initcov
	        mean = initmean
	        state = initpar
	        call targetf(state,pi)
	        j = 1
	        c = nburn
	        naccept = 0
	        write(*,*) '++ running main chain...'
	   		do i = 1, nchain, 1
	   			c = c + 1
	   			call proposal(cov, state, prop)
				! compute new target
				call targetf(prop,newpi)
				! compute acceptance probability
				accept = acceptf(newpi,pi)
				! update chain if proposal accepted
				call random_number(unif)
				if (unif .lt. accept) then
					state = prop
					pi = newpi
					naccept = naccept + 1
				end if
				newmean = 1.0d0/dble(c+1) * (state + dble(c)*mean)
				

				x = reshape(mean,(/npar,1/))
				y = reshape(newmean,(/npar,1/))
				z = reshape(state,(/npar,1/))
				newcov = (dble(c-1)/dble(c)) * cov 
				newcov = newcov + sigma*matmul(x,transpose(x))
				newcov = newcov - (sigma/dble(c)*dble(c+1))*matmul(y,transpose(y))
				newcov = newcov + (sigma/dble(c))*matmul(z,transpose(z))
				newcov = newcov + (sigma/dble(c)*eps)*eye

				mean = newmean
				cov = newcov

				!write(*,*) i, pi, state
				if (i .ge. nchain - nkeep + 1) then
					chain(:,j) = state
					j = j + 1
				end if
	   		end do

	   		write(*,*) 'Completed...'
	   		write(*,*) 'acceptance rate', dble(naccept)/dble(nchain)
	   		write(*,*) 'Std of proposal distribution'
	   		do i = 1, npar, 1
	   			write(*,*) i, sqrt(cov(i,i))
	   		end do
	   		

	   end subroutine dochain



	   subroutine dostats(targetf, chain, imean, istd, nq, quant)
	   		double precision chain(npar,nkeep)
	   		double precision, allocatable :: stats(:,:)
	   		logical imean, istd, iquant
	   		integer nstats, i, c, j, nq
	   		double precision quant(nq), x(nkeep), pi
	   		integer idquant(nq)
	   		external targetf
	   		nstats = 0
	   		if (imean) then
	   			nstats = nstats + 1
	   		end if	
	   		if (istd) then
	   			nstats = nstats + 1
	   		end if
	   		if (nq>0) then
	   			iquant = .true.
	   			nstats = nstats + nq
	   			do i = 1, nq, 1
					idquant(i) = ceiling(quant(i)*dble(nkeep))
				end do
	   		end if
	   		allocate(stats(npar,nstats))

	   		do i = 1, npar, 1
	   			j = 1
	   			if (imean) then
	   				stats(i,j) = sum(chain(i,:))/dble(nkeep)
	   				j = j + 1
	   			end if
	   			if (istd) then
	   				stats(i,j) = 0.0d0
	   				do c = 1, nkeep, 1
	   					stats(i,j) = stats(i,j) + chain(i,c)**2  / dble(nkeep)
	   				end do
	   				stats(i,j) = dsqrt(stats(i,j) - stats(i,1)**2)
	   				j = j + 1
	   			end if
	   			if (iquant) then
	   				x = chain(i,:)
	   				call sort(nkeep,x)
	   				do c = 1, nq, 1
	   					stats(i,j) = x(idquant(c))
	   					j = j + 1
	   				end do
	   			end if
	   		end do
			
	   		write(*,*) '** Statistics from running chain '
	   		do i = 1, npar, 1
	   			write(*,*) i, stats(i,:)
	   		end do
	   		call targetf(stats(:,1),pi)
	   		write(*,*) '- Function value at posterior mean = ', pi


	   end subroutine dostats

		subroutine sort(n,x)
		  double precision x(n), a
		  integer j, i, n
		  do j=2, n
		    a=x(j)
		    do i=j-1,1,-1
		      if (x(i)<=a) goto 10
		      x(i+1)=x(i)
		    end do
			i=0
		10  x(i+1)=a
		  end do
		end subroutine

	   subroutine doburn(targetf, initmean, initpar)
	   		double precision initmean(npar), burn(npar,nburn)
	   		double precision state(npar), initpar(npar)
	   		double precision prop(npar), newpi, pi, accept, unif, x(npar,1)
	   		external targetf
	   		integer select, i
	   		write(*,*) ''
	   		write(*,*) '++ burning-in chain...'
	   		initmean(:) = 0.0d0
	   		! function at initial parameter
	   		call targetf(initpar, pi)
	   		! draw id of element of burn chain used as start point
	   		call random_number(unif)
	   		select = ceiling(unif*dble(nburn))
	   		do i = 1, nburn, 1
	   			call proposal(initcov, state, prop)
				! compute new target
				call targetf(prop,newpi)
				! compute acceptance probability
				accept = acceptf(newpi,pi)
				! update chain if proposal accepted
				call random_number(unif)
				if (unif .lt. accept) then
					state = prop
					pi = newpi
				end if
				initmean = initmean + state / dble(nburn)
				burn(:,i) = state
				if (i .eq. select) then
					initpar = state
				end if	
				!write(*,*) i, pi, state
	   		end do

	   		! update initial covariance using what is obtained here
	   		initcov(:,:) = 0.0d0
	   		do i = 1, nburn, 1
	   			x = reshape((burn(:,i) - initmean),(/npar,1/))
	   			initcov = initcov + matmul(x,transpose(x))/dble(nburn-1)
	   		end do
	   		write(*,*) '- finished burn in (draw, mean, std) : '
	   		do i = 1, npar, 1
	   			write(*,*) i, initpar(i), initmean(i), sqrt(initcov(i,i))
	   		end do
	   		write(*,*) ' '

	   end subroutine doburn

		subroutine proposal(cov, state, prop)
			integer i
			double precision state(npar), prop(npar), draw(npar), cov(npar,npar)
			double precision  chol(npar,npar), unif, shock(npar)
			! take random gaussian draws
			do i = 1, npar, 1
				call random_number(unif)
				draw(i) = quann(unif)
			end do
			! apply covariance
			chol = cov
			call cholesky(chol,npar)
			shock = matmul(chol,draw)
			prop = state + shock
		end subroutine proposal

		double precision function acceptf(targ,oldtarg)
			double precision targ, oldtarg
			acceptf = dexp(targ - oldtarg) 
			if (isnan(acceptf)) then
				acceptf = 0.0d0
			else	
				if (acceptf.gt. 1.0d0) then
					acceptf = 1.0d0
				end if
			end if
		end function acceptf

		! wrapper for normal CDF from dcdflib.a
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

		subroutine cholesky(A,n)
		  ! formal vars
		  integer :: n      ! number of rows/cols in matrix
		  double precision    :: A(n,n) ! matrix to be decomposed
		  ! local vars
		  integer :: j      ! iteration counter
		  ! begin loop
		  do j = 1,n
		    ! perform diagonal component
		    A(j,j) = sqrt(A(j,j) - dot_product(A(j,1:j-1),A(j,1:j-1)))
		    ! perform off-diagonal component
		    if (j < n) A(j+1:n,j) = (A(j+1:n,j) - matmul(A(j+1:n,1:j-1),A(j,1:j-1))) / &
		   	&           A(j,j)

		  end do
		end subroutine cholesky




end module mcmc