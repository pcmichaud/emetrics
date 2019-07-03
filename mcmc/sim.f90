program main
	implicit none
	integer n, k, c, i, j, keep, start_keep, start_release, q5, q95
	double precision, allocatable :: chain(:,:), par(:), x(:,:), y(:), xb(:), prop(:)
	double precision q, draw, quann, newpi, accept, targetf, acceptf, pi, eps, tune
	double precision, allocatable :: cov(:,:), cov0(:,:), mean(:,:), eye(:,:), oldmean(:,:), state(:,:), postmean(:), postci(:,:)
	logical free

	! simulate data for a probit model
	n = 10000
	k = 4
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
	write(*,*) 'P(y=1) = ' , sum(y)/dble(n)
	
	! start MCMC
	c = 10000						! length of chain
	keep = 1000						! draws to keep at end to compute stats
	start_keep = c - keep + 1		! point in chain when starts keeping
	start_release = 100				! burn in period for covariance matrix estimation
	allocate(chain(k,keep))			! array to store chain that is kept
	allocate(prop(k), state(k,1))	! allocations
	state(:,1) = 0.0d0 				! state of chain to start is set to zero
	pi = targetf(y,x,n,k,state(:,1))	! value of likelihood at current state
	allocate(cov(k,k),cov0(k,k),mean(k,1), oldmean(k,1), eye(k,k))	
	oldmean(:,1) = par
	cov0 = 0.0d0 					! initial covariance matrix for proposal distribution
	do i = 1, k, 1
		cov0(i,i) = 0.1d0
	end do
	cov = cov0
	eps = 1.0d-6					! so that covariance does not degenerate
	tune = (2.4d0**2) / dble(k) 	! tuning parameter (Gelman)
	eye = 0.0d0 					! an identity matrix
	do i = 1, k, 1
		eye(k,k) = 1.0d0
	end do
	free = .false.  				! will start with initial covariance matrix and release later
	do i = 2, c, 1 					! start chain
		! make a proposal 
		if (i.ge.start_release) then
			free = .true.
			call proposal(cov,state(:,1), prop, k, free)
		else
			call proposal(cov0,state(:,1), prop, k, free)	
		end if	
		! compute new target
		newpi = targetf(y,x,n,k,prop)
		! compute acceptance probability
		accept = acceptf(newpi,pi)
		! update chain if proposal accepted
		call random_number(draw)
		if (draw .lt. accept) then
			state(:,1) = prop
			pi = newpi
		end if		
		! update mean and covariance of chain
		mean = (dble(i-1)*oldmean + state)/dble(i)
		cov = (dble(i-1)/dble(i))*cov + tune/dble(i) * (dble(i)*matmul(oldmean,transpose(oldmean)) &
				- dble(i+1)*matmul(mean,transpose(mean)) + matmul(state,transpose(state)) + eps*eye)
		oldmean = mean
		if (i .ge. start_keep) then
			chain(:, i - start_keep + 1) = state(:,1)
		end if	
	end do 
	allocate(postmean(k), postci(k,2))
	q5 =  floor(0.05*keep)
	q95 = floor(0.95*keep)
	do i = 1, k, 1
		postmean(i) = sum(chain(i,:))/dble(keep)
		call sort(keep,chain(i,:))
		postci(i,1) = chain(i,q5)
		postci(i,2) = chain(i,q95)		
	end do
	write(*,*) 'True, posterior mean, 5th percentile, 95th percentile '
	do i = 1, k, 1
		write(*,*) i, par(i), postmean(i), postci(i,1), postci(i,2)
	end do

end program main

double precision function targetf(y, x, n, k, par)
	integer n, k
	double precision y(n), x(n,k), xb(n), q(n), prob(n), par(k)
	double precision probn
	xb = matmul(x,par)
	q = 2.0d0*y - 1.0d0
	do i = 1, n, 1
		prob(i) = probn(q(i)*xb(i))
	end do
	targetf = sum(dlog(prob))
end function 


subroutine proposal(cov, par, newpar, k, free)
	integer k, i
	double precision par(k), newpar(k), draw(k), quann, cov(k,k), d, chol(k,k)
	logical free
	call random_number(d)
	! take random gaussian draws
	do i = 1, k, 1
		call random_number(d)
		draw(i) = quann(d)
	end do
	! apply covariance
	chol = cov
	call cholesky(chol,k)
	draw = matmul(chol,draw)
	newpar = par + draw
end subroutine proposal

double precision function acceptf(targ,oldtarg)
	double precision targ, oldtarg
	acceptf = dexp(targ - oldtarg) 
	if (isnan(accepf)) then
		acceptf = 0.0d0
	else	
		if (acceptf.gt. 1.0d0) then
			acceptf = 1.0d0
		end if
	end if
end function 

! wrapper for normal CDF from dcdflib.a
double precision function probn(value)
	double precision value,  q, bound
	integer ifail, status
	!double precision g01eaf
	!external g01eaf
	!ifail = 1
	call cdfnor(1,probn,q, value, 0.0d0, 1.0d0, status, bound)
	! function that gives the cumulative probability from normal
	!probn = g01eaf('L',value,ifail)
	!write(*,*) prob
end function

! wrapper for inverse normal cdf from dcdflib.a
double precision function quann(prob)
	double precision prob, bound
	integer ifail, status
	!double precision g01faf
	!external g01faf
	ifail = 1
	call cdfnor(2,prob,1.0d0-prob, quann, 0.0d0, 1.0d0, status, bound)
	! function that gives the cumulative probability from normal
	!quann = g01faf('L',prob,ifail)
	!write(*,*) prob
end function

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

! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/sort1_f90.txt
subroutine sort(n,x)
  double precision x(n), a
  integer j, i
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

