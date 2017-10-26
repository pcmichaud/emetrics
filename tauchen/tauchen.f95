! Simple example of how to discretize a continuous AR(1) process
! Pierre-Carl Michaud, 2017
! e(t) = rho*e(t-1) + eps(t), eps(t)~N(0,sigmaeps)
! compile first into a library the probability functions:
! cd ~/dcdflib.f/src
! rm *.o
! gfortran -c *.f
! ar r ../../dcdflib.a *.o
! Then compile this program
! gfortran tauchen.f95 dcdflib.a -o test
! execute using 
! ./test

program main
implicit none			
integer ne, i
double precision, allocatable :: pte(:), eprob(:,:), pstate(:)
double precision rho, sigmaeps, m, mean, var, eqvar

! transition probabilities shocks
ne = 10
allocate(pte(ne))
allocate(eprob(ne,ne))

! m is the for the bounds of the points (as multiple of stationary sd)
m = 2.0d0
! parameters of process above
rho = 0.97d0
sigmaeps = 0.1d0
eqvar = sigmaeps/(1.0d0-rho**2)
call tauchen(ne,0.0d0,rho,sigmaeps, m, pte,eprob)

write(*,*) 'results in the following discrete points : '
do i = 1, ne, 1
	write(*,*) i, pte(i)
end do
write(*,*) ''
write(*,*) 'with transition matrix (origin rows to destination columns) :'
do i = 1, ne, 1
	write(*,*) eprob(i,:)
end do
write(*,*) 'stationary distribution :'
allocate(pstate(ne))
pstate(:) = 1.0d0/dble(ne)
do i = 1, 100000, 1
	pstate = matmul(transpose(eprob),pstate)
end do
do i = 1, ne, 1
	write(*,*) i,pte(i), pstate(i)
end do
mean = 0.0d0
do i = 1, ne, 1
	mean = mean + pstate(i)*pte(i)
end do
write(*,*) 'mean (simulated, true) ', mean, 0.0d0
var = 0.0d0
do i = 1, ne, 1
	var = var + pstate(i)*(pte(i)**2)
end do
var = var - mean**2
write(*,*) 'variance (simulated, true) ', var, eqvar

end program main

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

! tauchen discretization of AR(1) process (adaption from other programs out there)
subroutine tauchen(n,mean,rho,sigmaeps, m, z,zprob)
	integer n,k, j, i
	double precision mean,rho,sigmaeps,m,z(n),zprob(n,n),zstep, probn
	if (n==1) then
		z= mean
		zprob= 1.0d0
	else
		z(n)=m*sqrt(sigmaeps/(1.0d0-rho**2))
		z(1)=-z(n)
		zstep=(z(n)-z(1))/dble(n-1)

		do i= 2, n-1, 1 
			z(i)=z(1)+zstep*dble(i-1)
		end do

		z= z + mean/(1.0d0-rho)

		do j=1, n, 1
			do k=1, n, 1
				if (k==1) then
					zprob(j,k)=probn((z(1)-mean-rho*z(j)+0.5d0*zstep)/dsqrt(sigmaeps))
				end if
				if (k==n) then
					zprob(j,k)=1.0d0-probn((z(n)-mean-rho*z(j) &
						-0.5d0*zstep)/dsqrt(sigmaeps))
				end if
				if (k .gt. 1 .and. k .lt. n) then
					zprob(j,k)=probn((z(k)-mean-rho*z(j)+0.5d0*zstep)/dsqrt(sigmaeps)) &
						- probn((z(k)-mean-rho*z(j)-0.5d0*zstep)/dsqrt(sigmaeps))
				end if
	  		end do
		end do
	end if
end subroutine tauchen