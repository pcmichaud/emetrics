program main

	implicit none

	integer i
	! age
	integer nages, agemin, agemax, a
	integer, allocatable :: ages(:)
	double precision, allocatable :: probm(:)
	double precision mortpar(2)
	! earnings shocks
	integer ne, retage, e
	double precision rho, sigmaeps, reprate, basewage, wagepar(2), eqsigmaeps
	double precision, allocatable :: gride(:), probe(:,:), cumprobe(:,:), wages(:)
	! wealth
	integer nw, w
	double precision wmin, wmax, gapw, rate
	double precision, allocatable :: gridw(:)

	! preferences
	double precision sigma, beta, phi, kappa

	! solutions
	double precision, allocatable :: value(:,:,:), cons(:,:,:), cash(:,:,:)
	double precision ax, bx, cx, x0, x1, x2, x3, f1, f2, f0, f3, xmin, fmin
	integer ee, w0, w1
	double precision nextw, wu 
	logical retired
	double precision r, c, tol, cump
	parameter (tol = 1.0d-6, r= 0.61803399d0, c= 1.0d0 - r)
	
	! simulating from these rules
	double precision draw,  wagesim, wsim, csim, cashsim
	integer agesim, esim, nsim

	! ages
	agemin = 20
	agemax = 110
	nages = agemax - agemin + 1
	allocate(ages(nages))
	do i = 1, nages, 1
		ages(i) = agemin + i - 1
	end do
	! earnings shocks
	ne = 10
	allocate(gride(ne), probe(ne,ne), cumprobe(ne,ne))
	! relatively typical parameters
	rho = 0.97d0
	sigmaeps = 0.03d0
	eqsigmaeps = sigmaeps/(1.0d0 - rho**2)
	call tauchen(ne,0.0d0,rho,sigmaeps,2.0d0,gride,probe)
	do i = 1, ne, 1
		gride(i) = dexp(gride(i) - 0.5d0*eqsigmaeps)
		cump = 0.0d0
		do ee = 1, ne, 1
			cump = cump + probe(i,ee)
			cumprobe(i,ee) = cump
		end do
	end do
	allocate(wages(nages))
	basewage = 40.0d3
	retage = 65
	reprate = 0.7d0
	wagepar(1) = 0.05d0
	wagepar(2) = -0.001d0
	write(*,*) 'wages'
	do i = 1, nages, 1
		if (ages(i) .lt. retage) then
			wages(i) = basewage*dexp(wagepar(1)*dble(i) + wagepar(2)*dble(i)**2)
		else if (ages(i) .eq. retage) then
			wages(i) = reprate*wages(i-1)
		else
			wages(i) = wages(i-1)	
		end if

		write(*,*) ages(i), wages(i)
	end do
	! wealth
	nw = 64
	wmin = 0.0d0
	wmax = 10*maxval(wages)
	gapw = (dsqrt(wmax) - dsqrt(wmin))/dble(nw-1)
	allocate(gridw(nw))
	do i = 1, nw, 1
		gridw(i) = dsqrt(wmin) + dble(i-1)*gapw
	end do

	! mortality risk (see mortality.do to calibrate gompertz process)
	mortpar(1) = dexp(-10.21d0)
	mortpar(2) =  0.0886d0
	write(*,*) 'mortality probabilities'
	allocate(probm(nages))
	do i = 1, nages, 1
		probm(i) = mortpar(1)*dexp(mortpar(2)*dble(ages(i)))
		write(*,*) ages(i), probm(i)
	end do

	! preference parameter
	sigma = 1.6d0
	beta = 0.98d0
	phi = 6.0d0
	kappa = 100.0d3
	rate = 0.03d0

	! allocate memory for solution
	allocate(value(nages,ne,nw), cons(nages,ne,nw), cash(nages,ne,nw))

	! work solution
	! last year
	do e = 1, ne, 1
		do w = 1, nw, 1
			! cash on hand
			cash(nages,e,w) = gridw(w)**2 + wages(nages)*gride(e)
			cons(nages,e,w) = 0.0d0
			! last year dies for sure, so compute bequest
			value(nages,e,w) = phi * ((cash(nages,e,w) + kappa)**(1.0d0 - sigma))/(1.0d0-sigma)
		end do
	end do 

	! other years
	do a = nages-1, 1, -1
		if (ages(a).ge.retage) then
			retired = .true.
		else
			retired = .false.
		end if	
		do e = 1, ne, 1
			do w = 1, nw, 1

				! compute cash on hand
				cash(a,e,w) = gridw(w)**2 + wages(a)*gride(e)
				! solve for optimal consumption using golden section search
				! use golden
				ax = 1.0d0
				bx = 0.25d0*cash(a,e,w)
				cx = cash(a,e,w)

				x0 = ax
				x3 = cx
				if (abs(cx-bx).gt.abs(bx-ax)) then
					x1 = bx
					x2 = bx + c*(cx-bx)
				else
					x2 = bx
					x1 = bx - c*(bx-ax)
				end if
				call getvalue(retired, cash(a,e,w), x1, gridw, nw, value(a+1,:,:), probe(e,:), e, ne, &
									probm(a), sigma, beta, phi, kappa, rate, f1)
				f1 = -f1


				call getvalue(retired, cash(a,e,w), x2, gridw, nw, value(a+1,:,:), probe(e,:), e, ne, &
									probm(a), sigma, beta, phi, kappa, rate, f2)
				f2 = -f2

				do while (abs(x3-x0).gt.(tol*(abs(x1)+abs(x2))))
					if (f2.lt.f1) then
						x0 = x1
						x1 = x2
						x2 = r*x1 + c*x3
						f0 = f1
						f1 = f2
						call getvalue(retired, cash(a,e,w), x2, gridw, nw, value(a+1,:,:), probe(e,:), e, ne, &
											probm(a), sigma, beta, phi, kappa, rate, f2)
						f2 = -f2
					else
						x3 = x2
						x2 = x1
						x1 = r*x2 + c*x0
						f3 = f2
						f2 = f1
						call getvalue(retired, cash(a,e,w), x1, gridw, nw, value(a+1,:,:), probe(e,:), e, ne, &
											probm(a), sigma, beta, phi, kappa, rate, f1)
						f1 = -f1
					end if
				end do
				if (f1.lt.f2) then
					xmin = x1
					fmin = f1
				else
					xmin = x2
					fmin = f2
				end if	

				value(a,e,w) = -fmin
				cons(a,e,w) = xmin
			end do
		end do 

	end do


	! output rules for analyis at various ages
	open(1, file='rules.dat')
	do i = 1, nages, 1
		do e = 1, ne, 1
			do w = 1, nw, 1
				write(1,*) ages(i), e, wages(i)*gride(e), cash(i,e,w), cons(i,e,w), value(i,e,w)
			end do
		end do 
	end do
	close(1)



	nsim = 10000
	open(1, file='simulated.dat')
	do i = 1, nsim, 1

		! initial conditions
		agesim = ages(1)
		! shock
		esim = floor(0.5d0*dble(ne))
		! initial wealth
		wsim = 0.0d0
		! retirement status
		retired = .false.
		do a = 1, nages, 1
			if (ages(a) .ge. retage) then
				retired = .true.
			end if	
			! get position
			call scale(dsqrt(wsim),dsqrt(wmin),dsqrt(wmax),nw,gridw,w0,w1,wu) 
			! find optimal consumption
			csim = cons(a,esim,w0) + wu*(cons(a,esim,w1) - cons(a,esim,w0))

			! compute income and cash
			wagesim = wages(a)*gride(esim)
			cashsim = wsim + wagesim

			write(1,*) i, ages(a), esim, wagesim, wsim, cashsim, csim
			call random_number(draw)
			if (draw .lt. probm(a)) then
				exit
			else	
				! update wealth
				wsim = (1.0d0 + rate)*(cashsim - csim)
				if (wsim .lt. wmin) then
					wsim = wmin
				end if
				if (wsim .gt. wmax) then
					wsim = wmax
				end if	

				! update shock if necessary	
				if (.not. retired) then
					call random_number(draw)
					do e = 1, ne, 1
						if (draw .lt. cumprobe(esim,e)) then
							esim = e
							exit
						end if	
					end do
				end if
			end if

		end do

	end do
	close(1)




end program main


subroutine getvalue(retired, cash, cons, gridw, nw, value, probe, e, ne, &
									probm, sigma, beta, phi, kappa, rate, func)
		integer ne, nw, w0, w1, ee, e
		double precision cash, cons, gridw(nw), value(ne,nw), probe(ne)
		double precision probm, sigma, beta, phi, kappa, rate, func		
		double precision nextw, pv, wmin, wmax, wu
		logical retired
		nextw = (1.0d0 + rate)*(cash - cons)
		wmin = minval(gridw(:))
		wmax = maxval(gridw(:))
		call scale(dsqrt(nextw),dsqrt(wmin),dsqrt(wmax),nw,gridw,w0,w1,wu)
		func = cons**(1.0d0-sigma)/(1.0d0-sigma)
		if (retired) then
			ee = e
			pv = value(ee,w0) + wu*(value(ee,w1) - value(ee,w0))
			func = func + beta*(1.0d0 - probm)*pv
		else
			do ee = 1, ne, 1
				pv = value(ee,w0) + wu*(value(ee,w1) - value(ee,w0))
				func = func + beta*probe(ee)*(1.0d0 - probm)*pv
			end do
		end if
		func = func + beta*probm*(phi*((kappa + nextw)**(1.0d0 - sigma))/(1.0d0-sigma))
end subroutine



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

subroutine scale(z,zmin,zmax,nz,gridz,z0,z1,zu)
	integer nz,z0,z1
	double precision z, zmin, zmax, gridz(nz),gapz, zu
	gapz = (zmax - zmin)/dble(nz-1)
	z0 = floor((z - zmin)/gapz) + 1
	if (z0.lt.1) then
		z = zmin
		z0 = 1
		z1 = 2
	else if (z0.ge.nz) then
		z1 = nz
		z0 = nz-1
	else 
		z1 = z0 + 1
	end if
	zu = (z - gridz(z0))/(gridz(z1)-gridz(z0))

end subroutine scale

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

