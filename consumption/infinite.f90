program main

	implicit none
	integer nw, ne, i,  e, w
	double precision wmin, wmax, gapw, wage
	double precision rho, sigmaeps, eqsigmaeps
	double precision, allocatable :: gridw(:), gride(:), probe(:,:), cumprobe(:,:)
	double precision beta, rate, sigma
	double precision r, c, tol
	parameter (tol = 1.0d-5, r= 0.61803399d0, c= 1.0d0 - r)
	double precision ax,bx,cx,x0,x1,x2,x3
	double precision f0,f1,f2,f3,xmin,fmin
	double precision nextw, wu, pv, criterion
	integer w0, w1, ee
	double precision, allocatable :: value(:,:), tvalue(:,:), cons(:,:), cash(:,:)
	
	nw = 100
	wage = 20.0d3
	wmin = 0.0d0
	wmax = 10*wage
	gapw = (dsqrt(wmax) - dsqrt(wmin))/dble(nw-1)
	allocate(gridw(nw))
	do i = 1, nw, 1
		gridw(i) = dsqrt(wmin) + dble(i-1)*gapw
	end do
	
	ne = 10
	allocate(gride(ne), probe(ne,ne), cumprobe(ne,ne))
	rho = 0.97d0
	sigmaeps = 0.1d0
	eqsigmaeps = sigmaeps/(1.0d0 - rho**2)
	call tauchen(ne,0.0d0,rho,sigmaeps,2.0d0,gride,probe)
	do i = 1, ne, 1
		gride(i) = dexp(gride(i) - 0.5d0*eqsigmaeps)
	end do

	! parameters
	beta = 0.94d0
	rate = 0.03d0
	sigma = 5.0d0


	allocate(value(ne,nw), cons(ne,nw), tvalue(ne,nw), cash(ne,nw))

	value(:,:) = 0.0d0
	criterion = 1.0d5
	do while (criterion .gt. 1.0d-12) 
		do e = 1, ne, 1
			do w = 1, nw, 1

				! cash on hand
				cash(e,w) = gridw(w)**2 + wage*gride(e)

				! use golden
				ax = 1.0d0
				bx = 0.25d0*cash(e,w)
				cx = cash(e,w)

				x0 = ax
				x3 = cx
				if (abs(cx-bx).gt.abs(bx-ax)) then
					x1 = bx
					x2 = bx + c*(cx-bx)
				else
					x2 = bx
					x1 = bx - c*(bx-ax)
				end if


				f1 = x1**(1.0d0-sigma)/(1.0d0-sigma)
				nextw = (1.0d0 + rate)*(cash(e,w) - x1)
				call scale(dsqrt(nextw),dsqrt(wmin),dsqrt(wmax),nw,gridw,w0,w1,wu)
				do ee = 1, ne, 1
					pv = value(ee,w0) + wu*(value(ee,w1) - value(ee,w0))
					f1 = f1 + beta*probe(e,ee)*pv
				end do
				f1 = -f1
				f2 = x2**(1.0d0-sigma)/(1.0d0-sigma)
				nextw = (1.0d0 + rate)*(cash(e,w) - x2)
				call scale(dsqrt(nextw),dsqrt(wmin),dsqrt(wmax),nw,gridw,w0,w1,wu)
				do ee = 1, ne, 1
					pv = value(ee,w0) + wu*(value(ee,w1) - value(ee,w0))
					f2 = f2 + beta*probe(e,ee)*pv
				end do
				f2 = -f2

				do while (abs(x3-x0).gt.(tol*(abs(x1)+abs(x2))))
					if (f2.lt.f1) then
						x0 = x1
						x1 = x2
						x2 = r*x1 + c*x3
						f0 = f1
						f1 = f2
						f2 = x2**(1.0d0-sigma)/(1.0d0-sigma)
						nextw = (1.0d0 + rate)*(cash(e,w) - x2)
						call scale(dsqrt(nextw),dsqrt(wmin),dsqrt(wmax),nw,gridw,w0,w1,wu)
						do ee = 1, ne, 1
							pv = value(ee,w0) + wu*(value(ee,w1) - value(ee,w0))
							f2 = f2 + beta*probe(e,ee)*pv
						end do
						f2 = -f2
					else
						x3 = x2
						x2 = x1
						x1 = r*x2 + c*x0
						f3 = f2
						f2 = f1
						f1 = x1**(1.0d0-sigma)/(1.0d0-sigma)
						nextw = (1.0d0 + rate)*(cash(e,w) - x1)
						call scale(dsqrt(nextw),dsqrt(wmin),dsqrt(wmax),nw,gridw,w0,w1,wu)
						do ee = 1, ne, 1
							pv = value(ee,w0) + wu*(value(ee,w1) - value(ee,w0))
							f1 = f1 + beta*probe(e,ee)*pv
						end do
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

				tvalue(e,w) = -fmin
				cons(e,w) = xmin

			end do
		end do	

		criterion = maxval(dabs(tvalue - value))
		value = tvalue
	end do

	! consumption function
	write(*,*) 'solution is (wealth, lowest shock, mean, highest shock)'
	open(1,file='rules.dat')
	do i = 1, nw, 1
		write(*,*) cons(1,i)/cash(1,i), cons(5,i)/cash(5,i), cons(ne,i)/cash(ne,i)
		write(1,*) cash(1,i), cons(1,i), cash(5,i), cons(5,i), cash(ne,i), cons(ne,i)
	end do
	close(1)


end program main


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

