
subroutine golden(ax,bx,cx,xmin,fmin)
	! assumes existence of a function, f() returns dble
	double precision r, c, tol
	parameter (tol = 1.0d-3, r= 0.61803399d0, c= 1.0d0 - r)
	double precision ax,bx,cx,x0,x1,x2,x3
	double precision f0,f1,f2,f3,xmin,fmin

	x0 = ax
	x3 = cx
	if (abs(cx-bx).gt.abs(bx-ax)) then
		x1 = bx
		x2 = bx + c*(cx-bx)
	else
		x2 = bx
		x1 = bx - c*(bx-ax)
	end if

	call distance(x1,f1) 
	call distance(x2,f2) 

	do while (abs(x3-x0).gt.(tol*(abs(x1)+abs(x2))))
		if (f2.lt.f1) then
			x0 = x1
			x1 = x2
			x2 = r*x1 + c*x3
			f0 = f1
			f1 = f2
			call distance(x2, nagg, taxrate, f2) 
		else
			x3 = x2
			x2 = x1
			x1 = r*x2 + c*x0
			f3 = f2
			f2 = f1
			call distance(x1, nagg, taxrate, f1) 
		end if
	end do
	if (f1.lt.f2) then
		xmin = x1
		fmin = f1
	else
		xmin = x2
		fmin = f2
	end if	

end subroutine golden