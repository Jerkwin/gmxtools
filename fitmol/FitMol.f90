Program FitMol
	integer, parameter:: n = 2000, Lxyz = 1, Lout = 2

	integer i, Ierr, Narg, Lsta, Natm, Icnt, Nfit, Nmol, Idone
	real*8  RMSD, RMSDmin, XYZref(3,n), XYZmol(3,n), XYZrot(3,n), Wref(n), Xcnt(3), Ycnt(3), q(4), u(3, 3)
	logical	YesBat
	character*8 Sref(n), Smol(n)
	character*256	txt, Tips
	common /Lsta/ Lsta

	integer, allocatable:: Ifit(:), idx(:), pidx(:)
	real*8,  allocatable:: w(:), x(:,:), y(:,:), z(:,:)

	Icnt = 0
	YesBat = .false.

	Narg = Iargc()
	if(Narg>0) then
		YesBat = .true.
		call getarg(1, txt)
		i = index(txt, '.xyz', .true.)-1; if(i>0) txt = txt(1:i)
		open(unit=Lxyz, file=trim(txt)//'.xyz', status='old', IOstat=Ierr)
		open(unit=Lout, file=trim(txt)//'~Fit.xyz', IOstat=Ierr)
		if(Narg>1) then
			call getarg(2, txt)
			read(txt, *) Icnt
		end if
		goto 100
	end if

	print*, '>>>>>>>>>>>>>>>>Program FitMol Running <<<<<<<<<<<<<<<<'
	print*, '>>>>>>>>>>>>>>>>        Jicun LI       <<<<<<<<<<<<<<<<'
	print*, '>> 2019-09-23: No. of atoms can be different for each mol.'
	print*, '>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	print*

	print*, '>>Type XYZ File Name (*.xyz/exit):'
	read*, txt
	if(index(txt, 'exit')/=0) stop '>>Normal Exit of Program FitMol.'

	i = index(txt, '.xyz', .true.)-1; if(i>0) txt = txt(1:i)
	open(unit=Lxyz, file=trim(adjustl(txt))//'.xyz', status='old', IOstat=Ierr)
	do while(Ierr/=0)
		print*, '>>XYZ Input File  .\'//trim(adjustl(txt))//'.xyz  NOT Exist !'
		print*, '>>Revise Input File Name (*.xyz/exit):'
		read*, txt
		if (index(txt, 'exit')/=0) stop '>>Normal Exit of Program FitMol.'
		i = index(txt, '.xyz')-1; if(i>0) txt = txt(1:i)
		open(unit=Lxyz, file=trim(adjustl(txt))//'.xyz', status='old', IOstat=Ierr)
	end do
	open(unit=Lout, file=trim(adjustl(txt))//'~Fit.xyz', IOstat=Ierr)

100	print*, '>>XYZ Input File:  .\', trim(adjustl(txt))//'.xyz'

	read(Lxyz, *, IOstat=Ierr) Natm
	read(Lxyz, '(A256)', IOstat=Ierr) txt
	if(Icnt==0) then
		read(txt, *, IOstat=Ierr) Icnt
		if(Ierr/=0) Icnt=0
	end if

	print*, '>>Input File Readin.'
	print*
	write(*, '(A, I4)') ' >>>>#Atoms of Ref-Mol      ', Natm
	if(Icnt==0) write(*, '(A, I4)') ' >>>>Rotation Origin is the weighted centriod.'
	if(Icnt/=0) write(*, '(A, I4)') ' >>>>Rotation Origin Fixed at Atom     ', Icnt

	if (.NOT.YesBat) then
		Lsta = 6
		print*
		print*, '>>FitMol Calculation Will Run.'
		print*
		print*, '>>Run Now ? Y/N'
		read*, txt
		if(txt=='N' .or. txt=='n') stop '>>Normal Exit of Program FitMol.'
	end if

	print*
	print*, '>>FitMol Running...'
	print*

	rewind(Lxyz)
	read(Lxyz, *) Natm
	read(Lxyz, '(A256)') Tips
	read(Lxyz, '(A256)') txt; backspace(Lxyz)
	read(txt, *, IOstat=Ierr) Sref(1), Wref(1), Wref(1), Wref(1), Wref(1)
	!print*, ierr

	Wref=1.d0
	do i=1, Natm
		if(Ierr/=0) read(Lxyz, *) Sref(i), XYZref(1,i), XYZref(2,i), XYZref(3,i)
		if(Ierr==0) read(Lxyz, *) Sref(i), XYZref(1,i), XYZref(2,i), XYZref(3,i), Wref(i)
	end do

	call center(Icnt, Natm, XYZref, Wref, Xcnt)
	!if(.NOT. YesBat) write(Lsta, '(A, 3F14.9)') 'Center of Ref-Mol: ', Xcnt

	write(Lout, *) Natm
	write(Lout, '(A)') trim(Tips)//' Ref-Mol'
	write(Lout, '(A4, 3F14.9, F8.1)') (Sref(i), XYZref(:, i), Wref(i), i=1, Natm)

	Nfit=0
	do i=1, Natm
		if(Wref(i)>0.d0) Nfit=Nfit+1
	end do

	allocate(Ifit(Nfit), w(Nfit), idx(0:Nfit), pidx(0:Nfit), &
	&		x(3, Nfit), y(3, Nfit), z(3, Nfit))

	Nfit=0
	do i=1, Natm
		if(Wref(i)>0.d0) then
			Nfit=Nfit+1
			Ifit(Nfit)=i
			w(Nfit)=Wref(i)
			x(1:3, Nfit)=XYZref(1:3,i)
		end if
	end do

	if(.NOT. YesBat) write(Lsta, '(A, I4, /)') '#Atoms to fit: ', Nfit

	Nmol=0; Ierr=0
	do
		read(Lxyz, *, IOstat=Ierr) Natm; if(Ierr/=0) exit
		read(Lxyz, '(A256)', IOstat=Ierr) Tips
		do i=1, Natm
			read(Lxyz, *) Smol(i), XYZmol(1, i), XYZmol(2, i), XYZmol(3, i)
		end do

		call center(Icnt, Natm, XYZmol, Wref, Ycnt)

		Nmol=Nmol+1
		RMSDmin=1.E9

		Idone=0; Imol=0
		do while(Idone<Nfit)
			call NextQuickPerm(Idone, Nfit, idx, pidx)

			imat=1
			do i=1, Nfit
				ii=Ifit(i); ip=Ifit(idx(i-1))
				if(Smol(ip) /= Sref(ii)) then; imat=0; exit; end if
			end do
			if(imat==0) cycle

			Imol=Imol+1

			do i=1, Nfit
				y(1:3, i)=XYZmol(1:3, Ifit(idx(i-1)))
			end do

			u = 0.d0; u(1,1) = 1.d0; u(2,2) = 1.d0; u(3,3) = 1.d0

			if (.NOT. YesBat)  then
				write(Lsta, *) '>> Mol:', Nmol, '    permutation:', Imol
				!call analyz(Nfit, x, y, w, u)
				write(Lsta, *) ' #Atom    Xmol        Ymol        Zmol      ' &
				&			// '#Atom    Xref        Yref        Zref      Weight'
				write(Lsta, '(2(I4,A3,3F12.6), 2F8.3)') &
				&	(Ifit(idx(i-1)), Smol(i), y(:,i), Ifit(i), Sref(i), x(:,i), &
				&	w(i), norm2(y(:,i)-x(:,i)), i=1, Nfit)
				write(Lsta, *)
			end if

			call qtrfit(Nfit, y, x, w, q, u, Ierr)
			call rotmol(Nfit, y, z, u)
			!if (.NOT.YesBat) call analyz(Nfit, z, x, w, u)

			RMSD = 0.d0
			do	i = 1, Nfit
				RMSD = RMSD + w(i)*norm2(x(:,i)-z(:,i))**2
			end do
			RMSD = sqrt(RMSD/dble(Nfit))
			RMSDmin=min(RMSDmin, RMSD)

			call rotmol(Natm, XYZmol, XYZrot, u)

			if (.NOT.YesBat) then
				write(Lsta, *) ' #Atom    Xfit        Yfit        Zfit      ' &
				&			// '#Atom    Xref        Yref        Zref    |Mol-Ref|  Weight'
				write(Lsta, '(2(I4,A3,3F12.6), 2F8.3)') (i, Smol(i), XYZrot(:,i), &
				&			i, Sref(i), XYZref(:,i), norm2(XYZrot(:,i)-XYZref(:,i)), Wref(i), i=1, Natm)
				write(Lsta, *)
			end if

			write(Lout, *) Natm
			write(Lout, '(A, 2(A, F12.6))') trim(Tips), ' RMSD(Fit)=', RMSD, ' minRMSD=', RMSDmin
			do i = 1, Natm
				write(Lout, '(A4, 3F14.9)') Smol(i), XYZrot(:,i)
			end do
		end do

	end do

	close(Lout)
	close(Lxyz)

	print*
	print*, '>>FitMol Achieved.'

	if (.NOT.YesBat) then
		print*, '>>Any Key to Exit...'
		read*
	end if

	print*
	stop '>>Normal Termination of Program FitMol.'

End Program FitMol
!
!
Subroutine center(Icnt, n, x, w, Xcnt)
! center a molecule about its weighted centroid or other origin
	integer i, n, Icnt
	real*8  x(3, *), w(*), Xcnt(3), wnorm

	if(Icnt/=0) then
		Xcnt(1:3) = x(1:3, Icnt)
	else
		wnorm = 1.d0/sum(w(1:n))
		Xcnt(1) = dot_product( x(1,1:n), w(1:n) )*wnorm
		Xcnt(2) = dot_product( x(2,1:n), w(1:n) )*wnorm
		Xcnt(3) = dot_product( x(3,1:n), w(1:n) )*wnorm
	end if

	do	i = 1, n
		x(:, i) = x(:, i) - Xcnt(:)
	end do
End
!
!
Subroutine analyz(n, x, y, w, u)
! analyze the x to y fit and the fitting matrix

	integer i, j, n, Lsta
	real*8  x(3, n), y(3, n), w(n), u(3, 3)
	real*8  Rerr, wnorm, urnorm (3), ucnorm (3)
	common /Lsta/ Lsta

! find the mean sum of squares error
	Rerr = 0.0D0
	wnorm = 0.0D0
	do	i = 1, n
		Rerr = Rerr + w (i) * ( (y (1, i) - x (1, i))**2 + &
	&							(y (2, i) - x (2, i))**2 + &
	&							(y (3, i) - x (3, i))**2 )
		wnorm = wnorm + w (i)
	end do
	Rerr = SQRT (Rerr / wnorm)

! find the row and column norms of u
	ucnorm (1) = u (1, 1)**2 + u (2, 1)**2 + u (3, 1)**2
	ucnorm (2) = u (1, 2)**2 + u (2, 2)**2 + u (3, 2)**2
	ucnorm (3) = u (1, 3)**2 + u (2, 3)**2 + u (3, 3)**2
	urnorm (1) = u (1, 1)**2 + u (1, 2)**2 + u (1, 3)**2
	urnorm (2) = u (2, 1)**2 + u (2, 2)**2 + u (2, 3)**2
	urnorm (3) = u (3, 1)**2 + u (3, 2)**2 + u (3, 3)**2

! write the error and u norms
	write (Lsta, *) '       [Rot Matrix]                                    |Row|'
	do	i = 1, 3
		write (Lsta, '(5X, 3F14.9, F14.3)') (u (i, j), j = 1, 3), urnorm (i)
	end do
	write (Lsta, '(8X, A, F6.3, 2F14.3, F14.9, A)') '|Col|', (ucnorm (j), j = 1, 3), Rerr, '=RMSD'
End
!
!
Subroutine qtrfit(n, x, y, w, q, u, nr)
! Find the quaternion, q, [and left rotation matrix, u] that minimizes
!   |qTXq - Y| ^ 2    [|uX - Y| ^ 2]
!
! This is equivalent to maximizing Re (qTXTqY).
!
! This is equivalent to finding the largest eigenvalue and corresponding
! eigenvector of the matrix
!   [A2   AUx  AUy  AUz ]
!   [AUx  Ux2  UxUy UzUx]
!   [AUy  UxUy Uy2  UyUz]
!   [AUz  UzUx UyUz Uz2 ]
!
! where
!   A2   = Xx Yx + Xy Yy + Xz Yz
!   Ux2  = Xx Yx - Xy Yy - Xz Yz
!   Uy2  = Xy Yy - Xz Yz - Xx Yx
!   Uz2  = Xz Yz - Xx Yx - Xy Yy
!   AUx  = Xz Yy - Xy Yz
!   AUy  = Xx Yz - Xz Yx
!   AUz  = Xy Yx - Xx Yy
!   UxUy = Xx Yy + Xy Yx
!   UyUz = Xy Yz + Xz Yy
!   UzUx = Xz Yx + Xx Yz
!
! The left rotation matrix, u, is obtained from q by
!   u = qT1q
!
! INPUT
!   n      - number of points
!   x      - test vector
!   y      - reference vector
!   w      - weight vector
!
! OUTPUT
!   q      - the best-fit quaternion
!   u      - the best-fit left rotation matrix
!   nr     - number of jacobi sweeps required
!
	integer i, n, nr
	real*8  xxyx, xxyy, xxyz, xyyx, xyyy, xyyz, xzyx, xzyy, xzyz, &
	&		x(3, n), y(3, n), w(n), q(0 : 3), u(3, 3), &
	&		c(0 : 3, 0 : 3), v(0 : 3, 0 : 3), d(0 : 3)

! generate the upper triangle of the quadratic form matrix
	xxyx = 0.0D0
	xxyy = 0.0D0
	xxyz = 0.0D0
	xyyx = 0.0D0
	xyyy = 0.0D0
	xyyz = 0.0D0
	xzyx = 0.0D0
	xzyy = 0.0D0
	xzyz = 0.0D0
	do	i = 1, n
		xxyx = xxyx + x (1, i) * y (1, i) * w (i)
		xxyy = xxyy + x (1, i) * y (2, i) * w (i)
		xxyz = xxyz + x (1, i) * y (3, i) * w (i)
		xyyx = xyyx + x (2, i) * y (1, i) * w (i)
		xyyy = xyyy + x (2, i) * y (2, i) * w (i)
		xyyz = xyyz + x (2, i) * y (3, i) * w (i)
		xzyx = xzyx + x (3, i) * y (1, i) * w (i)
		xzyy = xzyy + x (3, i) * y (2, i) * w (i)
		xzyz = xzyz + x (3, i) * y (3, i) * w (i)
	end do

	c(0, 0) = xxyx + xyyy + xzyz

	c(0, 1) = xzyy - xyyz
	c(1, 1) = xxyx - xyyy - xzyz

	c(0, 2) = xxyz - xzyx
	c(1, 2) = xxyy + xyyx
	c(2, 2) = xyyy - xzyz - xxyx

	c(0, 3) = xyyx - xxyy
	c(1, 3) = xzyx + xxyz
	c(2, 3) = xyyz + xzyy
	c(3, 3) = xzyz - xxyx - xyyy

! diagonalize c
	nr = 16
	call jacobi (c, 4, 4, d, v, nr)

! extract the desired quaternion
	q (0) = v (0, 3)
	q (1) = v (1, 3)
	q (2) = v (2, 3)
	q (3) = v (3, 3)

! Generate a left rotation matrix from a normalized quaternion
	u (1, 1) = q (0) ** 2 + q (1) ** 2 - q (2) ** 2 - q (3) ** 2
	u (2, 1) = 2.0D0 * (q (1) * q (2) - q (0) * q (3))
	u (3, 1) = 2.0D0 * (q (1) * q (3) + q (0) * q (2))

	u (1, 2) = 2.0D0 * (q (2) * q (1) + q (0) * q (3))
	u (2, 2) = q (0) ** 2 - q (1) ** 2 + q (2) ** 2 - q (3) ** 2
	u (3, 2) = 2.0D0 * (q (2) * q (3) - q (0) * q (1))

	u (1, 3) = 2.0D0 * (q (3) * q (1) - q (0) * q (2))
	u (2, 3) = 2.0D0 * (q (3) * q (2) + q (0) * q (1))
	u (3, 3) = q (0) ** 2 - q (1) ** 2 - q (2) ** 2 + q (3) ** 2
End
!
!
Subroutine rotmol (n, x, y, u)
! rotate a molecule
	integer i, n
	real*8  x(3, n), y(3, n), u(3, 3)

	do	i = 1, n
		y (1, i) = u(1, 1) * x(1, i) + u(1, 2) * x(2, i) + u(1, 3) * x(3, i)
		y (2, i) = u(2, 1) * x(1, i) + u(2, 2) * x(2, i) + u(2, 3) * x(3, i)
		y (3, i) = u(3, 1) * x(1, i) + u(3, 2) * x(2, i) + u(3, 3) * x(3, i)
	end do
End
!
!
Subroutine jacobi (a, n, np, d, v, nrot)
! Jacobi diagonalizer with sorted output.  Same calling sequence as
! EISPACK routine, but must specify nrot!

	integer i, j, k, l, n, np, nrot
	real*8  a(np, n), d(n), v(np, n), onorm, dnorm, b, dma, q, t, c, s, atemp, vtemp, dtemp

	do	j = 1, n
		do	i = 1, n
			v (i, j) = 0.0D0
		end do
		v (j, j) = 1.0D0
		d (j) = a (j, j)
	end do

	do	l = 1, nrot
		dnorm = 0.0D0
		onorm = 0.0D0
		do	j = 1, n
			dnorm = dnorm + ABS (d (j))
			do	i = 1, j - 1
				onorm = onorm + ABS (a (i, j))
			end do
		end do

		if (onorm / dnorm .LE. 0.0D0) goto 19999

		do	j = 2, n
			do	i = 1, j - 1
				b = a (i, j)
				if (ABS (b) .GT. 0.0D0) then
					dma = d (j) - d (i)
					if (ABS (dma) + ABS (b) .LE. ABS (dma)) then
						t = b / dma
					else
						q = 0.5D0 * dma / b
						t = SIGN (1.0D0 / (ABS (q) + SQRT (1.0D0 + q * q)), q)
					end if
					c = 1.0D0 / SQRT (t * t + 1.0D0)
					s = t * c
					a (i, j) = 0.0D0
					do	k = 1, i - 1
						atemp    = c * a (k, i) - s * a (k, j)
						a (k, j) = s * a (k, i) + c * a (k, j)
						a (k, i) = atemp
					end do
					do	k = i + 1, j - 1
						atemp    = c * a (i, k) - s * a (k, j)
						a (k, j) = s * a (i, k) + c * a (k, j)
						a (i, k) = atemp
					end do
					do	k = j + 1, n
						atemp    = c * a (i, k) - s * a (j, k)
						a (j, k) = s * a (i, k) + c * a (j, k)
						a (i, k) = atemp
					end do
					do	k = 1, n
						vtemp    = c * v (k, i) - s * v (k, j)
						v (k, j) = s * v (k, i) + c * v (k, j)
						v (k, i) = vtemp
					end do
					dtemp = c * c * d (i) + s * s * d (j) - 2.0D0 * c * s * b
					d (j) = s * s * d (i) + c * c * d (j) + 2.0D0 * c * s * b
					d (i) = dtemp
				end if
			end do
		end do
	end do

19999 continue

	nrot = l

	do	j = 1, n - 1
		k = j
		dtemp = d (k)
		do	i = j + 1, n
			if (d (i) .LT. dtemp) then
				k = i
				dtemp = d (k)
			end if
		end do
		if (k .GT. j) then
			d (k) = d (j)
			d (j) = dtemp
			do	i = 1, n
				dtemp    = v (i, k)
				v (i, k) = v (i, j)
				v (i, j) = dtemp
			end do
		end if
	end do
End

subroutine NextPermut(Ndat, P, YesDone)
	integer i, j, tmp, Ibeg, Iend, Ndat, YesDone, P(*)

	if (Ndat==0) then; YesDone=1; return; end if

	if(YesDone==-1) then
		YesDone=0
		do i=1, Ndat; P(i)=i; end do
		return
	end if

	!从后向前查找，看是否有后面的数大于前面的数的情况P(i-1)<Pi
	!若有则停在后一个数的位置。
	!若没有后面的数大于前面的数的情况，说明已经到了最后一个排列，返回
	do Iend=Ndat, 1, -1
		if ( P(Iend-1) < P(Iend) ) exit
	end do

	if (Iend==1) then; YesDone=1; return; end if

	!从后查到Iend，查找大于P(Iend-1)的最小的数，记入Ibeg
	Ibeg=Iend
	do i=Ndat, Iend, -1
		if (P(Iend-1)<P(i) .and. P(i)<P(Ibeg)) Ibeg=i
	end do

	!交换p(Ibeg)和p(Iend-1)
	tmp=P(Ibeg); P(Ibeg)=P(Iend-1); P(Iend-1)=tmp

	!倒置p(Iend)到p(Ndat)
	j=Ndat
	do i=Iend, Ndat
		tmp=P(j); P(j)=P(i); P(i)=tmp
		j=j-1
		if(i>=j) exit
	end do
end

subroutine qp()
	integer, parameter:: N=6
	integer idx(0:N), pidx(0:N)
	integer num, Idone

   !print*, "\nEXAMPLE_01 of head permutations using a linear array of %d objects.\n",N
   !display(a);          // remove comment to display linear array a[] (optional)

	call NextQuickPerm(0, N, idx, pidx)

	num=1
	Idone = 1;   ! setup first swap points to be 1 and 0 respectively (i & j)
	do while(Idone < N)
		num=num+1
		call NextQuickPerm(Idone, N, idx, pidx)
		print*, idx
	end do
	print*,  num

end

subroutine NextQuickPerm(i, n, a, p)
	integer i, j, t, a(0:*), p(0:*)

	if(i==0) then
		do j=0, n      ! initialize arrays; a[N] can be any type
			p(j) = j   ! p[N] > 0 controls iteration and the index boundary for i
			a(j) = j+1 ! a[i] value is not revealed and can be arbitrary
		end do
		i=1
		return
	end if

	p(i) = p(i)-1                   ! decrease index "weight" for i by one
	j = mod(i,2) * p(i)             ! IF i is odd then j = p[i] otherwise j = 0
	t = a(j); a(j) = a(i); a(i) = t ! swap(a[i], a[j])

	i = 1                           ! reset index i to 1 (assumed)
	do while(p(i) == 0)
		p(i) = i                    ! reset p[i] zero value
		i = i+1                     ! set new index value for i (increase by one)
	end do
end
