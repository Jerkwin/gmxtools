Program FitMol
	integer, parameter:: n = 2000, Lxyz = 1, Lout = 2

	integer i, ierr, Narg, Lsta, Natm, Icnt
	real*8  RMSD, x(3, n), y(3, n), z(3, n), w(n), Xcnt(3), q(4), u(3, 3)
	logical	YesBat, YesOrg
	character*8 Satm(n)
	character*256	txt, Finp, Tite
	common /Lsta/ Lsta

	Icnt = 0
	YesBat = .false.
	YesOrg = .false.

	Narg = Iargc()-1
	if(Narg>0) then
		YesBat = .true.
		call getarg(1, Finp)
		i = index(Finp, '.xyz')-1
		if(i>0) Finp = Finp(1:i)
		open(unit=Lxyz, file=trim(Finp)//'.xyz', status='old', IOstat=Ierr)
		open(unit=Lout, file=trim(Finp)//'~Fit.xyz', IOstat=Ierr)
		read(Lxyz, *, IOstat=Ierr) Natm
		read(Lxyz, '(A256)', IOstat=Ierr) Tite
		if(Narg==2) then
			call getarg(2, Finp)
			read(Finp, *) Icnt
			YesOrg = .true.
		end if
		goto 100
	end if

	print*, '>>>>>>>>>>>>>>>>Program FitMol Running <<<<<<<<<<<<<<<<'
	print*,	'>>>>>>>>>>>>>>>>    Jicun LI           <<<<<<<<<<<<<<<<'
	print*, '>>>>>>>>>>>>>>>>  2019-08-21 14:44:52  <<<<<<<<<<<<<<<<'

	print*, '>>Type XYZ File Name (*.xyz/Exit):'
	read*, Finp
	if (index(Finp, 'exit')/=0 .or. index(Finp, 'Exit')/=0) &
	&	stop '>>Normal Exit of Program FitMol.'

	i = index(Finp, '.xyz')-1
	if(i>0) Finp = Finp(1:i)
	open(unit=Lxyz, file=trim(adjustl(Finp))//'.xyz', status='old', IOstat=Ierr)
	do  while(Ierr/=0)
		print*, '>>XYZ Input File  .\'//trim(adjustl(Finp))//'.xyz  NOT Exist !'
		print*, '>>Revise Input File Name (*.xyz/Exit):'
		read*, Finp
		if (index(Finp, 'exit')/=0 .or. index(Finp, 'Exit')/=0) &
	&		stop '>>Normal Exit of Program FitMol.'
		i = index(Finp, '.xyz')-1
		if(i>0) Finp = Finp(1:i)
		open(unit=Lxyz, file=trim(adjustl(Finp))//'.xyz', status='old', IOstat=Ierr)
	end do
	open(unit=Lout, file=trim(adjustl(Finp))//'~Fit.xyz', IOstat=Ierr)

	read(Lxyz, *, IOstat=Ierr) Natm
	read(Lxyz, '(A256)', IOstat=Ierr) Tite
	read(Tite, *, IOstat=Ierr) Icnt
	if (Ierr==0) YesOrg = .true.

	print*, '>>Input File is  .\', trim(adjustl(Finp))//'.xyz'
100	print*, '>>Input File Readin.'
	print*
	write(*, '(A, I4)') '>>>>No. of Atoms         ', Natm
	if(YesOrg) write(*, '(A, I4)') '>>>>New Origin Fixed at Atom     ', Icnt

	if (.NOT.YesBat) then
		Lsta = 6
		print*
		print*, '>>FitMol Calculation Will Run.'
		print*
		print*, '>>Run File Now ? Y/N'
		read*, txt
		if (txt.EQ.'N' .or. txt.EQ.'n') stop '>>Normal Exit of Program FitMol.'
	end if

	print*
	print*, '>>FitMol Calculation Running...'

	read(Lxyz, '(A256)', IOstat=Ierr) Tite
	read(Tite, *, IOstat=Ierr) txt, txt, txt, txt, w(1)
	rewind(Lxyz)

	read(Lxyz, *)
	read(Lxyz, *) Tite
	if(Ierr==0) then
		read(Lxyz, *) (Satm(i), x(1, i), x(2, i), x(3, i), w(i), i=1, Natm)
	else
		w = 1.D0
		read(Lxyz, *) (Satm(i), x(1, i), x(2, i), x(3, i), i=1, Natm)
	end if

	if (YesOrg) Xcnt(1:3) = x(1:3, Icnt)
	call center(Natm, x, w, YesOrg, Xcnt)

	if (.NOT. YesBat) write(Lsta, '(A, /, 3F14.9/)') 'Center of Reference Molecule X', Xcnt

	write(Lout, *) Natm
	write(Lout, '(A)') trim(Tite)//' Ref Mol'
	write(Lout, '(A4, 3F14.9, F8.1)') (Satm(i), x(1, i), x(2, i), x(3, i), w(i), i=1, Natm)

	Ierr=0
	do
		read(Lxyz, *, IOstat=Ierr) Natm
		if (Ierr/=0) exit
		read(Lxyz, '(A256)', IOstat=Ierr) Tite
		read(Lxyz, *) (txt, y(1, i), y(2, i), y(3, i), i=1, Natm)

		if (YesOrg) Xcnt(1:3) = y(1:3, Icnt)
		call center(Natm, y, w, YesOrg, Xcnt)

		if (.NOT. YesBat)  then
			write(Lsta, '(/, A, /, 3F14.9/)') '>>>>Center of Fitting Molecule Y', Xcnt
			write(Lsta, '(A)') 'Atom     Xref          Yref          Zref' &
	&						//'          Xfit          Yfit          Zfit'
			do	i=1, Natm
				write(Lsta, '(X, A4, 6F14.9)') &
	&				Satm(i), x(1,i), x(2,i), x(3,i), y(1,i), y(2,i), y(3,i)
			end do
			write(Lsta, *)
		end if

		u = 0.d0; u(1,1) = 1.d0; u(2,2) = 1.d0; u(3,3) = 1.d0
		if(.NOT.YesBat) then
			call analyz (n, x, y, w, u)
			write (Lsta, '(/, A)') 'Fitting    X = RY'
		end if

		call qtrfit (Natm, y, x, w, q, u, ierr)
		call rotmol (Natm, y, z, u)
		if (.NOT.YesBat) call analyz (Natm, z, x, w, u)

		if (YesOrg) Xcnt(1:3) = z(1:3, Icnt)
		call center(Natm, z, w, YesOrg, Xcnt)

		if (.NOT.YesBat) then
			write(Lsta, *)
			write(Lsta, '(A)') 'Atom      RYx           RYy           RYz' &
	&						// '           Xx            Xy            Xz            |RY-X|'
			do	i=1, Natm
				write(Lsta, '(X, A4, 7F14.9)') Satm(i), z(1,i), z(2,i), z(3,i), &
			&				x(1,i), x(2,i), x(3,i), sqrt( (x(1,i)-z(1,i))**2 &
			&				+ (x(2,i)-z(2,i))**2 + (x(3,i)-z(3,i))**2)
			end do
		end if

		RMSD = 0.D0
		do	i = 1, Natm
			RMSD = RMSD + w(i)* ( (x(1,i)-z(1,i))**2 + (x(2,i)-z(2,i))**2 &
	&							+ (x(3,i)-z(3,i))**2 )
		end do
		RMSD = sqrt(RMSD/dble(Natm))

		write(Lout, *) Natm
		write(Lout, '(A, A, F12.6)') trim(Tite), ' RMSD(Fit)=', RMSD
		do	i = 1, Natm
			write(Lout, '(A4, 3F14.9)') Satm(i), z(1,i), z(2,i), z(3,i)
		end do
	end do

	close(Lout)
	close(Lxyz)

	print*
	print*, '>>FitMol Calculation Achieved.'

	if (.NOT.YesBat) then
		print*, '>>Any Key to Exit...'
		read*
	end if

	print*
	stop '>>Normal Termination of Program FitMol.'

End Program FitMol
!
!
Subroutine analyz (n, x, y, w, u)
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
	write (Lsta, *) '   Rot Matrix                                  Row Norm'
	do	i = 1, 3
		write (Lsta, '(3F14.9, 4X, F10.4)') (u (i, j), j = 1, 3), urnorm (i)
	end do
	write (Lsta, *) '   Col Norm                                        RMSD'
	write (Lsta, '(3F14.4, F14.9)') (ucnorm (j), j = 1, 3), Rerr

End
!
!
Subroutine center (n, x, w, YesOrg, Xcnt)
! center a molecule about its weighted centroid or other origin
	integer i, n
	real*8  x(3, n), w(n), Xcnt(3), wnorm
	logical YesOrg

	if(.NOT.YesOrg) then
		Xcnt = 0.0D0
		wnorm = 0.0D0
		do	i = 1, n
			Xcnt(1) = Xcnt(1) + x(1, i) * w(i)
			Xcnt(2) = Xcnt(2) + x(2, i) * w(i)
			Xcnt(3) = Xcnt(3) + x(3, i) * w(i)
			wnorm = wnorm + w (i)
		end do
		Xcnt(1) = Xcnt(1)/wnorm
		Xcnt(2) = Xcnt(2)/wnorm
		Xcnt(3) = Xcnt(3)/wnorm
	end if

	do	i = 1, n
		x(1, i) = x(1, i) - Xcnt(1)
		x(2, i) = x(2, i) - Xcnt(2)
		x(3, i) = x(3, i) - Xcnt(3)
	end do
End
!
!
Subroutine qtrfit (n, x, y, w, q, u, nr)
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

	c (0, 0) = xxyx + xyyy + xzyz

	c (0, 1) = xzyy - xyyz
	c (1, 1) = xxyx - xyyy - xzyz

	c (0, 2) = xxyz - xzyx
	c (1, 2) = xxyy + xyyx
	c (2, 2) = xyyy - xzyz - xxyx

	c (0, 3) = xyyx - xxyy
	c (1, 3) = xzyx + xxyz
	c (2, 3) = xyyz + xzyy
	c (3, 3) = xzyz - xxyx - xyyy

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
