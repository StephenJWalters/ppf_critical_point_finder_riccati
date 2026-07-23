! Riccati shooting method for Orr-Sommerfeld, fully analytic version, apart from Jacobian curvature terms.
!
! R, S^Re, S^alpha, S^c are all propagated together in ONE integration pass.
! The inner Newton solve for c uses the analytic dF/dc = f(R,S^c) -- no
! finite difference in c anywhere. Once c has converged, dci_dRe and f2
! (=dci_dalpha) fall out of the SAME pass via S^Re, S^alpha -- no separate
! "extended" call needed. Only the alpha-curvature terms (d2ci_dalpha2,
! d2ci_dalpha_dRe) still use a finite difference, of already-analytic f2/
! dci_dRe values across a small step h_alpha.
!
! compile: gfortran critical_point_finder_riccati4.f90 -O3 -o riccati_finder

module riccati_mod
    implicit none
    integer, parameter :: qp = selected_real_kind(33)
	integer, parameter :: rkstages=35
	real(qp), dimension(rkstages,rkstages) 	:: Aij
	real(qp), dimension(rkstages) 			:: Bj, Ci
	
contains

subroutine rkcoeffs(method,steps,a,b)
	implicit none
	integer,intent(in)						:: steps
	character(*),intent(in)					:: method
	real,intent(out),dimension(steps)		:: b
	real,intent(out),dimension(steps,steps)	:: a
	integer::i,j
	real::line
	a=0.
	open (unit=99, file='rkcoeffs/'//method//'Aij.txt', status='old', action='read')
	i=2
	j=1
	do
		read(99, *,END=200) line
		a(i,j)=line
		if (i .eq. j+1) then
			i=i+1
			j=1
		else
			j=j+1
		endif
	enddo
	200 continue
	close(99)
	open (unit=99, file='rkcoeffs/'//method//'B.txt', status='old', action='read')
	i=1
	do
		read(99, *,END=201) line
		!write(*,*) line
		b(i)=line
		i=i+1
	enddo
	201 continue
	close(99)
end subroutine

! ------------------------------------------------------------------
! d(y), e(y) and ALL first partials (Re, alpha, c). d_c, e_c are
! constants (no y-dependence) but returned here for uniformity.
! ------------------------------------------------------------------
subroutine coeffs_full(alpha, Re, c, y, d, e, d_Re, e_Re, d_alpha, e_alpha, d_c, e_c)
	real(qp),    intent(in)  :: alpha, Re, y
	complex(qp), intent(in)  :: c
	complex(qp), intent(out) :: d, e, d_Re, e_Re, d_alpha, e_alpha, d_c, e_c
	complex(qp) :: ii, Uc
	real(qp)    :: U

	ii = cmplx(0._qp, 1._qp, qp)
	U  = 1._qp - y**2
	Uc = cmplx(U, 0._qp, qp) - c

	e = 2._qp*alpha**2 + ii*alpha*Re*Uc
	d = -alpha**4 - ii*alpha*Re*(alpha**2*Uc - 2._qp)

	e_Re = ii*alpha*Uc
	d_Re = ii*alpha*(2._qp - alpha**2*Uc)

	e_alpha = 4._qp*alpha + ii*Re*Uc
	d_alpha = -4._qp*alpha**3 - 3._qp*ii*alpha**2*Re*Uc + 2._qp*ii*Re

	e_c = -ii*alpha*Re
	d_c =  ii*alpha**3*Re
end subroutine coeffs_full

! ------------------------------------------------------------------
! RHS for R plus all three sensitivities S^Re, S^alpha, S^c together
! ------------------------------------------------------------------
subroutine riccati_rhs_full(alpha, Re, c, y, R, SRe, Salpha, Sc, &
							 Rdot, SRedot, Salphadot, Scdot)
	real(qp),    intent(in)  :: alpha, Re, y
	complex(qp), intent(in)  :: c
	complex(qp), intent(in)  :: R(2,2), SRe(2,2), Salpha(2,2), Sc(2,2)
	complex(qp), intent(out) :: Rdot(2,2), SRedot(2,2), Salphadot(2,2), Scdot(2,2)
	complex(qp) :: d, e, d_Re, e_Re, d_alpha, e_alpha, d_c, e_c
	complex(qp) :: g, g_Re, g_alpha, g_c

	call coeffs_full(alpha, Re, c, y, d, e, d_Re, e_Re, d_alpha, e_alpha, d_c, e_c)

	g = d*R(1,1) + e
	Rdot(1,1) = R(2,1) - R(1,2)*g
	Rdot(1,2) = R(2,2) - R(1,1) - d*R(1,2)**2
	Rdot(2,1) = cmplx(1._qp,0._qp,qp) - R(2,2)*g
	Rdot(2,2) = -R(2,1) - d*R(1,2)*R(2,2)

	g_Re = d_Re*R(1,1) + d*SRe(1,1) + e_Re
	SRedot(1,1) = SRe(2,1) - SRe(1,2)*g - R(1,2)*g_Re
	SRedot(1,2) = SRe(2,2) - SRe(1,1) - d_Re*R(1,2)**2 - 2._qp*d*R(1,2)*SRe(1,2)
	SRedot(2,1) = -SRe(2,2)*g - R(2,2)*g_Re
	SRedot(2,2) = -SRe(2,1) - d_Re*R(1,2)*R(2,2) - d*(SRe(1,2)*R(2,2) + R(1,2)*SRe(2,2))

	g_alpha = d_alpha*R(1,1) + d*Salpha(1,1) + e_alpha
	Salphadot(1,1) = Salpha(2,1) - Salpha(1,2)*g - R(1,2)*g_alpha
	Salphadot(1,2) = Salpha(2,2) - Salpha(1,1) - d_alpha*R(1,2)**2 - 2._qp*d*R(1,2)*Salpha(1,2)
	Salphadot(2,1) = -Salpha(2,2)*g - R(2,2)*g_alpha
	Salphadot(2,2) = -Salpha(2,1) - d_alpha*R(1,2)*R(2,2) &
					  - d*(Salpha(1,2)*R(2,2) + R(1,2)*Salpha(2,2))

	g_c = d_c*R(1,1) + d*Sc(1,1) + e_c
	Scdot(1,1) = Sc(2,1) - Sc(1,2)*g - R(1,2)*g_c
	Scdot(1,2) = Sc(2,2) - Sc(1,1) - d_c*R(1,2)**2 - 2._qp*d*R(1,2)*Sc(1,2)
	Scdot(2,1) = -Sc(2,2)*g - R(2,2)*g_c
	Scdot(2,2) = -Sc(2,1) - d_c*R(1,2)*R(2,2) - d*(Sc(1,2)*R(2,2) + R(1,2)*Sc(2,2))
end subroutine riccati_rhs_full

! ------------------------------------------------------------------
! integrate R, S^Re, S^alpha, S^c together, general explicit RK tableau
! (every "...i" work array is reset to its base value every stage --
!  this is the fix for the earlier NaN bug, extended to all 4 matrices)
! ------------------------------------------------------------------
subroutine integrate_riccati_full(alpha, Re, c, nsteps, Rfinal, SRefinal, Salphafinal, Scfinal)
	real(qp),    intent(in)  :: alpha, Re
	complex(qp), intent(in)  :: c
	integer,     intent(in)  :: nsteps
	complex(qp), intent(out), dimension(2,2) :: Rfinal, SRefinal, Salphafinal, Scfinal

	complex(qp), dimension(2,2) :: R, Ri, SRe, SRei, Salpha, Salphai, Sc, Sci
	complex(qp), dimension(2,2,rkstages) :: R_prime, SR_prime, SA_prime, SC_prime
	real(qp)    :: y, h, yi
	integer     :: istep, s, k

	h = 2._qp / real(nsteps, qp)
	R = cmplx(0._qp,0._qp,qp); SRe = cmplx(0._qp,0._qp,qp)
	Salpha = cmplx(0._qp,0._qp,qp); Sc = cmplx(0._qp,0._qp,qp)
	y = -1._qp

	do istep = 1, nsteps
		do s = 1, rkstages
			yi      = y + ci(s)*h
			Ri      = R
			SRei    = SRe
			Salphai = Salpha
			Sci     = Sc
			do k = 1, s-1
				Ri      = Ri      + Aij(s,k)*R_prime(:,:,k)*h
				SRei    = SRei    + Aij(s,k)*SR_prime(:,:,k)*h
				Salphai = Salphai + Aij(s,k)*SA_prime(:,:,k)*h
				Sci     = Sci     + Aij(s,k)*SC_prime(:,:,k)*h
			end do
			call riccati_rhs_full(alpha, Re, c, yi, Ri, SRei, Salphai, Sci, &
								   R_prime(:,:,s), SR_prime(:,:,s), SA_prime(:,:,s), SC_prime(:,:,s))
		end do
		do s = 1, rkstages
			R      = R      + Bj(s)*R_prime(:,:,s)*h
			SRe    = SRe    + Bj(s)*SR_prime(:,:,s)*h
			Salpha = Salpha + Bj(s)*SA_prime(:,:,s)*h
			Sc     = Sc     + Bj(s)*SC_prime(:,:,s)*h
		end do
		y = y + h
	end do

	Rfinal = R; SRefinal = SRe; Salphafinal = Salpha; Scfinal = Sc
end subroutine integrate_riccati_full

! det(R) and its three analytic parameter-derivatives, product rule
subroutine dF_all(R, SRe, Salpha, Sc, F, dFdRe, dFdalpha, dFdc)
	complex(qp), intent(in)  :: R(2,2), SRe(2,2), Salpha(2,2), Sc(2,2)
	complex(qp), intent(out) :: F, dFdRe, dFdalpha, dFdc

	F        = R(1,1)*R(2,2) - R(1,2)*R(2,1)
	dFdRe    = R(2,2)*SRe(1,1)    - R(2,1)*SRe(1,2)    - R(1,2)*SRe(2,1)    + R(1,1)*SRe(2,2)
	dFdalpha = R(2,2)*Salpha(1,1) - R(2,1)*Salpha(1,2) - R(1,2)*Salpha(2,1) + R(1,1)*Salpha(2,2)
	dFdc     = R(2,2)*Sc(1,1)     - R(2,1)*Sc(1,2)     - R(1,2)*Sc(2,1)     + R(1,1)*Sc(2,2)
end subroutine dF_all

! ------------------------------------------------------------------
! Unified analytic Newton for c at fixed (alpha,Re). No FD anywhere.
! Returns the converged c AND the sensitivities at that point, so the
! caller gets dci_dRe / f2 "for free" from the same integration pass.
! ------------------------------------------------------------------
subroutine find_c_analytic(alpha, Re, nsteps, c_guess, tol, maxit, c, dci_dRe, f2)
	real(qp),    intent(in)  :: alpha, Re, tol
	integer,     intent(in)  :: nsteps, maxit
	complex(qp), intent(in)  :: c_guess
	complex(qp), intent(out) :: c
	real(qp),    intent(out) :: dci_dRe, f2

	complex(qp), dimension(2,2) :: R, SRe, Salpha, Sc
	complex(qp) :: F, dFdRe, dFdalpha, dFdc
	integer     :: iter

	c = c_guess
	do iter = 1, maxit
		call integrate_riccati_full(alpha, Re, c, nsteps, R, SRe, Salpha, Sc)
		call dF_all(R, SRe, Salpha, Sc, F, dFdRe, dFdalpha, dFdc)
		if (abs(F) < tol) exit
		c = c - F/dFdc
	end do

	dci_dRe = aimag(-dFdRe/dFdc)
	f2      = aimag(-dFdalpha/dFdc)
end subroutine find_c_analytic

! ------------------------------------------------------------------
! residual + Jacobian for the outer 2x2 Newton loop.
! f1, f2, dci_dRe: fully analytic. Curvature terms: FD of f2/dci_dRe
! across h_alpha (only remaining finite difference in the pipeline).
! ------------------------------------------------------------------
subroutine ci_and_jac_riccati(alpha, Re, nsteps, tol_c, h_alpha, c_inout, &
							   f1, f2, dci_dRe, d2ci_dalpha2, d2ci_dalpha_dRe)
	real(qp),    intent(in)    :: alpha, Re, tol_c, h_alpha
	integer,     intent(in)    :: nsteps
	complex(qp), intent(inout) :: c_inout
	real(qp),    intent(out)   :: f1, f2, dci_dRe, d2ci_dalpha2, d2ci_dalpha_dRe

	complex(qp) :: c0, c_p, c_m
	real(qp)    :: f2_p, f2_m, dci_dRe_p, dci_dRe_m

	call find_c_analytic(alpha, Re, nsteps, c_inout, tol_c, 30, c0, dci_dRe, f2)
	f1 = aimag(c0)

	call find_c_analytic(alpha+h_alpha, Re, nsteps, c0, tol_c, 30, c_p, dci_dRe_p, f2_p)
	call find_c_analytic(alpha-h_alpha, Re, nsteps, c0, tol_c, 30, c_m, dci_dRe_m, f2_m)

	d2ci_dalpha2    = (f2_p - f2_m) / (2._qp*h_alpha)
	d2ci_dalpha_dRe = (dci_dRe_p - dci_dRe_m) / (2._qp*h_alpha)

	c_inout = c0
end subroutine ci_and_jac_riccati

subroutine solve2x2(J, rhs, dx)
	real(qp), intent(in)  :: J(2,2), rhs(2)
	real(qp), intent(out) :: dx(2)
	real(qp) :: det
	det   = J(1,1)*J(2,2) - J(1,2)*J(2,1)
	dx(1) = (rhs(1)*J(2,2) - rhs(2)*J(1,2)) / det
	dx(2) = (J(1,1)*rhs(2) - J(2,1)*rhs(1)) / det
end subroutine solve2x2

end module riccati_mod


program riccati_finder
    use riccati_mod
    implicit none
    real(qp)    :: Re, alpha, Re_c, alpha_c, Re_c_ref, alpha_c_ref, c_r_ref
    real(qp)    :: f1, f2, dci_dRe, d2ci_dalpha2, d2ci_dalpha_dRe
    real(qp)    :: J(2,2), rhs(2), dx(2), tol, h_alpha
    complex(qp) :: c
    integer     :: nsteps, iter, max_iter, istage, gridsize, i, first_diff
	real(8)		:: t1,t2
	character(len=45) :: str_ref, str_cur
	
	call cpu_time(t1)
	! --- 1. Initialization ---
    call rkcoeffs('feagin14', rkstages, Aij, Bj)
    do istage = 1, rkstages
        ci(istage) = sum(Aij(istage,:))
    end do
	
	! converged triple from the N=256 NSC spectral run being validated
	Re_c_ref    = 5772.22181620969823713375856640255277_qp
	alpha_c_ref	= 1.02054744928533350624946400664133281_qp
	c_r_ref 	= 0.264000260478798443339340265061735948_qp
	
	open(unit=10, file='convergence_table_riccati.tex', status='replace')
	write(10,'(a)') '\begin{table}[h]'
	write(10,'(a)') '\centering'
	write(10,'(a)') '\begin{tabular}{cl}'
	write(10,'(a)') '\hline'
	write(10,'(a)') '$N$ & $Re_c$ \\'
	write(10,'(a)') '\hline'
	write(str_ref, '(f45.35)') Re_c_ref
	
	! --- 2. Newton-Raphson Iteration ---
	do gridsize=1,5
		! start from Orszag's triple each time
		Re    = 5772.22_qp
		alpha = 1.02056_qp
		c     = cmplx(0.264_qp, 0._qp, qp)

		nsteps   = 4000*2**gridsize
		h_alpha  = 1.0e-9_qp
		tol      = 1.0e-27_qp
		max_iter = 15
		
		do iter = 1, max_iter
			call ci_and_jac_riccati(alpha, Re, nsteps, 1.0e-29_qp, h_alpha, c, &
									 f1, f2, dci_dRe, d2ci_dalpha2, d2ci_dalpha_dRe)
			Re_c = Re; alpha_c = alpha
			write(*,'(a,i3,a,e20.10,a,e20.10)') 'iter', iter, '  f1=', f1, '  f2=', f2
			if (abs(f1) < tol .and. abs(f2) < tol) exit

			J(1,1) = dci_dRe;         J(1,2) = f2
			J(2,1) = d2ci_dalpha_dRe; J(2,2) = d2ci_dalpha2
			rhs(1) = -f1; rhs(2) = -f2
			call solve2x2(J, rhs, dx)

			Re    = Re    + dx(1)
			alpha = alpha + dx(2)
			if (iter.eq.max_iter) then
				print*,'max newton iterations reached'
			endif
		end do
		! --- 3. Convergence Output ---
		print*,'nsteps',nsteps
		write(*,'(a,(1x, g0))') 'Re_c (Riccati)    = ', Re_c
		write(*,'(a,(1x, g0))') 'alpha_c (Riccati) = ', alpha_c
		write(*,'(a,(1x, g0))') 'c_r (Riccati)     = ', real(c,qp)
		
		write(*,'(a,e20.10)') 'Re_c difference             = ', Re_c_ref - Re_c
		write(*,'(a,e20.10)') 'alpha_c difference             = ', alpha_c_ref - alpha_c
		write(*,'(a,e20.10)') 'c_r difference             = ', real(c,qp) - c_r_ref
		

		write(str_cur, '(1x, g0)') Re_c
		write(*,'(i0, 3(1x, g0))') nsteps, Re_c, alpha_c, real(c,qp)
		do i = 1, len_trim(str_cur)
			if (str_cur(i:i) /= str_ref(i:i)) then
				first_diff = i
				exit
			end if
		end do
		! write latex row: common prefix normal, differing digits represented by dots
		write(10,'(a,i7,a,a,a,a,a)') ' ', nsteps, ' & $', trim(str_cur(1:first_diff-1)), '\ldots$ \\'

		call cpu_time(t2)
		print*,'runtime', t2-t1
	enddo
	
	write(10,'(a)') '\hline'
	write(10,'(a)') '\end{tabular}'
	write(10,'(a)') '\caption{Spectral convergence of critical point values with number of modes $N$}'
	write(10,'(a)') '\label{tab:convergence}'
	write(10,'(a)') '\end{table}'
	close(10)
end program riccati_finder