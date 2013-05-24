C=======================================================================
C This contains implementations of several functions provided by
C Rogue Wave's IMSL (C) Numerical Library used in TIMCOAT, a
C research code developed at MIT.  Since IMSL is commercial, and
C TIMCOAT uses just a handful of the functions, it makes sense to
C use built-in functions (unless performance is demonstrably better
C otherwise).
C
C J. Roberts, 5/20/2013
C=======================================================================

C=======================================================================
C DLSARG solves the linear system Ax = b
C=======================================================================
      subroutine dlsarg(n, A, lda, b, ipath, x)

        implicit none
        ! input/output variables
        integer :: n,                    ! size of system
     &             lda                   ! unused
        double precision :: A(n, n),     ! coefficient matrix
     &                      b(n),        ! right hand side
     &                      ipath,       ! unused
     &                      x(n)         ! solution
        ! local variables
        double precision :: Acopy(n, n), ! copy to hold upper triangle
     &                      bcopy(n)     ! copy to overwrite

        ! copy A and B, since the original subroutine specifies that
        ! the input data is not overwritten
        Acopy = A
        bcopy = b

        ! solve the system
        call gausselim(n, Acopy, bcopy, x)

      end subroutine dlsarg

C=======================================================================
C GAUSSELIM solves the linear system Ax = b via gaussian elimination
C=======================================================================
      subroutine gausselim(n, A, b, x)    

        implicit none
        ! input/output variables
        integer :: n                     ! size of system
        double precision :: A(n, n),     ! coefficient matrix
     &                      b(n),        ! right hand side
     &                      x(n)         ! solution
        ! local variables
        double precision :: temp_row(n)  ! swapping row
        integer, dimension(1) :: i_max   ! store pivot index
        integer :: i, j, k               ! row/column indices
        double precision :: temp,        ! swapping value
     &                      m            ! row multiplier

        ! loop over all rows
        do k = 1, n-1         

           ! find row with largest value of |a(i, j)|, i = k, ..., n
           i_max = maxloc(abs(A(k:n, k)))
           i = i_max(1) + k - 1

           ! check whether largest magnitude in the column is near zero
           if (abs(a(i, k)) .le. 1E-16) then             
              stop "Zero pivot.  Solve fails!"   
           end if         

           ! swap rows i and k of A and the associated elements of b
           if (k .ne. i) then            
              temp_row = A(i, :)            
              A(i, :)  = A(k, :)             
              A(k, :)  = temp_row                  
              temp   = b(i)            
              b(i)   = b(k)            
              b(k)   = temp         
           end if         

           ! do all rows below the pivot.
           do i = k+1, n    
              if (abs(A(k, k)) .le. 1D-16) stop "zero on diagonal"
              m = A(i, k) / A(k, k)            
              A(i, k+1:n) = A(i, k+1:n) - A(k, k+1:n) * m        
              b(i) = b(i) - m * b(k)      
              ! lower triangle is filled with zeros
              A(i, k) = 0.0_8  
           end do

        end do   
       
        ! backsubstitute to get solution
        do i = n, 1, -1     
           temp = b(i)         
           do j = i+1, n            
              temp = temp - A(i, j) * b(j)         
           end do         
           b(i) = temp / A(i, i)      
        end do
        x = b
  
      end subroutine gausselim

C=======================================================================
C DRCURV fits data to a polynomial using least squares
C=======================================================================
      subroutine drcurv(nobs, xdata, ydata, ndeg, b, sspoly, stat)

        implicit none
        ! input/output
        integer :: nobs,                         ! number of points
     &             ndeg                          ! degree of polynomial
        double precision :: xdata(nobs),         ! x data
     &                      ydata(nobs),         ! y data (at x points)
     &                      b(ndeg+1),           ! coefficients
     &                      sspoly(NDEG+1),      ! sum of squares / deg
     &                      stat(10)             ! vector of statistics
        ! local data
        double precision :: A(nobs, NDEG+1),     ! matrix A
     &                      At(NDEG+1, nobs),    ! its transpose, A'
     &                      AtA(NDEG+1, NDEG+1), ! vandermonde, A'A
     &                      AtY(NDEG+1),         ! A' * y
     &                      yappx(nobs)          ! approximate y
        integer :: i, j

        ! Given
        !         [x1, x2, ..., xN]
        ! and
        !         [y1, y2, ..., yN
        ! we get a Mth degree polynomial fit via least squares
        ! minimization.
        ! 
        ! Let     |1, x1, x1^2, ... x1^k|   |c0|   |y1|
        !         |1, x2, ...           |   |c1|   |y2|
        !         |  ...                | x |..| = |..|
        !         |1, xN, ...           |   |cN|   |yN|
        ! We can't solve this in general, since k+1 != N.  We get
        ! around this by noting
        !         A'*A*c = A'*y
        ! which yields a square system of (k+1)*(k+1).  We obtain c via 
        !         c = inv(A'*A)*A'*y
        ! a simple linear solve.
        
        ! Fill matrices and vectors
        do i = 0, NDEG
          do j = 1, nobs
            A(j, i+1) = xdata(j)**i
          end do
        end do
        At  = transpose(A)
        AtA = matmul(At, A)
        AtY = matmul(At, ydata)
  
        ! Solve the least squares system
        b = 0.0
        call gausselim(ndeg+1, AtA, AtY, b)

        ! Approximate y data
        yappx = 0_8
        do i = 0, ndeg
          yappx = yappx + b(i+1)*ydata**i
          sspoly(i+1) = sum((yappx - sum(ydata)/nobs)**2)
        end do

        ! Statistics --- I'm assuming definitions for some of these
        stat(1)  = sum(xdata)/nobs                       ! x mean
        stat(2)  = sum(ydata)/nobs                       ! y mean
        stat(3)  = (sum(xdata**2)/nobs-stat(1))/(nobs-1) ! x sample var
        stat(4)  = (sum(ydata**2)/nobs-stat(2))/(nobs-1) ! y sample var
        stat(6)  = NDEG+1                                ! # dof's
        stat(7)  = sum((stat(2)-yappx)**2)   ! regression sum of squares
        stat(8)  = NDEG+1                    ! # dof's for error (???)
        stat(9)  = sum((ydata-yappx)**2)     ! error sum of squares
        stat(10) = 0_8                       ! # NaN's (???)
        stat(5)  = stat(7)/(stat(7)+stat(9)) ! R^2
        
      end subroutine drcurv

C=======================================================================
C DZPLRC drives Laguerre's method to find all roots of a polynomial
C
C Note, this implementation is just a slight modification of ZROOTS
C from Numerical Recipes in Fortran 77.
C=======================================================================
      subroutine dzplrc(m, a_in, roots)

        implicit none
        ! input/output
        integer :: m                  ! number of roots
        double precision :: a_in(m+1) ! polynomial coefficients
        double complex :: roots(m)    ! roots
        ! local
        integer :: MAXM
        double precision :: EPS
        parameter (EPS = 1.0D-12, MAXM = 101)
        double complex :: a(m+1)
        logical :: polish = .true.
        integer :: i, j, jj, its
        double complex :: ad(MAXM), x, b, c

        a = dcmplx(a_in)
        do j = 1, m+1
          ad(j) = a(j)
        end do
        do j = m, 1, -1
          x = dcmplx(0.0_8, 0.0_8)
          call laguer(ad, j, x, its)
          if(abs(aimag(x)) .le. 2.0_8 * EPS**2 * abs(dble(x))) then
            x = dcmplx(abs(x), 0.0_8)
          end if
          roots(j) = x
          b = ad(j+1)
          do jj = j, 1, -1
            c = ad(jj)
            ad(jj) = b
            b = x*b + c
          end do
        end do
        if (polish) then
          do j = 1, m
            call laguer(a, m, roots(j), its)
          end do
        end if
        do j = 2, m
          x = roots(j)
          do i = j-1, 1, -1
            if(dble(roots(i)) .le. dble(x)) exit
            roots(i+1) = roots(i)
          end do
          i = 0
          roots(i+1) = x
        end do

      end subroutine dzplrc

C=======================================================================
C LAGUERRE finds a single root of a polynomial given an initial guess
C
C Note, this implementation is just a slight modification of LAGUER
C from Numerical Recipes in Fortran 77.
C=======================================================================
      subroutine laguer(a, m, x, its)

        implicit none
        ! input/output
        integer :: m,              ! number of roots
     &             its             ! total iterations
        double complex :: a(m+1),  ! polynomial coefficients
     &                    x        ! roots
        ! local
        integer :: MAXIT,MR,MT
        double precision :: EPSS
        parameter (EPSS=1.0D-12, MR=8, MT=15, MAXIT = MT*MR)
        integer :: iter, j
        double precision :: abx,abp,abm,err,frac(MR)
        double complex :: dx,x1,b,d,f,g,h,sq,gp,gm,g2
        save frac
        data frac /.5_8,.25_8,.75_8,.13_8,.38_8,.62_8,.88_8,1.0_8/

        do iter = 1, MAXIT
          its = iter
          b = a(m+1)
          err = abs(b)
          d = dcmplx(0.0_8, 0.0_8)
          f = dcmplx(0.0_8, 0.0_8)
          abx = abs(x)
          do j = m, 1, -1
            f = x*f + d
            d = x*d + b
            b = x*b + a(j)
            err = abs(b) + abx*err
          end do
          err = EPSS*err
          if(abs(b) .le. err) then
            return
          else
            g = d/b
            g2 = g*g
            h = g2 - 2.0_8*f/b
            sq = sqrt((m-1)*(m*h-g2))
            gp = g + sq
            gm = g - sq
            abp = abs(gp)
            abm = abs(gm)
            if(abp .lt. abm) gp=gm
            if (max(abp, abm) .gt. 0.0_8) then
              dx = m / gp
            else
              dx = exp(dcmplx(log(1.0_8+abx), dble(iter)))
            end if
          end if
          x1 = x - dx
          if(x .eq. x1) return
          if (mod(iter, MT) .ne. 0) then
            x = x1
          else
            x = x - dx*frac(iter/MT)
          end if
        end do
        !print *, 'Warning: too many iterations in laguer.'

      end subroutine laguer
