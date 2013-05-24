      program test_math_functions

        external gausselim
        external drcurv

        integer :: i, j, n, ndeg
        parameter (n = 10, ndeg = 2)
        double precision :: A(n, n), x(n), b(n), y(n), 
     &                      sspoly(10), stat(10), c(ndeg+1)
        double complex   :: roots(ndeg)


        !----------------------------------------------------------------------!
        ! linear solver                                                        !
        !----------------------------------------------------------------------! 
        !
        ! test case:
        !                 | 2 -1  0 ...
        !             A = |-1  2 -1 ...
        !                 | 0 -1  2 -1 ..
        !
        ! and         
        !             b = [1,1,...]'
        !
        ! solution:   x = [5,9,12,14,15,15,14,12,9,5]'
   
        A = 0.0
        x = 0.0
        b = 1.0
        do i = 1, n
          A(i, i) = 2.0
          if (i .gt. 1) then
            A(i, i-1) = -1.0
          end if
          if (i .lt. n) then
            A(i, i+1) = -1.0
          end if
        end do
    
        call gausselim(n, A, b, x)
        if (abs(x(1) - 5_8) > 1D-12) stop "linear solver failed"

        !----------------------------------------------------------------------!
        ! polynomial fit                                                       !
        !----------------------------------------------------------------------! 
        !
        ! ydata = [5,9,12,14,15,15,14,12,9,5].^2
        ! xdata = [1,2,...,10]
        ! solution: [-85.8, 110,-10]
        y = x**2
        do i = 1, n
          x(i) = dble(i)
        end do
 
        call drcurv(n, x, y, ndeg, c, sspoly, stat)
        if (abs(c(1)+85.8_8) > 1.0e-10) stop "Error in polynomials fit!"
        if (abs(c(2)-110_8)  > 1.0e-10) stop "Error in polynomials fit!"
        if (abs(c(3)+10_8)   > 1.0e-10) stop "Error in polynomials fit!"

        !----------------------------------------------------------------------!
        ! polynomial roots                                                     !
        !----------------------------------------------------------------------! 
        !
        ! p(x) = -85.8 + 110*x - 10*x^2
        call dzplrc(3, c, roots)
        if ((abs(dble(roots(1)) - 0.8448952750770239_8) > 1D-12)  .or.
     &      (abs(dble(roots(2)) - 10.155104724922966_8) > 1D-12)) then
          stop "root finder failed"
        end if

        print *, "tests successful"
      end program test_math_functions
