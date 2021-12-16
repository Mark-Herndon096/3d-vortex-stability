!=================================================================================
! Written by Mark A. Herndon
! Lehigh University, Department of Mechanical Engineering and Mechanics
! Subroutines for 5th order accurate Runge-Kutta integration for systems of
! ordinary differential equations
! Thomas algorithm for solution to tri-diagonal matrix equation
!=================================================================================
MODULE mod_numerical_routines
    IMPLICIT NONE
    ABSTRACT INTERFACE
        FUNCTION DERIVATIVE(x_0, m, h, ch, z_opts, num_opts)
        IMPLICIT NONE
        INTEGER,                          INTENT(IN)    :: m
        INTEGER, OPTIONAL,                INTENT(INOUT) :: num_opts
        INTEGER, OPTIONAL, DIMENSION(:), INTENT(INOUT)  :: z_opts
        REAL(KIND=8),                     INTENT(IN)    :: h, ch
        REAL(KIND=8), DIMENSION(m),       INTENT(IN)    :: x_0
        REAL(KIND=8), DIMENSION(m)                      :: DERIVATIVE
        END FUNCTION
    END INTERFACE
    ABSTRACT INTERFACE
        FUNCTION root_function(x,z_opt)
            IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: x
            REAL(KIND=8), INTENT(IN), OPTIONAL :: z_opt
            REAL(KIND=8)             :: root_function
        END FUNCTION
    END INTERFACE  
    ! Butchers Tableau for Dorman-Prince Runge-Kutta method
    !_____________________________________________________
    ! c1 |
    ! c2 |  a21,                                           
    ! c3 |  a31, a32,                                      
    ! c4 |  a41, a42, a43,                                 
    ! c5 |  a51, a52, a53, a54,                            
    ! c6 |  61,  a62, a63, a64, a65,                       
    ! c7 |  a71, a72, a73, a74, a75, a76,                 
    !-----------------------------------------------------  
    !    |  b1,  b2,  b3,  b4,  b5,  b6,  b7,                
    !    |  b1s, b2s, b3s, b4s, b5s, b6s, b7s
    REAL(KIND=8), PARAMETER :: a21 =  1.d0/5.d0,          &
                               a31 =  3.d0/40.d0,         &
                               a32 =  9.d0/40.d0,         &
                               a41 =  44.d0/45.d0,        &
                               a42 = -56.d0/15.d0,        &
                               a43 =  32.d0/9.d0,         &
                               a51 =  19372.d0/6561.d0,   &
                               a52 = -25360.d0/2187.d0,   &
                               a53 =  64448.d0/6561.d0,   &
                               a54 = -212.d0/729.d0,      &
                               a61 =  9017.d0/3168.d0,    &
                               a62 = -355.d0/33.d0,       &
                               a63 =  46732.d0/5247.d0,   &
                               a64 =  49.d0/176.d0,       &
                               a65 = -5103.d0/18656.d0,   &
                               a71 =  35.d0/384.d0,       &
                               a72 =  0.d0,               &
                               a73 =  500.d0/1113.d0,     &
                               a74 =  125.d0/192.d0,      &
                               a75 = -2187.d0/6784.d0,    &
                               a76 =  11.d0/84.d0,        &
                               b1  =  35.d0/384.d0,       &
                               b2  =  0.d0,               &
                               b3  =  500.d0/1113.d0,     &
                               b4  =  125.d0/192.d0,      &
                               b5  = -2187.d0/6784.d0,    &
                               b6  =  11.d0/84.d0,        &
                               b7  =  0.d0,               &
                               b1s =  5179.d0/57600.d0,   &
                               b2s =  0.d0,               &
                               b3s =  7571.d0/16695.d0,   &
                               b4s =  393.d0/640.d0,      &
                               b5s = -92097.d0/339200.d0, &
                               b6s =  187.d0/2100.d0,     &
                               b7s =  1.d0/40.d0,         &
                               c1  =  0.d0,               &
                               c2  =  1.d0/5.d0,          &
                               c3  =  3.d0/10.d0,         &
                               c4  =  4.d0/5.d0,          &
                               c5  =  8.d0/9.d0,          &
                               c6  =  1.d0,               &
                               c7  =  1.d0 
CONTAINS
!======================================================================
! ode_dprk --> Dormand-Prince Runge-Kutta method (without adaptive step)
! INPUTS: 
!          nt    !< number of time steps
!          m     !< dimension of input vector
!          h     !< time step (dt)
!          DERIV !< derivative function from interface DERIVATIVE
! OUTPUTS: 
!          x  !< return solution vector x(1:m,1:nt)
! Place initial condition x_0 in x(:,1)
!======================================================================
SUBROUTINE ode_dprk(x,m,nt,h,DERIV, z_opts, num_opts)
    IMPLICIT NONE
    PROCEDURE(DERIVATIVE)                        :: DERIV
    INTEGER,      INTENT(IN)                     :: m, nt
    INTEGER,      INTENT(INOUT), OPTIONAL           :: num_opts
    INTEGER,      INTENT(INOUT), DIMENSION(:), OPTIONAL :: z_opts
    REAL(KIND=8), INTENT(INOUT)                  :: h
    REAL(KIND=8), INTENT(INOUT), DIMENSION(m,nt) :: x
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE      :: k1, & !< RK stage variables
                                                    k2, &
                                                    k3, &
                                                    k4, &
                                                    k5, &
                                                    k6, &
                                                    k7  
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE      :: y1, & !< RK stage variables
                                                    y2, &
                                                    y3, &
                                                    y4, &
                                                    y5, &
                                                    y6
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE      :: x_0, x_new
    INTEGER                                      :: n 
    REAL(KIND=8)                                 :: ch

    ALLOCATE(k1(m))
    ALLOCATE(k2(m))
    ALLOCATE(k3(m))
    ALLOCATE(k4(m))
    ALLOCATE(k5(m))
    ALLOCATE(k6(m))
    ALLOCATE(k7(m))
    ALLOCATE(y1(m))
    ALLOCATE(y2(m))
    ALLOCATE(y3(m))
    ALLOCATE(y4(m))
    ALLOCATE(y5(m))
    ALLOCATE(y6(m))
    ALLOCATE(x_0(m))
    ALLOCATE(x_new(m))
    
    x_0 = x(:,1);
    IF (PRESENT(z_opts)) THEN 
        ! Primary loop over each time series
        DO n = 1, nt-1
            z_opts(1) = n
            ch = c1*h
            k1 = DERIV(x_0, m, h, ch, z_opts, num_opts)
            y1 = x_0 + h*a21*k1
            ch = c2*h
            k2 = DERIV(y1, m, h, ch, z_opts, num_opts)
            y2 = x_0 + h*(a31*k1 + a32*k2)
            ch = c3*h
            k3 = DERIV(y2, m, h, ch, z_opts, num_opts)
            y3 = x_0 + h*(a41*k1 + a42*k2 + a43*k3)
            ch = c4*h
            k4 = DERIV(y3, m, h, ch, z_opts, num_opts)
            y4 = x_0 + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4)
            ch = c5*h
            k5 = DERIV(y4, m, h, ch, z_opts, num_opts)
            y5 = x_0 + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)
            ch = c6*h
            k6 = DERIV(y5, m, h, ch, z_opts, num_opts)
            y6 = x_0 + h*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 +a76*k6)
            ch = c7*h
            k7 = DERIV(y6, m, h, ch, z_opts, num_opts)
            
            x_new  = x_0 + h*(b1s*k1 + b2s*k2 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6 + b7s*k7)

            x(:,n+1) = x_new 
            x_0      = x_new 
        END DO
    ELSE
    ! Primary loop over each time series
        DO n = 1, nt-1
            ch = c1*h
            k1 = DERIV(x_0,m,h,ch)
            y1 = x_0 + h*a21*k1
            ch = c2*h
            k2 = DERIV(y1,m,h,ch)
            y2 = x_0 + h*(a31*k1 + a32*k2)
            ch = c3*h
            k3 = DERIV(y2,m,h,ch)
            y3 = x_0 + h*(a41*k1 + a42*k2 + a43*k3)
            ch = c4*h
            k4 = DERIV(y3,m,h,ch)
            y4 = x_0 + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4)
            ch = c5*h
            k5 = DERIV(y4,m,h,ch)
            y5 = x_0 + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)
            ch = c6*h
            k6 = DERIV(y5,m,h,ch)
            y6 = x_0 + h*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 +a76*k6)
            ch = c7*h
            k7 = DERIV(y6,m,h,ch)
            
            x_new  = x_0 + h*(b1s*k1 + b2s*k2 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6 + b7s*k7)

            x(:,n+1) = x_new 
            x_0      = x_new 
        END DO
    END IF

    DEALLOCATE(k1)
    DEALLOCATE(k2)
    DEALLOCATE(k3)
    DEALLOCATE(k4)
    DEALLOCATE(k5)
    DEALLOCATE(k6)
    DEALLOCATE(k7)
    DEALLOCATE(y1)
    DEALLOCATE(y2)
    DEALLOCATE(y3)
    DEALLOCATE(y4)
    DEALLOCATE(y5)
    DEALLOCATE(y6)
    DEALLOCATE(x_0)
    DEALLOCATE(x_new)
    
END SUBROUTINE ode_dprk
!=================================================================================
! RK5 is a 5th order accurate Runge-Kutta scheme based on the Dormand
! Prince method. Adaptive step sizes will be implemented in the future.
! Input require is initial position x_0 as a vector or order m and step size h
! DERIV function is required as an external function corresponding to the form
! of ABSTRACT INTERFACE DERIVATIVE
!
! Adaptize step size goal: take initial step size and calulate error. If solution
! is within ~8-10% of x_0 then grow step size. If error between O(h^4) and O(h^5)
! is beyond tolerance --> shrink step size and advance to time n+1 and return
! for larger step sizes than provided h, interpolate calculated points at nodes
!=================================================================================
SUBROUTINE RK5(x_0,x_new,h,m,DERIV)
    IMPLICIT NONE
    PROCEDURE(DERIVATIVE) :: DERIV   ! ABSTRACT INTERFACE DERIVATIVE => DERIV
    INTEGER :: n
    INTEGER, INTENT(INOUT) :: m      ! Order of equations
    REAL(KIND=8), INTENT(INOUT) :: h ! Initial step size
    REAL(KIND=8), DIMENSION(m), INTENT(INOUT)  :: x_0 ! Starting point
    REAL(KIND=8), DIMENSION(m), INTENT(OUT) :: x_new  ! value at end of integration
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_O4   ! 4th order accurate solution
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y1, y2, y3, y4, y5, y6
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: k1, k2, k3, k4, k5, k6, k7
    REAL(KIND=8) :: error, tolerance, ch  ! ch represents fractional step size
    ! Butchers Tableau for Dormand-Prince embedded RK5(4) method
    REAL(KIND=8) :: a21,                                            &
                    a31, a32,                                       &
                    a41, a42, a43,                                  &
                    a51, a52, a53, a54,                             &
                    a61, a62, a63, a64, a65,                        &
                    a71, a72, a73, a74, a75, a76,                   &
                    b1,  b2,  b3,  b4,  b5,  b6, b7,                &
                    b1_2,  b2_2,  b3_2,  b4_2,  b5_2,  b6_2, b7_2,  &
                    c1,  c2,  c3,  c4,  c5,  c6, c7

    c1  = 0.d0
    c2  = 1.d0/5.d0
    c3  = 3.d0/10.d0
    c4  = 4.d0/5.d0
    c5  = 8.d0/9.d0
    c6  = 1.d0
    c7  = 1.d0

    a21 = 1.d0/5.d0
    a31 = 3.d0/40.d0
    a32 = 9.d0/40.d0
    a41 = 44.d0/45.d0
    a42 = -56.d0/15.d0
    a43 = 32.d0/9.d0
    a51 = 19372.d0/6561.d0
    a52 = -25360.d0/2187.d0
    a53 = 64448.d0/6561.d0
    a54 = -212.d0/729.d0
    a61 = 9017.d0/3168.d0
    a62 = -355.d0/33.d0
    a63 = 46732.d0/5247.d0
    a64 = 49.d0/176.d0
    a65 = -5103.d0/18656.d0
    a71 = 35.d0/384.d0
    a72 = 0.d0
    a73 = 500.d0/1113.d0
    a74 = 125.d0/192.d0
    a75 = -2187/6784.d0
    a76 = 11.d0/84.d0

    b1 = a71
    b2 = a72
    b3 = a73
    b4 = a74
    b5 = a75
    b6 = a76
    b7 = 0.d0

    b1_2 = 5179.d0/57600.d0
    b2_2 = 0.d0
    b3_2 = 7571.d0/16695.d0
    b4_2 = 393.d0/640.d0
    b5_2 = -92097.d0/339200.d0
    b6_2 = 187.d0/2100.d0
    b7_2 = 1.d0/40.d0

    ALLOCATE(k1(m))
    ALLOCATE(k2(m))
    ALLOCATE(k3(m))
    ALLOCATE(k4(m))
    ALLOCATE(k5(m))
    ALLOCATE(k6(m))
    ALLOCATE(k7(m))
    ALLOCATE(y1(m))
    ALLOCATE(y2(m))
    ALLOCATE(y3(m))
    ALLOCATE(y4(m))
    ALLOCATE(y5(m))
    ALLOCATE(y6(m))
    ALLOCATE(x_O4(m))
          ch = c1*h
          k1 = DERIV(x_0,m,h,ch)
          y1 = x_0 + h*a21*k1
          ch = c2*h
          k2 = DERIV(y1,m,h,ch)
          y2 = x_0 + h*(a31*k1 + a32*k2)
          ch = c3*h
          k3 = DERIV(y2,m,h,ch)
          y3 = x_0 + h*(a41*k1 + a42*k2 + a43*k3)
          ch = c4*h
          k4 = DERIV(y3,m,h,ch)
          y4 = x_0 + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4)
          ch = c5*h
          k5 = DERIV(y4,m,h,ch)
          y5 = x_0 + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)
          ch = c6*h
          k6 = DERIV(y5,m,h,ch)
          y6 = x_0 + h*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 +a76*k6)
          ch = c7*h
          k7 = DERIV(y6,m,h,ch)
          ! 5th order accurate solution
    !      x_new = x_0 + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
          ! 4th order accurate solution
          x_new  = x_0 + h*(b1_2*k1 + b2_2*k2 + b3_2*k3 + b4_2*k4 + b5_2*k5 + b6_2*k6 + b7_2*k7)

          !error = MAXVAL(x_new - x_O4)
!          WRITE(*,*) 'error : ', error

    DEALLOCATE(k1)
    DEALLOCATE(k2)
    DEALLOCATE(k3)
    DEALLOCATE(k4)
    DEALLOCATE(k5)
    DEALLOCATE(k6)
    DEALLOCATE(k7)
    DEALLOCATE(y1)
    DEALLOCATE(y2)
    DEALLOCATE(y3)
    DEALLOCATE(y4)
    DEALLOCATE(y5)
    DEALLOCATE(y6)
    DEALLOCATE(x_O4)
END SUBROUTINE RK5
!=================================================================================
SUBROUTINE TSOLVE (IL, IU, BB, DD, AA, CC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IL, IU
    REAL(KIND=8), DIMENSION(IL:IU), INTENT(IN) :: AA, BB
    REAL(KIND=8), DIMENSION(IL:IU), INTENT(INOUT) :: CC, DD

    INTEGER :: LP, I, J
    REAL(KIND=8) :: R

    LP = IL + 1
    ! LU DECOMPOSITION AND FORWARD SUBSTITUTION FOR CC(I)
    DO I = LP, IU
        R = BB(I)/DD(I-1)
        DD(I) = DD(I) - R*AA(I-1)
        CC(I) = CC(I) - R*CC(I-1)
    END DO
    ! BACKWARD SUBSTITUTION
    CC(IU) = CC(IU)/DD(IU)
    DO I = LP, IU
        J = IU - I + IL
        CC(J) = (CC(J) - AA(J)*CC(J+1))/DD(J)
    END DO
    ! SOLUTION STORED IN CC
END SUBROUTINE TSOLVE
!=================================================================================
SUBROUTINE bisection_method(fun, x, a, b, tol, z_opt)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT)        :: x
    REAL(KIND=8), INTENT(IN)           :: a
    REAL(KIND=8), INTENT(IN)           :: b
    REAL(KIND=8), INTENT(IN)           :: tol
    REAL(KIND=8), INTENT(IN), OPTIONAL :: z_opt
    PROCEDURE(root_function)           :: fun

    INTEGER                  :: i, j
    REAL(KIND=8)             :: dx, mid, x1, x2, c, prod

    x1 = a; x2 = b;
    
    IF (PRESENT(z_opt)) THEN 
    
    DO i = 1, 60
        mid  = (x1 + x2)/2.d0;
        prod = fun(x1,z_opt)*fun(mid,z_opt);
        IF ( prod .LE. 0.d0 ) THEN
            x2  = mid
            mid = (x1 + x2)/2.d0
            dx  = ABS(x2 - x1);
        ELSE IF ( prod .GT. 0.d0 ) THEN
            x1  = mid;
            mid = (x1 + x2)/2.d0;
            dx  = ABS(x1 - x2);
        END IF
        
        IF ( dx .LE. tol ) THEN
            EXIT
        END IF
    END DO

    x = mid;

    ELSE

    DO i = 1, 60
        mid  = (x1 + x2)/2.d0;
        prod = fun(x1)*fun(mid);
        IF ( prod .LE. 0.d0 ) THEN
            x2  = mid
            mid = (x1 + x2)/2.d0
            dx  = ABS(x2 - x1);
        ELSE IF ( prod .GT. 0.d0 ) THEN
            x1  = mid;
            mid = (x1 + x2)/2.d0;
            dx  = ABS(x1 - x2);
        END IF
        
        IF ( dx .LE. tol ) THEN
            EXIT
        END IF
    END DO
    
    x = mid;
    END IF

        
END SUBROUTINE bisection_method
!=================================================================================
END MODULE mod_numerical_routines
