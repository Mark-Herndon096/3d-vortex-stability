!======================================================================
!	Author: Mark Herndon
!	Date: 12/15/21
!	Description: Main driver for vortex stability analysis
!======================================================================
PROGRAM main
    ! Add ONLY : statements to module USE statements
    USE mod_io
    USE mod_global
    USE omp_lib
    USE LAPACK95, ONLY : GESVD
    IMPLICIT NONE
    INTERFACE
        FUNCTION self_induction(kappa)
            IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: kappa
            REAL(KIND=8)             :: self_induction
        END FUNCTION self_induction
    END INTERFACE
    INTEGER :: k
    CHARACTER(LEN=100) :: input_fname    
    REAL(KIND=8) :: test_val

    CALL GET_COMMAND_ARGUMENT(1,input_fname)
    CALL read_input_data(TRIM(input_fname))
    IF ( GE == .TRUE. ) THEN
        CALL SET_GROUND_EFFECT
    END IF
    
    CALL generate_ka_array
    
    !$OMP PARALLEL
    !$OMP DO    
    DO k = 1, nk
        omega_array(k) = self_induction(ka(k))
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
    OPEN(1,FILE='omega.x',ACTION='WRITE',STATUS='REPLACE',ACCESS='STREAM',FORM='UNFORMATTED')
    WRITE(1) nk
    WRITE(1) ka
    WRITE(1) omega_array
    CLOSE(1)
    
END PROGRAM main
!======================================================================
!======================================================================
FUNCTION self_induction(kappa)
    USE mod_numerical_routines, ONLY : bisection_method, root_function
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: kappa
    PROCEDURE(root_function) :: bessel_root
    PROCEDURE(root_function) :: dispersion
    REAL(KIND=8) :: x !< return value for bisection root finding method
    REAL(KIND=8) :: a !< Left intial end point for bisection method
    REAL(KIND=8) :: b !< Right intial end point for bisection method
    REAL(KIND=8) :: tol = 1E-16 !< Tolerance for root approximation
    REAL(KIND=8) :: eps = 1E-10; !< Adjust bessl_root_val by eps
    REAL(KIND=8) :: bessel_root_val !< clearly indicate root of bessel
    REAL(KIND=8) :: self_induction
    
    a = 0.01d0; b = 5.d0; ! First root of bessel J in this interval
    CALL bisection_method(bessel_root, x, a, b, tol)
    bessel_root_val = x;
    b = bessel_root_val - eps;
    CALL bisection_method(dispersion, x, a, b, tol, kappa)
    self_induction = ((2*kappa/SQRT(kappa**2 + x**2)) - 1.d0);
END FUNCTION self_induction
!======================================================================
!======================================================================
FUNCTION dispersion(beta, kappa)
    USE special_function_interface, ONLY : BESSELJ0, BESSELJ1, BESSELJN, &
                                           BESSELK0, BESSELK1, BESSELKN
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: beta, kappa
    REAL(KIND=8)             :: dispersion
    REAL(KIND=8)             :: J0, J1, J2, J1_p 
    REAL(KIND=8)             :: K0, K1, K2, K1_p 

    J0 = BESSELJ0(beta);  J1 = BESSELJ1(beta);  J2 = BESSELJN(2,beta);
    K0 = BESSELK0(kappa); K1 = BESSELK1(kappa); K2 = BESSELKN(2,kappa);

    J1_p =  (J0 - J2)/2.d0;
    K1_p = -(K0 + K2)/2.d0;

    dispersion = (1.d0/beta)*(J1_p/J1) + K1_p/(kappa*K1) + &
                 SQRT(beta**2 + (kappa)**2)/(kappa*beta**2) 
END FUNCTION dispersion
!=======================================================================
!=======================================================================
FUNCTION bessel_root(x)
    USE special_function_interface, ONLY : BESSELJ1
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: x
    REAL(KIND=8)             :: bessel_root
    REAL(KIND=8)             :: y
    
    y = BESSELJ1(x)
    bessel_root = y

END FUNCTION bessel_root
!=======================================================================
!=======================================================================
