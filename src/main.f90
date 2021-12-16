!======================================================================
!	Author: Mark Herndon
!	Date: 12/15/21
!	Description: Main driver for vortex stability analysis
!======================================================================
PROGRAM main
    ! Add ONLY : statements to module USE statements
    USE mod_io
    USE mod_global
    USE mod_numerical_routines
    USE omp_lib
    USE LAPACK95, ONLY : GESVD
    IMPLICIT NONE
    PROCEDURE(DERIVATIVE) :: vortex_t_deriv
    PROCEDURE(DERIVATIVE) :: vortex_deriv
    PROCEDURE(self_induction) :: omega_func   !< Special functions
    INTEGER :: k, num_opts, jj
    INTEGER, DIMENSION(:), ALLOCATABLE :: opts
    CHARACTER(LEN=100) :: input_fname    
    REAL(KIND=8) :: test_val, xx, tstart, tf, tend

    CALL GET_COMMAND_ARGUMENT(1,input_fname)
    CALL read_input_data(TRIM(input_fname))
    IF ( GE == .TRUE. ) THEN
        CALL SET_GROUND_EFFECT
    END IF
    
    tau(1) = 0.d0
    DO n = 1, nt-1
        tau(n+1) = dtau*(REAL(n,KIND=8))
    END DO    

    CALL generate_kb_array
    CALL set_init
    CALL ode_dprk(yz, m, nt, dtau, vortex_t_deriv)
    CALL write_trajectories(nk)
    
    num_opts = 2
    ALLOCATE(opts(num_opts))
    
    DO k = 1, nk
        DO jj = 1, nv
            omega(jj,k) = omega_func(kb(k)*a(jj))
        END DO
    END DO


    tstart = OMP_get_wtime()
    DO k = 1, nk
        opts(1) = k
        DO jj = 1, m
            CALL ode_dprk(phi(:,jj,:,k), m, nt, dtau, vortex_deriv, opts, num_opts)
        END DO
        tend = OMP_get_wtime()
        tf   = tend - tstart
        WRITE(*,'(A,I4,X,A,X,I4,3X,A,3X,F12.6,3X,A)') 'COMPLETED ITERATION ', k, '/',nk,'. . .', tf, 'SECONDS ELAPSED'
    END DO    
     
    !!$OMP PARALLEL
    !!$OMP DO    
    !DO k = 1, nk
    !    omega_array(k) = self_induction(ka(k))
    !END DO
    !!$OMP END DO
    !!$OMP END PARALLEL
    !
    !OPEN(1,FILE='omega.x',ACTION='WRITE',STATUS='REPLACE',ACCESS='STREAM',FORM='UNFORMATTED')
    !WRITE(1) nk
    !WRITE(1) ka
    !WRITE(1) omega_array
    !CLOSE(1)
    
END PROGRAM main
!=======================================================================
!=======================================================================
FUNCTION vortex_deriv(x_0, m, h, ch, z_opts, num_opts)
    USE mod_global, ONLY : ge, pi, nv, nvt, gam, mutual_induction, ka, &
                           kb, a, yz, omega
    IMPLICIT NONE
    INTEGER,                    INTENT(IN)    :: m
    INTEGER, OPTIONAL,          INTENT(IN)    :: num_opts
    INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN) :: z_opts
    REAL(KIND=8),               INTENT(IN)    :: h, ch
    REAL(KIND=8), DIMENSION(m), INTENT(IN)    :: x_0
    REAL(KIND=8), DIMENSION(m)                :: vortex_deriv
    ! FUNCTION SPECIFIC VARIABLES AND PARAMETERS
    PROCEDURE(mutual_induction)               :: psi, phi    !< Special functions
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_temp      !< Temporary y array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_temp      !< Temporary z array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: eta_temp    !< Temporary eta array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: zeta_temp   !< Temporary zeta array

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_deriv     !< y derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_deriv     !< z derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: eta_deriv   !< eta derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: zeta_deriv  !< zeta derivative array
    INTEGER                                   :: i, j        !< Loop index integers
    INTEGER                                   :: n, k
    REAL(KIND=8)                              :: y_mn, z_mn  !< Relative y and z coordinates
    REAL(KIND=8)                              :: r_mn        !< Relative radius
    REAL(KIND=8)                              :: sum_y       !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_z       !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_eta     !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_zeta    !< Sum holder for y eq
    REAL(KIND=8)                              :: V1_mn       !< First term for eta equation
    REAL(KIND=8)                              :: V2_mn       !< Second term for eta equation
    REAL(KIND=8)                              :: V3_mn       !< Third term for eta equation
    REAL(KIND=8)                              :: V4_mn       !< Fourth term for eta equation
    REAL(KIND=8)                              :: W1_mn       !< First term for zeta equation
    REAL(KIND=8)                              :: W2_mn       !< Second term for zeta equation
    REAL(KIND=8)                              :: W3_mn       !< Third term for zeta equation
    REAL(KIND=8)                              :: W4_mn       !< Fourth term for zeta equation

    ALLOCATE(eta_temp(nvt))
    ALLOCATE(zeta_temp(nvt))
    ALLOCATE(y_temp(nvt))
    ALLOCATE(z_temp(nvt))
    ALLOCATE(eta_deriv(nv))
    ALLOCATE(zeta_deriv(nv))
    
    n = z_opts(1); k = z_opts(2);

    ! Initialize sum holder values to 0.0
    sum_eta = 0.0; sum_zeta = 0.0;

    ! Map x_0 values into temporary arrays
    DO i = 1, nv
        y_temp(i)    = yz(i,n)    + ch*(yz(i,n+1)-yz(i,n))/h 
        z_temp(i)    = yz(i+nv,n) + ch*(yz(i+nv,n+1)-yz(i+nv,n))/h 
        eta_temp(i)  = x_0(i)
        zeta_temp(i) = x_0(i+nv)
    END DO

    IF ( GE == .TRUE. ) THEN
     DO i = nv+1, nvt
         y_temp(i)     =  y_temp(i-nv)
         z_temp(i)     = -z_temp(i-nv)
         eta_temp(i)   =  x_0(i-nv)
         zeta_temp(i)  = -x_0(i)
         omega(i,k)      = omega(i-nv,k)
     END DO
    END IF

    DO i = 1, nv
        DO j = 1, nvt
            IF ( i .EQ. j ) THEN
                CYCLE
            ELSEIF ( i .NE. j ) THEN
                z_mn  = z_temp(j) - z_temp(i)
                y_mn  = y_temp(j) - y_temp(i)
                r_mn  = SQRT(y_mn**2 + z_mn**2)

                V1_mn =   gam(j)*2.d0*y_mn*z_mn/(r_mn**4)
                V2_mn = -(gam(j)*2.d0*y_mn*z_mn/(r_mn**4))*PHI(kb(k)*r_mn) 
                V3_mn = -(gam(j)/(r_mn**2))*(1.d0 - (2.d0*z_mn**2/r_mn**2))
                V4_mn =  (gam(j)/(r_mn**2))*(PSI(kb(k)*r_mn) - ((2.d0*z_mn**2)/r_mn**2)*PHI(kb(k)*r_mn))

                W1_mn =  -gam(j)*2.d0*y_mn*z_mn/(r_mn**4)
                W2_mn =  (gam(j)*2.d0*y_mn*z_mn/(r_mn**4))*PHI(kb(k)*r_mn) 
                W3_mn =  (gam(j)/(r_mn**2))*(1.d0 - (2.d0*y_mn**2/r_mn**2))
                W4_mn = -(gam(j)/(r_mn**2))*(PSI(kb(k)*r_mn) - ((2.d0*y_mn**2)/r_mn**2)*PHI(kb(k)*r_mn))

                sum_eta  = sum_eta  + V1_mn*eta_temp(i)  + V2_mn*eta_temp(j)  + V3_mn*zeta_temp(i) + V4_mn*zeta_temp(j)
                sum_zeta = sum_zeta + W1_mn*zeta_temp(i) + W2_mn*zeta_temp(j) + W3_mn*eta_temp(i)  + W4_mn*eta_temp(j)
            END IF
        END DO
        eta_deriv(i)  = sum_eta  + (gam(i)/(a(i)**2))*omega(i,k)*zeta_temp(i) 
        zeta_deriv(i) = sum_zeta - (gam(i)/(a(i)**2))*omega(i,k)*eta_temp(i) 
        sum_eta       = 0.d0
        sum_zeta      = 0.d0
    END DO

    DO i = 1, nv
        vortex_deriv(i) = eta_deriv(i)
        vortex_deriv(i+nv) = zeta_deriv(i)
    END DO
    

    DEALLOCATE(y_temp)
    DEALLOCATE(z_temp)
    DEALLOCATE(eta_temp)
    DEALLOCATE(zeta_temp)
    DEALLOCATE(eta_deriv)
    DEALLOCATE(zeta_deriv)
END FUNCTION vortex_deriv
!=======================================================================
!=======================================================================
FUNCTION vortex_t_deriv(x_0, m, h, ch)
    USE mod_global, ONLY : ge, pi, nv, nvt, gam, omega
    IMPLICIT NONE
    INTEGER,                    INTENT(IN)    :: m
    REAL(KIND=8),               INTENT(IN)    :: h, ch
    REAL(KIND=8), DIMENSION(m), INTENT(IN)    :: x_0
    REAL(KIND=8), DIMENSION(m)                :: vortex_t_deriv
    ! FUNCTION SPECIFIC VARIABLES AND PARAMETERS
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_temp   !< Temporary y array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_temp   !< Temporary z array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_deriv  !< y derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_deriv  !< z derivative array
    INTEGER                                   :: i, j     !< Loop index integers
    REAL(KIND=8)                              :: y_mn     !< Relative y coordinates
    REAL(KIND=8)                              :: z_mn     !< Relative z coordinates
    REAL(KIND=8)                              :: r_mn     !< Relative radius
    REAL(KIND=8)                              :: sum_y    !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_z    !< Sum holder for y eq

    ALLOCATE(y_temp(nvt))
    ALLOCATE(z_temp(nvt))

    ALLOCATE(y_deriv(nv))
    ALLOCATE(z_deriv(nv))

    ! Initialize sum holder values to 0.0
    sum_y = 0.0; sum_z = 0.0;

    ! Map x_0 values into temporary arrays
    DO i = 1, nv
        y_temp(i)    = x_0(i)
        z_temp(i)    = x_0(i+nv)
    END DO
    IF ( GE == .TRUE. ) THEN
     DO i = nv+1, nvt
         y_temp(i)    =  x_0(i-nv)
         z_temp(i)    = -x_0(i)
     END DO
    END IF

    DO i = 1, nv
        DO j = 1, nvt
            IF ( i .EQ. j ) THEN
                CYCLE
            ELSEIF ( i .NE. j ) THEN
                z_mn  = z_temp(j) - z_temp(i)
                y_mn  = y_temp(j) - y_temp(i)
                r_mn  = SQRT(y_mn**2 + z_mn**2)
                
                sum_y = sum_y + GAM(j)*z_mn/r_mn**2
                sum_z = sum_z - GAM(j)*y_mn/r_mn**2
            END IF
        END DO
        y_deriv(i)    = sum_y
        z_deriv(i)    = sum_z
        sum_y         = 0.d0
        sum_z         = 0.d0
    END DO
    DO i = 1, nv
        vortex_t_deriv(i)    = y_deriv(i)
        vortex_t_deriv(i+nv) = z_deriv(i)
    END DO
    
    DEALLOCATE(y_temp)
    DEALLOCATE(z_temp)
    DEALLOCATE(y_deriv)
    DEALLOCATE(z_deriv)

END FUNCTION vortex_t_deriv
!=======================================================================
!======================================================================
FUNCTION omega_func(kappa)
    USE mod_numerical_routines, ONLY : bisection_method, root_function
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: kappa
    PROCEDURE(root_function) :: bessel_root
    PROCEDURE(root_function) :: dispersion
    REAL(KIND=8) :: x !< return value for bisection root finding method
    REAL(KIND=8) :: a !< Left intial end point for bisection method
    REAL(KIND=8) :: b !< Right intial end point for bisection method
    REAL(KIND=8) :: tol = 1E-16 !< Tolerance for root approximation
    REAL(KIND=8) :: eps = 1E-16; !< Adjust bessl_root_val by eps
    REAL(KIND=8) :: bessel_root_val !< clearly indicate root of bessel
    REAL(KIND=8) :: omega_func
    REAL(KIND=8) :: tiny_eps
    
    tiny_eps = eps/2.d0;
    a = tiny_eps; b = 5.d0; ! First root of bessel J in this interval
    CALL bisection_method(bessel_root, x, a, b, tol)
    bessel_root_val = x;
    b = bessel_root_val - eps;
    CALL bisection_method(dispersion, x, a, b, tol, kappa)
    omega_func = ((2*kappa/SQRT(kappa**2 + x**2)) - 1.d0);
END FUNCTION omega_func
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
FUNCTION psi(beta)
    USE special_function_interface, ONLY : BESSELK0, BESSELK1
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: beta
    REAL(KIND=8)             :: psi
    
    IF ( beta == 0.d0 ) THEN
        psi = 1.d0
    ELSE IF ( beta .GE. 600 ) THEN
        psi = 0.d0
    ELSE
        psi = (beta**2)*BESSELK0(ABS(beta)) + & 
              ABS(beta)*BESSELK1(ABS(beta))
    END IF
END FUNCTION psi
!=======================================================================
!=======================================================================
FUNCTION phi(beta)
    USE special_function_interface, ONLY : BESSELKN
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: beta
    REAL(KIND=8)             :: phi

    IF ( BETA == 0.d0 ) THEN
        phi = 1.d0
    ELSE IF ( BETA .GE. 600.d0 ) THEN
        phi = 0.d0
    ELSE 
        phi = (1.d0/2.d0)*(beta**2)*BESSELKN(2,ABS(beta))
    END IF
END FUNCTION phi
!=======================================================================
!=======================================================================
