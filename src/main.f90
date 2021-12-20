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
    IF ( ge == .TRUE. ) THEN
        CALL set_ground_effect
    END IF

    CALL generate_kb_array
    CALL set_init
    CALL ode_dprk(yz, m, nt, dtau, vortex_t_deriv)
    CALL write_trajectories(nk)
    CALL calculate_omega 
    num_opts = 2
    ALLOCATE(opts(num_opts))


    tstart = OMP_get_wtime()
    DO k = 1, nk
        opts(2) = k
        DO jj = 1, m
            CALL ode_dprk(phi(:,jj,:,k), m, nt, dtau, vortex_deriv, opts, num_opts)
        END DO
            CALL calculate_singular_values(k)
        tend = OMP_get_wtime()
        tf   = tend - tstart
        WRITE(*,'(A,I4,X,A,X,I4,3X,A,3X,F12.6,3X,A)') 'COMPLETED ITERATION ', k, '/',nk,'. . .', tf, 'SECONDS ELAPSED'
    END DO
    
    CALL set_optimal_init
    opts(2) = 244
        
    WRITE(*,*) 'yz_perturb(1,1) = ', yz_perturb(1,1)
    WRITE(*,*) 'yz_perturb(2,1) = ', yz_perturb(2,1)
    WRITE(*,*) 'yz_perturb(3,1) = ', yz_perturb(3,1)
    WRITE(*,*) 'yz_perturb(4,1) = ', yz_perturb(4,1)

    CALL ode_dprk(yz_perturb, m, nt, dtau, vortex_deriv, opts, num_opts)
    CALL write_solution_file(nk)
     
    OPEN(1,FILE='omega.x',ACTION='WRITE',STATUS='REPLACE',ACCESS='STREAM',FORM='UNFORMATTED')
    WRITE(1) nv, nk
    WRITE(1) kb
    WRITE(1) omega
    CLOSE(1)
    
END PROGRAM main
!=======================================================================
!=======================================================================
FUNCTION vortex_deriv(x_0, m, h, ch, z_opts, num_opts)
    USE mod_global, ONLY : ge, pi, nv, nvt, gam, mutual_induction, ka, &
                           kb, a, yz, omega, nt
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
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: omega_temp   !< Temporary zeta array

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
    REAL(KIND=8)                              :: v1_mn       !< First term for eta equation
    REAL(KIND=8)                              :: v2_mn       !< Second term for eta equation
    REAL(KIND=8)                              :: v3_mn       !< Third term for eta equation
    REAL(KIND=8)                              :: v4_mn       !< Fourth term for eta equation
    REAL(KIND=8)                              :: w1_mn       !< First term for zeta equation
    REAL(KIND=8)                              :: w2_mn       !< Second term for zeta equation
    REAL(KIND=8)                              :: w3_mn       !< Third term for zeta equation
    REAL(KIND=8)                              :: w4_mn       !< Fourth term for zeta equation

    ALLOCATE(eta_temp(nvt))
    ALLOCATE(zeta_temp(nvt))
    ALLOCATE(y_temp(nvt))
    ALLOCATE(z_temp(nvt))
    ALLOCATE(omega_temp(nvt))
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
    
    !IF ( n .EQ. 1 ) THEN
    !    WRITE(*,*) '=============================================='
    !    WRITE(*,*) 'x_0(1) = ', x_0(1)
    !    WRITE(*,*) 'x_0(2) = ', x_0(2)
    !    WRITE(*,*) 'x_0(3) = ', x_0(3)
    !    WRITE(*,*) 'x_0(4) = ', x_0(4)
    !    WRITE(*,*) '============================================='
    !END IF
   
    IF ( GE == .TRUE. ) THEN
     DO i = nv+1, nvt
         y_temp(i)     =  y_temp(i-nv)
         z_temp(i)     = -z_temp(i-nv)
         eta_temp(i)   =  x_0(i-nv)
         zeta_temp(i)  = -x_0(i)
         !omega(i,k)    =  omega(i-nv,k)
         omega_temp(i)    =  omega(i-nv,k)
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

                v1_mn =   gam(j)*2.d0*y_mn*z_mn/(r_mn**4)
                v2_mn = -(gam(j)*2.d0*y_mn*z_mn/(r_mn**4))*phi(kb(k)*r_mn) 
                v3_mn = -(gam(j)/(r_mn**2))*(1.d0 - (2.d0*z_mn**2/r_mn**2))
                v4_mn =  (gam(j)/(r_mn**2))*(psi(kb(k)*r_mn) - ((2.d0*z_mn**2)/r_mn**2)*phi(kb(k)*r_mn))

                w1_mn =  -gam(j)*2.d0*y_mn*z_mn/(r_mn**4)
                w2_mn =  (gam(j)*2.d0*y_mn*z_mn/(r_mn**4))*phi(kb(k)*r_mn) 
                w3_mn =  (gam(j)/(r_mn**2))*(1.d0 - (2.d0*y_mn**2/r_mn**2))
                w4_mn = -(gam(j)/(r_mn**2))*(psi(kb(k)*r_mn) - ((2.d0*y_mn**2)/r_mn**2)*phi(kb(k)*r_mn))

                sum_eta  = sum_eta  + v1_mn*eta_temp(i)  + v2_mn*eta_temp(j)  + v3_mn*zeta_temp(i) + v4_mn*zeta_temp(j)
                sum_zeta = sum_zeta + w1_mn*zeta_temp(i) + w2_mn*zeta_temp(j) + w3_mn*eta_temp(i)  + w4_mn*eta_temp(j)
            END IF
        END DO
        eta_deriv(i)  = sum_eta  + (gam(i)/(a(i)**2))*omega_temp(i)*zeta_temp(i) 
        zeta_deriv(i) = sum_zeta - (gam(i)/(a(i)**2))*omega_temp(i)*eta_temp(i) 
        sum_eta       = 0.d0
        sum_zeta      = 0.d0
    END DO

    DO i = 1, nv
        vortex_deriv(i) = eta_deriv(i)
        vortex_deriv(i+nv) = zeta_deriv(i)
    END DO
    
    !WRITE(*,*) '--------------------------------------'
    !WRITE(*,*) 'k          = ', k
    !WRITE(*,*) 'kb(k)      = ', kb(k) 
    !WRITE(*,*) 'a(1)       = ', a(1) 
    !WRITE(*,*) 'a(2)       = ', a(2) 
    !WRITE(*,*) 'gam(1)     = ', gam(1) 
    !WRITE(*,*) 'gam(2)     = ', gam(2) 
    !WRITE(*,*) 'omega(1,k) = ', omega(1,k) 
    !WRITE(*,*) 'omega(2,k) = ', omega(2,k) 
    !WRITE(*,*) '--------------------------------------'
    !WRITE(*,*) '--------------------------------------'

    DEALLOCATE(y_temp)
    DEALLOCATE(z_temp)
    DEALLOCATE(eta_temp)
    DEALLOCATE(zeta_temp)
    DEALLOCATE(omega_temp)
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
SUBROUTINE calculate_singular_values(ii)
    USE mod_global, ONLY : m, phi, s, nt, V
    USE LAPACK95, ONLY : GESVD
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    INTEGER :: i, j, ni, nj, nn
    INTEGER(KIND=8) :: INFO
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: a_tmp, u, vt
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: s_array
    CHARACTER(LEN=1) :: JOB = 'U'
    ALLOCATE(a_tmp(m,m))
    ALLOCATE(s_array(m))
    ALLOCATE(vt(m,m))
    
    DO nn = 1, nt
        a_tmp(:,:) = phi(:,:,nn,ii)
!        WRITE(*,'("Matrix A"/(<m>F8.4))'), ((A_tmp(i,j), i = 1, m), j = 1, m)
        CALL GESVD(A=A_tmp,S=s_array,VT=VT,JOB=JOB,INFO=INFO)
    !    WRITE(*,'("Matrix V"/(<m>F8.4))'), ((VT(i,j), i = 1, m), j = 1, m)
!        WRITE(*,'("Singular Values"/(1F8.4))'), (s_array(i), i = 1, m)
        s(ii,nn) = s_array(1)
        DO j = 1, m
            V(j,nn,ii) = VT(1,j)
        !    WRITE(*,*) VT(j,1)
        END DO
    END DO
    
    DEALLOCATE(A_tmp); DEALLOCATE(s_array); DEALLOCATE(VT)
END SUBROUTINE calculate_singular_values
!======================================================================
!======================================================================
SUBROUTINE calculate_omega 
    USE mod_global, ONLY : omega, nk, nv, kb, a, self_induction
    IMPLICIT NONE
    PROCEDURE(self_induction) :: omega_func   !< Special functions
    INTEGER :: k, jj    
    DO k = 1, nk
        DO jj = 1, nv
            omega(jj,k) = omega_func(kb(k)*a(jj))
        END DO
    END DO

END SUBROUTINE calculate_omega 
!======================================================================
!======================================================================
