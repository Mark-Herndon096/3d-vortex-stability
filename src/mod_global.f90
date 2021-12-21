!======================================================================
!	Author: Mark Herndon
!	Date: 12/15/21
!	Description: Global variables and memory allocation routines
!======================================================================
MODULE mod_global
    IMPLICIT NONE
    ! CONSTANTS
    REAL(KIND=8), PARAMETER :: pi = 4.0*ATAN(1.0)

    ! USER-SPECIFIED PARAMETERS
    INTEGER      :: nt    !< # of time steps
    INTEGER      :: nv    !< # of vortices in real plane
    INTEGER      :: nk    !< # of wavenumbers
    INTEGER      :: nvt   !< Total # of vortices in ground-image system
    REAL(KIND=8) :: kspan !< Wavenumber span (kb)
    LOGICAL      :: ge    !< In Ground Effect Logical

    REAL(KIND=8) :: dtau  !< Time step
    REAL(KIND=8) :: b_0   !< Initial vortex separation

    ! USER-SPECIFIED INITIAL CONDITIONS
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: y_0    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: z_0    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: eta_0  
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: zeta_0 


    ! VORTEX POSITION AND PERTURBATION AMPLITUDE ARRAYS
    ! DIMENSION(nvt,nt) --> Ex. y(vortex index, time index) == position
    ! vortex (vortex index) at t = (time index)
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: y    
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: z    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: yz
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: eta  
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: zeta 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: yz_perturb 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: omega !< Self induced rotation frequency 

    ! VORTEX CIRCULATION STRENGTH AND ORIENTATION
    ! DIMENSION(nvt) --> Ex. gam(vortex index 1) ...
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: gam   
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: a
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: ka
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: kb 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: omega_array

    ! DERIVED PARAMETERS
    REAL(KIND=8) :: b     
    REAL(KIND=8) :: h     
    INTEGER      :: m     !< Dimension of VORT array
    INTEGER      :: k_ind, t_ind
    ! GLOBAL VARIABLES
    INTEGER :: n    !< Time integration indexing integer
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: vc_0
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: vc_new
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: vp_0
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: vp_new
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tau

    ! PROPAGATOR MATRIX
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)     :: eye
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)     :: s
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: phi
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)   :: v

    ABSTRACT INTERFACE
        FUNCTION mutual_induction(beta)
            REAL(KIND=8), INTENT(IN) :: beta
            REAl(KIND=8)             :: mutual_induction
        END FUNCTION
    END INTERFACE
    ABSTRACT INTERFACE
        FUNCTION self_induction(kappa)
            IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: kappa
            REAL(KIND=8)             :: self_induction
        END FUNCTION self_induction
    END INTERFACE

CONTAINS

!! SUBROUTINES AND FUNCTIONS
!======================================================================
SUBROUTINE allocate_Variables
    IMPLICIT NONE
    
    m = nv*2
    
    ALLOCATE(y_0(nv)) 
    ALLOCATE(z_0(nv)) 
    ALLOCATE(eta_0(nv)) 
    ALLOCATE(zeta_0(nv)) 

    !ALLOCATE(y(nv,nt))
    !ALLOCATE(z(nv,nt))
    ALLOCATE(yz(m,nt))
    ALLOCATE(yz_perturb(m,nt))
    ALLOCATE(eta(nv,nt))
    ALLOCATE(zeta(nv,nt))
    
    ALLOCATE(gam(nvt))
    ALLOCATE(a(nvt))
    ALLOCATE(omega(nvt,nk))

    ALLOCATE(vc_0(m))
    ALLOCATE(vc_new(m))
    ALLOCATE(vp_0(m))
    ALLOCATE(vp_new(m))
    ALLOCATE(tau(nt))
    
    ALLOCATE(phi(m,m,nt,nk))
    
    ! omega / wavenumber arrays
    ALLOCATE(omega_array(nk))
    ALLOCATE(ka(nk))
    ALLOCATE(kb(nk))
    ALLOCATE(s(nk,nt))
    ALLOCATE(v(m,nt,nk))
    ALLOCATE(eye(m,m))

END SUBROUTINE allocate_variables    
!======================================================================
!=======================================================================
SUBROUTINE set_ground_effect
    IMPLICIT NONE
    INTEGER :: j

    DO j = nv+1, nvt
        gam(j) = -gam(j-nv)
        a(j)   =  a(j-nv)
    END DO

END SUBROUTINE set_ground_effect
!=======================================================================
!=======================================================================
SUBROUTINE generate_kb_array
    IMPLICIT NONE
    INTEGER      :: k
    REAL(KIND=8) :: dk !< Wave number step size
    
    dk = kspan/REAL(nk,KIND=8)
    DO k = 1, nk
        kb(k) = REAL(k,KIND=8)*dk
    END DO

END SUBROUTINE generate_kb_array
!=======================================================================
!=======================================================================
SUBROUTINE set_init
    IMPLICIT NONE
    INTEGER :: i, j


    DO i = 1, nv
        yz(i,1)    = y_0(i)
        yz(i+nv,1) = z_0(i)
    END DO
    
    eye(:,:) = 0.d0

    DO j = 1, m
        eye(j,j) = 1.d0
    END DO
    
    DO j = 1, nk
        DO i = 1, m
            phi(:,i,1,j) = eye(:,i)
        END DO
    END DO
    
    tau(1) = 0.d0

    DO i = 1, nt-1
        tau(i+1) = dtau*(REAL(i,KIND=8))
    END DO    
    
END SUBROUTINE set_init
!=======================================================================
!=======================================================================
SUBROUTINE set_optimal_init
    IMPLICIT NONE
    INTEGER :: i, j

    t_ind = INT(4.0d0/0.001d0);
    !WRITE(*,*) 'Enter t index: '
    !READ(*,*) t_ind
    !WRITE(*,*) 'Enter k index: '
    !READ(*,*) k_ind
    k_ind = MAXLOC(s(:,t_ind),1)    
    WRITE(*,*) ' k = ', kb(k_ind), ' at k_ind = ', k_ind, 't_ind = ', t_ind

    DO i = 1, nv
        yz_perturb(i,1)    = V(i,t_ind,k_ind)
        yz_perturb(i+nv,1) = V(i+nv,t_ind,k_ind)
    END DO

    DO i = 1, 2*nv
        WRITE(*,*) yz_perturb(i,1)
    END DO
    

END SUBROUTINE set_optimal_init
!=======================================================================
!=======================================================================
END MODULE mod_global
!======================================================================
