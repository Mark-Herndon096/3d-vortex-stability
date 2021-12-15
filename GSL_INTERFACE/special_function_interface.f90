!================================================================
!   FORTRAN interface module for C wrapper to GNU GSL special 
!   functions -- currently include Bessel functions
!   
!   Created by: Mark A. Herndon
!   Lehigh University, Department of Mechanical Engineering
!   and Mechanics 
!================================================================
MODULE special_function_interface
    IMPLICIT NONE

    INTERFACE BESSELJ0
        FUNCTION bessel_j0_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_j0_wrapper
        END FUNCTION 
    END INTERFACE

    INTERFACE BESSELJ1
        FUNCTION bessel_j1_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_j1_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELJN
        FUNCTION bessel_jn_wrapper(n, x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE)             :: bessel_jn_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELY0
        FUNCTION bessel_y0_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_y0_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELY1
        FUNCTION bessel_y1_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_y1_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELYN
        FUNCTION bessel_yn_wrapper(n, x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE)             :: bessel_yn_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELI0
        FUNCTION bessel_i0_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_i0_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELI1
        FUNCTION bessel_i1_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_i1_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELIN
        FUNCTION bessel_in_wrapper(n, x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE)             :: bessel_in_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELK0
        FUNCTION bessel_k0_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_k0_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELK1
        FUNCTION bessel_k1_wrapper(x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_k1_wrapper
        END FUNCTION
    END INTERFACE

    INTERFACE BESSELKN
        FUNCTION bessel_kn_wrapper(n, x) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(C_DOUBLE), INTENT(IN) :: x
            REAL(C_DOUBLE) :: bessel_kn_wrapper
        END FUNCTION
    END INTERFACE

END MODULE special_function_interface


