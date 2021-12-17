!======================================================================
!	Author: Mark Herndon
!	Date: 12/15/21
!	Description: Module for file read and writes
!======================================================================
MODULE mod_io
    IMPLICIT NONE

    !!TODO: VARIABLE AND INTERFACE DECLARATIONS

CONTAINS

!! SUBROUTINES AND FUNCTIONS
!======================================================================
SUBROUTINE read_input_data(input_fname)
    USE, INTRINSIC :: iso_fortran_env, ONLY : stderr => error_unit
    USE mod_global
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: input_fname
    INTEGER :: i, f_unit, io_stat
    LOGICAL :: ex_stat

    NAMELIST /CODE_DATA/ nt, dtau, nv, nvt, nk, kspan, ge
    NAMELIST /VORTEX_DATA/ y_0, z_0, eta_0, zeta_0, gam, a

    INQUIRE (FILE=TRIM(input_fname), EXIST=ex_stat)
    IF (ex_stat .EQ. 0) THEN
        WRITE(stderr, '(3a)') 'Error: file "', trim(input_fname), '" not found'
    END IF
    
    OPEN(FILE=TRIM(input_fname), newunit=f_unit, ACTION='READ',STATUS='OLD')
    READ(NML=CODE_DATA, UNIT=f_unit, iostat=io_stat) 
    IF (io_stat /= 0) THEN
        WRITE(stderr, '(3a)') 'Error reading namelist CODE_DATA in "', trim(input_fname),'"'
    END IF
    
    CALL allocate_variables
    
    READ(NML=VORTEX_DATA, UNIT=f_unit, iostat=io_stat) 
    IF (io_stat /= 0) THEN
        WRITE(stderr, '(3a)') 'Error reading namelist VORTEX_DATA in "', trim(input_fname),'"'
    END IF
    CLOSE(f_unit)
    

END SUBROUTINE read_input_data
!=======================================================================
!=======================================================================
SUBROUTINE WRITE_SOLUTION_FILE(ii)
    USE mod_global, ONLY : nv, nvt, nt, tau, ge, phi, nk, s, kb, V, yz
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    INTEGER :: zz
    CHARACTER(LEN=60) :: fname_1, fname_2
    zz = INT(yz(1+nv,1));
    IF ( GE == .FALSE. ) THEN
        WRITE(fname_2,'("DATA/perturbations-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        WRITE(1) phi
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) kb
        WRITE(1) V
        CLOSE(1)

        WRITE(fname_2,'("DATA/perturbations_2-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        !WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) kb
        WRITE(1) V
        CLOSE(1)
    ELSE
        WRITE(fname_2,'("DATA/perturbations-GE-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) kb
        WRITE(1) V
        CLOSE(1)

        WRITE(fname_2,'("DATA/perturbations_2-GE-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        !WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) kb
        WRITE(1) V
        CLOSE(1)
    END IF
END SUBROUTINE WRITE_SOLUTION_FILE
!=======================================================================
!=======================================================================
!======================================================================
!=======================================================================
SUBROUTINE write_trajectories(ii)
    USE mod_global, ONLY : nv, nt, yz, tau, ge
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    INTEGER :: zz
    CHARACTER(LEN=60) :: fname_1

    zz = INT(yz(nv+1,1));
    IF ( GE == .FALSE. ) THEN
        WRITE(fname_1,'("DATA/vortex_trajectories-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_1,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt
        WRITE(1) yz
        WRITE(1) tau
        CLOSE(1)
    ELSE
        WRITE(fname_1,'("DATA/vortex_trajectories-GE-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_1,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt
        WRITE(1) yz
        WRITE(1) tau
        CLOSE(1)
    END IF
    
END SUBROUTINE write_trajectories
!======================================================================
!=======================================================================
END MODULE mod_io
!======================================================================
