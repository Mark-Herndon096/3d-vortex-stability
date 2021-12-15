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
    INTEGER :: f_unit, io_stat
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
!======================================================================
END MODULE mod_io
!======================================================================
