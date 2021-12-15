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
    INTEGER :: k
    CHARACTER(LEN=100) :: input_fname    

    CALL GET_COMMAND_ARGUMENT(1,input_fname)
    CALL read_input_data(TRIM(input_fname))
    IF ( GE == .TRUE. ) THEN
        CALL SET_GROUND_EFFECT
    END IF
    
    CALL generate_ka_array

END PROGRAM main
!======================================================================
