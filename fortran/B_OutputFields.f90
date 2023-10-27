Subroutine B_OutputFields(IO, NI, NJ, X, Y, P, GradP, GradPError, V, DivV, DivVError)
    Real, Dimension(NI, NJ) :: X, Y
    Real, Dimension(0:NI, 0:NJ) :: P, DivV, DivVError
    Real, Dimension(0:NI, 0:NJ, 2) :: GradP, GradPError, V

    Write(IO, *) 'VARIABLES = "X", "Y", "P", "GradP_X", "GradP_Y", "GradPError_X", "GradPError_Y", "V_x", "V_y", "DivV", "DivVError"'
    Write(IO, *) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    Write(IO, '(100F14.7)') X(1:NI, 1:NJ)
    Write(IO, '(100F14.7)') Y(1:NI, 1:NJ)
    Write(IO, '(100F14.7)') P(1:NI - 1, 1:NJ - 1)
    Write(IO, '(100F25.17)') GradP(1:NI - 1, 1:NJ - 1, 1)
    Write(IO, '(100F25.17)') GradP(1:NI - 1, 1:NJ - 1, 2)
    Write(IO, '(100F25.17)') GradPError(1:NI - 1, 1:NJ - 1, 1)
    Write(IO, '(100F25.17)') GradPError(1:NI - 1, 1:NJ - 1, 2)
    Write(IO, '(100F25.17)') V(1:NI - 1, 1:NJ - 1, 1)
    Write(IO, '(100F25.17)') V(1:NI - 1, 1:NJ - 1, 2)
    Write(IO, '(100F25.17)') DivV(1:NI - 1, 1:NJ - 1)
    Write(IO, '(100F25.17)') DivVError(1:NI - 1, 1:NJ - 1)
End Subroutine 
