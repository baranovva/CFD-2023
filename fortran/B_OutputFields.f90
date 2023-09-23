Subroutine B_OutputFields(IO, NI, NJ, X, Y, P, GradP)
    Real, Dimension(NI, NJ) :: X, Y
    Real, Dimension(0:NI, 0:NJ) :: P
    Real, Dimension(0:NI, 0:NJ, 2) :: GradP

    Write(IO, *) 'VARIABLES = "X", "Y", "P", "GradPX", "GradPY"'
    Write(IO, *) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    Write(IO, '(100F14.7)') X(1:NI, 1:NJ)
    Write(IO, '(100F14.7)') Y(1:NI, 1:NJ)
    Write(IO, '(100F14.7)') P(1:NI - 1, 1:NJ - 1)
    Write(IO, '(100F25.17)') GradP(1:NI - 1, 1:NJ - 1, 1)
    Write(IO, '(100F25.17)') GradP(1:NI - 1, 1:NJ - 1, 2)

End Subroutine 
