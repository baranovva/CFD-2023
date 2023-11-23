Function Pressure(X,Y)
    !Pressure = x + y
    Pressure = x**2+y**2
    End Function

Function RLinearInterp(D1,D2,X1,X2)
    RLinearInterp = (X1*D2+X2*D1)/(D1+D2)
END FUNCTION

SUBROUTINE Calc_GradP_Exact(X,Y,GPE)
REAL:: X,Y
REAL :: GPE(2)

GPE(1) = 1.0
GPE(2) = 1.0
!GPE(1) = 2.0*X
!GPE(2) = 2.0*Y

    END SUBROUTINE
    
    SUBROUTINE Velocity(x,y,V)
    REAL :: x, y, V(2)
    
    V(1) = 1.0 + X
    V(2) = 1.0 + Y
    
    END SUBROUTINE
    
    FUNCTION DivVelocityExact(X,Y)
    REAL :: X, Y
    DivVelocityExact = 2.0
    END FUNCTION
    
    FUNCTION DivVelocityPExact(X,Y)
    REAL :: X, Y
    DivVelocityPExact = 2.0+3.0*X+3.0*Y
    END FUNCTION
    
    FUNCTION RLaplacianPExact(X,Y)
    REAL :: X, Y
    RLaplacianPExact = 4.0
    END FUNCTION