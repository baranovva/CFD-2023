Function Pressure(X,Y)
    !Pressure = x + y
    Pressure = x**3+y**3+x**2
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
    
    !V(1) = 1.0 + X
    !V(2) = 1.0 + Y
    V(1) = 5.0 + 5.0*Y + Y**2
    V(2) = 5.0 + 10.0*X + X**2
    !V(1) = -Y
    !V(2) = X
    
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
    RLaplacianPExact = 6.0*X+6.0*Y+2.0
    END FUNCTION
    
    FUNCTION RotVelocityExact(X,Y)
    REAL :: X, Y
    RotVelocityExact = 2.0*(X-Y)+5.0
    END FUNCTION