Function Pressure(X,Y)
    !Pressure = exp(x)+sin(y)
    Pressure = x**2+y**2+x+y+1.0
    !Pressure = exp(x)+sin(y)+x+y
    !Pressure = x**3+y**3+x**2
    End Function

Function RLinearInterp(D1,D2,X1,X2)
    RLinearInterp = (X1*D2+X2*D1)/(D1+D2)
END FUNCTION

SUBROUTINE Calc_GradP_Exact(X,Y,GPE)
REAL:: X,Y
REAL :: GPE(2)

!GPE(1) = 1.0
!GPE(2) = 1.0
!GPE(1) = exp(x)+1.0
!GPE(2) = cos(y)+1.0
GPE(1) = 2.0*X+1.0
GPE(2) = 2.0*Y+1.0

    END SUBROUTINE
    
    SUBROUTINE Velocity(x,y,V)
    REAL :: x, y, V(2)
    
    V(1) = 1.0 + X
    V(2) = 1.0 + Y
    !V(1) = 5.0 + 5.0*Y + Y**2
    !V(2) = 5.0 + 10.0*X + X**2
    !V(1) = -(Y**3+y)
    !V(2) = X**3+x
    V(1) = cos(y)
    V(2) = exp(x)
    
    END SUBROUTINE
    
    FUNCTION DivVelocityExact(X,Y)
    REAL :: X, Y
    DivVelocityExact = 2.0
    END FUNCTION
    
    FUNCTION DivVelocityPExact(X,Y)
    REAL :: X, Y
    !DivVelocityPExact = 2.0
    !DivVelocityPExact = 4.0+12.0*X**2+12.0*Y**2
    !DivVelocityPExact = 5.0*X+5.0*Y+4.0*X**2+4.0*Y**2+2.0
    !DivVelocityPExact = 5.0*X**3+5.0*Y**3+3.0*x**2+3.0*Y**2+3.0*X+3.0*Y+2.0
    DivVelocityPExact = (3.0+x)*exp(x)+(1.0+y)*cos(y)+2.0*sin(y)
    END FUNCTION
    
    FUNCTION RLaplacianPExact(X,Y)
    REAL :: X, Y
    RLaplacianPExact = exp(x)-sin(y)
    !RLaplacianPExact = 4.0
    !RLaplacianPExact = 12.0*X**2+12.0*Y**2+4.0
    END FUNCTION
    
    FUNCTION RotVelocityExact(X,Y)
    REAL :: X, Y
    !RotVelocityExact = 3.0*(X**2+Y**2)+2.0
    RotVelocityExact = exp(x)+sin(y)
    !RotVelocityExact = 2.0*(X-Y)+5.0
    END FUNCTION