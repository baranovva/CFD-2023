Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradPError,V,DivV,DivVError,LapP,LapPError,RotV,RotVError)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P, DivV, DivVError, LapP, LapPError, RotV, RotVError!, VelCD, T, TError
  Real,Dimension(0:NI,0:NJ,2)::GradP, GradPError, V

  Write(IO,*) 'VARIABLES = "X", "Y", "P", "GradP_X", "GradP_Y", "GradPError_X", "GradPError_Y", "V_x", "V_y", "DivV", "DivVError", "LapP", "LapPError", "RotV", "RotVError"' 
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.17)') GradP(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.17)') GradPError(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.17)') GradPError(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.17)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.17)') V(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.17)') DivV(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') DivVError(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') LapP(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') LapPError(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') RotV(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') RotVError(1:NI-1,1:NJ-1)

End Subroutine 

Subroutine B_OutputCD(IO,NI,NJ,X,Y,VelCD,T,TError)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ):: T, TError
  Real,Dimension(0:NI,0:NJ,2):: VelCD
  
  Write(IO,*) 'VARIABLES = "X", "Y", "VelC-D_X", "VelC-D_Y", "T", "T_Error"' 
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F25.17)') VelCD(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F25.17)') VelCD(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F25.17)') T(1:NI-1,1:NJ-1)
  Write(IO,'(100F25.17)') TError(1:NI-1,1:NJ-1)
  
End Subroutine