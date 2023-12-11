Program Main

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt',InputMesh='InputMesh.txt',Temper='Temperature.txt',Veloc='VelocityField.txt',OutputCD='OutputCD.plt'! names of input and output files
  character MeshFile*30        ! name of file with computational mesh
  integer, parameter:: IO = 12, II = 13, III=14, IIII=15 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume,X1,Y1,CellVolume1,T_et, T,TError ! scalar arrays
  real,allocatable,dimension(:,:,:):: GradP, GradPExact, GradPError, V, VelCD
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays
  real,allocatable,dimension(:,:,:):: CellCenter1,IFaceCenter1,IFaceVector1,JFaceCenter1,JFaceVector1
  real,allocatable,dimension(:,:):: DivV, DivVExact, DivVError
  real,allocatable,dimension(:,:):: LapP, LapPExact, LapPError
  real,allocatable,dimension(:,:):: RotV, RotVExact, RotVError
  INTEGER :: IT,NIter,MODE
  REAL :: CFL, Re, Pr

!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  !READ(IO,*) INT(NIter)
  NIter = 10
  MODE = 3
  CFL = 0.01
  Re = 100.0
  Pr = 1.0
  CLOSE(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile, 'and',InputMesh
  OPEN(IO,FILE = MeshFile)
  OPEN(II,FILE = InputMesh)
  READ(IO,*) NI,NJ
  READ(II,*) NK,NL
  WRITE(*,*) 'For spatial NI, NJ = ',NI,NJ
  WRITE(*,*) 'For C-D NK, NL = ',NK,NL

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(GradP(0:NI,0:NJ,2))
  allocate(GradPExact(0:NI,0:NJ,2))
  allocate(GradPError(0:NI,0:NJ,2))
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces
  allocate(V(0:NI,0:NJ,2))
  allocate(DivV(0:NI,0:NJ))
  allocate(DivVExact(0:NI,0:NJ))
  allocate(DivVError(0:NI,0:NJ))
  allocate(LapP(0:NI,0:NJ))
  allocate(LapPExact(0:NI,0:NJ))
  allocate(LapPError(0:NI,0:NJ))
  allocate(RotV(0:NI,0:NJ))
  allocate(RotVExact(0:NI,0:NJ))
  allocate(RotVError(0:NI,0:NJ))
  allocate(X1(NK,NL)) ! mesh nodes X-coordinates
  allocate(Y1(NK,NL)) ! mesh nodes Y-coordinates
  allocate(VelCD(0:NK,0:NL,2))
  allocate(T_et(0:NK,0:NL))
  allocate(T(0:NK,0:NL))
  allocate(TError(0:NK,0:NL))
  allocate(CellVolume1(NK-1,NL-1))
  allocate(CellCenter1(0:NK,0:NL,2)) ! Cell Centers
  allocate(IFaceCenter1( NK,NL-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector1( NK,NL-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter1( NK-1,NL,2)) ! Face Centers for J-faces
  allocate(JFaceVector1( NK-1,NL,2)) ! Face Vectors for I-faces
  
!=== READ FIELDS T AND V ===
  WRITE(*,*) 'Read fields from file: ',InputMesh
  OPEN(III,FILE = Veloc)
  OPEN(IIII,FILE = Temper)
  READ(III,*) ((VelCD(I,J,1),VelCD(I,J,2),I=0,NK),J=0,NL)
  READ(IIII,*) ((T_et(I,J),I=0,NK),J=0,NL)
  CLOSE(III)
  CLOSE(IIII)
  
!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile,'and',InputMesh
  READ(II,*) ((X1(I,J),Y1(I,J),I=1,NK),J=1,NL)
  READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  CLOSE(IO)
  CLOSE(II)
  
!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
  Call B_CalcMetric(NK,NL,X1,Y1,CellCenter1,CellVolume1,IFaceCenter1,IFaceVector1,JFaceCenter1,JFaceVector1) 
  
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
      Call Calc_GradP_Exact(CellCenter(I,J,1),CellCenter(I,J,2),GradPExact(I,J,:))
      CALL Velocity(CellCenter(I,J,1), CellCenter(i,j,2), V(I,J,:))
      DivVExact(I,J) = DivVelocityPExact(CellCenter(I,J,1), CellCenter(i,j,2))
      LapPExact(I,J) = RLaplacianPExact(CellCenter(I,J,1), CellCenter(i,j,2))
      RotVExact(I,J) = RotVelocityExact(CellCenter(I,J,1), CellCenter(i,j,2))
    ENDDO
  ENDDO

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives'       
  DO IT=1,NIter
  
  Call B_CalcGradient(NI,NJ,P,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
  GradPError=ABS(GradPExact-GradP)/GradPExact
  WRITE(*,*) 'Iteration: ', IT, 'MAXIMUM ERROR: ', MAXVAL(GradPError(1:NI-1,1:NJ-1,:))
 
  ENDDO
  write(*,*) 'Maximum Gradx-error:', maxval(GradPError(1:NI-1,1:NJ-1,1))
  write(*,*) 'Maximum Grady-error:', maxval(GradPError(1:NI-1,1:NJ-1,2))
 
  !=== CALCULATE DIVERGENCE ===
  WRITE(*,*) 'Calculate divergence'
  DivV = 0.0
  CALL B_CalcDivergence(NI,NJ,V,P,DivV,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,MODE)
  DivVError = ABS(DivVExact-DivV)/DivVExact
  WRITE(*,*) 'Maximum error of divergence: ', MAXVAL(DivVError(1:NI-1,1:NJ-1))
  
  !=== CALCULATE LAPLACIAN ===
  WRITE(*,*) 'Calculate laplacian'
  LapP = 0.0
  CALL B_CalcLaplacian(NI,NJ,P,GradP,LapP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
  LapPError = ABS(LapPExact-LapP)/LapPExact
  WRITE(*,*) 'Maximum error of laplacian: ', MAXVAL(LapPError(1:NI-1,1:NJ-1))
  
  !=== CALCULATE ROTTER ===
  WRITE(*,*) 'Calculate rotter'
  RotV = 0.0
  CALL B_CalcRotter(NI,NJ,V,RotV,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
  RotVError = ABS(RotVExact-RotV)/RotVExact
  WRITE(*,*) 'Maximum error of rotter: ', MAXVAL(RotVError(1:NI-1,1:NJ-1))
  
  !=== CALCULATE C-D ===
  WRITE(*,*) 'Calculate C-D'
  T = 1.0
  CALL B_CalcCD(NK,NL,VelCD,T,CellVolume1,CellCenter1,IFaceCenter1,IFaceVector1,JFaceCenter1,JFaceVector1,CFL,Re,Pr)
  TError = ABS(T_et-T)/T_et
  WRITE(*,*) 'Maximum error of C-D: ', MAXVAL(TError(1:NI-1,1:NJ-1)) !Рассмотреть флос с архи-чоу и без
  !OPEN(2, FILE='Residuals_FIELD.plt')
  !   Write(2,*) 'VARIABLES = "X", "Y", "RES"' 
  !Write(2,*) 'ZONE I=',NK,', J=',NL,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
  !Write(2,'(100F14.7)') X1(1:NI,1:NJ) 
  !Write(2,'(100F14.7)') Y1(1:NI,1:NJ)
  !Write(2,'(100F14.7)') RES(1:NI-1,1:NJ-1)
  !CLOSE(2)
    CLOSE(IO)
  
  !=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradPError,V,DivV,DivVError,LapP,LapPError,RotV,RotVError)
  Close(IO)
  OPEN(II,FILE=OutputCD)
  CALL B_OutputCD(II,NK,NL,X1,Y1,VelCD,T,TError)
  CLOSE(II)

END PROGRAM Main  
