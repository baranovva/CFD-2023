Program Main

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt' ! names of input and output files
  character MeshFile*30        ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume ! scalar arrays
  real,allocatable,dimension(:,:,:):: GradP, GradPExact, GradPError
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays
  INTEGER :: IT,NIter

!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  !READ(IO,*) INT(NIter)
  NIter = 10
  CLOSE(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ

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

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  CLOSE(IO)

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
      Call Calc_GradP_Exact(CellCenter(I,J,1),CellCenter(I,J,2),GradPExact(I,J,:))
    ENDDO
  ENDDO

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives'       
  DO IT=1,NIter
  
  Call B_CalcGradient(NI,NJ,P,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
  GradPError=ABS(GradPExact-GradP)!/GradPExact
  WRITE(*,*) 'Iteration: ', IT, 'MAXIMUM ERROR: ', MAXVAL(GradPError(1:NI-1,1:NJ-1,:))
 
  ENDDO
  write(*,*) 'Maximum Gradx-error:', maxval(GradPError(1:NI-1,1:NJ-1,1))
  write(*,*) 'Maximum Grady-error:', maxval(GradPError(1:NI-1,1:NJ-1,2))

!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradPError)
  Close(IO)

END PROGRAM Main  
