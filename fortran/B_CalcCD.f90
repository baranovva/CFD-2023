Subroutine B_CalcCD(NI,NJ,VelCD,T,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,CFL,Re,Pr,RES)

    INTEGER :: NI,NJ
    integer, parameter:: IO = 20
    REAL :: VelCD(0:NI,0:NJ,2), T(0:NI,0:NJ), CFL, Re, Pr
    REAL :: CellVolume(NI-1,NJ-1)
    REAL :: CellCenter(0:NI,0:NJ,2), IFaceCenter( NI,NJ-1,2), IFaceVector( NI,NJ-1,2), JFaceCenter( NI-1,NJ,2), JFaceVector( NI-1,NJ,2)
    
    REAL :: VOL, RF(2), SF(2), RC(2), RN(2), VF(2), DC, DN, DNC, NF(2), TC, TN, TF, RES(0:NI,0:NJ), RESM, EPS, T_EX(0:NI,0:NJ), dt, dx, dx_min, dx_max
    INTEGER :: IFace,I_N,J_N,ITER,MAXITER
    
    RES = 0.0
    MAXITER = 100000
    EPS = 1.0E-8
    ITER=1
    T(0,:) = 1.0
    T(NI,:) = 2.0
    !T(:,0) = ((4*T(:,1)-T(:,2))/3.0 + T(:,1))/2.0
    !T(:,NJ) = ((4*T(:,NJ-1)-T(:,NJ-2))/3.0 + T(:,NJ-1))/2.0
    T_EX=T
    OPEN(IO, FILE='Residuals.plt')
    WRITE(IO,'(a13)') 'Iteration RES'
    
    DO WHILE (ITER<=MAXITER) !.AND. MAXVAL(RES(1:NI-1,1:NJ-1))>=EPS)
    RES = 0.0
        DO I=1,NI-1
        DO J=1,NJ-1
            dx = 0.0
            dx_min = 100000.0
            dx_max = 0.0
            DO IFace=1,4
                SELECT CASE(IFace)
                
                CASE(1)
                    I_N=I-1 !»ндексы соседних €чеек
                    J_N=J
                    RF(:) = IFaceCenter(I,J,:) !–адиус-вектор центра грани
                    SF(:) = -IFaceVector(I,J,:) !¬ектор площади грани
                CASE(2)
                    I_N=I+1 
                    J_N=J
                    RF(:) = IFaceCenter(I+1,J,:) 
                    SF(:) = IFaceVector(I+1,J,:) 
                CASE(3)
                    I_N=I 
                    J_N=J-1
                    RF(:) = JFaceCenter(I,J,:) 
                    SF(:) = -JFaceVector(I,J,:) 
                CASE(4)
                    I_N=I 
                    J_N=J+1
                    RF(:) = JFaceCenter(I,J+1,:) 
                    SF(:) = JFaceVector(I,J+1,:)
                END SELECT
                
                VOL = CellVolume(I,J)
                RC(:) = CellCenter(I,J,:) ! центр текущей €чейки 
                RN(:) = CellCenter(I_N,J_N,:) !центр соседней €чейки
                
                DC = Norm2(RF(:)-RC(:)) !норм2 - модуль вектора, рассто€ние до грани
                DN = Norm2(RF(:)-RN(:)) ! рассто€ние от центра соседней до грани
                DNC = Norm2(CellCenter(I_N,J_N,:)-CellCenter(I,J,:)) !–ассто€ние между центрами €чеек
                NF(:) = SF(:)/Norm2(SF(:)) !ѕерпендикул€рный вектор грани
                IF (2*DC<dx_min) dx_min=2*DC
                IF (2*DC>dx_max) dx_max=2*DC
                
                VF(1) = RLinearInterp(DC,DN,VelCD(I,J,1),VelCD(I_N,J_N,1))
                VF(2) = RLinearInterp(DC,DN,VelCD(I,J,2),VelCD(I_N,J_N,2))
                TC = T_EX(I,J)
                TN = T_EX(I_N,J_N)
                TF = RLinearInterp(DC,DN,TC,TN)
    
                IF (J_N==NJ) THEN
                    RES(I,J) = RES(I,J)
                ELSEIF (J_N==0) THEN
                    RES(I,J) = RES(I,J) + (DOT_PRODUCT(NF(:),VF(:))*TF)*NORM2(SF(:))
                ELSE
                    RES(I,J) = RES(I,J) + (DOT_PRODUCT(NF(:),VF(:))*TF-(TN-TC)/(DNC*Re*Pr))*NORM2(SF(:))
                ENDIF
                 
            ENDDO

            dx = (dx_min + dx_max)/2.0
            dt = CFL*dx/1.0 
            T(I,J) = T_EX(I,J) - dt*RES(I,J)/VOL
            
        ENDDO
        ENDDO

    !T(:,0) = ((4*T(:,1)-T(:,2))/3.0 + T(:,1))/2.0
    !T(:,NJ) = ((4*T(:,NJ-1)-T(:,NJ-2))/3.0 + T(:,NJ-1))/2.0
    IF (ITER==1) RESM = MAXVAL(RES(1:NI-1,1:NJ-1))
    WRITE(IO,*) ITER, MAXVAL(RES(1:NI-1,1:NJ-1))/RESM

    ITER=ITER+1
    T_EX=T
    
    IF (MAXVAL(ABS(RES(1:NI-1,1:NJ-1)))<EPS) EXIT
    ENDDO
End Subroutine 
