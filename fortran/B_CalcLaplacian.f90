Subroutine B_CalcLaplacian(NI,NJ,P,GradP,LapP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

    INTEGER :: NI,NJ
    REAL :: P(0:NI,0:NJ), GradP(0:NI,0:NJ,2), LapP(0:NI,0:NJ)
    REAL :: CellVolume(NI-1,NJ-1)
    REAL :: CellCenter(0:NI,0:NJ,2), IFaceCenter( NI,NJ-1,2), IFaceVector( NI,NJ-1,2), JFaceCenter( NI-1,NJ,2), JFaceVector( NI-1,NJ,2)
    
    REAL :: GP(2),VOL, RF(2), SF(2), RC(2), RN(2), DNC, NF(2), RNC(2), GF(2), dpdn, dpdn_c
    INTEGER :: IFace,I_N,J_N
    
    DO I=1,NI-1
        DO J=1,NJ-1
            GP = 0.0
            DO IFace=1,4
                SELECT CASE(IFace)
                
                CASE(1)
                    I_N=I-1 !Индексы соседних ячеек
                    J_N=J
                    RF(:) = IFaceCenter(I,J,:) !Радиус-вектор центра грани
                    SF(:) = -IFaceVector(I,J,:) !Вектор площади грани
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
                RC(:) = CellCenter(I,J,:) ! центр текущей ячейки 
                RN(:) = CellCenter(I_N,J_N,:) !центр соседней ячейки
                
                DC = Norm2(RF(:)-RC(:)) !норм2 - модуль вектора, расстояние до грани
                DN = Norm2(RF(:)-RN(:)) ! расстояние от центра соседней до грани
                DNC = Norm2(CellCenter(I_N,J_N,:)-CellCenter(I,J,:)) !Расстояние между центрами ячеек
                NF(:) = SF(:)/Norm2(SF(:)) !Перпендикулярный вектор грани
                
                !for skew correction
                
                RNC(:)=(CellCenter(I_N,J_N,:)-CellCenter(I,J,:))/DNC !Вектор по линии центров
                GF(1) = RLinearInterp(DC,DN,GradP(I,J,1),GradP(I_N,J_N,1))
				GF(2) = RLinearInterp(DC,DN,GradP(I,J,2),GradP(I_N,J_N,2))
                
                dpdn = (P(I_N,J_N)-P(I,J))/DNC
                
                IF (DN<1E-7) THEN
                    dpdn_c=DOT_PRODUCT(GradP(I,J,:),NF(:))
                    dpdn=5.0/3.0*dpdn-2.0/3.0*dpdn_c !2-ORDER
                    GF(:)=GradP(I,J,:)
                ENDIF
                			
				dpdn = dpdn + DOT_PRODUCT((NF(:)-RNC(:)),GF(:)) !skew correction
                LapP(I,J)=LapP(I,J)+dpdn*Norm2(SF(:))
                
            ENDDO
            LapP(I,J) = LapP(I,J)/VOL
            
        ENDDO
    ENDDO
    
End Subroutine 
