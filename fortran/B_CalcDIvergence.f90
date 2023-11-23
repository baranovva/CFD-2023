Subroutine B_CalcDivergence(NI,NJ,V,P,DivV,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,MODE)

    INTEGER :: NI,NJ
    REAL :: V(0:NI,0:NJ,2), DivV(0:NI,0:NJ), P(0:NI,0:NJ), GradP(0:NI,0:NJ,2)
    REAL :: CellVolume(NI-1,NJ-1)
    REAL :: CellCenter(0:NI,0:NJ,2), IFaceCenter( NI,NJ-1,2), IFaceVector( NI,NJ-1,2), JFaceCenter( NI-1,NJ,2), JFaceVector( NI-1,NJ,2)
    
    REAL :: VOL, RF(2), SF(2), RC(2), RN(2), VF(2), DC, DN, PF, GFLUX, PN, GC, GB, GN
    INTEGER :: IFace,I_N,J_N,MODE
    !INTEGER :: I,J,I_N,J_N,IFace,NCELL(4,2)
    !REAL :: VOL, RC(2), RF(2), SF(2), GP(2)
    !REAL :: DC,DN,RN(2),PF
    
    DO I=1,NI-1
        DO J=1,NJ-1
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
                
                VF(1) = RLinearInterp(DC,DN,V(I,J,1),V(I_N,J_N,1))
                VF(2) = RLinearInterp(DC,DN,V(I,J,2),V(I_N,J_N,2))
                
                GFLUX = DOT_PRODUCT(SF(:),VF(:)) !Для того, чтобы понять куда течет: из или в
                SELECT CASE(MODE)
                CASE(1)!Центральная
                    
                    PF = RLinearInterp(DC,DN,P(I,J),P(I_N,J_N))
                    DivV(I,J) = DivV(I,J) + DOT_PRODUCT(PF*VF,SF)
                    
                CASE(2)!Противопоточная 1-го порядка
                    
                    IF (GFLUX>=0.0) THEN
                        PF=P(I,J)
                    ELSE
                        PF=P(I_N,J_N)
                        IF (DN<1E-6) PF=2*P(I_N,J_N)-P(I,J) !Экстраполяция для заграничных
                    ENDIF
                    DivV(I,J) = DivV(I,J) + DOT_PRODUCT(PF*VF,SF)
                    
                CASE(3)!Противопоточная 2-го порядка
                    
                    IF (GFLUX>=0.0) THEN
                        PF=P(I,J)+DOT_PRODUCT(RF(:)-CellCenter(I,J,:),GradP(I,J,:))
                    ELSE
                        PF=P(I_N,J_N)+DOT_PRODUCT(RF(:)-CellCenter(I_N,J_N,:),GradP(I_N,J_N,:))
                        IF (DN<1E-6) THEN
                            PN=2*P(I_N,J_N)-P(I,J)
                            GC=DOT_PRODUCT(GradP(I,J,:),CellCenter(I,J,:)-RF(:))
                            GB=(P(I,J)-P(I_N,J_N))
                            GN=4*GB-3*GC
                            PF=PN+GN
                        ENDIF
                    ENDIF
                    DivV(I,J) = DivV(I,J) + DOT_PRODUCT(PF*VF,SF)
                    
                ENDSELECT                  
                
                
            ENDDO
            DivV(I,J) = DivV(I,J)/VOL
            
        ENDDO
    ENDDO
    
End Subroutine 
