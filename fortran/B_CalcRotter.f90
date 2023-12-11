Subroutine B_CalcRotter(NI,NJ,V,RotV,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

    INTEGER :: NI,NJ
    REAL :: V(0:NI,0:NJ,2), RotV(0:NI,0:NJ)
    REAL :: CellVolume(NI-1,NJ-1)
    REAL :: CellCenter(0:NI,0:NJ,2), IFaceCenter( NI,NJ-1,2), IFaceVector( NI,NJ-1,2), JFaceCenter( NI-1,NJ,2), JFaceVector( NI-1,NJ,2)
    
    REAL :: VOL, RF(2), SF(2), RC(2), RN(2), VF(2), DC, DN
    INTEGER :: IFace,I_N,J_N
    
    DO I=1,NI-1
        DO J=1,NJ-1
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
                
                VF(1) = RLinearInterp(DC,DN,V(I,J,1),V(I_N,J_N,1))
                VF(2) = RLinearInterp(DC,DN,V(I,J,2),V(I_N,J_N,2))
    
                RotV(I,J) = RotV(I,J) + (SF(1)*VF(2)-SF(2)*VF(1))          
                                
            ENDDO
            RotV(I,J) = RotV(I,J)/VOL
            
        ENDDO
    ENDDO
    
End Subroutine 
