Subroutine B_CalcGradient(NI,NJ,P,GradP,CellVolume,CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

    INTEGER :: NI,NJ
    REAL :: P(0:NI,0:NJ), GradP(0:NI,0:NJ,2)
    REAL :: CellVolume(NI-1,NJ-1)
    REAL :: CellCenter(0:NI,0:NJ,2), IFaceCenter( NI,NJ-1,2), IFaceVector( NI,NJ-1,2), JFaceCenter( NI-1,NJ,2), JFaceVector( NI-1,NJ,2)
    
    REAL :: GP(2),VOL, RF(2), SF(2), RC(2), RN(2), RE(2), GPE(2),PF,PE
    INTEGER :: IFace,I_N,J_N
    !INTEGER :: I,J,I_N,J_N,IFace,NCELL(4,2)
    !REAL :: VOL, RC(2), RF(2), SF(2), GP(2)
    !REAL :: DC,DN,RN(2),PF
    
    DO I=1,NI-1
        DO J=1,NJ-1
            GP = 0.0
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
                
                PE = RLinearInterp(DC,DN,P(I,J),P(I_N,J_N))
				
				RE(1) = RLinearInterp(DC,DN,CellCenter(I,J,1),CellCenter(I_N,J_N,1))
				RE(2) = RLinearInterp(DC,DN,CellCenter(I,J,2),CellCenter(I_N,J_N,2))
				
				GPE(1) = RLinearInterp(DC,DN,GradP(I,J,1),GradP(I_N,J_N,1))
				GPE(2) = RLinearInterp(DC,DN,GradP(I,J,2),GradP(I_N,J_N,2))
				
				PF = PE + DOT_PRODUCT((RF(:)-RE(:)),GPE(:)) 
                GP(:) = GP(:)+PF*SF(:)
                
            ENDDO
            Gradp(I,J,:) = GP(:)/VOL
            
        ENDDO
    ENDDO
    
End Subroutine 
