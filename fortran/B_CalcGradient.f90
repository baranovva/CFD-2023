Subroutine B_CalcGradient(NI, NJ, P, GradP)

integer :: NI, NJ
real P(0:NI,0:NJ), GradP(0:NI,0:NJ,2)

do  i = 1, NI-1
    do j = 1, NJ-1

    GradP(i,j,:) = P(i,j)

    end do
end do

End Subroutine 
