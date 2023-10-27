PROGRAM OPENMP_3D_ARRAY
    INCLUDE 'omp_lib.h'
    INTEGER I, J, K, L, NI, NJ, NK, INFILE, OUTFILE
    DOUBLE PRECISION T1, T2
    DOUBLE PRECISION, ALLOCATABLE :: X(:, :, :), Y(:, :, :), Z(:, :, :), V(:, :, :)

    OPEN(INFILE, FILE = "CUBE.MSH")
    READ(INFILE, *) NI, NJ, NK
    ALLOCATE (X(NI, NJ, NK), Y(NI, NJ, NK), Z(NI, NJ, NK), V(NI - 1, NJ - 1, NK - 1))
    DO K = 1, NK
        DO J = 1, NJ
            DO I = 1, NI
                READ(INFILE, *) X(I, J, K), Y(I, J, K), Z(I, J, K)
            END DO
        END DO
    END DO
    CLOSE(INFILE)

    V = 0.0
    T1 = OMP_GET_WTIME()
    !$OMP PARALLEL PRIVATE(I,J,K) !���������� ������������ ������� ��������� � ���������� �����������
    DO L = 1, 30
        !$OMP DO SCHEDULE(static) COLLAPSE(3) !��������� ������� ������������� �������� ����� ����� ������
        DO K = 1, NK - 1
            DO J = 1, NJ - 1
                DO I = 1, NI - 1
                    V(I, J, K) = (X(I + 1, J, K) - X(I, J, K)) * (Y(I, J + 1, K) - Y(I, J, K)) * (Z(I, J, K + 1) - Z(I, J, K))
                END DO
            END DO
        END DO
        !$OMP END DO
    END DO
    !$OMP END PARALLEL
    T2 = OMP_GET_WTIME()
    PRINT *, 'computational time: ', T2 - T1, 's' !���������� ������� ������� ���������� ��������� (�.�. ��������� ������� ����� � 01.01.1970)
    WRITE(*, *) 'Max volume: ', MAXVAL(V) !���������� ������

    OPEN(OUTFILE, FILE = "CUBE_NEW.MSH")
    WRITE(OUTFILE, *) NI, NJ, NK
    DO K = 1, NK
        DO J = 1, NJ
            DO I = 1, NI
                WRITE(OUTFILE, *) X(I, J, K), Y(I, J, K), Z(I, J, K)
            END DO
        END DO
    END DO
    CLOSE(OUTFILE)
    DEALLOCATE(X, Y, Z, V)

END
