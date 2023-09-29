from numpy import savetxt


def output_fields(
        file_name: str,
        ni: int,
        nj: int,
        x: object,
        y: object,
        p: object,
        grad_p: object,
        grad_p_error: object
) -> None:
    with open(file_name, 'w') as f:
        f.write(
                'VARIABLES = "X", "Y",'
                ' "P", "GradPX", "GradPY",'
                ' GradPErrorX","GradPErrorY"'
                ' \n'
        )
        f.write(
                f'ZONE I={ni}, J={nj},'
                f' DATAPACKING=BLOCK,'
                f' VARLOCATION=([3-20]=CELLCENTERED)'
                f'\n'
        )

        savetxt(f, x, fmt='%14.7f')
        savetxt(f, y, fmt='%14.7f')
        savetxt(
                f,
                p[:ni - 1, :nj - 1],
                fmt='%14.7f'
        )

        # grad_p_x and grad_p_y
        savetxt(
                f,
                grad_p[:ni - 1, :nj - 1, 0],
                fmt='%25.17f'
        )
        savetxt(
                f,
                grad_p[:ni - 1, :nj - 1, 1],
                fmt='%25.17f'
        )

        # grad_p_error_x and grad_p_error_y
        savetxt(
                f,
                grad_p_error[:ni - 1, :nj - 1, 0],
                fmt='%25.17f'
        )
        savetxt(
                f,
                grad_p_error[:ni - 1, :nj - 1, 1],
                fmt='%25.17f'
        )
