from numpy import savetxt


def output_fields(file_name: str, ni: int, nj: int, x: object,
                  y: object, p: object, v: object, grad_p: object,
                  grad_p_error: object, div_v=None, div_v_error=None) -> None:
    with open(file_name, 'w') as f:
        f.write(f'VARIABLES = "X", "Y", "P", "VX", "VY", "GradPX", "GradPY", "GradPErrorX", "GradPErrorY"\n')

        '''  , "DivVx", "DivVy",'
                ' "DivVxError", "DivVyError"  '''
        f.write(f'ZONE I={ni}, J={nj},'
                f' DATAPACKING=BLOCK,'
                f' VARLOCATION=([3-20]=CELLCENTERED)'
                f'\n')

        savetxt(f, x[:, :], fmt='%14.7f')
        savetxt(f, y[:, :], fmt='%14.7f')

        savetxt(f, p[1:ni, 1:nj], fmt='%14.7f')

        # v_x and v_y
        savetxt(f, v[1:ni, 1:nj, 0], fmt='%25.17f')
        savetxt(f, v[1:ni, 1:nj, 1], fmt='%25.17f')

        # grad_p_x and grad_p_y
        savetxt(f, grad_p[1:ni, 1:nj, 0], fmt='%25.17f')
        savetxt(f, grad_p[1:ni, 1:nj, 1], fmt='%25.17f')

        # grad_p_error_x and grad_p_error_y
        savetxt(f, grad_p_error[1:ni, 1:nj, 0], fmt='%25.17f')
        savetxt(f, grad_p_error[1:ni, 1:nj, 1], fmt='%25.17f')
        '''
        # div_v_error_x and div_v_error_y
        savetxt(f, div_v_error[1:ni, 1:nj, 0], fmt='%25.17f')
        savetxt(f, div_v_error[1:ni, 1:nj, 1], fmt='%25.17f')
'''
