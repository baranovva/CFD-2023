from numpy import savetxt


def output_fields(file_name: str, ni: int, nj: int, x: object, y: object,
                  p=None, v=None, grad_p=None, grad_p_error=None,
                  div_v_p=None, div_v_p_error=None, t=None, t_error=None,
                  lap_p=None, lap_p_error=None, rot_v=None, rot_v_error=None) -> None:
    with open(file_name, 'w') as f:
        header = 'VARIABLES = "x", "y"'
        if p is not None:
            header += ', "p"'
        if v is not None:
            header += ', "v_x", "v_y"'
        if grad_p is not None:
            header += ', "grad(p)_x", "grad(p)_y"'
        if grad_p_error is not None:
            header += ', "grad(p)_x_error", "grad(p)_y_error"'
        if div_v_p is not None:
            header += ', "div(vp)"'
        if div_v_p_error is not None:
            header += ', "div(vp)_error"'
        if t is not None:
            header += ', "T"'
        if t_error is not None:
            header += ', "T_error"'
        if lap_p is not None:
            header += ', "Lap(p)"'
        if lap_p_error is not None:
            header += ', "Lap(p)_error"'
        if rot_v is not None:
            header += ', "rot(v)"'
        if rot_v_error is not None:
            header += ', "rot(v)_error"'
        f.write(f'{header}\n')
        f.write(f'ZONE I={ni}, J={nj}, DATAPACKING=BLOCK,VARLOCATION=([3-20]=CELLCENTERED)\n')

        savetxt(f, x[1:, 1:], fmt='%14.7f')
        savetxt(f, y[1:, 1:], fmt='%14.7f')
        if p is not None:
            savetxt(f, p[1:ni, 1:nj], fmt='%14.7f')
        if v is not None:
            savetxt(f, v[1:ni, 1:nj, 0], fmt='%25.17f')
            savetxt(f, v[1:ni, 1:nj, 1], fmt='%25.17f')
        if grad_p is not None:
            savetxt(f, grad_p[1:ni, 1:nj, 0], fmt='%25.17f')
            savetxt(f, grad_p[1:ni, 1:nj, 1], fmt='%25.17f')
        if grad_p_error is not None:
            savetxt(f, grad_p_error[1:ni, 1:nj, 0], fmt='%25.17f')
            savetxt(f, grad_p_error[1:ni, 1:nj, 1], fmt='%25.17f')
        if div_v_p is not None:
            savetxt(f, div_v_p[1:ni, 1:nj], fmt='%25.17f')
        if div_v_p_error is not None:
            savetxt(f, div_v_p_error[1:ni, 1:nj], fmt='%25.17f')
        if t is not None:
            savetxt(f, t[1:ni, 1:nj], fmt='%25.17f')
        if t_error is not None:
            savetxt(f, t_error[1:ni, 1:nj], fmt='%25.17f')
        if lap_p is not None:
            savetxt(f, lap_p[1:ni, 1:nj], fmt='%25.17f')
        if lap_p_error is not None:
            savetxt(f, lap_p_error[1:ni, 1:nj], fmt='%25.17f')
        if rot_v is not None:
            savetxt(f, rot_v[1:ni, 1:nj], fmt='%25.17f')
        if rot_v_error is not None:
            savetxt(f, rot_v_error[1:ni, 1:nj], fmt='%25.17f')
