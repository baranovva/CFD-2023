from numpy import savetxt


def B_OutputFields(
        OutputFile, NI, NJ,
        X, Y, P, grad_p
):
    with open(OutputFile, 'w') as f:
        f.write(
                'VARIABLES = "X", "Y",'
                ' "P", "GradPX", "GradPY"\n'
        )
        f.write(
                f'ZONE I={NI}, J={NJ},'
                f' DATAPACKING=BLOCK,'
                f' VARLOCATION=([3-20]=CELLCENTERED)\n'
        )
        savetxt(f, X, fmt='%14.7f')
        savetxt(f, Y, fmt='%14.7f')
        savetxt(
                f, P[:NI - 1, :NJ - 1],
                fmt='%14.7f'
        )
        savetxt(
                f, grad_p[:NI - 1, :NJ - 1, 0],
                fmt='%25.17f'
        )
        savetxt(
                f, grad_p[:NI - 1, :NJ - 1, 1],
                fmt='%25.17f'
        )
