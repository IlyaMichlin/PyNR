"""3GPP 38.211 5.2 V.16.0.0

"""

import numpy as np


def generate_pseudo_random_sequence(M_pn, c_init):
    """3GPP 38.211 5.2.1 V.16.0.0

    :param M_pn: pseudo random sequence length
    :param c_init: x2 initializarion
    :return: pseudo random sequence
    """
    # constant
    N_c = 1600

    # initialize x1
    x1 = np.zeros((M_pn+N_c+31,), dtype='int')
    x1[:2] = [1, 0]

    # generate x1
    for n in range(len(x1)-31):
        x1[n+31] = np.remainder(x1[n+3] + x1[n], 2)

    # initialize x2
    x2 = np.zeros((M_pn+N_c+31,), dtype='int')
    x2_init = [int(x) for x in bin(c_init)[2:]]
    x2[:len(x2_init)] = x2_init

    # generate x2
    for n in range(len(x2)-31):
        x2[n+31] = np.remainder(x2[n+3] + x2[n+2] + x2[n+1] + x2[n], 2)

    # generate c sequence
    c = np.zeros((M_pn,), dtype='int')
    for n in range(M_pn):
        c[n] = np.remainder(x1[n+N_c] + x2[n+N_c], 2)

    return c
