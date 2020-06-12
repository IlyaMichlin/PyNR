import numpy as np


def _sub_block_interleaving(d, N, E, K, Q_0, n_pc):
    """"""
    P = np.array([0,1,2,4,3,5,6,7,8,16,9,17,10,18,11,19,12,20,13,21,14,22,15,23,24,25,26,28,27,29,30,31])

    J = np.zeros((N,), dtype='int32')
    y = np.zeros((N,), dtype='int32')
    for n in range(N):
        i = (32*n)//N
        J[n] = P[i] * N//32 + n%(N//32)
        y[n] = d[J[n]]

    Q_roof_f_tmp = []
    if E < N:
        if K/E <= 7/16:
            Q_roof_f_tmp = np.intersect1d(Q_roof_f_tmp, J[:N-E])

            if E >= 3*N/4:
                Q_roof_f_tmp = np.intersect1d(Q_roof_f_tmp, np.arange(int(np.ceil(3*N/4 - E/2))))
            else:
                Q_roof_f_tmp = np.intersect1d(Q_roof_f_tmp, np.arange(int(np.ceil(9*N/16 - E/4))))
        else:
            Q_roof_f_tmp = np.intersect1d(Q_roof_f_tmp, J[E:N])

    Q_roof_i_tmp = Q_0 / Q_roof_f_tmp

    Q_roof_i = np.flip(np.argsort(P[Q_roof_i_tmp]))[:K+n_pc]

    Q_roof_f = Q_0 / Q_roof_i

    return Q_roof_f, Q_roof_i


def _bit_selection(y):
    """"""
    return 0


def _coded_bits_interleaving(e):
    """"""
    return 0


def polar_code_rate_matching(d, E, I_bil):
    """"""
    return 0
