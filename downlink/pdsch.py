import numpy as np

from downlink.generic_functions import gold_sequence


def scrambling(b, n_rnti, q, n_id):
    """"""
    M_pn = len(b)
    c_init = n_rnti*2**15+q*2**14+n_id
    c = gold_sequence(M_pn, c_init)
    b_tilda = (b+c)%2

    return 0
