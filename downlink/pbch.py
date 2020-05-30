import numpy as np

from downlink.general_procedures.channel_coding import polar_coding
from downlink.general_procedures.rate_matching import polar_code_rate_matching
from downlink.generic_functions import gold_sequence
from downlink.general_procedures.crc import crc_calculation


def payload_generation(a_roof, sfn, half_frame, L_roof_max, k_ssb, sspbch_index):
    """3GPP 38.212 7.1.1 V.16.0.0

    :param a_roof:transport block
    :param sfn: SFN
    :param half_frame: half frame bit
    :param L_roof_max: number of candidates
    :param k_ssb: SS/PBCH frequency difference
    :param sspbch_index: SS/PBCH block index
    :return: PBCH payload
    """
    A_roof = list(a_roof)
    [a_roof.append(sfn&n) for n in range(4)]
    a_roof.append(half_frame)
    if L_roof_max == 10:
        a_roof.append(list(bin(k_ssb)[2:])[-1])
        a_roof.append(0)
        a_roof.append(list(bin(sspbch_index)[2:])[-1])
        sspbch_indexs = [A_roof+7]
    elif L_roof_max == 20:
        a_roof.append(list(bin(k_ssb)[2:])[-1])
        [a_roof.append(sspbch_index&n) for n in [5,4]]
        sspbch_indexs = [A_roof+6,A_roof+7]
    if L_roof_max == 64:
        [a_roof.append(sspbch_index&n) for n in [6,5,4]]
        sspbch_indexs = [A_roof+5,A_roof+6,A_roof+7]
    else:
        a_roof.append(list(bin(k_ssb)[2:])[-1])
        [a_roof.append(0) for _ in range(2)]
        sspbch_indexs = []
    a_roof = np.array(a_roof)

    A = A_roof + 8
    j_sfn = 0
    j_hrf = 10
    j_ssb = 11
    j_other = 14
    G = np.array([16,23,18,17,8,30,10,6,24,7,0,5,3,2,1,4,9,11,12,13,14,15,19,20,21,22,25,26,27,28,29,31])

    a = np.zeros((A,))
    for i, a_roof_i in enumerate(a_roof):
        if A <= i <= A+3:
            a[G[j_sfn]] = a_roof_i
            j_sfn += 1
        elif i == A+4:
            a[G[j_hrf]] = a_roof_i
        elif A_roof+5 <= i <= A_roof+7:
            a[G[j_ssb]] = a_roof_i
            j_ssb += 1
        else:
            a[G[j_other]] = a_roof_i
            j_other += 1

    return a, sspbch_indexs


def scrambling(a, L_max, sfn, N_cell_id, sspbch_index):
    """3GPP 38.212 7.1.2 V.16.0.0

    :param a: bit sequence
    :param L_max: number of candidates
    :param sfn: SFN
    :param N_cell_id: cell ID
    :param sspbch_index: indexes of SS/BPCH indexs
    :return: scrambled bits
    """
    A = len(a)
    v = sfn&6>>1

    if L_max == 4 or L_max == 8:
        M = A-3
    elif L_max == 64:
        M = A-6
    else:
        raise Exception('Error: Incorrect L_max value of {}'.format(L_max))

    c_init = N_cell_id
    M_pn = v*M+A
    c = gold_sequence(M_pn, c_init)

    i = 0
    j = 0
    s = np.array((A,))
    while i < A:
        if A-7 <= i <= A-6 or i == A-4 or i in sspbch_index:
            s[i] = 0
        else:
            s[i] = c[j+v*M]
            j += 1
        i += 1

    a_prime = (a+s)%2

    return a_prime


def crc(a_prime):
    """3GPP 38.212 7.1.3 V.16.0.0

    :param a_prime: bit sequence
    :return: bits with CRC
    """
    b = crc_calculation(a_prime, 'crc24c')

    return b


def channel_coding(c):
    """"""
    n_max = 9
    I_il = 1
    n_pc = 0
    n_pc_wm = 0

    # TODO: implement according to 5.3.1
    d = polar_coding()

    return d


def rate_matching(d):
    """"""
    E = 864
    I_bil = 0

    # TODO: implement according to 5.4.1
    f = polar_code_rate_matching()
