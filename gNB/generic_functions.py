import numpy as np


def pi_2_bpsk(b):
    """3GPP 38.211 5.1.1 V.16.0.0

    :param b: bits
    :return: modulation symbols
    """
    return np.exp(1j*np.pi*np.remainder(list(range(len(b))), 2)/2) * ((1-2*b)+1j*(1-2*b)) / np.sqrt(2)


def bpsk(b):
    """3GPP 38.211 5.1.2 V.16.0.0

    :param b: bits
    :return: modulation symbols
    """
    return ((1-2*b)+1j*(1-2*b)) / np.sqrt(2)


def qpsk(b):
    """3GPP 38.211 5.1.3 V.16.0.0

    :param b: bits
    :return: modulation symbols
    """
    return ((1-2*b[::2])+1j*(1-2*b[1::2])) / np.sqrt(2)


def qam16(b):
    """3GPP 38.211 5.1.4 V.16.0.0

    :param b: bits
    :return: modulation symbols
    """
    return ((1-2*b[::4])*(2-(1-2*b[2::4]))+1j*(1-2*b[1::4])*(2-(1-2*b[3::4]))) / np.sqrt(10)


def qam64(b):
    """3GPP 38.211 5.1.5 V.16.0.0

    :param b: bits
    :return: modulation symbols
    """
    return ((1-2*b[::6])*(4-(1-2*b[2::6])*(2-(1-2*b[4::6])))+1j*(1-2*b[1::6])*(4-(1-2*b[3::6])*(2-(1-2*b[5::6])))) / np.sqrt(42)


def qam256(b):
    """3GPP 38.211 5.1.6 V.16.0.0

    :param b: bits
    :return: modulation symbols
    """
    return ((1-2*b[::8])*(8-(1-2*b[2::8])*(4-(1-2*b[4::8])*(2-(1-2*b[6::8]))))+1j*(1-2*b[1::8])*(8-(1-2*b[3::8])*(4-(1-2*b[5::8])*(2-(1-2*b[7::8]))))) / np.sqrt(170)


def gold_sequence(M_pn, c_init):
    """3GPP 38.211 5.2.1 V.16.0.0

    :param M_pn: pseudo random sequence length
    :param c_init: x2 initialization
    :return: pseudo random sequence
    """
    # constant
    N_c = 1600

    # initialize x1
    x1 = np.zeros((M_pn+N_c+31,), dtype='int8')
    x1[0] = 1

    # generate x1
    for n in range(len(x1)-31):
        x1[n+31] = np.remainder(x1[n+3] + x1[n], 2)

    # initialize x2
    x2 = np.zeros((M_pn+N_c+31,), dtype='int8')
    x2_init = list(bin(c_init)[2:])
    x2[:len(x2_init)] = x2_init

    # generate x2
    for n in range(len(x2)-31):
        x2[n+31] = np.remainder(x2[n+3] + x2[n+2] + x2[n+1] + x2[n], 2)

    # generate c sequence
    c = np.remainder(x1[N_c:N_c+M_pn] + x2[N_c:N_c+M_pn], 2)

    return c


def ofdm_type1(s, mu):
    """3GPP 38.211 5.3.1 V.16.0.0

    :param s: frequency signal
    :param mu: spacing configuration
    :return: time-continuous signal
    """
    return 0
