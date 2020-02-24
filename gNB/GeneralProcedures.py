"""3GPP 38.212 5 V.16.0.0

Data and control streams from/to MAC layer are encoded /decoded to offer transport and control services over the radio
transmission link. Channel coding scheme is a combination of error detection, error correcting, rate matching,
interleaving and transport channel or control information mapping onto/splitting from physical channels.
"""

import numpy as np

from gNB.crc import crc_encoder


def crc_calculation(a, poly):
    """3GPP 38.212 5.1 V.16.0.0

    :param a: input bits
    :param poly: cyclic generator polynomials name
    :return: input bits with CRC parity bits
    """
    if poly == 'CRC24A':
        L = 24
        g_crc = np.zeros((L+1,))
        g_crc[[24,23,18,17,14,11,10,7,6,5,4,3,1,0]] = 1
    elif poly == 'CRC24B':
        L = 24
        g_crc = np.zeros((L+1,))
        g_crc[[24,23,6,5,1,0]] = 1
    elif poly == 'CRC24C':
        L = 24
        g_crc = np.zeros((L+1,))
        g_crc[[24,23,21,20,17,15,13,12,8,4,2,1,0]] = 1
    elif poly == 'CRC16':
        L = 16
        g_crc = np.zeros((L+1,))
        g_crc[[16,12,5,1]] = 1
    elif poly == 'CRC11':
        L = 11
        g_crc = np.zeros((L+1,))
        g_crc[[11,10,9,5,0]] = 1
    elif poly == 'CRC6':
        L = 6
        g_crc = np.zeros((L+1,))
        g_crc[[6,5,0]] = 1

    b = crc_encoder(a, g_crc)

    return b

