"""3GPP 38.211 7.4.3 V.16.0.0

Generate and map SS/BPCH block on the downlink grid
"""

import numpy as np

from downlink.general_procedures.crc import crc_calculation
from downlink.generic_functions import gold_sequence
from NRConstants import FR1_HIGH, FR2_LOW, FR2_HIGH


def generate_cell_id(N_id_1, N_id_2):
    """3GPP 38.211 7.4.2.1 V.16.0.0

    :param N_id_1: PCI group, range(0, 336)
    :param N_id_2: physical-layer identity, range(0, 3)
    :return: PCI (Physical Cell ID)
    """
    return 3*N_id_1+N_id_2


def _generate_PSS(N_id_2):
    """3GPP 38.211 7.4.2.2.1 V.16.0.0

    :param N_id_2: physical-layer identity, range(0, 3)
    :return: PSS sequence
    """
    # generate x
    x = np.zeros((121,), dtype='int')
    x[:7] = [0, 1, 1, 0, 1, 1, 1]
    for n in range(120):
        x[n+7] = np.remainder(x[n+4]+x[n], 2)

    # generate PSS sequence
    d_pss = np.zeros((127,), dtype='int')
    for n in range(127):
        m = np.remainder(n+43*N_id_2, 127)
        d_pss[n] = 1 - 2*x[m]

    return d_pss


def _generate_SSS(N_id_1, N_id_2):
    """3GPP 38.211 7.4.2.3.1 V.16.0.0

    :param N_id_1: PCI group, range(0, 336)
    :param N_id_2: physical-layer identity, range(0, 3)
    :return: SSS sequence
    """
    # generate x
    x0 = np.zeros((127,), dtype='int')
    x1 = np.zeros((127,), dtype='int')
    x0[:7] = [1, 0, 0, 0, 0, 0, 0]
    x1[:7] = [1, 0, 0, 0, 0, 0, 0]
    for n in range(120):
        x0[n+7] = np.remainder(x0[n+4] + x0[n], 2)
        x1[n+7] = np.remainder(x1[n+1] + x1[n], 2)

    # generate SSS sequence
    d_sss = np.zeros((127,), dtype='int')
    for n in range(127):
        m0 = 15 * int(np.floor(N_id_1/112)) + 5*N_id_2
        m1 = np.remainder(N_id_1, 112)
        d_sss[n] = (1 - 2*x0[np.remainder(n+m0, 127)]) * (1 - 2*x1[np.remainder(n+m1, 127)])

    return d_sss


def _generate_pbch_dmrs(N_cell_id, n_hf, L_max_hat):
    """3GPP 38.211 7.4.1.4.1 V.16.0.0

    :param N_cell_id: PCI (Physical Cell ID)
    :param n_hf: 0 for the first half-frame in the frame and 1 for the second half-frame. 3GPP 38.213 4.1 V.16.0.0
    :param L_max_hat: maximum half frame index 3GPP. 38.213 4.1 V.16.0.0
    :return: DMRS PBCH sequence
    """
    # constant
    M_np = 288

    if L_max_hat <= 4:
        i_ssb = L_max_hat & 3
        i_ssb_hat = i_ssb + 4*n_hf
    elif L_max_hat > 4:
        i_ssb = L_max_hat & 7
        i_ssb_hat = i_ssb
    else:
        return False
    c_init = 2**11*(i_ssb_hat+1)*(np.floor(N_cell_id/4)+1) + 2**6*(i_ssb+1)*(np.remainder(N_cell_id, 4))

    c = gold_sequence(M_np, c_init)

    r = np.zeros((144,), dtype='compex')
    for m in range(144):
        r[m] = 1/np.sqrt(2)*(1-2*c[2*m]) + 1j/np.sqrt(2)*(1-2*c[2*m+1])

    return r


def _generate_sspbch_block(N_id_1, N_id_2, d_pbch, n_hf, L_max_hat, beta_pss=0, beta_pbch=0):
    """3GPP 38.211 7.4.3.1 V.16.0.0

    :param N_id_1: PCI group, range(0, 336)
    :param N_id_2: physical-layer identity, range(0, 3)
    :param d_pbch: PBCH data length of 275 symbols
    :param n_hf: 0 for the first half-frame in the frame and 1 for the second half-frame. 3GPP 38.213 4.1 V.16.0.0
    :param L_max_hat: maximum half frame index 3GPP. 38.213 4.1 V.16.0.0
    :param beta_pss: PSS power boost compared to SSS power. Can be 0dB or 3dB. 3GPP 38.213 4.1 V.16.0.0
    :param beta_pbch: PDCCH DMRS EPRE to SSS EPRE is within -8 dB and 8 dB. 3GPP 38.213 4.1 V.16.0.0
    :return: SS/PBCH block
    """
    # check betas
    if beta_pss != 0 | beta_pss != 3:
        raise Exception('Incorrect PSS power boost')
    if beta_pbch < -8 | beta_pbch > 8:
        raise Exception('Incorrect PDCCH DMRS EPRE power')

    # generate SS/PBCH block
    sspbch_block = np.zeros((240,4), dtype='complex')

    # generate PCI
    N_cell_id = generate_cell_id(N_id_1, N_id_2)

    # map PSS
    pss = _generate_PSS(N_id_2)
    pss *= 10**(beta_pss/10)
    sspbch_block[56:183, 0] = pss

    # map SSS
    sss = _generate_SSS(N_id_1, N_id_2)
    sspbch_block[56:183, 2] = sss

    # map PBCH DMRS
    pbch_dmrs = _generate_pbch_dmrs(N_cell_id, n_hf, L_max_hat)
    v = np.remainder(N_cell_id, 4)
    map_idxs = np.arange(0, 237, 4) + v

    sspbch_block[map_idxs, 1] = pbch_dmrs[:60]
    sspbch_block[map_idxs[:12], 2] = pbch_dmrs[60:72]
    sspbch_block[map_idxs[-12:], 2] = pbch_dmrs[72:-60]
    sspbch_block[map_idxs, 3] = pbch_dmrs[-60:]

    # map PBCH
    d_pbch *= 10**(beta_pbch/10)
    idx = 0
    for n in range(240):
        # continue if index is DMRS index
        if np.isin(n, map_idxs):
            continue

        sspbch_block[n, 1] = d_pbch[idx]
        idx += 1

    for n in range(48):
        # continue if index is DMRS index
        if np.isin(n, map_idxs):
            continue

        sspbch_block[n, 2] = d_pbch[idx]
        idx += 1

    for n in range(192, 240):
        # continue if index is DMRS index
        if np.isin(n, map_idxs):
            continue

        sspbch_block[n, 2] = d_pbch[idx]
        idx += 1

    for n in range(240):
        # continue if index is DMRS index
        if np.isin(n, map_idxs):
            continue

        sspbch_block[n, 3] = d_pbch[idx]
        idx += 1

    return sspbch_block


def _pbch_rate_matching(d):
    # TODO: implement according to 38.212 7.1.5
    E = 864
    I_BIL = 0

    # implement according to 38.212 5.4.1
    f = d

    return f


def _pbch_channel_coding(c):
    # TODO: implement channel coding according to 38.212 7.1.4
    n_max = 9
    I_IL = 1
    n_pc = 0
    n_PC_wm = 0

    # TODO: implement channel coding according to 38.212 5.3.1
    d = c

    return d


def _pbch_crc(a_prime):
    """generate and concatenate CRC using CRC24C polynomial

    :param a_prime: input list of bits
    :return: input list of bits with CRC
    """
    # polynomial
    poly = 'crc24c'

    # generate parity bits using the polynomial g_CRC24C(D) generator
    c = crc_calculation(a_prime, poly)

    return c


def _pbch_scrambling(a, A, L_max, SFN, N_cell_id):
    """"""
    # TODO: implement _pbch_scrambling from 38.212 7.1.2
    # define M
    if L_max == 4 or L_max == 8:
        M = A-3
    elif L_max == 64:
        M = A-6
    else:
        # return 0 if L_max is incorrect
        return 0

    # define v
    v = SFN & 3

    # generate c
    M_pn = A+v*M
    c_init = N_cell_id
    c = gold_sequence(M_pn, c_init)

    # generate s
    s = np.zeros((A,))
    j = 0
    for i in range(A):
        # TODO: fix the if statement
        if a[i]:
            s[i] = 0
        else:
            s[i] = c[j+v*M]
            j += 1

    # generate a_prime
    a_prime = np.remainder(a+s, 2)

    return a_prime


def _generate_pbch_payload(SFN):
    """"""
    # TODO: implement _generate_pbch_block from 38.212 7 and 38.321 6.1.1
    pbch = np.zeros((275,))

    return pbch


def _map_ss_pbch_block(dl_half_frame, sspbch_block, scs, carrier_frequency, shared_spectrum, paired_spectrum):
    """3GPP 38.213 4.1 V16.0.0

    :param dl_half_frame: DL half frame
    :param sspbch_block: SS/PBCH block to map
    :param scs: subcarrier spacing
    :param carrier_frequency: frame carrier frequency
    :param shared_spectrum: shared spectrum indicator
    :param paired_spectrum: paired/unpaired spectrum indicator
    :return: mapped DL half frame with SS/PBCH block
    """
    # calculate SS/BPCH frequency indexes
    freq_start = int(dl_half_frame.shape[0]/2) - 120
    freq_end = freq_start + 240

    # case A
    if scs == 15e3:  # hz
        # SS/PBCH indexes
        idxs = [2, 8]
        # SS/PBCH time shift between n
        shift = 14
        # get n
        if not shared_spectrum:
            if carrier_frequency <= 3e9:
                n = [0, 1]
            elif 3e9 < carrier_frequency < FR1_HIGH:
                n = [0, 1, 2, 3]
            else:
                raise Exception('Incorrect configuration for SS/PBCH configuration in case A')
        else:
            n = [0, 1, 2, 3, 4]
    # case B
    elif scs == 30e3 and (869e6 < carrier_frequency < 894e6 or 2110e6 < carrier_frequency < 2200e6):  # operating band n5 or n66
        # SS/PBCH indexes
        idxs = [4, 8, 16, 20]
        # SS/PBCH time shift between n
        shift = 28
        # get n
        if carrier_frequency <= 3e9:
            n = [0]
        elif 3e9 < carrier_frequency < FR1_HIGH:
            n = [0, 1]
        else:
            raise Exception('Incorrect configuration for SS/PBCH configuration in case B')
    # case C
    elif scs == 30e3:  # Hz
        # SS/PBCH indexes
        idxs = [2, 8]
        # SS/PBCH time shift between n
        shift = 14
        # get n
        if not shared_spectrum:
            if paired_spectrum:
                if 3e9 < carrier_frequency < FR1_HIGH:
                    n = [0, 1, 2, 3]
                else:
                    raise Exception('Incorrect configuration for SS/PBCH configuration in case C')
            else:  # unpaired spectrum
                if carrier_frequency <= 2.4e9:
                    n = [0, 1]
                elif 2.4e9 < carrier_frequency < FR1_HIGH:
                    n = [0, 1, 2, 3]
                else:
                    raise Exception('Incorrect configuration for SS/PBCH configuration in case C')
        else:  # shared spectrum
            n = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    # case D
    elif scs == 120e3:  # Hz
        # SS/PBCH indexes
        idxs = [4, 8, 16, 20]
        # SS/PBCH time shift between n
        shift = 28
        # get n
        if FR2_LOW < carrier_frequency < FR2_HIGH:
            n = [0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18]
        else:
            raise Exception('Incorrect configuration for SS/PBCH configuration in case D')
    # case E
    elif scs == 240e3:  # Hz
        # SS/PBCH indexes
        idxs = [8, 12, 16, 20, 32, 36, 40, 44]
        # SS/PBCH time shift between n
        shift = 56
        # get n
        if FR2_LOW < carrier_frequency < FR2_HIGH:
            n = [0, 1, 2, 3, 5, 6, 7, 8]
        else:
            raise Exception('Incorrect configuration for SS/PBCH configuration in case E')
    else:
        raise Exception('Incorrect configuration for SS/PBCH configuration')

    # map SS/PBCH block
    for idx in idxs:
        for nn in n:
            ss_pbch_idx = idx + shift * nn
            dl_half_frame[freq_start:freq_end, ss_pbch_idx:ss_pbch_idx+4] = sspbch_block

    return dl_half_frame


def ss_pbch_block(dl_frame, ssbSubcarrierSpacing, carrier_frequency, shared_spectrum, paired_spectrum, N_id_1, N_id_2, L_max_hat, beta_pss, beta_pbch, sfn):
    """generate and map SS/PBCH block on one DL frame

    :param dl_frame: DL frame
    :param ssbSubcarrierSpacing: subcarrier spacing
    :param carrier_frequency: frame carrier frequency
    :param shared_spectrum: shared spectrum indicator
    :param paired_spectrum: paired/unpaired spectrum indicator
    :param N_id_1: PCI group, range(0, 336)
    :param N_id_2: physical-layer identity, range(0, 3)
    :param L_max_hat: maximum half frame index 3GPP. 38.213 4.1 V.16.0.0
    :param beta_pss: PSS power boost compared to SSS power. Can be 0dB or 3dB. 3GPP 38.213 4.1 V.16.0.0
    :param beta_pbch: PDCCH DMRS EPRE to SSS EPRE is within -8 dB and 8 dB. 3GPP 38.213 4.1 V.16.0.0
    :return: DL frame with mapped SS/PBCH blocks
    """
    # generate PBCH block
    d_pbch = _generate_pbch_payload(sfn)

    # loop over the two half frames
    for n_hf in range(2):
        # generate SS/PBCH block
        sspbch_block = _generate_sspbch_block(N_id_1, N_id_2, d_pbch, n_hf, L_max_hat, beta_pss, beta_pbch)
        # map SS/PBCH block on the half frame
        dl_frame[:, int(dl_frame.shape[1]/2*n_hf):int(dl_frame.shape[1]/2*(n_hf+1))] =\
            _map_ss_pbch_block(
                dl_frame[:, int(dl_frame.shape[1]/2*n_hf):int(dl_frame.shape[1]/2*(n_hf+1))],
                sspbch_block,
                ssbSubcarrierSpacing,
                carrier_frequency,
                shared_spectrum,
                paired_spectrum)

    return dl_frame
