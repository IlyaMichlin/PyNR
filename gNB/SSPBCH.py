"""3GPP 38.211 7.4.3 V.16.0.0

Generate and map SS/BPCH block on the downlink grid
"""

import numpy as np

from gNB.SequenceGeneration import generate_pseudo_random_sequence


def _generate_cell_id(N_id_1, N_id_2):
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

    c = generate_pseudo_random_sequence(M_np, c_init)

    r = np.zeros((144,), dtype='compex')
    for m in range(144):
        r[m] = 1/np.sqrt(2)*(1-2*c[2*m]) + 1j/np.sqrt(2)*(1-2*c[2*m+1])

    return r


def _generate_sspbch_block(N_id_1, N_id_2, d_pbch, n_hf, L_max_hat, beta_pss=0, beta_pbch=1):
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
    # generate SS/PBCH block
    sspbch_block = np.zeros((240,4), dtype='complex')

    # generate PCI
    N_cell_id = _generate_cell_id(N_id_1, N_id_2)

    # map PSS
    pss = _generate_PSS(N_id_2)
    pss *= beta_pss
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
    d_pbch *= beta_pbch
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


def _map_sspbch(dl_grid, sspbch_block, N_cell_id, k_sbb):
    # TODO: implement _map_sspbch
    """3GPP 38.211 7.4.3 V.16.0.0
    :return:
    """


def ss_pbch_block(dl_grid, scs, shared_spectrum=False, large_operat_freq=False):
    # TODO: implement ss_pbch_block
    """

    :param dl_grid:
    :param scs:
    :param shared_spectrum:
    :param large_operat_freq:
    :return:
    """
    if scs == 15:  #khz
        if not shared_spectrum:
            if not large_operat_freq:
                n = [0, 1]
            else:
                n = [0, 1, 2, 3]
        else:
            n = [0, 1, 2, 3, 4]
        idxs =
        
    elif scs == 30:  #khz
