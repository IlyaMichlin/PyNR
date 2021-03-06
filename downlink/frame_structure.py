import numpy as np


# constants
delta_f_max = 480e3  # Hz. Maximum subcarrier spacing
N_f = 4096
T_c = 1/(delta_f_max*N_f)  # basic time unit for NR

delta_f_ref = 15e3  # Hz
N_f_ref = 2048
T_s = 1/(delta_f_ref*N_f_ref)  # basic time unit for LTE

k = T_s/T_c  # 64. The ratio between T_s and T_c

T_f = (delta_f_max * N_f / 100) * T_c  # 10ms. Radio frame duration
T_sf = (delta_f_max * N_f / 1000) * T_c  # 1ms. Subframe duration

N_rb_sc = 12  # Number of subcarriers per resource block


def get_delta_f(mu, cyclicPrefix):
    """3GPP 38.211 4.2 V.16.0.0

    :param mu: subcarrierSpacing
    :param cyclicPrefix: cyclic prefix
    :return: subcarrier spacing in Hz
    """
    if (cyclicPrefix == 'extended') & (mu == 2):
        return 2**mu * 15e3
    elif (cyclicPrefix == 'normal') & (0 <= mu <= 4):
        return 2**mu * 15e3
    else:
        raise Exception('Error: Incorrect cyclicPrefix and mu combination: cyclicPrefix={}, mu={}'.format(cyclicPrefix, mu))


def timing_advance(N_ta, N_ta_offset, msg):
    """3GPP 38.211 4.3.1 V.16.0.0

    :param N_ta: timing advance between downlink and uplink
    :param N_ta_offset: a fixed offset used to calculate the timing advance [5, TS 38.213]
    :param msg: message type
    :return: timing advance between downlink and uplink
    """
    if msg == 'msgA':
        return 0

    return (N_ta+N_ta_offset)*T_c


def get_Ns(mu, cyclicPrefix):
    """3GPP 38.211 4.3.2 V.16.0.0

    :param mu: subcarrierSpacing
    :param cyclicPrefix: cyclic prefix
    :return: N_symb_slot, N_slot_frame, N_slot_subframe
    """
    if (cyclicPrefix == 'extended') & (mu == 2):
        N_symb_slot = 12
        N_slot_frame = 40
        N_slot_subframe = 4
    elif cyclicPrefix == 'normal':
        N_symb_slot = 14
        if mu == 0:
            N_slot_frame = 10
            N_slot_subframe = 1
        elif mu == 1:
            N_slot_frame = 20
            N_slot_subframe = 2
        elif mu == 2:
            N_slot_frame = 40
            N_slot_subframe = 4
        elif mu == 3:
            N_slot_frame = 80
            N_slot_subframe = 8
        elif mu == 4:
            N_slot_frame = 160
            N_slot_subframe = 16
        else:
            raise Exception('Error: Incorrect mu for normal prefix')
    else:
        raise Exception('Error: Incorrect mu and prefix combination')

    return N_symb_slot, N_slot_frame, N_slot_subframe
