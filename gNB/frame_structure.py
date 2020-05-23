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


def get_delta_f(mu):
    """3GPP 38.211 4.2 V.16.0.0

    :param mu: subcarrierSpacing
    :return: subcarrier spacing in Hz
    """
    return 2**mu * 15e3


def get_Ns(mu, prefix):
    """3GPP 38.211 4.3.2 V.16.0.0

    :param mu: subcarrierSpacing
    :param prefix: cyclicPrefix
    :return:
    """
    if (prefix == 'extended') & (mu == 2):
        N_symb_slot = 12
        N_slot_frame = 40
        N_slot_subframe = 4
    elif prefix == 'normal':
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


# def symbols_in_subframe(mu):
#     """"""


def timing_advance(N_ta, N_ta_offset, msg):
    """

    :param N_ta: timing advance between downlink and uplink
    :param N_ta_offset: a fixed offset used to calculate the timing advance [5, TS 38.213]
    :param msg: message type
    :return: timing advance between downlink and uplink
    """
    if msg == 'msgA':
        return 0

    return (N_ta+N_ta_offset)*T_c
