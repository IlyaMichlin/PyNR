import numpy as np


def power(x, in_dB=True):
    """calculate signal mean power

    :param x: input signal vector
    :return: power in dBm
    """
    s_power = np.mean(np.abs(x)**2)

    if in_dB:
        return 10*np.log10(s_power)

    return s_power
