import numpy as np


def _xor(a, b):
    # initialize result
    result = []

    # Traverse all bits, if bits are
    # same, then XOR is 0, else 1
    for i in range(1, len(b)):
        if a[i] == b[i]:
            result.append(0)
        else:
            result.append(1)

    return result


def _mod2div(divident, divisor):
    # Performs Modulo-2 division
    # Number of bits to be XORed at a time.
    pick = len(divisor)

    # Slicing the divident to appropriate
    # length for particular step
    tmp = divident[:pick]

    while pick < len(divident):
        if tmp[0] == 1:
            # replace the divident by the result
            # of XOR and pull 1 bit down
            tmp = _xor(divisor, tmp)
        else:  # If leftmost bit is '0'
            # If the leftmost bit of the dividend (or the
            # part used in each step) is 0, the step cannot
            # use the regular divisor; we need to use an
            # all-0s divisor.
            tmp = _xor([0] * pick, tmp)

        tmp.append(divident[pick])

        # increment pick to move further
        pick += 1

    # For the last n bits, we have to carry it out
    # normally as increased value of pick will cause
    # Index Out of Bounds.
    if tmp[0] == 1:
        tmp = _xor(divisor, tmp)
    else:
        tmp = _xor([0] * pick, tmp)

    checkword = tmp
    return checkword


def _crc_encoder(data, key):
    # Function used at the sender side to encode
    # data by appending remainder of modular division
    # at the end of data.
    # https://www.geeksforgeeks.org/cyclic-redundancy-check-python/
    l_key = len(key)

    # Appends n-1 zeroes at end of data
    appended_data = data + [0] * (l_key - 1)
    remainder = _mod2div(appended_data, key)

    # Append remainder in the original data
    codeword = data + remainder
    return codeword


def crc_calculation(a, poly):
    """3GPP 38.212 5.1 V.16.0.0

    :param a: input bits
    :param poly: cyclic generator polynomials name
    :return: input bits with CRC parity bits
    """
    if poly == 'crc24a':
        L = 24
        g_crc = np.zeros((L+1,), dtype='int8')
        g_crc[[24,23,18,17,14,11,10,7,6,5,4,3,1,0]] = 1
    elif poly == 'crc24b':
        L = 24
        g_crc = np.zeros((L+1,), dtype='int8')
        g_crc[[24,23,6,5,1,0]] = 1
    elif poly == 'crc24c':
        L = 24
        g_crc = np.zeros((L+1,), dtype='int8')
        g_crc[[24,23,21,20,17,15,13,12,8,4,2,1,0]] = 1
    elif poly == 'crc16':
        L = 16
        g_crc = np.zeros((L+1,), dtype='int8')
        g_crc[[16,12,5,1]] = 1
    elif poly == 'crc11':
        L = 11
        g_crc = np.zeros((L+1,), dtype='int8')
        g_crc[[11,10,9,5,0]] = 1
    elif poly == 'crc6':
        L = 6
        g_crc = np.zeros((L+1,), dtype='int8')
        g_crc[[6,5,0]] = 1
    else:
        raise Exception('Error: Incorrect CRC polynomial ({})'.format(poly))

    b = _crc_encoder(a, g_crc)

    return b
