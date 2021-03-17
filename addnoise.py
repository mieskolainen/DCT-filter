# Create noisy image file
#
# Mikael Mieskolainen

import matplotlib.pyplot as plt
import numpy as np

MAXVAL = 255

def readimg(filename, M, N):
    """ Read binary file
    """
    with open(filename, 'rb') as f:
        #M, N = ...  # parse; advance file-pointer to data segment
        data  = np.fromfile(f, dtype='<f8', count=M*N)
        array = np.reshape(data, [M, N], order='C')
    return array

def writeimg(img, filename):

    with open(filename, 'wb') as f:
        img.tofile(f, sep='', format='%s')
    return True


A_o = readimg('lena.bin', M=512, N=512)

# Add Gaussian noise
sigma = 20
A_n = A_o + sigma*np.random.normal(loc=0, scale=1, size=A_o.shape)

# Truncate and re-scale (note that these violate naive image statistics)
A_n[A_n < 0] = 0
A_n = A_n / np.max(A_n[:])
A_n *= (np.max(A_o[:]))

writeimg(A_n, 'lena_noisy.bin')
