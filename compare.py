# Open the denoised file
#
# Mikael Mieskolainen

import matplotlib.pyplot as plt
import numpy as np
from skimage import exposure
from PIL import Image


def readimg(filename, M, N):
  """ Read binary file
  """
  with open(filename, 'rb') as f:
    #M, N = ...  # parse; advance file-pointer to data segment
    data  = np.fromfile(f, dtype='<f8', count=M*N)
    array = np.reshape(data, [M,N], order='C')
 
  ind = ~np.isfinite(array)
  if np.sum(ind) > 0:
  	print(f'{np.sum(ind)} NaN found in {filename}')
 
  array[ind] = 0

  return array


M=512
N=512

A_o = readimg('lena.bin', M=M, N=N)
A_n = readimg('lena_noisy.bin', M=M, N=N)
A_d = readimg('lena_denoised.bin', M=M, N=N)

# -----------------------------------------------------

## Normalize the scales
MAXVAL = 2**8-1

A_o = A_o / np.max(A_o[:]) * MAXVAL
A_n = A_n / np.max(A_n[:]) * MAXVAL
A_d = A_d / np.max(A_d[:]) * MAXVAL


## Save images
im = Image.fromarray(A_o).convert("L")
im.save('./img/original.png')

im = Image.fromarray(A_n).convert("L")
im.save('./img/noisy.png')

im = Image.fromarray(A_d).convert("L")
im.save('./img/denoised.png')


## Compute metrics
MSE_n  = np.mean((A_n - A_o)**2)
MSE_d  = np.mean((A_d - A_o)**2)

pchi2_n = np.sum((A_n - A_o)**2 / A_o)
pchi2_d = np.sum((A_d - A_o)**2 / A_o)


## Plot
fig,ax = plt.subplots(ncols=3, nrows=1, figsize=(25,35))
cmap = 'viridis'

ax[0].imshow(A_o, cmap=plt.get_cmap(cmap))
ax[1].imshow(A_n, cmap=plt.get_cmap(cmap))
ax[2].imshow(A_d, cmap=plt.get_cmap(cmap))

ax[0].set_title('Original')
ax[1].set_title(f'Noisy MSE = {MSE_n:0.1f}')
ax[2].set_title(f'Denoised MSE = {MSE_d:0.1f}')

print(f'Noisy    MSE = {MSE_n:0.1f}, pseudo-chi2 = {pchi2_n:0.1f}')
print(f'Denoised MSE = {MSE_d:0.1f}, pseudo-chi2 = {pchi2_d:0.1f}')

#plt.show()
plt.savefig('./img/collage.png', bbox_inches='tight')
plt.close()

print('Figures saved under ./img/')
