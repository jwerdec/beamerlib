import numpy as np

# Lineshape functions
gaussian = lambda x, x0, A, w, offset: A * np.exp(- (x-x0)**2/w**2) + offset
asymw = lambda x, x0, w0, alpha: 2*w0/(1+np.exp(alpha*(x-x0)))
asymgauss = lambda x, A, x0, w0, alpha, offset : \
    gaussian(x, x0, A, asymw(x, x0, w0, alpha), offset)
