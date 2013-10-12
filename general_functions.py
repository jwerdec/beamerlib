import numpy as np

gaussian = lambda x, x0, A, w, offset: A * np.exp(- (x-x0)**2/w**2) + offset
