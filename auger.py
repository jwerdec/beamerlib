import numpy as np
import matplotlib.pylab as plt
import scipy.signal as sig

class Auger(object):
    def __init__(self, filename, plot=True):
        self.__filename = filename
        self.__x, self.__y = np.loadtxt(filename, unpack=True)
        if plot:
            self.plot()

    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y

    def plot(self):
        plt.figure(1, figsize=(8,6), dpi=100)
        plt.plot(self.__x, self.__y, 'k-')
        plt.title(r'Auger electron spectrum (%s)' % self.__filename)
        plt.xlabel(r'$\mathrm{Kinetic\ Energy\ /\ eV}$')
        plt.ylabel(r'$\mathrm{d}N/\mathrm{d}E \ / \ \mathrm{a.u.}$')
        plt.show()

if __name__ == '__main__':
    TestAuger = Auger('testdata/Auger1.dat')
