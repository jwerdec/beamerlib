from __future__ import division 
from lmfit import LMFit
from polar_plot import SemiPolarPlot
from numpy import arange, cos, radians, max, min, array, linspace, savetxt
from helper import normalize
from geometry import REMPISetup
import matplotlib.pyplot as plt
from scipy.integrate import quad

class AngularDistribution(object):
    def __init__(self, Setup, Integrals, Z=arange(1,9.1, 0.5),
                 ICE=None, fit=True, plot=True,
		 Normalized=True, exclude=[]):
        self.__ICE = ICE
	self.__Angles = Setup(Z)
	self.__excluded = []
        if ICE != None:
            for i in range(len(Integrals)):
                Integrals[i] = Integrals[i]/ICE(Setup.Z0 - Z[i])
        self.__Integrals = array(Integrals)   
	if Normalized:
	    self.__Integrals /= self.__Integrals.max()
	self.__fit = None
	if fit:
	    self.fit(verbose=False, plot=False, exclude=exclude)
	if plot:
	    self.plot()

    def testfunc(self, theta, A, m, theta0):
        return A * cos(radians(theta-theta0))**m

    def fit(self, func=None, p0={'A':1, 'theta0': 3, 'm': 3}, verbose=True, 
	    plot=False, exclude=[]):
        self.__excluded = exclude
        Integrals = self.__Integrals
        if func == None:
            func = self.testfunc
        ex = set(exclude)
        idx = set(range(len(Integrals)))
        idx = list(idx.difference(ex))
        Integrals = Integrals[idx]
        Angles = self.__Angles[idx]
        self.__fit = LMFit(func, Angles, Integrals, p0, verbose=verbose, 
                           plot=plot)

    def plot(self):
	y = self.__Integrals / self.Fit.P['A']*0.9
	fig = plt.figure(1, figsize=(8, 8))
	fig.clf()
	ax1, polar_ax = SemiPolarPlot(fig)
        exclude = self.__excluded
        idx = list(set(range(len(y))).difference(set(exclude)))
	polar_ax.plot(self.__Angles[idx], y[idx], 'bo',
		      label=r'$\mathrm{Experimental\ Data}$')
	theta = linspace(-90, 90, 300)
        if self.__fit is not None:
            label = r'$%1.2f\ \cos^{%2.1f}(\theta-%1.2f^\circ)$' %\
                (self.Fit.P['A'], self.Fit.P['m'], self.Fit.P['theta0'])
            polar_ax.plot(theta, 0.9*self.Fit(theta)/self.Fit.P['A'],
                          'b-', label=label)
            label = r'$\cos(\theta)$' % self.Fit.P['theta0'] 
            polar_ax.plot(theta, cos(radians(theta)),
                          'k--', label=label)
	if len(self.__excluded) > 0:
	    polar_ax.plot(self.__Angles[self.__excluded],
			  y[self.__excluded], 'ro',
			  label=r'$\mathrm{Excluded\ Data\ Points}$')
	ax1.legend(loc='best', numpoints=1)
	ax1.set_title(r"""Angular Distribution""")
	self.fig = fig
	self.axis1 = ax1
	self.polar_axis = polar_ax
	self.savefig = fig.savefig
        self.__y = y
	plt.show(fig)

    def savetxt(self, filenameprefix):
        b = self.Fit.P['theta0']
        savetxt(filenameprefix + '_raw.dat',
                array([self.__Angles, self.__y]).T)
        savetxt(filenameprefix + '_fit.dat',
                array([linspace(-90+b, 90-b, 300),
                       0.9*self(linspace(-90+b, 90-b, 300))/self.Fit.P['A']]).T)

    def __call__(self, theta):
        return 0.9*self.Fit(theta)/self.__fit.P['A']

    @property
    def Angles(self):
	return self.__Angles

    @property
    def Integrals(self):
	return 0.9*self.__Integrals/self.__fit.P['A']

    @property
    def Fit(self):
	return self.__fit

    @property
    def Area(self):
        b = self.Fit.P['theta0']
        return quad(self, -90+b, 90+b)[0]

    @property
    def FitParams(self):
        return self.Fit.P['m'], self.Fit.P['theta0']


# TESTCODE   
if __name__ == '__main__':
    ZOffset = 0.394666
    Ys0 = 29.232596443386047
    Setup = REMPISetup(Positions={'REMPIx': 0, 'Z0': 5.1, 'Z': 5,
                                  'Z Offset': ZOffset, 'Ys': 35,
                                  'Ys0': Ys0}, visualize=True)
    angles = Setup(arange(1, 9.1, 0.5))
    signal = cos(radians(angles-8))**12
    _ = AngularDistribution(Setup, signal, exclude=[3,4])
