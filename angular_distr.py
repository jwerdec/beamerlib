from __future__ import division 
from lmfit import LMFit
from polar_plot import SemiPolarPlot
from numpy import arange, cos, radians, max, min, array, linspace
from sys import stderr
from helper import normalize

class AngularDistribution(object):
    def __init__(self, Integrals, ICE=None, fit=True, plot=True,
		 Normalized=True, exclude=[]):
	if len(REMPISetups) != len(Integrals):
	    stderr.write('Error: Dimensions of Setups and ' + 
			     'Integrals do not match!')
            return
        self.__ICE = ICE
	self.__Setups = REMPISetups
	self.__Normalized = Normalized
	self.__Angles = array([setup.Angle for setup in self.__Setups])
	self.__excluded = []
        if ICE != None:
            for i in range(len(Integrals)):
                Integrals[i] = Integrals[i]/ICE(self.__Setups[i].ZR)
        self.__Integrals = array(Integrals)   
	self.__NormIntegrals = normalize(Integrals)
	self.__fit = None
	if fit:
	    self.fit(verbose=False, plot=False, exclude=exclude)
	if plot:
	    self.plot()

    def testfunc(self, theta, A, m, theta0):
        return A * cos(radians(theta-theta0))**m

    def fit(self, func=None, p0={'A':1, 'theta0': 3, 'm': 3}, verbose=True, 
	    plot=False, exclude=[], Normalized=True):
        self.__Normalized = Normalized
        self.__excluded = exclude
        if func == None:
            func = self.testfunc
        if Normalized:
            Integrals = self.__NormIntegrals
        else:
            Integrals = self.__Integrals
        ex = set(exclude)
        idx = set(range(len(Integrals)))
        idx = list(idx.difference(ex))
        Integrals = Integrals[idx]
        Angles = self.__Angles[idx]
        self.__fit = LMFit(func, Angles, Integrals, p0, verbose=verbose, 
                           plot=plot)

    def plot(self):
	if self.__Normalized:
	    y = self.__NormIntegrals
	else:
	    y = self.__Integrals
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
            polar_ax.plot(theta, self.Fit(theta), 'b-', label=label)
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
	plt.show(fig)

    @property
    def Angles(self):
	return self.__Angles

    @property
    def Setups(self):
	return self.__Setups

    @property
    def Integrals(self):
	return self.__Integrals

    @property
    def NormIntegrals(self):
	return self.__NormIntegrals

    @property
    def Fit(self):
	return self.__fit

# TESTCODE   
if __name__ == '__main__':
    from tof_analysis import SurfacePosMeas
    angles = linspace(-45, 45, 17) 
    signal = -cos(radians(angles+8))**12
    pos = arange(1,9.5,0.5)
    surfpos = arange(30,33,0.2)
    power = [0.3, 0.3, 0.285, 0.275, 0.255, 0.24, 0.23, 0.2, .179, .155,
             .129, .107, .082, .066, .047]
    SPM1 = SurfacePosMeas(surfpos, power, IRPos=7)
    Setups = [REMPISetup(SPM1, {'REMPI': 0, 'Center ZRM': 5.1, 'ZRM': zrm,
                                'Surface Y': 35}) for zrm in pos]
    AngDistr = AngularDistribution(Setups, signal)
