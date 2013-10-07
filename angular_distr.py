from __future__ import division 
from lmfit import LMFit
from polar_plot import SemiPolarPlot
from numpy import arange, cos, radians, max, min, array, linspace
from sys import stderr

class AngularDistribution(object):
    def __init__(self, REMPISetups, Integrals, ICE=None, fit=True, plot=True,
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
	self.__NormIntegrals = self.__norm(Integrals)
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

    def __norm(self, x):
	    return x/max(x)

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

from matplotlib.path import Path
import matplotlib.patches as patches 
from numpy import arctan, sqrt, degrees
import matplotlib.pylab as plt
import matplotlib.gridspec as gs

class REMPISetup(object):
    def __init__(self, SurfacePosMeas,
		 Positions={'REMPI': 0, 'Center ZRM': 5.1, 'ZRM': 3.5,
			    'Surface Y': None},
		 visualize=False):
        self.__SurfacePosMeas = SurfacePosMeas
        self.__Pos = Positions
        self.__l, self.__angle = self.calc(Positions)
        if visualize == True:
            self.visualize()
    
    # PUBLIC ATTRIBUTES

    @property
    def ZR(self):
        return self.__Pos['Center ZRM'] - self.__Pos['ZRM']

    @property
    def Positions(self):
        return self.__Pos
    
    @property
    def SurfacePosMeas(self):
        return self.__SurfacePosMeas
        
    @property
    def FlightLength(self):
        return self.__l
        
    @property
    def Angle(self):
        return self.__angle
    
    def calc(self, Pos):
        mmperdiv = 2
        mmperManipulatorDiv = 0.5
        REMPIx = Pos['REMPI']
        Center = Pos['Center ZRM']
        Surf = Pos['Surface Y'] 
        REMPIy = (Center - Pos['ZRM'])*2.54/mmperdiv
        self.__Pos['REMPI'] = (REMPIx, REMPIy)
        SurfRef = self.__SurfacePosMeas.Y0
        IRRef = self.__SurfacePosMeas.IRPos
        SurfPos = (SurfRef - Surf)*mmperManipulatorDiv / mmperdiv + IRRef
        self.__Surface = SurfPos
        self.__Pos['Surface Beamtool'] = SurfPos
        l = sqrt((SurfPos - REMPIx)**2 + REMPIy**2)*mmperdiv
        angle = degrees(arctan(REMPIy /(SurfPos - REMPIx)))
        self.__REMPIx = REMPIx
        self.__REMPIy = REMPIy
        return l, angle
        
    def SetRempi(self, horizontal, ZRM):
        self.__Pos['REMPI'] = horizontal
        self.__Pos['ZRM'] = ZRM
        self.__distances = self.calc(self.__Pos)
        
    def visualize(self):
        """
        ...
        """
        fig, ax1, ax2 = self.__plot()
        self.savefig = fig.savefig
        plt.show(fig)
        
    def __plot(self):
        _xlim = (10,-6)
        _ylim = (-4,5)
        _xticks = range(_xlim[0], _xlim[1], -1)
        _yticks = range(_ylim[0], _ylim[1], 1)
        
        fig = plt.figure(figsize=(12,5), dpi=100)
        grdspc = gs.GridSpec(1, 2, width_ratios=[2,1])
        ax1 = fig.add_subplot(grdspc[0])
        ax2 = fig.add_subplot(grdspc[1])
        newxax = ax1.twiny()
        newyax = ax1.twinx()
        
        ax1.set_title('Beam Position')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.xaxis.set_ticks_position('bottom')
        ax1.spines['bottom'].set_position(('data',0))
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['left'].set_position(('data',0))
        ax1.set_xlim(_xlim)
        ax1.set_ylim(_ylim)
        ax1.xaxis.set_ticks(_xticks)
        ax1.xaxis.set_ticklabels([])
        ax1.yaxis.set_ticks(_yticks)
        ax1.yaxis.set_ticklabels([])
        
        ax2.set_frame_on(False)
        for spine in ax2.spines.values():
            spine.set_visible(False)
        ax2.patch.set_visible(False)
        ax2.xaxis.set_ticks([])
        ax2.yaxis.set_ticks([])
        
        for ax in [newxax, newyax]:
            ax.set_frame_on(True)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.patch.set_visible(False)
            
        newxax.spines['left'].set_visible(False)
        newxax.xaxis.set_ticks_position('bottom')
        newxax.xaxis.set_label_position('bottom')
        newxax.spines['bottom'].set_position(('outward', 10))
        newxax.set_xlim(_xlim)
        newxax.xaxis.set_ticks(_xticks)
        newxax.set_xlabel('Position on Beamtool / div')
        
        newyax.spines['bottom'].set_visible(False)
        newyax.yaxis.set_ticks_position('left')
        newyax.yaxis.set_label_position('left')
        newyax.spines['left'].set_position(('outward', 10))
        newyax.set_ylim(_ylim)
        newyax.yaxis.set_ticks(_yticks)
        newyax.set_ylabel('Position on Beamtool / div')
        
        width = 0.2
        height = 2
        verts = [
             (self.__Surface, -height/2),      # right bottom
             (self.__Surface, height/2),       # right top
             (self.__Surface-width, height/2), # left top
             (self.__Surface-width, -height/2),# left bottom
             (0., 0.) #ignored
             ]
        codes = [
             Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY
             ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='gold', lw=0)
        ax1.add_patch(patch)
        ax1.plot(self.__REMPIx, self.__REMPIy, 'bo', markersize=12)
        string = r"""
        Positions:
$(x,y)_\mathrm{UV} =\, (%2.2f,\ %2.2f)\ \mathrm{div}$
$x_\mathrm{S} =\, %2.2f \ \mathrm{mm}$

Angle:
$\theta_\mathrm{f} = \, %2.2f ^\circ$
""" %(self.__REMPIx, self.__REMPIy, self.__Surface, self.__angle)
        ax2.annotate(string, xy=(0,0.5), xycoords='axes fraction')
        return fig, ax1, ax2

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
