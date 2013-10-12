from __future__ import division
from matplotlib.path import Path
import matplotlib.patches as patches 
from numpy import arctan, sqrt, degrees
import matplotlib.pylab as plt
import matplotlib.gridspec as gs
from matplotlib.ticker import MultipleLocator, NullFormatter
from scipy.special import erfc
from lmfit import LMFit
from helper import normalize
import numpy as np
from general_functions import gaussian

class IncidentBeamPos(object):
    def __init__(self, Z, Integrals, Z0, p0=None, plot=False):
        Integrals = Integrals/np.array(Integrals).max()
        Z = np.array(Z)
        if p0 == None:
            p0 = {'A':1, 'x0': 0.2, 'w':0.1}
        self.__fit = LMFit(self.__func, Z0 - Z, Integrals,
                           p0=p0, verbose=False)
        self.__Z = Z
        self.__Z0 = Z0
        self.__Integrals = Integrals
        self.__offset = self.__fit.P['x0']
        if plot:
            self.plot()

    @property
    def offset(self):
        return self.__offset

    def __func(self, x, x0, A, w):
        return gaussian(x, x0, A, w, 0)

    def __repr__(self):
        return "IncidentBeamPos(%s, %s, %f)" %\
            (list(self.__Z), list(self.__Integrals), self.__Z0)

    def plot(self):
        plt.plot(self.__Z0 - self.__Z, self.__Integrals, 'bo')
        xspace = np.linspace(self.__Z0 - self.__Z[-1] - 1,
                             self.__Z0 - self.__Z[0] + 1, 200)
        plt.plot(xspace, self.__fit(xspace), 'b-',
                 label=r'$A \, \exp \left( - \frac{(x - %.3f)^2}{%.3f^2} \right)$' % (self.offset, self.__fit.P['w'])
                 )
        plt.xlabel(r'Vetical Offset $Z_0 - Z\ / \ 0.1$"')
        plt.ylabel(r'REMPI Signal / a.u.')
        plt.legend(loc='best')
        plt.show()


class SurfacePosMeas(object):
    def __init__(self, Y, Power, BeamPos=7):
        self.__Power = np.array(Power)
        self.Irel = self.__Power/self.__Power.max()*100
        self.Y = np.array(Y)
        self.__BeamPos = BeamPos
        self.__xlim = (self.Y.min()-0.01*self.Y.min(),
                       self.Y.max()+0.01*self.Y.max())
        self.__ylim = (0, self.Irel.max()+0.05*self.Irel.max())
        self.fit()
        self.plot()
        
    def __func(self, y, a, Imax, y0):
        return Imax * erfc(a*(y-y0))
    
    @property
    def Ys0(self):
        return self.__Ys0

    @property
    def Fit(self):
        return self.__fit
    
    def fit(self, p0={}):
        if p0=={}:
            p0={'a': 1, 'Imax': 50, 'y0': self.Y[len(self.Y)-1]}
        self.__fit = LMFit(self.__func, self.Y, self.Irel, p0=p0, verbose=False,
                    plot=False)
        self.__Ys0 = self.__fit.P['y0']*0.5 + self.__BeamPos*2
        self.__a = self.__fit.P['a']

    def plot(self):
        x = np.linspace(self.__xlim[0], self.__xlim[1], 100)
        fig, axes = plt.subplots(figsize=(8,5), dpi=100)
        axes.set_title(r'IR Power vs. Surface Position')
        axes.plot(self.Y, self.Irel, 'o', color='w')
        axes.plot(x, self.__fit(x), 'b-',
                  label=(r'$\mathrm{erfc}[ %.2f \cdot \left(Y_\mathrm{S} - %.2f \right)]$') % (self.__a, self.__Ys0)
                  )
        axes.set_ylabel(r'$I/I_\mathrm{max}\ [\%]$')
        axes.set_xlim(self.__xlim)
        axes.set_ylim(self.__ylim)
        axes.xaxis.set_major_locator(MultipleLocator(0.5))
        axes.xaxis.set_minor_locator(MultipleLocator(0.1))
        axes.yaxis.set_major_locator(MultipleLocator(10))
        axes.grid(which='minor', axis='x')
        axes.grid(which='major', axis='x')
        axes.grid(which='major', axis='y')
        axes.legend(loc='best')
        plt.show(fig)


class REMPISetup(object):
    def __init__(self, Positions={'REMPIx': 0, 'Z0': 5.1, 'Z': 3.5,
                                  'Z Offset': 0, 'Ys': 33, 'Ys0': 30},
		 visualize=False):
        self.__Pos = Positions
        self.__Surf = self.__surf()
        self.__REMPIx = self.__tomm(self.__Pos['REMPIx'])
        self.__REMPIy = self.__vert(self.__Pos['Z'])
        self.__angle = self(self.__Pos['Z'])
        self.__l = sqrt((self.__hor())**2 +
                        (self.__REMPIy - self.__Pos['Z Offset']*2.54)**2)
        if visualize == True:
            self.visualize()

    def __call__(self, Z):
        return self.calc_angle(Z)

    @property
    def Z0(self):
        return self.__Pos['Z0']

    @property
    def Positions(self):
        return self.__Pos
        
    @property
    def FlightDistance(self):
        return self.__l
        
    @property
    def Angle(self):
        return self.__angle

    def __tomm(self, div):
        return div*2

    def __todiv(self, mm):
        return mm/2

    def __vert(self, Z):
        return (self.__Pos['Z0'] - Z)*2.54

    def __surf(self):
        return (self.__Pos['Ys0'] - self.__Pos['Ys']*0.5)

    def __hor(self):
        return self.__Surf - self.__REMPIx

    def calc_angle(self, Z):
        return degrees(arctan((self.__vert(Z) -
                               self.__Pos['Z Offset']*2.54)/self.__hor()))
        
    def visualize(self):
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
        soffset = 0.5
        Surf = self.__todiv(self.__Surf)
        verts = [
             (Surf, -height/2+soffset),      # right bottom
             (Surf, height/2+soffset),       # right top
             (Surf+width, height/2+soffset), # left top
             (Surf+width, -height/2+soffset),# left bottom
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
        beamspace = np.linspace(Surf, _xlim[1], 200)
        ax1.plot(beamspace, [self.__Pos['Z Offset']*2.54/2 for _ in beamspace],
                 'g-')
        ax1.plot([self.__todiv(self.__Surf), self.__todiv(self.__REMPIx)],
                 [self.__Pos['Z Offset']*2.54/2, self.__todiv(self.__REMPIy)],
                 'g-')
        ax1.plot(self.__todiv(self.__REMPIx),
                 self.__todiv(self.__REMPIy), 'bo', markersize=12)
        string = r"""
        Positions:
$(x,y)_\mathrm{UV} =\, (%2.2f,\ %2.2f)\ \mathrm{mm}$
$x_\mathrm{S} =\, %2.2f \ \mathrm{mm}$

Angle:
$\theta_\mathrm{f} = \, %2.2f ^\circ$

Flight distance:
$L = \, %2.2f \ \mathrm{mm}$
""" %(self.__REMPIx, self.__REMPIy, self.__Surf, self.__angle, self.__l)
        ax2.annotate(string, xy=(0,0.5), xycoords='axes fraction')
        return fig, ax1, ax2


if __name__ == '__main__':
    obj = IncidentBeamPos(np.arange(4, 5.6, 0.2),
                          [0.00036014073169634458, 0.0065634184421959766,
                           0.0095181994254066316, 0.01066251499384108,
                           0.0089403275416902097, 0.0072330786125905013, 
                           0.0010042068897813349, 0.0001348028515445096],
                          5.0, plot=True)
    string = repr(obj)
    print string
    newobj = eval(string)
    print "The new offset is %f." % newobj.offset

    surfpos = np.arange(29,31.2,0.2)
    power = [1.04, 1.04, 1.03, 1.02, 0.99, 0.98, 0.97,0.64, 0.24, 0.036, 0.021]
    SPM = SurfacePosMeas(surfpos, power, BeamPos=7)

    Setup = REMPISetup(Positions={'REMPIx': 0, 'Z0': 5.1, 'Z': 5,
                                  'Z Offset': newobj.offset, 'Ys': 35,
                                  'Ys0': SPM.Ys0},
                       visualize=True)
