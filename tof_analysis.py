from __future__ import division
from matplotlib.pyplot import plot
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, cos, sin, log, exp, linspace
from matplotlib.ticker import MultipleLocator, NullFormatter
from lmfit import LMFit
from matplotlib.path import Path
import matplotlib.patches as patches 
import matplotlib.gridspec as gs
import scipy.constants as const
from scipy.special import erfc
from scipy.integrate import quad
from scipy.optimize import leastsq
import helper

class SurfacePosMeas(object):
    
    def __init__(self, Y=[], Power=[], IRPos=7):
        """
        ...
        """
        self.__Power = np.array(Power)
        self.Irel = self.__Power/self.__Power.max()*100
        self.Y = np.array(Y)
        self.__IRPos = IRPos
        self.__xlim = (self.Y.min()-0.01*self.Y.min(),
                       self.Y.max()+0.01*self.Y.max())
        self.__ylim = (0, self.Irel.max()+0.05*self.Irel.max())
        self.fit()
        
    def __func(self, y, a, Imax, y0):
        return Imax * erfc(a*(y-y0))
    
    @property
    def Y0(self):
        return self.__Y0
        
    @property
    def IRPos(self):
        return self.__IRPos
    
    def fit(self, p0={}):
        if p0=={}:
            p0={'a': 1, 'Imax': 50, 'y0': self.Y[len(self.Y)-1]}
        fit = LMFit(self.__func, self.Y, self.Irel, p0=p0, verbose=False,
                    plot=False)
        self.fit = fit
        self.__Y0 = fit.P['y0']
        self.__a = fit.P['a']
        x = linspace(self.__xlim[0], self.__xlim[1], 100)
        fig, axes = plt.subplots(figsize=(8,5), dpi=100)
        axes.set_title(r'IR Power vs. Surface Position')
        axes.plot(self.Y, self.Irel, 'o', color='w')
        axes.plot(x, fit(x), 'b-',
                  label=(r'$\mathrm{erfc}[ %.2f \cdot \left(Y_\mathrm{S} - %.2f \right)]$') % (self.__a, self.__Y0)
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

class TaggingSetup(object):
    """
    ...
    """
    
    # CONSTRUCTOR
    def __init__(self, Positions={'IR': (7,0), 'REMPI': 0, 'Center ZRM': 5.1,
                                  'ZRM': 3.5, 'Surface Y': None},
                 SurfacePosMeas=None, visualize=False):
        """
        ...
        """
        self.__SurfacePosMeas = SurfacePosMeas
        self.__Pos = Positions
        self.__distances = self.calc(Positions)
        if visualize == True:
            self.visualize()
            
    # PUBLIC ATTRIBUTES
    @property
    def Positions(self):
        return self.__Pos
    
    @property
    def SurfacePosMeas(self):
        return self.__SurfacePosMeas
    
    @property
    def distances(self):
        return self.__distances
        
    @property
    def FlightLength(self):
        return self.__l
    
    def calc(self, Pos):
        mmperdiv = 2
        mmperYs = 0.5
        IRx, IRy = Pos['IR']
        REMPIx = Pos['REMPI']
        Center = Pos['Center ZRM']
        REMPIy = (Center - Pos['ZRM'])*2.54/mmperdiv
        self.__Pos['REMPI'] = (REMPIx, REMPIy)
        if self.__SurfacePosMeas == None:
            self.__plot = self.__plotIRMPI
            distances = {'IR-MPI': sqrt((IRx - REMPIx)**2+(IRy - REMPIy)**2) *
                         mmperdiv}
            self.__l = distances['IR-MPI']
        else:
            self.__plot = self.__plotIRSurfMPI
            Surf = Pos['Surface Y']
            SurfRef = self.__SurfacePosMeas.Y0
            IRRef = self.__SurfacePosMeas.IRPos
            SurfPos = (SurfRef - Surf) * mmperYs / mmperdiv + IRRef
            self.__Surface = SurfPos
            self.__Pos['Surface Beamtool'] = SurfPos
            distances = {'IR-S': 0, 'MPI-S':0, 'IR-S-MPI': 0}
            distances['IR-S'] = sqrt((SurfPos - IRx)**2 + IRy**2) * mmperdiv
            distances['MPI-S'] = sqrt((SurfPos - REMPIx)**2 + REMPIy**2
                                      ) * mmperdiv
            self.__l = distances['IR-S'] + distances['MPI-S']
            distances['IR-S-MPI'] = self.__l
        self.__IRx = IRx
        self.__IRy = IRy
        self.__REMPIx = REMPIx
        self.__REMPIy = REMPIy
        return distances
        
    def SetRempi(self, horizontal, ZRM):
        self.__Pos['REMPI'] = horizontal
        self.__Pos['ZRM'] = ZRM
        self.__distances = self.calc(self.__Pos)
        
    def SetIR(self, horizontal, vertical):
        self.__Pos['IR'] = (horizontal, vertical)
        self.__distances = self.calc(self.__Pos)
        
    def show(self):
        string = ''
        for item in self.__distances:
            string += "\n%s = %f" % (item, self.__distances[item])
        print """
Positions (div on beamtool):
REMPI Beam: %f, %f
IR Beam:    %f, %f

Distances (mm): %s
""" %(self.__REMPIx, self.__REMPIy,\
     self.__IRx, self.__IRy, string)
    
    def visualize(self):
        """
        ...
        """
        fig, ax1, ax2 = self.__plotbasic()
        self.savefig = fig.savefig
        self.__plot(fig, ax1, ax2)
        plt.show(fig)
    
    # PRIVATE ATTRIBUTES
            
    def __plotbasic(self):
        _xlim = (11,-12)
        _ylim = (-6,7)
        _xticks = range(_xlim[0], _xlim[1], -2)
        _yticks = range(_ylim[0], _ylim[1], 2)
        
        fig = plt.figure(figsize=(12,5), dpi=100)
        grdspc = gs.GridSpec(1, 2, width_ratios=[2,1])
        ax1 = fig.add_subplot(grdspc[0])
        ax2 = fig.add_subplot(grdspc[1])
        newxax = ax1.twiny()
        newyax = ax1.twinx()
        
        ax1.set_title('Beam Positions')
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
    
        ax1.plot(self.__IRx, self.__IRy, 'ro', markersize=12)
        ax1.plot(self.__REMPIx, self.__REMPIy, 'bo', markersize=12)
        return fig, ax1, ax2
        
    def __plot(self, fig, ax1, ax2):
        pass
    
    def __plotIRSurfMPI(self, fig, ax1, ax2):
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
        
        string = r"""
Positions:
$(x,y)_\mathrm{UV} =\, (%2.2f,\ %2.2f)\ \mathrm{div}$
$(x,y)_\mathrm{IR} =\, (%2.2f,\ %2.2f)\ \mathrm{div}$
$x_\mathrm{S} =\, %2.2f \ \mathrm{div}$

Distances:
$d(\mathrm{IR-S}) = \, %2.2f \ \mathrm{mm}$
$d(\mathrm{S-UV}) = \, %2.2f \ \mathrm{mm}$
$d(\mathrm{IR-S-UV}) = \, %2.2f \ \mathrm{mm}$
""" %(self.__REMPIx, self.__REMPIy,\
          self.__IRx, self.__IRy,\
          self.__Surface, \
          self.__distances['IR-S'],\
          self.__distances['MPI-S'],\
          self.__distances['IR-S-MPI'])
        ax2.annotate(string, xy=(0,0), xycoords='axes fraction')
        
    def __plotIRMPI(self, fig, ax1, ax2):
        ax1.plot([self.__IRx, self.__REMPIx], [self.__IRy, self.__REMPIy],
                 color='black', linestyle='--')
        
        ax1.annotate(r'%2.2f mm' % self.__distances['IR-MPI'],
                     (0.5*(self.__IRx-abs(self.__REMPIx))+1,
                      0.5*(self.__REMPIy-abs(self.__IRy)) + 0.6),
                     rotation=np.degrees(np.arctan((self.__REMPIy-self.__IRy)/
                                                   (self.__IRx-self.__REMPIx))))
        string = r"""
Positions:
$(x,y)_\mathrm{UV} =\, (%2.2f,\ %2.2f)\ \mathrm{div}$
$(x,y)_\mathrm{IR} =\, (%2.2f,\ %2.2f)\ \mathrm{div}$

Distance:
$d(\mathrm{IR-UV}) = \, %2.2f \ \mathrm{mm}$
""" %(self.__REMPIx, self.__REMPIy, self.__IRx, self.__IRy, \
        self.__distances['IR-MPI'])
        ax2.annotate(string, xy=(0,0), xycoords='axes fraction')

class RawTOFData(object):
    """
    ...
    """
    
    # CONSTRUCTOR

    def __init__(self, filename):
        self.__filename = filename
        self.__RawTable = np.loadtxt(filename)
        self.__numpoints = len(self.__RawTable)
        self.__numcols = len(self.__RawTable[0])

    @property
    def TOF(self):
        return self[0]

    @property
    def Signal(self):
        return self[1]

    @property
    def Baseline(self):
        return self[2]

    @property 
    def sorted(self):
        return self.__sort(self.__RawTable)

    @property
    def filename(self):
        return self.__filename

    @property
    def numcols(self):
        return self.__numcols

    @property
    def numpoints(self):
        return self.__numpoints

    def __getitem__(self, k):
        return self.__RawTable[:, k]

    def __sort(self, datatable):
        return np.transpose(datatable[datatable[:,0].argsort()])

    def plot(self):
        """
        ...
        """
        fig, axes = plt.subplots(figsize=(8,6), dpi=100)
        axes.set_title(self.__filename)
        axes.set_xlabel(r'$\mathrm{Time\ of\ Flight\, /\, ns }$')
        axes.set_ylabel(r'$\mathrm{Intensity}$')
        for i in range(1, self.__numcols, 1):
            plot(self[0], self[i], '-', label='Column %d' % (i+1))
        axes.legend(loc='best')
        plt.show(fig)

class RempiTOA(object):
    """
    ...
    """
    def __init__(self, filename, numbasecorr=20, Normalize=True, plot=True):
        self.__filename = filename
        self.__numbasecorr = numbasecorr
        RawData = RawTOFData(filename)
        self.__RawData = RawData
        idx = RawData.TOF.argsort()
        self.__TOF = RawData.TOF[idx]/1000
        Signal = (RawData.Baseline - RawData.Signal)[idx]
        self.__Signal = helper.subtract_baseline(Signal, left=self.__numbasecorr)
        if Normalize:
            self.__Signal = helper.normalize(self.__Signal)
            self.__ylabel = r'$\mathrm{Normalized\ REMPI\ Signal}$'
        else: 
            self.__ylabel = r'$\mathrm{REMPI\ Signal}$'
        if plot:
            self.plot()
    
    def plot(self):
        fig, axes = plt.subplots(figsize=(8,5), dpi=100)
        axes.set_title(r'$\mathrm{%s}$' % self.__filename)
        axes.plot(self.__TOF, self.__Signal, 'bo')
        axes.set_ylabel(self.__ylabel)
        axes.set_xlim((self.__TOF[0], self.__TOF[-1]))
        axes.set_xlabel(r'$\mathrm{Time\ of\ Arrival\ /\ \mu s}$')
        self.fig = fig
        self.axes = axes
        self.savefig = self.fig.savefig
        plt.show(fig)
        
    # PUBLIC ATTRIBUTES   

    @property
    def RawData(self):
        return self.__RawData

    @property
    def TOF(self):
        return self.__TOF

    @property
    def Setup(self):
        return self.__setup

class TaggingTOF(object):
    """
    ...
    """
    # CONSTRUCTOR
    def __init__(self, Setup, filename, IRDelay,
                 mode='Flux_vs_TOF', mass=28, 
                 verbose=False, plot=True, fit=True, numbasecorr=20,
                 Normalize=True):
        """
        ...
        """
        self.__mass = mass
        self.__filename = filename
        self.__IRDelay = IRDelay
        self.__setup = Setup
        self.__l = Setup.FlightLength
        self.__offset = 400
        self.__Normalize = Normalize
        self.__numbasecorr = numbasecorr
        self.__RawData = RawTOFData(filename)
        TOF, Signal = self.__treat_data(self.__RawData)
        v, E = self.__invert(TOF)
        Flux = self.__density_to_flux(Signal, v)
        if Normalize:
            Flux = [helper.normalize(flux, 5) for flux in Flux]
        self.__TOF = (TOF, Flux[0])
        self.__v = (v, Flux[1])
        self.__E = (E, Flux[2])
        self.set_mode(mode)
        self.__fit = None
        self.__moments = None
        if fit:
            self.fit(verbose=verbose)
        if plot:
            self.plot()

    # PRIVATE ATTRIBUTES    

    def __treat_data(self, RawData):     
        idx = RawData.TOF.argsort()
        TOF = RawData.TOF[idx]/1000 -self.__IRDelay + self.__offset
        Signal = (RawData.Baseline - RawData.Signal)[idx]
        Signal = helper.subtract_baseline(Signal, left=self.__numbasecorr)
        idx, = np.where(TOF > 0)
        numomitt = len(TOF) - len(idx)
        if numomitt > 0:
            print 'Warning: Omitting the first %d datapoints because TOF <= 0!' % numomitt
        return TOF[idx], Signal[idx]

    def __invert(self, TOF):
        v = self.__l*1000 / TOF
        E = 0.5 * self.__mass * v**2 * const.u/const.eV 
        return v, E
        
    def __density_to_flux(self, Signal, v):
        TOFFlux = Signal*v
        vFlux = Signal*self.__l/v
        EFlux = vFlux/(v*self.__mass)
        return TOFFlux, vFlux, EFlux
        
    def __FluxTOFfit(self, x, F0, alpha, x0):
        return F0 * self.__l**4/x**4 * exp(- (self.__l/alpha)**2 *
                                            (1/x - 1/x0)**2)
    
    def __FluxVfit(self, x, F0, alpha, x0):
        return F0 * x**3 * exp(-( (x-x0) /alpha)**2)
    
    def __FluxEfit(self, x, F0, alpha, x0):
        return F0 * x/self.__mass**2  * exp(- 2/(self.__mass*alpha**2) * 
                                              (sqrt(x)- sqrt(x0))**2 )
    
    def __DensityVfit(self, x, F0, alpha, x0):
        return F0 * x**2 * exp(- ((x-x0)/alpha)**2)
    
    def __estimate_p0(self, x, y):
        maxpos = np.where(y == y.max())[0][0]
        x0 = x[maxpos]
        F0 = y[maxpos]
        i = maxpos
        while y[i] > F0/2:
            i += 1
        right = i
        i = maxpos
        while y[i] > F0/2:
            i -= 1
        left = i
        alpha = x[right] - x[left]
        if self.__mode == 'Flux_vs_E':
            F0 = y[maxpos]*self.__mass**2/x[maxpos]
        return {'F0':F0, 'x0':x0, 'alpha':alpha}
        
    def __get_moments(self):
        if self.__moments == None:
            x, _ = self.__fitargs['data']
            zeroth, first, second = helper.Moments(self.__fit, n=[0,1,2],
                                                   limits=(0, x[-1]*10))
            self.__moments = ([zeroth[0], first[0], second[0]],
                              first[0]/zeroth[0], second[0]/zeroth[0])
        return self.__moments
    
    def __call__(self, x):
        return self.__fit(x)
              
    # PUBLIC ATTRIBUTES   

    @property
    def Mass(self):
        return self.__mass

    @property
    def RawData(self):
        return self.__RawData

    @property
    def TOF(self):
        return self.__TOF
        
    @property 
    def v(self):
        return self.__v
        
    @property 
    def E(self):
        return self.__E

    @property
    def Fit(self):
        return self.__fit

    @property
    def Setup(self):
        return self.__setup
        
    @property
    def Mean(self):
        return self.__get_moments()[1]
      
    @property
    def StdDev(self):
        return self.__get_moments()[2]
        
    @property
    def Moments(self):
        return self.__get_moments()[0]

    def get_fitfunc(self):
        return self.__fitargs['func']

    def get_values(self, limits=()):
        """
        Returns 2d array with the x and y data in the current mode 
        Limits: Tuple of x values definding a range of values that are returned
        """
        x, y = self.__fitargs['data']
        tab = np.array([x, y])
        return helper.extract_from_Table(tab, limits=limits)
        
    def set_mode(self, mode):
        pltset = {
            'xlabel': r'$\mathrm{Time\ of\ Flight\, /\, \mu s }$',
            'ylabel': r'',
            }
        if self.__Normalize:
            pltset['ylabel'] = r'$\mathrm{Normalized}$'
        if mode == 'Flux_vs_TOF':
            self.__fitargs = {'func': self.__FluxTOFfit, 'data': self.__TOF}
            pltset['ylabel'] = pltset['ylabel'] + r' $\mathrm{Flux}$'
        elif mode == 'Flux_vs_v':
            self.__fitargs = {'func': self.__FluxVfit, 'data': self.__v}
            pltset['ylabel'] = pltset['ylabel'] + r' $\mathrm{Flux}$'
            pltset['xlabel'] = r'$v \, / \, \mathrm{m}\, \mathrm{s}^{-1}$'
        elif mode == 'Flux_vs_E':
            self.__fitargs = {'func': self.__FluxEfit, 'data': self.__E}
            pltset['ylabel'] = pltset['ylabel'] + r' $\mathrm{Flux}$'
            pltset['xlabel'] = r'$E \, / \, \mathrm{eV}$'
        else:
            raise Exception('Error: Mode %s unknown!' % mode)
            return None
        self.__mode = mode
        self._plotsettings = pltset

    def fit(self, func=None, p0={}, verbose=False, plot=False):
        self.__moments = None
        if func == None:
            func = self.__fitargs['func']
        x, y = self.__fitargs['data']
        if p0 == {}:
            p0 = self.__estimate_p0(x, y)
        self.__fit = LMFit(func, x, y, p0, plot=plot, verbose=verbose)
       
    def plot(self):
        fig, axes = plt.subplots(figsize=(8,5), dpi=100)
        axes.set_title(r'$\mathrm{%s}$' % self.__filename)
        x, y = self.__fitargs['data']
        pltset = self._plotsettings
        xlim = ( 0.928*x[0], 1.02*x[-1] )
        if xlim[0] > xlim[1]:
            xlim = (xlim[1], xlim[0])
        axes.plot(x, y, 'bo')
        xspace = linspace(xlim[0], xlim[1], 1000)
        axes.set_xlim(xlim)
        if self.__fit:
            axes.plot(xspace, self.Fit(xspace), 'r-')   
        axes.set_ylabel(pltset['ylabel'])
        axes.set_xlabel(pltset['xlabel'])
        self.fig = fig
        self.axes = axes
        self.savefig = self.fig.savefig
        plt.show(fig)
  
class MultiTOF(object):
    """
    ...
    """
         
    #CONSTRUCTOR
    
    def __init__(self, filenames=[], Descriptions=[], Setups=[],
                 IRDelays=[], mode='Flux_vs_TOF', plot=True, xlim=()):
        """
        ...
        """
        self.__n = len(filenames)
        self.__nrange = range(self.__n)
        self.__filenames = filenames
        if Descriptions == []:
            self.__descriptions = ['' for _ in self.__nrange]
        else:
            self.__descriptions = Descriptions
        self.__xlim = xlim
        if type(Setups) == type([]):
            self.__Setups = Setups
        else:
            self.__Setups = [Setups for _ in self.__nrange]
        if type(IRDelays) == type([]):
            self.__IRDelays = IRDelays
        else:
            self.__IRDelays = [IRDelays for _ in self.__nrange]
        self.__mode = mode
        self.__TOFData = [TaggingTOF(self.__Setups[i], self.__filenames[i],
                                      self.__IRDelays[i], mode=mode,
                                     plot=False) for i
                          in self.__nrange]    
        if plot:
            self.plot()

    def __getitem__(self, k):
        return self.__TOFData[k]
            
    def set_xlim(self, value=()):
        self.__xlim = value

    @property
    def Setups(self):
        return self.__Setups

    @property
    def RawData(self):
        return [i.RawData for i in self.__TOFData]

    @property
    def Fig(self):
        return self.__fig
        
    def plot(self, xlim=(), Descriptions=[]):
        """
        ...
        """
	if Descriptions != []:
	    self.__descriptions = Descriptions
        if xlim!=():
            self.__xlim = xlim
        numrows = (self.__n + 1)//2
        fig, axes = plt.subplots(numrows, 2, dpi=100,\
                                     figsize=(16, 5*numrows), sharey=True)
        ax = fig.add_subplot(111)
        ax.set_frame_on(False)
        ax.yaxis.set_major_formatter(NullFormatter())
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])
        fig.subplots_adjust(wspace=0, hspace=0)
        axes = [item for sublist in axes for item in sublist]

        x = [None for _ in self.__TOFData]
        y = [None for _ in self.__TOFData]
        for i in range(len(x)):
            x[i], y[i] = self.__TOFData[i].get_values()
        if self.__xlim == ():
            xmax = max(x[0])
            xmin = min(x[0])
            for xvec in x:
                if max(xvec) > xmax:
                    xmax = max(xvec)
                if min(xvec) < xmin:
                    xmin = min(xvec)
            xlim = (0.98*xmin , 1.02*xmax)
        else:
            xlim = self.__xlim

        for axis in axes[::2]:
            axis.set_ylabel(self.__TOFData[0]._plotsettings['ylabel'])
        for axis in axes[-2:]:
            axis.set_xlabel(self.__TOFData[0]._plotsettings['xlabel'])
        for axis in axes[:-2]:
            axis.xaxis.set_ticklabels([])

        xspace = linspace(xlim[0], xlim[1], 1000)
        
        for i in self.__nrange:
            axes[i].set_title(r'$\mathrm{%s}$' % self.__filenames[i])
            axes[i].set_ylim((-0.19,1.1))
            axes[i].set_xlim(xlim)
            axes[i].plot(x[i], y[i], 'bo')
            axes[i].plot(xspace, self.__TOFData[i](xspace), 'r-')
            axes[i].annotate(self.__descriptions[i], xy=(0.8, 0.8),\
                                 xycoords='axes fraction')
            
        self.__axes = axes
        self.__fig = fig
        plt.show(fig)

# TESTCODE
if __name__ == '__main__':
    # Determination of Surface Position
    Y = np.linspace(33, 35, 11)
    Power = [11.3, 11.0, 11.0, 10.7, 10.7, 10.6, 10.7, 10.3, 8.88, 5.16, 1.91]
    test = SurfacePosMeas(Y, Power)

    # Setup
    Positions={'IR': (7,0), 'REMPI': 0,
               'Center ZRM': 5.1, 'ZRM': 3.5, 'Surface Y': None}
    testsetup = TaggingSetup(Positions, visualize=True)
    Positions={'IR': (7,0), 'REMPI': 0, 'Center ZRM': 5.1,
               'ZRM': 3.5, 'Surface Y': 33.8}
    testsetup2 = TaggingSetup(Positions, SurfacePosMeas=test, visualize=True)

    # Rempi TOA
    TestTOA = RempiTOA('testdata/1208005.dat')

    # Tagging TOF
    TestTOF = TaggingTOF(testsetup2, 'testdata/1308033.dat', 606, verbose=True,
                         mode='Flux_vs_E')
    print TestTOF.Mean, TestTOF.StdDev, TestTOF.Moments
