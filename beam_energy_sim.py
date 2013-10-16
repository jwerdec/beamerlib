import scipy.constants as const
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl

class Gas(object):
    def __init__(self, Name, Mass, Cp):
        self.__name = Name
        self.__M = Mass
        self.__Cp = Cp
        
    @property
    def Name(self):
        return self.__name
        
    @property
    def Mass(self):
        return self.__M
        
    @property
    def Cp(self):
        return self.__Cp

CO = Gas(r'CO', 28, 7/2*const.R)
H2 = Gas(r'H$_2$', 2.016, 28.836)
NO = Gas(r'NO', 30, 7/2*const.R)
Ar = Gas(r'Ar', 40, 20.786)
Ne = Gas(r'Ne', 20.179, 20.786)
He = Gas(r'He', 4, 20.786)
N2 = Gas(r'N$_2$', 28.01, 29.124)
HCl = Gas(r'HCl', 36, 7/2)

def BeamEnergiesPlot(T, Gas, Carriers=[H2, He, Ne]):
    fig, axes = plt.subplots(figsize=(16,9), dpi=100)
    x = np.linspace(0, 1, 500)
    for carrier in Carriers:
        axes.plot(x, SeededBeam(T, Gas, carrier)(x), linewidth=2,
                  label=r'$E_\mathrm{kin}$ of %s in %s' % (Gas.Name, carrier.Name))
    axes.set_title(r'Kinetic energy of a %s molecular beam seeded in various gases at %d K' % (Gas.Name, T))
    axes.set_xlabel(r'$x$(%s)' % Gas.Name)
    axes.set_ylabel(r'$E_\mathrm{kin} \, / \, \mathrm{eV}$ (solid line)')
    ax2 = axes.twinx()
    for carrier in Carriers:
        ax2.plot(x, SeededBeam(T, Gas, carrier).calc_v(x), linewidth=2,
                 ls='-.', label=r'$v$ of %s in %s' % (Gas.Name, carrier.Name))
    ax2.set_ylabel(r'$v\ / \ \mathrm{m} \ \mathrm{s}^{-1}$ (dashed dotted line)')
    axes.set_xlim((0,.5))
    axes.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
    axes.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.01))
    axes.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05))
    axes.grid(which='minor', axis='x')
    axes.grid(which='major', axis='x')
    axes.grid(which='minor', axis='y')
    axes.grid(which='major', axis='y')
    axes.legend(loc='best')
    return fig, axes

class SeededBeam(object):
    def __init__(self, T, Gas, Carrier, plot=False):
        self.__T = T
        self.__Gas = Gas
        self.__Carrier = Carrier
        if plot:
            self.plot()
        
    def calc(self, x):
        T = self.__T
        Cp = self.__Gas.Cp
        CpCarrier = self.__Carrier.Cp
        M = self.__Gas.Mass
        MCarrier = self.__Carrier.Mass
        v2 = 2*T*(x*Cp + (1-x)*CpCarrier)/((1-x)*MCarrier + x*M)*1000
        v = np.sqrt(v2)
        Ekin = 0.5*M/const.Avogadro/const.eV*v2/1000
        return v, Ekin
    
    def calc_Ekin(self, x):
        return self.calc(x)[1]
    
    def calc_v(self, x):
        return self.calc(x)[0]
    
    def calc_x(self, Ekin):
        pass
        
    def __call__(self, x):
        return self.calc_Ekin(x)
    
    def plot(self, xlim=(0,1)):
        fig, axes = plt.subplots(figsize=(16,9), dpi=100)
        x = np.linspace(0, 1, 500)
        axes.plot(x, self(x), linewidth=2,
                  label=r'%s in %s' % (self.__Gas.Name, self.__Carrier.Name))
        axes.set_title(r'Kinetic energy of %s seeded in %s at %d K' %\
                           (self.__Gas.Name, self.__Carrier.Name, self.__T))
        axes.set_xlabel(r'$x$(%s)' % self.__Gas.Name)
        axes.set_ylabel(r'$E_\mathrm{kin} \, / \, \mathrm{eV}$')
        axes.set_xlim(xlim)
        axes.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.01))
        axes.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05))
        axes.grid(which='minor', axis='x')
        axes.grid(which='major', axis='x')
        axes.grid(which='minor', axis='y')
        axes.grid(which='major', axis='y')
        self.fig = fig
        self.axes = axes
        plt.show(fig)

if __name__ == '__main__':
    SeededBeam(298, CO, H2).plot()
    fig, axes = BeamEnergiesPlot(198, CO)
    plt.show(fig)
