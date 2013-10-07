import numpy as np
from scipy.integrate import quad
from scipy.signal import find_peaks_cwt
from math import floor, ceil

def normalize(x, numpoints=10):
    return x/x[x.argsort()[-numpoints:]].mean()

def subtract_baseline(x, left=20):
    return x - x[:left].mean()

def sort_table(table, col=0):
    return table[:, table[col, :].argsort()]

def find_peaks(x, y, widthrange, rel_threshold=0.1):
    """Peak-finding in a 2d dataset.
    
    Parameters
    ----------
    x, y : array_like
        Input arrays.
    widthrange : tuple
        Lower and upper limit of peak widths to find.
    rel_threshold : float, optional
        Peaks with a height lower than this value times the height of the
        maximum in 'y' are ignored.

    Returns
    -------
    list
        Array indices of where the peaks were found in 'y'.

    See Also
    --------
    scipy.signal.find_peaks_cwt : Peak-finding using a continous wavelet
    transform technique.
    """
    dx = abs(x[1] - x[0])
    minwidth, maxwidth = widthrange
    widths = np.arange(floor(minwidth/dx), ceil(maxwidth/dx))
    peakpos = find_peaks_cwt(y, widths)
    maxy = max(y)
    return [pos for pos in peakpos if y[pos] >= rel_threshold*maxy]

def Moments(func, n, limits=(0, np.inf), args={}):
    """
    Calculate the nth moment around 0 of a function func.
    Parameter args={} can be used to pass additional arguments to the function
    limits: tuple of limits of integration (default=(0, infinity))
    """
    if args=={}:
        _func = lambda x: func(x)
    else:
        _func = lambda x: func(x, **args)
    if type(n) == type([]):
        result = [i for i in n]
        for i in n:
            function = lambda x: _func(x)*x**i
            result[i] = quad(function, limits[0], limits[1])
    else: 
        function = lambda x: _func(x)*x**n
        result = quad(function, limits[0], limits[1])
    return result

def extract_from_Table(Table, limits=(),  col=0):
    """
    Returns a new 2d array with only those rows from the input table whos
    value in column col are within the specified limits or equal to them.
    The ordering is preserved.
    """
    if limits==():
        return Table
    low = limits[0]
    high = limits[1]
    idxlower = set(np.where(Table[col]<=high)[0])
    idxhigher = set(np.where(Table[col]>=low)[0])
    idx = idxhigher.intersection(idxlower)
    return Table[:, list(idx)]

def extract_from_Array(array, limits=()):
    """
    Returns a new 1d array with only those rows from the input array whos
    value are within the specified limits or equal to them.
    The ordering is preserved.
    """
    if limits==():
        return array
    low = limits[0]
    high = limits[1]
    idxlower = set(np.where(array<=high)[0])
    idxhigher = set(np.where(array>=low)[0])
    idx = idxhigher.intersection(idxlower)
    return array[list(idx)]

def construct_filename(prefix_str, num):
    if type(num) == type([]):
        return ['%s%.3d.dat' % (prefix_str, i) for i in num]
    else:
        return '%s%.3d.dat' % (prefix_str, num)

import Tkinter as Tk
import tkFileDialog
import os

def SelectFileDialog(relative=True):
    root = Tk.Tk()
    root.withdraw()
    current_path = os.getcwd()
    file_paths = tkFileDialog.askopenfilenames(
        initialdir=current_path,
        filetypes=[('data files', '.dat'), ('all files', '.*')])
    if relative:
        file_paths = [os.path.relpath(path, start=current_path)
                      for path in file_paths]
    return file_paths

import  mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import SubplotHost
from mpl_toolkits.axisartist import GridHelperCurveLinear
from mpl_toolkits.axisartist.grid_finder import MaxNLocator, FixedLocator
from numpy import pi

def SemiPolarPlot(fig):
    # see demo_curvelinear_grid.py for details
    tr = Affine2D().scale(pi/180., 1.) + PolarAxes.PolarTransform()
    extreme_finder = angle_helper.ExtremeFinderCycle(20, 20,
                                                     lon_cycle = 360,
                                                     lat_cycle = None,
                                                     lon_minmax = None,
                                                     lat_minmax = (0, 1),
                                                     )

    grid_locator1 = angle_helper.LocatorDMS(11)
    grid_locator2 = FixedLocator([0.25, 0.5, 1., 0.75])
    tick_formatter1 = angle_helper.FormatterDMS()
    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        grid_locator2=grid_locator2,
                                        tick_formatter1=tick_formatter1,
                                        tick_formatter2=None
                                        )


    ax1 = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    ax1.axis["top"].set_visible(False)
    ax1.axis["right"].set_visible(False)
    ax1.axis["bottom"].set_visible(False)
    ax1.axis["lon"] = ax1.new_floating_axis(1, 1)
    ax1.set_aspect(1)
    ax1.set_xlim(0, 2)
    ax1.set_ylim(-1., 1.)
    ax1.grid(True)
    
    curved_ax = ax1.get_aux_axes(tr)

    curved_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9
    
    return ax1, curved_ax

    
if __name__ == '__main__':
    import matplotlib.pylab as plt
    from numpy import linspace, cos, radians
    angles = linspace(-45, 45, 20)
    signal = cos(radians(angles-8))**12
    fig = plt.figure(1, figsize=(8, 8))
    fig.clf()
    ax1, curved_ax = SemiPolarPlot(fig)
    curved_ax.plot(angles, signal, 'bo',\
                       label=r'$\mathrm{Synthetic\ Data}$')
    theta = linspace(-90, 90, 100)
    label = r'$\cos(\theta)$' 
    curved_ax.plot(theta, cos(radians(theta)), 'k--', label=label)
    ax1.legend(loc='best', numpoints=1)
    ax1.set_title(r'Angular Distribution')
    plt.show(fig)
