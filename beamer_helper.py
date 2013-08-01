"""
Requirements:
    Python v2.7 or later (not compatible with Python 3)
    numpy, matplotlib
Tested with:
    Linux Mint 14 Cinnamon 64bit
    python 2.7.3
    numpy 1.6.2
    matplotlib 1.2.1

Version History:
"""

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


def construct_filename(prefix_str, num):
    if type(num) == type([]):
        return ['%s%.3d.dat' % (prefix_str, i) for i in num]
    else:
        return '%s%.3d.dat' % (prefix_str, num)


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

