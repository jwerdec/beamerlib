import numpy as np
from scipy.integrate import quad
from scipy.signal import find_peaks_cwt
from math import floor, ceil
import Tkinter as Tk
import tkFileDialog
import os

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

if __name__ == '__main__':
    pass
