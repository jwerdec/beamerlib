from __future__ import division 
from PyQt4 import Qt
from PyQt4.QtGui import QWidget, QMainWindow, QHBoxLayout, QVBoxLayout
from guiqwt.curve import CurvePlot
from guiqwt.plot import PlotManager
from guiqwt.builder import make
from guiqwt.signals import SIG_RANGE_CHANGED
from guiqwt.shapes import XRangeSelection
from guiqwt.styles import ShapeParam, LineStyleParam
from guiqwt.tools import SelectTool
from numpy import trapz, polyfit, polyval, array, where, loadtxt
from scipy.integrate import quad
from PyQt4.QtGui import QApplication
from helper import construct_filename

class iScopeWidget(QWidget):
    """
    Filter testing widget
    parent: parent widget (QWidget)
    x, y: NumPy arrays
    func: function object (the signal filter to be tested)
    """
    def __init__(self, parent, x, y):
        QWidget.__init__(self, parent)
        self.setMinimumSize(320, 200)
        self.x = x
        self.y = y
        #---guiqwt related attributes:
        self.plot = None
        self.curve_item = None
        #---
        
    def setup_widget(self):
        #---Create the plot widget:
        x = self.x
        y = self.y
        self.plot = CurvePlot(self)
        self.curve_item = make.curve([], [], color='b')
        self.plot.add_item(self.curve_item)
        self.plot.set_antialiasing(True)
        width = x[-1] - x[0]
        self.intrange = make.range(x[0]+0.4*width, x[-1]-0.4*width)
        self.plot.add_item(self.intrange)
        self.lbgrange = make.range(x[0]+0.3*width, x[-1]-0.65*width)
        self.plot.add_item(self.lbgrange)
        self.lbgrange.pen = Qt.QPen(Qt.QColor('blue'))
        self.lbgrange.brush = Qt.QBrush(Qt.QColor(0,0,120,100))
        self.rbgrange = make.range(x[0]+0.7*width, x[-1]-0.1*width)
        self.rbgrange.pen = Qt.QPen(Qt.QColor('blue'))
        self.rbgrange.brush = Qt.QBrush(Qt.QColor(0,0,120,100))
        self.label1 = make.label(r"", "TR", (0, 0), "TR")
        self.plot.add_item(self.rbgrange)
        self.bg_item = make.curve([], [], color='r')
        self.plot.add_item(self.bg_item)
        self.fit_bg()
        self.plot.add_item(self.label1)
        self.connect(self.plot, SIG_RANGE_CHANGED, self.fit_bg)
        #---
        vlayout = QVBoxLayout()
        vlayout.addWidget(self.plot)
        self.setLayout(vlayout)
        self.update_curve()
        
    def fit_bg(self):
        degree = 3
        table = array([self.x, self.y])
        low, high = self.lbgrange.get_range()
        idxlower = set(where(table[0]<=high)[0])
        idxhigher = set(where(table[0]>=low)[0])
        idx1 = list(idxhigher.intersection(idxlower))
        low, high = self.rbgrange.get_range()
        idxlower = set(where(table[0]<=high)[0])
        idxhigher = set(where(table[0]>=low)[0])
        idx2 = list(idxhigher.intersection(idxlower))
        idx = idx1 + idx2
        x, y = table[:, idx]
        self.coeff = polyfit(x, y, degree)
        
        left, right = self.intrange.get_range()
        #bg = abs(quad(lambda x: polyval(self.coeff, x), left, right)[0])
        idxlower = set(where(table[0]<=right)[0])
        idxhigher = set(where(table[0]>=left)[0])
        idx = list(idxhigher.intersection(idxlower))
        x, y = table[:, idx]
        self.int = abs(trapz(y-polyval(self.coeff, x), x=x))
        
        self.update_label()
        self.update_curve()
        
    def update_label(self):
        self.label1.set_text(u"""trapz(red) - int(blue) = %e""" % self.int)
            
    def update_curve(self):
        #---Update curve
        self.curve_item.set_data(self.x, self.y)
        self.plot.replot()
        y = polyval(self.coeff, self.x)
        self.bg_item.set_data(self.x, y)
    
class iScope(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle("iScope (TM)")
        
        hlayout = QHBoxLayout()
        central_widget = QWidget(self)
        central_widget.setLayout(hlayout)
        self.setCentralWidget(central_widget)
        #---guiqwt plot manager
        self.manager = PlotManager(self)
        #---
        
    def add_plot(self, x, y):
        widget = iScopeWidget(self, x, y)
        widget.setup_widget()
        self.centralWidget().layout().addWidget(widget)
        #---Register plot to manager
        self.manager.add_plot(widget.plot)
        #---
        
    def setup_window(self):
        #---Add toolbar and register manager tools
        toolbar = self.addToolBar("tools")
        self.manager.add_toolbar(toolbar, id(toolbar))
        self.manager.register_standard_tools()
        self.manager.tools[0].activate()

def IntegrateScopeTrace(fileprefix, filenumbers):
    intensities = []

    for dat in [loadtxt(filename, unpack=True) for filename in
		construct_filename(fileprefix, filenumbers)]:
        x, y = dat
        app = QApplication.instance()
        if not app:
            app = QApplication([])
        win = iScope()
        win.add_plot(x*1e6, y)
        #---Setup window
        win.setup_window()
        win.resize(800,600)
        #---
        win.show()
        try:
            app.exec_()
            intensities.append(win.centralWidget().layout().itemAt(0).widget().int)
        except:
            pass
    
    print 'Number of elements: %d' % len(intensities)
    print intensities
