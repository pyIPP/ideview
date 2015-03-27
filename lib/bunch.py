import os, sys, numpy as np, matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.gridspec as gridspec
from matplotlib.figure import Figure
from IPython import embed
from copy import copy,deepcopy
import json
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')

import dd_20140805 as dd
import kk_abock as kk
sys.path.append('/afs/ipp/home/g/git/python/repository/')

from tooltip import createToolTip
import fconf

import Tkinter as tk

class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class PlotBunch(Bunch):
    kind = 'profile'

class PlotInfoBunch(PlotBunch):
    kind = 'profile'
    parameter = 'pressure'
#    name = '%s-%s' % (kind, parameter)

#print PlotInfoBunch.name
