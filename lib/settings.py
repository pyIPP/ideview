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

class Settings(object):
    def __init__(self, filename='%s/.ideview' % os.environ['HOME']):
        self.filename = os.path.expandvars(filename)
        if os.path.isfile(self.filename):
            self._load()
        else:
            self._new()

    startDict = {'selectedPlots':['profile-pressure', 'profile-q', 'contour-rho','trace-Wmhd'],
                'lastDiag': 'IDE','lastShot': '31163','lastExp': 'ABOCK','lastEd': '0','windowSize': '800x600'}

    def _new(self):  # some default settings
        self.__dict__.update(self.startDict)
        self.save()

    def _load(self):
        self.__dict__.update(json.load(open(self.filename)))

    def save(self):
        d = copy(self.__dict__)
        d.pop('filename')
        json.dump(d, open(self.filename, 'w'))

    def __setattr__(self, name, value):
        self.__dict__[name] = value
        if name != 'filename':
            self.save()

    def __getattr__(self, name):
        if name not in self.__dict__:
            self.__dict__[name] = self.startDict[name]
        return self.__dict__[name]
