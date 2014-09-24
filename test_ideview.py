#!/usr/bin/env python

import imp
ideview = imp.load_source('ideview', 'ideview')
from IPython import embed
import matplotlib.pyplot as plt
import numpy as np

be = ideview.ShotfileBackend()

times = be.getAvailableTimes()

names = be.getAvailablePlotNames()

plots = be.getPlotsForTimePoint(names, times[len(times)/2])


for i, key in enumerate(plots):
    f = plt.figure(i)
    f.canvas.set_window_title(key)
    ps = plots[key]
    if not hasattr(ps, 'data'):
        continue
    if ps.kind == 'profile':
        for p in ps.data:
            plt.plot(p['x'], p['y'], p['ls'])
    elif ps.kind == 'contour':
        for p in ps.data:
            plt.contour(p['x'], p['y'], p['z'], levels=np.arange(0, 1.1, 0.1))
            ax = plt.gca()
            ax.set_aspect('equal')
    elif ps.kind == 'trace':
        for p in ps.data:
            plt.plot(p['x'], p['y'], p['ls'])

    if hasattr(ps, 'ylim'):
        plt.ylim(ps.ylim)

plt.show()

#embed()
