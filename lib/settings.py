import os
from copy import copy
import json

class Settings(object):
    def __init__(self, filename='%s/.ideview' % os.environ['HOME']):
        self.filename = os.path.expandvars(filename)
        if os.path.isfile(self.filename):
            self._load()
        else:
            self._new()

    startDict = {'selectedPlots':[], 'selectedGyrotrons': [],'lastDiag': 'IDE','lastShot': '31113','lastExp': 'AUGD','lastEd': '0','windowConditions': '800x600+5+5'}

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
