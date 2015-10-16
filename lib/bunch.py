

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
