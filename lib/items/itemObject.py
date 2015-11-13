

class itemObject(object):
    """docstring for itemObject"""
    
    cache = None

    availableAreabases = None

    def __init__(self, shotfiles):
        super(itemObject, self).__init__()
        self.shotfiles = shotfiles
    
    def isAvailable(self):
        return False

    def getAvailableAreabases(self):
        return self.availableAreabases

    def getData(self, form=None, coordinates=None, t_ind=None):
        return 





#q.getData('temporal_2D', 'rho_pol')