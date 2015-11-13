from itemObject import itemObject

class q(itemObject):
    """docstring for q"""

    availableAreabases = ['rho_pol']

    def __init__(self, shotfiles):
        super(q, self).__init__(shotfiles)
        self.idf = [shotfile for shotfile in shotfiles if 'IDF' in shotfile.getObjectNames().values()][0]

    def isAvailable(self):
        for shotfile in self.shotfiles:
            if 'q_sa' in shotfile.getObjectNames().values():
                return True
        return False

    def getData(self, form='radial_1D', coordinates='rho_pol', t_ind=None):
        if self.cache is None:
            self.cache = {}
            self.cache['q_sa'] = self.idf('q_sa')
            self.cache['q_sa_plu'] = self.idf('q_sa_plu')
            self.cache['q_sa_min'] = self.idf('q_sa_min')
            self.cache['rhopol_q'] = self.idf('rhopol_q')
        rp = self.cache['rhopol_q'] 
        q, q_plus, q_min = self.cache['q_sa'], self.cache['q_sa_plu'], self.cache['q_sa_min']

        if form == 'radial_1D':
            if coordinates == 'rho_pol':
                x = rp.data[t_ind]
            elif coordinates == 'rho_tor':
                pass
            toReturn = {'x':x, 'y':-q.data[t_ind], 'yp':-q_min.data[t_ind], 'ym':-q_plus.data[t_ind]}
        elif form == 'temporal_2D':
            if coordinates == 'rho_pol':
                y = rp.data[0]
            x = q.time
            toReturn = {'x':x, 'y':y, 'z':-q.data.T}

        return toReturn
