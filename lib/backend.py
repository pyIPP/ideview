import sys
import numpy as np
from IPython import embed
from copy import copy,deepcopy
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
from lib.bunch import PlotBunch
#from lib.settings import Settings
import dd_20140805 as dd
import kk_abock as kk
sys.path.append('/afs/ipp/home/g/git/python/repository/')


class Backend(object):
    def __init__(self):
        pass

    def getAvailableTimes(self):
        return None

    def getAvailablePlotNames(self):
        return None

    def getPlotsForTimePoint(self, time):
        return None

class ShotfileBackend(Backend):
    openedEditions = {}

    def __init__(self, diag='TRE',experiment='TODSTRCI', shot=30579, edition=0):
        super(ShotfileBackend, self).__init__()
        self.experiment = experiment
        self.shot = shot
        self.edition = edition
        self.equ = []
        self.openedEditions = {}
        #self.settings = Settings()
        #load shotfile with 2d quantities
        self.equ.append(dd.shotfile(diag, shot, experiment, edition))
        self.openedEditions[diag] = self.equ[-1].edition
        try:
            if diag == 'IDE': 
                idf_ed = self.equ[-1].getParameter('idefg_ed', 'idf_ed').data
                idg_ed = self.equ[-1].getParameter('idefg_ed', 'idg_ed').data
                self.equ.append(dd.shotfile('IDF', shot, experiment, idf_ed))
                self.openedEditions['IDF'] = self.equ[-1].edition
                self.equ.append(dd.shotfile('IDG', shot, experiment, idg_ed))
                self.openedEditions['IDG'] = self.equ[-1].edition
        except Exception as e:
            print 'WARNING! Please use the old edition.'
            self.equ.append(dd.shotfile('IDF', shot, experiment, edition))
            self.openedEditions['IDF'] = self.equ[-1].edition
            self.equ.append(dd.shotfile('IDG', shot, experiment, edition))
            self.openedEditions['IDG'] = self.equ[-1].edition

        #shotfiles with 1d quantities
        diags_1d =  {'EQE':'GQE','FPQ':'FPK','EQI':'GQI','EQH':'GQH',
                     'FPP':'GPI','EQR':'FPG'}
        
        if diag in diags_1d:
            self.equ.append(dd.shotfile(diags_1d[diag], shot, experiment, edition))
            self.openedEditions[diags_1d[diag]] = self.equ[-1].edition

        #shotfile with the errorbars
        
        if diag == 'MGS':
            raise Warning('MGS is not supported yet')

        self.eq = kk.kk()
        self.eq.Open(shot, experiment, diag, edition)
        self.times = self.equ[0]('time').data

    def __del__(self):
        for key in self._cache.keys():
            del self._cache[key]

    def getAvailableTimes(self):
        return self.times
 
    __plotNames = ['profile-pressure', 'profile-pcon', 'profile-pressure/pcon', 'profile-q', 'profile-ECcur_tot', 'profile-ECcur_gyr', 'profile-Jpol_tot', 'profile-Jpol_pla', 'profile-Itps', 'profile-Itps_av', 'profile-Itfs', 'profile-Itfs_av', 'contour-pfl', 'contour-rho', 'trace-Wmhd', 'timecontour-pressure', 'timecontour-q', 'timecontour-Dpsi', 'timecontour-Iext', 'timecontour-pcon', 'timecontour-pol', 'timecontour-mse','timecontour-I_tor', 'timecontour-I_torcon', 'timecontour-Bprob', 'profile-I_tor', 'profile-Dpsi', 'res(profile)-Dpsi', 'profile-Iext', 'res(profile)-Iext', 'res(profile)-pcon', 'profile-pol', 'profile-mse', 'profile-I_torcon', 'profile-Bprob', 'res(profile)-Bprob', 'trace-Bprob', 'trace-Dpsi', 'trace-Iext', 'trace-Rmag', 'trace-Zmag', 'trace-Rin', 'trace-Raus', 'trace-betapol', 'trace-betapol+li/2', 'trace-Itor', 'trace-Rxpu', 'trace-Zxpu', 'trace-ahor', 'trace-bver', 'trace-bver/ahor', 'trace-XPfdif', 'trace-delR', 'trace-delZ', 'trace-q0', 'trace-q25', 'trace-q50', 'trace-q75', 'trace-q95', 'trace-delR_oben', 'trace-delR_unten', 'trace-k_oben', 'trace-k_unten', 'trace-dRxP', 'trace-eccd_tot', 'trace-Itax', 'trace-ecrhmax(R)', 'trace-ecrhmax(z)', 'trace-ecrhmax(y)', 'trace-ecrhmax(rho)', 'trace-li']

    def getAvailablePlotNames(self):
        return self.__plotNames

    def getAvailableGyrotrons(self):
        if self.getData('ecrhpos') != None:
            return self.getData('ecrhpos').data.shape[3]

    def lookfordata(self, name, onlyMES=False):
        try:
            mes = '%s_mes' %name
            MES = self.getData(mes)
            FIT = None
            if MES is None:
                mes = '%smes' %name
                MES = self.getData(mes)
            if MES is None:
                mes = name
                MES = self.getData(mes)
            if not onlyMES:
                if name == 'betpol':
                    UNC = self.getData('betp_unc')
                else:
                    unc = '%s_unc' %name
                    UNC = self.getData(unc)
                    if UNC is None:
                        unc = '%sunc' %name
                        UNC = self.getData(unc)

                fit = '%s_fit' %name
                FIT = self.getData(fit)
                if FIT is None:
                    fit = '%sfit' %name
                    FIT = self.getData(fit)
            # workaround:
            if MES.area is None:
                lookforab = '%s_rp'%name
                MES.area = self.getData(lookforab)
                if name == 'q_sa':
                    MES.area = self.getData('rhopol_q')
                    
                if not onlyMES and FIT is not None and FIT.area is None:
                    FIT.area = self.getData(lookforab)
            if onlyMES:
                return MES
            else:
                return MES, FIT, UNC
        except Exception as e:
            print e, '\n Something went wrong while looking for the data'
            pass

    def profiledata(self, name, t, t_ind, unc=True, fit=True):
        try:
            MES, FIT, UNC = self.lookfordata(name)
            XMAX = MES.area.data.max() if (MES.area.data.max() < 1000) else 1
            XMIN = MES.area.data.min() if (MES.area.data.min() > -1000) else 0
            YMAX = np.nanmax(MES.data)
            YMIN = np.nanmin(MES.data)

            shape = MES.area.data.shape
            MES_Area_data = MES.area.data[0 if shape[0]==1 else t_ind]
            MES_Data = MES.data[t_ind]
            if name != 'q_sa' and name != 'eccurgyr':
                data=[{'x': MES_Area_data, 'y': MES_Data, 'marker':'', 'ls':'-'}]

            
            if (not UNC is None) and unc:
                UNC_Data = UNC.data[t_ind]
                if 'con' in name:
                    data[0].update({'c':'k'})
                    data.extend([{'x': MES_Area_data, 'y': MES_Data + UNC_Data,'marker':'_', 'ls':'', 'c':'k'},
                                 {'x': MES_Area_data, 'y': MES_Data-UNC_Data, 'marker':'_', 'ls':'', 'c':'k'}])
                else:
                    data.extend([{'x': MES_Area_data, 'y': MES_Data + UNC_Data, 'ls':'--'},
                                 {'x': MES_Area_data, 'y': MES_Data-UNC_Data, 'ls': '--'}])
            if (not FIT is None) and fit:
                FIT_Area_data = FIT.area.data[0 if shape[0]==1 else t_ind]
                FIT_Data = FIT.data[t_ind]
                data.extend([{'x': FIT_Area_data, 'y': FIT_Data, 'c':'r','ls':'-'}])

            if name == 'q_sa':
                qsa = MES
                YMIN = 0
                YMAX = 10
                rho = qsa.area.data[t_ind]
                rho[0] = np.nan  #diffifulties with elevated profiles and q=1
                q = abs(qsa.data)
                
                XIQ1  = np.interp(1,q[t_ind],rho,left=np.nan)
                XIQ2  = np.interp(2,q[t_ind],rho,left=np.nan)
                XIQ3_2= np.interp(1.5,q[t_ind],rho,left=np.nan)
                #XIQ3  = np.interp(3,q[ind],rho,left=np.nan)

                data = [{'x':  qsa.area.data[t_ind], 'y': q[t_ind]},{'x': XIQ1}, {'x': XIQ2},{'x':XIQ3_2},
                      {'y': 1,'c':'k','alpha':.3}, {'y': 2,'c':'k','alpha':.3},{'y': 1.5,'c':'k','alpha':.3}]
                
                qsap = self.getData('q_sa_plu')
                qsap.area = self.getData('rhopol_q')
                qsam = self.getData('q_sa_min')
                qsam.area = self.getData('rhopol_q')
                if not qsap is None and not qsam is None:
                    #qsap = qsap(tBegin=t, tEnd=t)
                    #qsam = qsam(tBegin=t, tEnd=t)
                    #embed()
                    data.extend([{'x': qsap.area.data[t_ind], 'y': np.abs(qsap.data)[t_ind], 'ls': '--'},
                                {'x': qsam.area.data[t_ind], 'y': np.abs(qsam.data)[t_ind], 'ls': '--'}])

            if name == 'eccurgyr':
                data = []
                colours = ['r','b', 'g', 'brown', 'cyan', 'magenta', 'purple', 'orange']
                MES.area = self.getData('eccd_rp')
                for i in range(0, int(MES.data.shape[2])):
                    colour = colours[i]
                    if MES.area.data[t_ind].max() > 1:                        
                        data.extend([{'x': np.array([0,1]), 'y': np.array([0,0]), 'ls': '--', 'c': colour, 'label':'gyr%i'%(i+1), 'exc' : True}])
                    else:
                        data.extend([{'x': MES.area.data[t_ind], 'y': MES.data[t_ind, :, i], 'ls': '--', 'c': colour, 'label':'gyr%i'%(i+1), 'exc' : True}])
                total = self.getData('eccurtot')
                if MES.area.data[t_ind].max() > 1:                        
                    data.extend([{'x': np.array([0,1]), 'y': np.array([0,0]), 'ls': '-', 'c': 'k', 'label':'total', 'exc' : True}])
                else:
                    data.extend([{'x': MES.area.data[t_ind], 'y': total.data[t_ind, :], 'ls': '-', 'c': 'k', 'label':'total', 'exc' : True}])

            return PlotBunch(data=data, setting={'xlim': (XMIN,XMAX), 'ylim':(YMIN,YMAX)})

        except Exception as e:
            print e,  '\n%s-profile is not available in that shot for t=%s'%(name, t)
            pass

    def resdata(self, name, t, t_ind):
        try:
            MES, FIT, UNC = self.lookfordata(name)
            RESarray = np.array([])
            RES = (MES.data[t_ind]-FIT.data[t_ind])/UNC.data[t_ind]
            tmparray = np.array([RES])
            RESarray = np.append(RESarray,tmparray)
     
            #data=[{'y': 0,'ls':'--','c':'k'}]
            # {'x': MES.area.data[0], 'y':RESarray,'marker':'+','ls':'', 'label': '%s residuum'%name, }
            data = []
            for i, y in enumerate(RESarray):
                if i == 0:
                    data.append({'x': [MES.area.data[0,i]]*2, 'y': [0, y], 'marker':'+', 'markersize':10, 
                        'label': '%s residuum'%name, 'exc':True}) # exc: exception for residua labels
                else:
                    data.append({'x': [MES.area.data[0,i]]*2, 'y': [0, y], 'marker':'+', 'markersize':10})
            return PlotBunch(data=data, setting={'ylim':(-3, 3)})

        except (TypeError,AttributeError) as e:
            print e, '\n%s-residuum is not available in that shot for t=%s'%(name, t)
            pass

    def tracedata(self, name, t, ecrh = None):
        try:
            MES, FIT, UNC = self.lookfordata(name)
            if not 'signalGroup' in str(type(MES)):
                qlist = ['q0', 'q25', 'q50', 'q75', 'q95']
                if name in qlist:
                    q = abs(MES.data)
                    uncplus = '%s_plus'%name
                    UNCplus = self.getData(uncplus)
                    uncminus = '%s_minu'%name
                    UNCminus = self.getData(uncminus)
                    data=[{'x': MES.time, 'y': q, 'ls': '-'},{'x': t,'c': 'k'},
                          {'x': MES.time, 'y': abs(UNCplus.data), 'ls': '--'}, {'x' : MES.time, 'y': abs(UNCminus.data), 'ls': '--'}]
                    return PlotBunch(kind='trace',data=data, setting={'ylim':(0,10)})
                else:
                    data=[{'x': MES.time, 'y': MES.data, 'ls': '-'},{'x': t,'c': 'k'}]
                    if not UNC is  None:
                        data.extend([{'x': MES.time, 'y': MES.data+UNC.data, 'ls': '--'},
                                     {'x': MES.time, 'y': np.maximum(0,MES.data-UNC.data), 'ls': '--'}])

                    return PlotBunch(kind='trace',data=data)

            elif name == 'ecrhmax':
                data = [{'x': t,'c': 'k'}]
                colours = ['r','b', 'g', 'brown', 'cyan', 'magenta', 'purple', 'orange']
                if ecrh == 'R':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,0,i], 'label':'gyro%i'%(i+1),  'c': colour, 'exc' : True})
                if ecrh == 'z':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,1,i], 'label':'gyro%i'%(i+1), 'c': colour, 'exc' : True})
                if ecrh == 'y':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,2,i], 'label':'gyro%i'%(i+1), 'c': colour, 'exc' : True})
                if ecrh == 'rho':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,3,i], 'label':'gyro%i'%(i+1), 'c': colour, 'exc' : True})

                return PlotBunch(kind='trace',data=data)

            else:
                tmp = ((MES.data-FIT.data)/UNC.data)**2
                points = MES.data.shape[-1]
                sum_ = np.sum(tmp, axis=1)
                TRACE = sum_/points
                data=[{'x':MES.time, 'y':TRACE, 'ls':'-'},{'x':t,'c':'k'}]
                
                return PlotBunch(kind='trace',data=data)

        except (TypeError,AttributeError)as e:
            print e, '\n%s-trace is not available in that shot for t=%s'%(name, t)
            pass

    def timecontourdata(self, name, t):
        try:
            MES = self.lookfordata(name, onlyMES=True)
            if name == 'q_sa':
                return PlotBunch(kind='timecontour', data=[{'x':MES.time, 'y':MES.area.data[0],
                                                            'z': np.abs(MES.data),'levels':[1, 1.5, 2, 3, 4, 5]}, {'x':t, 'c':'k'}], setting= {'ylim':(0,1)})
            else:
                return PlotBunch(kind='timecontour', data=[{'x':MES.time, 'y':MES.area.data[0],
                                                            'z':MES.data}, {'x':t, 'c':'k'}])
        except (TypeError,AttributeError) as e:
            print e, '\n%s-timecontour is not available in that shot for t=%s'%(name, t)
            pass

    _cache = {}

    def getData(self, name, workaround_time=None):

        if name == 'pfm' and workaround_time != None:
            return self.eq.get_pfm(workaround_time)

        #workaround for the loading of the data from the equilibrium shotfiles
        if name in ['Qpsi','Pres','Jpol'] and workaround_time!= None and not name in self._cache:
            
            tvec = self.getData('time').data

            t_index = slice(None,len(tvec))

            if name in ['Pres','Jpol']: 
                prof = self.getData(name)[t_index,::2] 
            else:
                prof = self.getData(name)[t_index,:]

            psiAx, psiSep = self.getData('PFxx').data[None,t_index, :2].T
            phi = self.getData('PFL')[t_index,:]
            rho = np.sqrt(np.abs((phi-psiAx).data/(psiSep-psiAx)))
            ind = np.argsort(rho,axis=1)
            
            for i,ii in enumerate( ind): prof.data[i] = prof.data[i][ii]
            for i,ii in enumerate( ind): rho[i] = rho[i][ii]

            prof.area  = Data()
            prof.area.data = rho
            prof.time = tvec
            
            self._cache[name] = prof

        if name in ['Qpsi','Pres','Jpol'] and workaround_time!= None and name in self._cache:
            if np.isfinite(workaround_time):
                prof = self._cache[name]
                t_index = np.argmin(np.abs(prof.time-workaround_time))
                prof_ = deepcopy(prof)
                prof_.data = prof.data[t_index]
                prof_.time = prof.time[t_index]
                prof_.area.data = prof.area.data[t_index]
                return prof_#__

        if name not in self._cache:
            self._cache[name] = None   #returned in the case that nothing was found
            #embed()
            for diag in self.equ:
                if name in diag.getObjectNames().values():
                    try:# solving Exception: Error by DDsinfo (8.1): status of signal doesn't permit access
                        self._cache[name] = copy(diag(name))
                    except Exception:
                        pass
                    break
  
        return self._cache[name]

    def getSinglePlotForTimePoint(self, name, t):
        t_ind = np.abs(self.getData('time').data - t).argmin()
        if 'timecontour' in name:

            if name == 'timecontour-pressure':
                return self.timecontourdata('pres', t)
            
            elif name == 'timecontour-q':
                return self.timecontourdata('q_sa', t)

            elif name == 'timecontour-Bprob':
                return self.timecontourdata('Bprob', t)

            elif name == 'timecontour-I_tor':
                return self.timecontourdata('It', t)

            elif name == 'timecontour-Dpsi':
                return self.timecontourdata('Dpsi', t)

            elif name == 'timecontour-Iext':
                return self.timecontourdata('Iext', t)

            elif name == 'timecontour-pcon':
                return self.timecontourdata('pcon', t)

            elif name == 'timecontour-pol':
                return self.timecontourdata('pol', t)

            elif name == 'timecontour-mse':
                return self.timecontourdata('mse', t)

            elif name == 'timecontour-I_torcon':
                return self.timecontourdata('jtcon', t)

        elif 'contour' in name:
            t_index = np.abs(self.times - t).argmin()
            #Ri = self.getData('Ri').data[t_index]
            #Zj = self.getData('Zj').data[t_index]
            #pfm = self.getData('PFM').data[t_index]
            #embed()
            
            tmp = self.getData('pfm', t)
            Ri = tmp['Ri']; Zj = tmp['zj']; pfm = tmp['pfm']

            psiAx, psiSep = self.getData('PFxx').data[t_index, :2]
            #pfm = np.sqrt(np.abs((pfm-psiAx)/(psiSep-psiAx)))
            #embed()
            
            #ikCAT = self.getData('ikCAT').data
            #PFxx = self.getData('PFxx').data
            #sep_ind = np.array([2,0,3,1])[ikCAT-1]
            #psiAx,psiSep = PFxx[0,:], PFxx[sep_ind,np.arange(len(mag))]        
            
            data = [{'x': Ri, 'y': Zj, 'z': pfm, 'psiSep': psiSep, 'psiAx': psiAx,'levels':np.arange(0,2,.1)}]

            if 'pfl' in name:
                pfm = pfm - psiSep
                lvls = np.sort(np.insert(np.arange(pfm.min(), pfm.max(), 0.05), 0, psiSep))
                data[0].update({'z': pfm,'levels':lvls})

            elif 'rho' in name:
                pfm = np.sqrt(np.abs((pfm-psiAx)/(psiSep-psiAx)))
                data[0].update({'z': pfm,'levels':np.arange(0,2,.1)})

                if self.getData('ecrhpos') != None and self.getData('ecrhmax') != None:
                    MESs = [self.getData('ecrhpos'), self.getData('ecrhposu'), self.getData('ecrhposl')]
                    MESmax = self.getData('ecrhmax')
                    colours = ['r','b', 'g', 'brown', 'cyan', 'magenta', 'purple', 'orange']
                    gyrind = []
                    for i in range(MESmax.data.shape[2]):
                        colour = colours[i]
                        for MES in MESs:
                            R = MES.data[t_index,:,0,i]
                            z = MES.data[t_index,:,1,i]
                            ind = np.where((R==0) + (z==0))
                            R = np.delete(R, ind)
                            z = np.delete(z, ind)
    
                            data.append({'x':R, 'y':z, 'ls':'-', 'c':colour, 'gyro_ind':i})
                        data.append({'x':MESmax.data[t_index,0,i], 'y':MESmax.data[t_index,1,i], 'ls':'', 'marker':'o', 'c':colour, 'gyro_ind':i})
                        gyrind.append(i)
                #embed()

            return PlotBunch(kind='contour', data=data)

        elif 'trace' in name:

            if name == 'trace-Wmhd':
                return self.tracedata('Wmhd',t)
            elif name == 'trace-Bprob':
                return self.tracedata('Bprob', t)
            elif name == 'trace-Dpsi':
                return self.tracedata('Dpsi', t)
            elif name == 'trace-Iext':
                return self.tracedata('Iext', t)
            elif name == 'trace-Rmag':
                return self.tracedata('Rmag', t)
            elif name == 'trace-Zmag':
                return self.tracedata('Zmag', t)
            elif name == 'trace-Rin':
                return self.tracedata('Rin', t)
            elif name == 'trace-Raus':
                return self.tracedata('Raus', t)
            elif name == 'trace-betapol':
                return self.tracedata('betpol', t)
            elif name == 'trace-betapol+li/2':
                return self.tracedata('bpli2', t)
            elif name == 'trace-li':
                return self.tracedata('li', t)
            elif name == 'trace-Itor':
                return self.tracedata('Itor', t)
            elif name == 'trace-Rxpu':
                return self.tracedata('Rxpu', t)
            elif name == 'trace-Zxpu':
                return self.tracedata('Zxpu', t)
            elif name == 'trace-ahor':
                return self.tracedata('ahor', t)
            elif name == 'trace-bver':
                return self.tracedata('bver', t)
            elif name == 'trace-bver/ahor':
                return self.tracedata('k', t)
            #elif name == 'trace-Vol':
                #return self.tracedata('Vol', t)
            elif name == 'trace-XPfdif':
                return self.tracedata('XPfdif', t)
            elif name == 'trace-delR':
                return self.tracedata('delR', t)
            elif name == 'trace-delZ':
                return self.tracedata('delZ', t)
            elif name == 'trace-q0':
                return self.tracedata('q0', t)
            elif name == 'trace-q25':
                return self.tracedata('q25', t)
            elif name == 'trace-q50':
                return self.tracedata('q50', t)
            elif name == 'trace-q50':
                return self.tracedata('q50', t)
            elif name == 'trace-q75':
                return self.tracedata('q75', t)
            elif name == 'trace-q95':
                return self.tracedata('q95', t)
            elif name == 'trace-delR_oben':
                return self.tracedata('delRoben', t)
            elif name == 'trace-delR_unten':
                return self.tracedata('delRuntn', t)
            elif name == 'trace-k_oben':
                return self.tracedata('koben', t)
            elif name == 'trace-k_unten':
                return self.tracedata('kuntn', t)
            elif name == 'trace-dRxP':
                return self.tracedata('dRXP', t)
            #elif name == 'trace-ikCAT':
                #return self.tracedata('ikCAT', t)
            elif name == 'trace-eccd_tot':
                return self.tracedata('eccd_tot', t)
            elif name == 'trace-ecrhmax(R)':
                return self.tracedata('ecrhmax', t, ecrh = 'R')
            elif name == 'trace-ecrhmax(z)':
                return self.tracedata('ecrhmax', t, ecrh = 'z')
            elif name == 'trace-ecrhmax(y)':
                return self.tracedata('ecrhmax', t, ecrh = 'y')
            elif name == 'trace-ecrhmax(rho)':
                return self.tracedata('ecrhmax', t, ecrh = 'rho')
            elif name == 'trace-Itax':
                return self.tracedata('Itax', t)

        elif 'profile' in name:
            if name == 'profile-q':
                return self.profiledata('q_sa', t, t_ind)
                '''qsa = self.getData('q_sa')
                if qsa is None:
                    qsa = self.getData('Qpsi',t)
                else:
                    qsa = qsa(tBegin=t, tEnd=t)

                ind = np.diff(np.abs(qsa.data)) > 0 #detect resonance surfaces only in the region of positive shear
                ind[-1] = True
                rho = qsa.area.data[ind]
                rho[0] = np.nan  #diffifulties with elevated profiles and q=1
                q = abs(qsa.data)
                
                XIQ1  = np.interp(1,q[ind],rho,left=np.nan)
                XIQ2  = np.interp(2,q[ind],rho,left=np.nan)
                XIQ3_2= np.interp(1.5,q[ind],rho,left=np.nan)
                #XIQ3  = np.interp(3,q[ind],rho,left=np.nan)

                data=[{'x':  qsa.area.data, 'y': q},{'x': XIQ1}, {'x': XIQ2},{'x':XIQ3_2},
                      {'y': 1,'c':'k','alpha':.3}, {'y': 2,'c':'k','alpha':.3},{'y': 1.5,'c':'k','alpha':.3}]
                
                qsap = self.getData('q_sa_plu')
                qsam = self.getData('q_sa_min')
                if not qsap is None and not qsam is None:
                    qsap = qsap(tBegin=t, tEnd=t)
                    qsam = qsam(tBegin=t, tEnd=t)
                    data.extend([{'x': qsap.area.data, 'y': np.abs(qsap.data), 'ls': '--'},
                                {'x': qsam.area.data, 'y': np.abs(qsam.data), 'ls': '--'}])
   
                return PlotBunch(data=data,setting={'ylim':(0, 10),'xlim':(0,1)})'''
                
            elif name == 'profile-pressure':
                return self.profiledata('pres', t, t_ind)

            elif name == 'profile-I_tor':
                return self.profiledata('It', t, t_ind)

            elif name == 'profile-Dpsi':
                return self.profiledata('Dpsi', t, t_ind, fit=False)

            elif name == 'res(profile)-Dpsi':
                return self.resdata('Dpsi', t, t_ind)

            elif name == 'profile-Iext':
                return self.profiledata('Iext', t, t_ind, fit=False)

            elif name == 'res(profile)-Iext':
                return self.resdata('Iext', t, t_ind)

            elif name == 'profile-pcon':
                return self.profiledata('pcon',t, t_ind)

            elif name == 'profile-pressure/pcon':
                p = self.getData('pres')
                if p is None:
                    p = self.getData('Pres',t)
                else:
                    p = p(tBegin=t, tEnd=t)
         
                punc = self.getData('pres_unc',t)
                pcon = self.getData('pcon_mes')
                pcon = pcon(tBegin=t, tEnd=t)
                maxpres = self.getData('Pres').data.max()
                pconunc = self.getData('pcon_unc')(tBegin=t, tEnd=t)
                
                data=[{'x': p.area.data, 'y': p.data, 'ls': '-', 'label':'pressure', 'exc': True},
                      {'x': pcon.area.data, 'y': pcon.data, 'ls':'', 'marker':'+', 'c':'k', 'label':'constraint', 'exc': True}]
                
                if not punc is None:
                    punc = punc(tBegin=t, tEnd=t)
                    data.extend([{'x': p.area.data, 'y': p.data+punc.data, 'ls': '--', 'label':'uncertainty', 'c':'b', 'exc': True}, {'x': p.area.data, 'y': p.data-punc.data, 'ls': '--', 'c':'b', 'exc': True}])
                if not pconunc is None:
                    data.extend([{'x': pcon.area.data, 'y': pcon.data+pconunc.data,'marker':'_', 'ls': '', 'label':'uncertainty', 'c':'k', 'exc': True}, {'x': pcon.area.data, 'y': pcon.data-pconunc.data,'marker':'_', 'ls': '', 'c':'k', 'exc': True}])
                
                return PlotBunch(data=data, setting={'ylim':(0, maxpres),'xlim':(0,1)})

            elif name == 'res(profile)-pcon':
                return self.resdata('pcon', t, t_ind)

            elif name == 'profile-pol':
                return self.profiledata('pol', t, t_ind)

            elif name == 'profile-mse':
                return self.profiledata('mse', t, t_ind)

            elif name == 'profile-I_torcon':
                return self.profiledata('jtcon', t, t_ind)

            elif name == 'profile-Bprob':
                return self.profiledata('Bprob', t, t_ind, unc=False)

            elif name == 'res(profile)-Bprob':
                return self.resdata('Bprob', t, t_ind)

            elif name == 'profile-ECcur_tot':
                return self.profiledata('eccurtot', t, t_ind)

            elif name == 'profile-ECcur_gyr':
                return self.profiledata('eccurgyr', t, t_ind)

            elif name == 'profile-Jpol_tot':
                return self.profiledata('Jpol_tot', t, t_ind)

            elif name == 'profile-Jpol_pla':
                return self.profiledata('Jpol_pla', t, t_ind)

            elif name == 'profile-Itps':
                return self.profiledata('Itps', t, t_ind)

            elif name == 'profile-Itps_av':
                return self.profiledata('Itps_av', t, t_ind)

            elif name == 'profile-Itfs':
                return self.profiledata('Itfs', t, t_ind)

            elif name == 'profile-Itfs_av':
                return self.profiledata('Itfs_av', t, t_ind)

            elif name == 'profile-eccur':
                return self.profiledata('eccurtot', t, t_ind)


    def getPlotsForTimePoint(self, names, t):
        
        toReturn = {}
        for name in names:
            if name in self.__plotNames:
                try:
                    toReturn[name] = self.getSinglePlotForTimePoint(name, t)
                except IndexError as e:
                    print e
                    
        return toReturn
