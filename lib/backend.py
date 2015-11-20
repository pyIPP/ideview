import sys
import numpy as np
from IPython import embed
from copy import copy,deepcopy
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
from lib.bunch import PlotBunch
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

    # List of plots sorted by kind
    # Respect the order!
    __plotNames = ['profile-pressure',
                   'profile-pcon',
                   'profile-pressure/pcon',
                   'profile-q',
                   'profile-ECcur_tot+j_BS+j_nbcd',
                   'profile-ECcur_gyr',
                   'profile-Jpol_tot',
                   'profile-Jpol_pla',
                   'profile-Itps',
                   'profile-Itps_av',
                   'profile-Itfs',
                   'profile-Itfs_av',
                   'profile-pol',
                   'profile-mse',
                   'profile-I_torcon',
                   'profile-Bprob',
                   'profile-Dpsi',
                   'profile-I_tor',
                   'profile-Iext',
                   'profile-T_e+T_i',
                   'profile-n_e+n_i',
                   'profile-sigma',
                   'profile-ncft',
                   'profile-Zeff']

    __plotNames += ['res(profile)-Dpsi',
                    'res(profile)-Iext',
                    'res(profile)-pcon',
                    'res(profile)-Bprob']

    __plotNames += ['trace-Wmhd',
                    'trace-Bprob',
                    'trace-Dpsi',
                    'trace-Iext',
                    'trace-Rmag',
                    'trace-Zmag',
                    'trace-Rin',                   
                    'trace-Raus',
                    'trace-betapol',
                    'trace-betapol+li/2',
                    'trace-Itor',
                    'trace-Rxpu',
                    'trace-Zxpu',
                    'trace-ahor',
                    'trace-bver',
                    'trace-bver/ahor',
                    'trace-XPfdif',
                    'trace-delR',
                    'trace-delZ',
                    'trace-q0',
                    'trace-q25',
                    'trace-q50',
                    'trace-q75',
                    'trace-q95',
                    'trace-delR_oben',
                    'trace-delR_unten',
                    'trace-k_oben',
                    'trace-k_unten',
                    'trace-dRxP',
                    'trace-eccd_tot',
                    'trace-Itax',
                    'trace-ecrhmax(R)',
                    'trace-ecrhmax(z)',
                    'trace-ecrhmax(y)',
                    'trace-ecrhmax(rho)',
                    'trace-li',
                    'trace-Vol',
                    'trace-pol',
                    'trace-mse',
                    'trace-Shear_q1']

    __plotNames += ['timecontour-pressure',
                    'timecontour-q',
                    'timecontour-Dpsi',
                    'timecontour-Iext',
                    'timecontour-pcon',
                    'timecontour-pol',
                    'timecontour-mse',
                    'timecontour-I_tor',
                    'timecontour-I_torcon',
                    'timecontour-Bprob']

    __plotNames += ['contour-pfl',
                    'contour-rho']


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
            if MES.area is None:
                lookforab = '%s_rp'%name
                MES.area = self.getData(lookforab)
                if name == 'q_sa':
                    MES.area = self.getData('rhopol_q')
                    
                if not onlyMES and FIT is not None and FIT.area is None:
                    FIT.area = self.getData(lookforab)

            if name == 'It':
                lookforab = 'R_Itor' # 'rhop_It'
                MES.area = self.getData(lookforab)


            if onlyMES:
                return MES
            else:
                return MES, FIT, UNC
        except Exception as e:
            print e, '\nSomething went wrong while looking for the data'
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
                data=[{'x': MES_Area_data, 'y': MES_Data, 'marker':'', 'ls':'-', 'c':'k'}]

            
            if (not UNC is None) and unc:
                UNC_Data = UNC.data[t_ind]
                if 'con' in name:
                    data.extend([{'x': MES_Area_data, 'y': MES_Data + UNC_Data, 'ls':'--', 'c':'k'},
                                 {'x': MES_Area_data, 'y': MES_Data-UNC_Data, 'ls':'--', 'c':'k'}])
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

                data = [{'x':  qsa.area.data[t_ind], 'y': q[t_ind], 'c':'k'},{'x': XIQ1}, {'x': XIQ2},{'x':XIQ3_2},
                      {'y': 1,'c':'k','alpha':.3}, {'y': 2,'c':'k','alpha':.3},{'y': 1.5,'c':'k','alpha':.3}]
                
                qsap = self.getData('q_sa_plu')
                qsap.area = self.getData('rhopol_q')
                qsam = self.getData('q_sa_min')
                qsam.area = self.getData('rhopol_q')
                if not qsap is None and not qsam is None:
                    #qsap = qsap(tBegin=t, tEnd=t)
                    #qsam = qsam(tBegin=t, tEnd=t)
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

            #elif name == 'It':
            #    embed()
            #    Valy = np.array([])
            #    Valx = np.array([])
            #    Unc = np.array([])
            #    MIN = 1
            #    data = [
                #for i in range(MES_Area_data.shape[0]):
                #    if MES_Area_data[i] < MIN:
                #        I = i
                #        MIN = MES_Area_data[i]
                #for i in range(I):
                #    if (I-i > 0) and (I+i < MES_Area_data.shape[0]):
                #        MEANy = (MES_Data[I+i]+MES_Data[I-i])/2
                #        MEANx = (MES_Area_data[I+i]+MES_Area_data[I-i])/2
                #        MEANunc = (UNC_Data[I+i]+UNC_Data[I-i])/2
                #        Valy = np.append(Valy,[MEANy])
                #        Valx = np.append(Valx,[MEANx])
                #        Unc = np.append(Unc,[MEANunc])
                #data = [{'x':Valx, 'y':Valy, 'marker':'', 'ls':'-', 'c':'k'}]
                #
                #data.extend([{'x': Valx, 'y': Valy + Unc, 'ls':'--'},
                #                 {'x': Valx, 'y': Valy- Unc, 'ls': '--'}])
                #XMIN = 0
                #XMAX = 1.2

            elif name == 'eccurtot':
                data = []
                eccurtot = self.getData('eccurtot')
                j_BS = self.getData('cde_bs')
                j_nbcd = self.getData('cde_nbcd')
                data.extend([{'x' : MES_Area_data, 'y': eccurtot.data[t_ind], 'marker':'', 'ls':'-', 'c':'k', 'label':'ECcur_tot', 'exc': True}])
                data.extend([{'x' : MES_Area_data, 'y': j_BS.data[t_ind], 'marker':'', 'ls':'-', 'c':'b', 'label':'j_BS', 'exc': True}])
                data.extend([{'x': MES_Area_data, 'y': j_nbcd.data[t_ind], 'marker':'', 'ls':'-', 'c':'r', 'label':'j_nbcd', 'exc': True}])
                #if (YMIN > np.nanmin(j_BS.data)):
                #    YMIN = np.nanmin(j_BS.data)
                #if (YMIN > np.nanmin(j_nbcd.data)):
                YMIN = np.nanmin(eccurtot.data+j_nbcd.data + j_BS.data)
                #if (YMAX < np.nanmax(j_BS.data)):
                #    YMAX = np.nanmax(j_BS.data)
                #if (YMAX < np.nanmax(j_nbcd.data)):
                YMAX = np.nanmax(eccurtot.data+j_nbcd.data + j_BS.data)
                if YMIN < -1.5e6: YMIN = -1.5e6
                if YMAX > 1.5e6: YMAX = 1.5e6
                #print YMIN, YMAX
                #embed()
                XMIN = 0
                XMAX = 1.2

            elif name == 'cde_te':
                data = []
                YMAX = 0
                YMIN = 0
                T_e = self.getData('cde_te')
                T_i = self.getData('cde_ti')
                Area_base = T_e.area.data[t_ind]
                DataT_e = T_e.data[t_ind]
                DataT_i = T_i.data[t_ind]
                data.extend([{'x':Area_base, 'y':DataT_e, 'marker':'', 'ls':'-', 'c':'k', 'label':'T_e', 'exc': True}])
                data.extend([{'x':Area_base, 'y':DataT_i, 'marker':'', 'ls':'-', 'c':'b', 'label':'T_i', 'exc': True}])
                XMAX = T_e.area.data.max() if (T_e.area.data.max() < 1000) else 1
                XMIN = T_e.area.data.min() if (T_e.area.data.min() > -1000) else 0
                Ted = np.ravel(T_e.data)
                Tid = np.ravel(T_i.data)
                Tedsorted = list(Ted[np.isfinite(Ted)][Ted[np.isfinite(Ted)].argsort()])
                Tidsorted = list(Tid[np.isfinite(Tid)][Tid[np.isfinite(Tid)].argsort()])
                while Tedsorted[-1] > 1e5:
                    del Tedsorted [-1]
                while Tidsorted[-1] > 1e5: # avoiding of broken values that are not inf or NaN
                    del Tidsorted [-1]
                YMIN = 0
                YMAX = max(Tedsorted[-1], Tidsorted[-1])

            elif name == 'cde_ne':
                data = []
                YMAX = 0
                YMIN = 0
                n_e = self.getData('cde_ne')
                n_i = self.getData('cde_ni')
                Area_base = n_e.area.data[t_ind]
                Datan_e = n_e.data[t_ind]
                Datan_i = n_i.data[t_ind]
                data.extend([{'x':Area_base, 'y':Datan_e, 'marker':'', 'ls':'-', 'c':'k', 'label':'n_e', 'exc': True}])
                data.extend([{'x':Area_base, 'y':Datan_i, 'marker':'', 'ls':'-', 'c':'b', 'label':'n_i', 'exc': True}])
                XMAX = n_e.area.data.max() if (n_e.area.data.max() < 1000) else 1
                XMIN = n_e.area.data.min() if (n_e.area.data.min() > -1000) else 0
                #if YMAX < np.nanmax(n_e.data):
                #    YMAX = np.nanmax(n_e.data)
                #if YMAX < np.nanmax(n_i.data):
                #    YMAX = np.nanmax(n_i.data)
                YMIN = 0
                ned = np.ravel(n_e.data)
                nid = np.ravel(n_i.data)
                YMAX = max(ned[np.isfinite(ned)][ned[np.isfinite(ned)].argsort()][-1],
                           nid[np.isfinite(nid)][nid[np.isfinite(nid)].argsort()][-1])

            elif name == 'cde_zeff':
                YMIN = 0
                YMAX = 10
                # Broken value (for shot RRF:32297:4 at t = 2.28) destroys scale otherwise. Same for ncft

            elif name == 'ncft':
                YMIN = 0
                YMAX = 1
                
            elif name == 'sig_snc':
                YMIN = 0
                a = np.ravel(MES.data)
                YMAX = a[np.isfinite(a)][a[np.isfinite(a)].argsort()][-1]

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
            if not 'signalGroup' in str(type(MES)): # SIGNALS in here
                if 'q' in name and name != 'shear_q1':
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
                elif name == 'Vol':
                    data=[{'x': MES.time, 'y': MES.data, 'ls': '-'},{'x':t, 'c':'k'}]
                    return PlotBunch(kind='trace',data=data)
                elif name == 'betpol':
                    data=[{'x': MES.time, 'y': MES.data, 'ls': '-'},{'x': t,'c': 'k'}]
                    if UNC is not None:
                        data.extend([{'x': MES.time, 'y': MES.data+UNC.data, 'ls': '--'},
                                     {'x': MES.time, 'y': np.maximum(0,MES.data-UNC.data), 'ls': '--'}])
                    return PlotBunch(kind='trace',data=data, 
                        setting={'ylim':(0,min(3, max(MES.data+UNC.data)))})
                else:
                    data=[{'x': MES.time, 'y': MES.data, 'ls': '-'},{'x': t,'c': 'k'}]
                    if UNC is not None:
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
                elif ecrh == 'z':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,1,i], 'label':'gyro%i'%(i+1), 'c': colour, 'exc' : True})
                elif ecrh == 'y':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,2,i], 'label':'gyro%i'%(i+1), 'c': colour, 'exc' : True})
                elif ecrh == 'rho':
                    for i in range(MES.data.shape[2]):
                        colour = colours[i]
                        data.append({'x':MES.time, 'y':MES.data[:,3,i], 'label':'gyro%i'%(i+1), 'c': colour, 'exc' : True})
                return PlotBunch(kind='trace',data=data)

            elif name == 'mse':
                data = [{'x':t,'c':'k'}]
                for j in range(MES.data.shape[1]):
                    Mesarray = np.array([])
                    Fitarray = np.array([])
                    for i in range(MES.data.shape[0]):
                        Mes = MES.data[i][j]
                        Fit = FIT.data[i][j]
                        Mesarray = np.append(Mesarray,[Mes])
                        Fitarray = np.append(Fitarray,[Fit])                        
                    
                    data.append({'x':MES.time, 'y':Mesarray, 'ls':'-'})
                    data.append({'x':MES.time, 'y':Fitarray, 'ls':'-', 'c':'r'})
                    
                return PlotBunch(kind='trace',data=data)

            elif name == 'pol':
                data = [{'x':t,'c':'k'}]
                #Mesarray = np.array([])
                #Fitarray = np.array([])
                #Tracearray = np.array([])
                #print 'pol', MES.data.shape
                #for i in range(MES.data.shape[0]):
                #    Mes = MES.data[i][1]
                #    Fit = FIT.data[i][1]
                #    Trace = (MES.data[i][1]-FIT.data[i][1])/UNC.data[i][1]
                #    Mesarray = np.append(Mesarray,[Mes])
                #    Fitarray = np.append(Fitarray,[Fit])
                #    Tracearray = np.append(Tracearray, [Trace])
                chans = [0,1]
                Mesarray = MES.data[:, chans]
                Fitarray = FIT.data[:, chans]
                Tracearray = (MES.data[:, chans] - FIT.data[:, chans])/UNC.data[:, chans]
                
                for i in range(len(chans)):                
                    data.append({'x':MES.time, 'y':Mesarray[:, i], 'ls':'-', 'c':'k', 'exc':True})
                    data.append({'x':MES.time, 'y':Fitarray[:, i], 'ls':'-', 'c':'r', 'exc':True})
                    data.append({'x':MES.time, 'y':Tracearray[:, i], 'ls':'-', 'c':'green', 'exc':True})
                    if i == 0:
                        data[-1]['label'] = 'res'
                        data[-2]['label'] = 'fit'
                    
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
            for diag in self.equ:
                if name == 'Vol' and 'IDG' in diag.getObjectNames().values():
                    try:
                        self._cache[name] = copy(diag(name))
                    except Exception:
                        pass
                    break

                else:
                    if name in diag.getObjectNames().values() and name != 'Vol':
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
            
            tmp = self.getData('pfm', t)
            Ri = tmp['Ri']; Zj = tmp['zj']; pfm = tmp['pfm']

            psiAx, psiSep = self.getData('PFxx').data[t_index, :2]        
            
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
            elif name == 'trace-Vol':
                return self.tracedata('Vol', t)
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
            elif name == 'trace-pol':
                return self.tracedata('pol', t)
            elif name == 'trace-mse':
                return self.tracedata('mse', t)
            elif name == 'trace-Shear_q1':
                return self.tracedata('shear_q1', t)

        elif 'profile' in name:
            if name == 'profile-q':
                return self.profiledata('q_sa', t, t_ind)
                
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
                maxpres = np.nanmax(p.data)
                pconunc = self.getData('pcon_unc')(tBegin=t, tEnd=t)
                
                data=[{'x': p.area.data, 'y': p.data, 'ls': '-', 'label':'pressure', 'exc': True},
                      {'x': pcon.area.data, 'y': pcon.data, 'ls':'', 'marker':'+', 'mew':2, 'ms':8, 'c':'k', 'label':'constraint', 'exc': True}]
                
                if not punc is None:
                    punc = punc(tBegin=t, tEnd=t)
                    data.extend([{'x': p.area.data, 'y': p.data+punc.data, 'ls': '--', 'label':'uncertainty', 'c':'b', 'exc': True}, {'x': p.area.data, 'y': p.data-punc.data, 'ls': '--', 'c':'b', 'exc': True}])
                if not pconunc is None:
                    data.extend([{'x': pcon.area.data, 'y': pcon.data+pconunc.data,'marker':'_', 'mew': 2, 'ms':8, 'ls': '', 'label':'uncertainty', 'c':'k', 'exc': True}, {'x': pcon.area.data, 'y': pcon.data-pconunc.data,'marker':'_', 'mew':2, 'ms':8, 'ls': '', 'c':'k', 'exc': True}])
                
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

            elif name == 'profile-ECcur_tot+j_BS+j_nbcd':
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

            elif name == 'profile-T_e+T_i':
                return self.profiledata('cde_te', t, t_ind)

            elif name == 'profile-n_e+n_i':
                return self.profiledata('cde_ne', t, t_ind)

            elif name == 'profile-sigma':
                return self.profiledata('sig_snc', t, t_ind)

            elif name == 'profile-ncft':
                return self.profiledata('ncft', t, t_ind)

            elif name == 'profile-Zeff':
                return self.profiledata('cde_zeff', t, t_ind)

    def getPlotsForTimePoint(self, names, t):
        
        toReturn = {}
        for name in names:
            if name in self.__plotNames:
                try:
                    toReturn[name] = self.getSinglePlotForTimePoint(name, t)
                except IndexError as e:
                    print e
                    
        return toReturn
