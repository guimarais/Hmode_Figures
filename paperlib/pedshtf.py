import dd
import numpy as np
import kk_abock
from getsig import getsig
from ddremoveELMData import *

class objview(object):
    def __init__(self, d):
        self.__dict__=d

def pedshtf(shotnr=30554, exper='guimas', edition=3, nr_diags=4, elm_exper='AUGD'):
    """Reads a PED shotfile and returns an object with electron density data"""
    
    ped = dd.shotfile('PED', shotnr, experiment='guimas', edition=edition)
    t1 = ped('t1')
    t2 = ped('t2')
    nedata = ped('neData')
    mskz = nedata.data!=0.0
    rhos = nedata.area.data[0,mskz]
    dens = nedata.data[mskz]
    diagind = ped('DiagIndx')
    nedpts = ped('neDPts')
    indexesr = nedpts.data[nedpts.data!=0]
    tstart = ped('tstart')
    tstop = ped('tstop')
    neRshift = ped('neRshift')
    ped.close()
    
    #Check the intervals
    indi = np.zeros(nr_diags).astype(int)
    indf = np.zeros(nr_diags).astype(int)
    #Initialize arrays
    indf[0] = np.array(indexesr[0]-1).astype(int)

    for i in range(1,nr_diags):
        indi[i] = np.sum(indexesr[0:i]).astype(int)
        indf[i] = np.sum(indexesr[0:i+1]).astype(int)-1
        
    #Gets the average separatrix position
    fpg = getsig(shotnr, 'FPG', 'Raus')
    fpgmsk = ddremoveELMData(shotnr, fpg.time, preft=0.002, suft=0.004, elm_exper=elm_exper)
    fpgind = (fpg.time>=t1.data)&(fpg.time<=t2.data)&fpgmsk
    fpgavg = np.mean(fpg.data[fpgind])
    
    import kk_abock
    eq = kk_abock.kk()
    eq.Open(shotnr, diag='FPP')
    radius = eq.rhopol_to_Rz((np.float(t1.data)+np.float(t2.data))/2.0, rhos, 0.0)
    rad = radius['R']
    eq.Close()
    
    return objview({'t1':np.array(t1.data),
                    't2':np.array(t2.data),
                    'rhos':np.array(rhos),
                    'rad':np.array(rad),
                    'dens':np.array(dens),
                    'indi':np.array(indi),
                    'indf':np.array(indf),
                    'neRshift':np.array(neRshift.data),
                    'fpgavg':np.array(fpgavg)})