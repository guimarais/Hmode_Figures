import dd
import numpy as np
import kk_abock
from getsig import getsig
from ddremoveELMData import *

class objview(object):
    def __init__(self, d):
        self.__dict__=d

def neAugped(shotnr=30554, exper='guimas', edition=3, nr_diags=4, elm_exper='AUGD', magdiag='FPP'):
    """Reads the density information of a PED shotfile and returns an object with electron density data
    
    Parameters
    -----------
    shotnr: int
        Number of the shot
    exper: str
        AUGPED only saves private shotfiles, so you must explicit the name of the author of the private shotfile
    edition: int
        Edition of the private shotfile
    nr_diags: int
        How many diagnostics are stored in the private shotfile: still under automation.
    elm_exper: str
        Where to get the ELM data. This is still needed to set 
    
    Returns
    -----------
    The return is a single object with the following entries:
    
    t1: float
        Initial time for the analysis
    t2: float
        Final time for the analysis
    rhos: np.array(float)
    rad: np.array(float)
    dens: np.array(float)
    indi: np.array(int)
    indf: np.array(int)
    neRshift: np.array(float)
    fpgavg: float    
    """
    
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
    if magdiag=='FPP':
        sep_shotfile = 'FPG'
    elif magdiag=='EQH':
        sep_shotfile = 'GQH'
    elif magdiag=='EQI':
        sep_shotfile = 'GQI'
    elif magdiag=='IDE':
        sep_shotfile=='IDG'
    else:#Use FPP by default
        sep_shotfile = 'FPG'
        magdiag = 'FPP'
        print('Unknown magnectic reconstruction shotfile, reverting to FPP!')
        
    fpg = getsig(shotnr, sep_shotfile, 'Raus')
    fpgmsk = ddremoveELMData(shotnr, fpg.time, preft=0.002, suft=0.004, elm_exper=elm_exper)
    fpgind = (fpg.time>=t1.data)&(fpg.time<=t2.data)&fpgmsk
    fpgavg = np.mean(fpg.data[fpgind])
        
    eq = kk_abock.kk()
    eq.Open(shotnr, diag=magdiag)
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