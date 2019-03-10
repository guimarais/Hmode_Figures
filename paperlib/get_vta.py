import dd
import numpy as np
from ipfnpytools.trz_to_rhop import trz_to_rhop

class objview(object):
    def __init__(self, d):
        self.__dict__=d

def get_vta(shotnr, tBegin=0.0, tEnd=10.0, magdiag='EQH', verbose=False):
    """Reads the VTA(Vertical Thomson Array) shotfile and maps data to rho pol. The Edge and Core systems have different time bases so they are treated as independent systems.
    
    Parameters
    -----------
    shotnr: int
        Number of the shot
    tBegin: float
        Beginning of the time window for the data block
    tEnd: float
        End of the time window for the data block
    magdiag: string ('EQH', 'EQI', 'FPP', 'IDE')
        String for the magnetic equilibria to be used in mapping data to rho poloidal
    verbose: bool
        Flag to pass to trz_to_rhop to print the rho_pol progess
    
    Returns
    -----------
    vta: object
        Object containing data from the core and edge VTA. Data is already sorted by ascending rho_pol values.
        For each system, Core (c) and Edge (e), the vta object has the following elements:
        time_#: Timebase of the e/c system.
        rho_#:  Rho_pol of the e/c system.
        ne_#:   Density of the e/c system.
        Te_#:   Temperature of the e/c system.        
        
    Example
    -----------
    vta = get_vta(30733, tBegin=1.0, tEnd=1.5, magdiag='FPP')
    """

    #Reads data from the VTA shotfile
    vta = dd.shotfile('VTA', shotnr)
    ne_c = vta('Ne_c', tBegin=tBegin, tEnd=tEnd)
    te_c = vta('Te_c', tBegin=tBegin, tEnd=tEnd)
    r_c  = vta('R_core', tBegin=tBegin, tEnd=tEnd)
    z_c  = vta('Z_core', tBegin=tBegin, tEnd=tEnd)
    
    ne_e = vta('Ne_e', tBegin=tBegin, tEnd=tEnd)
    te_e = vta('Te_e', tBegin=tBegin, tEnd=tEnd)
    r_e  = vta('R_edge', tBegin=tBegin, tEnd=tEnd)
    z_e  = vta('Z_edge', tBegin=tBegin, tEnd=tEnd)
    vta.close()
    
    #Adjusts R and Z dimensions to use afterwards in 'trz_to_rhop'
    ## Edge
    rmap_e = np.tile(r_e.data, [len(z_e.data), 1]).T
    zmap_e = np.tile(z_e.data, [len(ne_e.data),1])
    
    ## Core
    zmap_c = np.tile(z_c.data, [len(ne_c.data),1])
    rmap_c = np.tile(r_c.data, [len(z_c.data), 1]).T

    #Converts R and Z to rho_pol
    rho_e = trz_to_rhop(ne_e.time, rmap_e, zmap_e, shot=shotnr, eq='FPP', squeeze=True, verbose=verbose)
    rho_c = trz_to_rhop(ne_c.time, rmap_c, zmap_c, shot=shotnr, eq='FPP', squeeze=True, verbose=verbose)

    #Sort by ascending rho
    #Edge
    rr_e = np.zeros_like(rho_e)
    nn_e = np.zeros_like(rho_e)
    tt_e = np.zeros_like(rho_e)

    for i in range(len(rho_e)):
        dum = np.argsort(rho_e[i])
        rr_e[i,:] = rho_e[i, dum]
        nn_e[i,:] = ne_e.data[i,dum]
        tt_e[i,:] = te_e.data[i,dum]

    #Core
    rr_c = np.zeros_like(rho_c)
    nn_c = np.zeros_like(rho_c)
    tt_c = np.zeros_like(rho_c)

    for i in range(len(rho_c)):
        dum = np.argsort(rho_c[i])
        rr_c[i,:] = rho_c[i, dum]
        nn_c[i,:] = ne_c.data[i,dum]
        tt_c[i,:] = te_c.data[i,dum]
    
    return objview({'time_e': np.array(ne_e.time),
                    'rho_e': np.array(rr_e),
                    'ne_e': np.array(nn_e),
                    'Te_e': np.array(tt_e),
                    'time_c': np.array(ne_c.time),
                    'rho_c': np.array(rr_c),
                    'ne_c': np.array(nn_c),
                    'Te_c': np.array(tt_c)})
                    