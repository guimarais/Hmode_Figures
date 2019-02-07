import dd
import kk_abock
import numpy as np
from ipfnpytools.trz_to_rhop import trz_to_rhop

class objview(object):
    def __init__(self, d):
        self.__dict__=d

def get_vta(shotnr, tBegin=0.0, tEnd=10.0, magdiag='EQH'):
    """Reads the VTA shotfile and maps data to rho pol"""
    eq = kk_abock.kk()
    eq.Open(shotnr, diag=magdiag)
    
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
    
    zmap_e = np.tile(z_e.data, [len(ne_e.data),1])
    rmap_e = np.tile(r_e.data, [len(z_e.data), 1]).T

    zmap_c = np.tile(z_c.data, [len(ne_c.data),1])
    rmap_c = np.tile(r_c.data, [len(z_c.data), 1]).T

    rho_e = trz_to_rhop(ne_e.time, rmap_e, zmap_e, shot=shotnr, eq='FPP', squeeze=True)
    rho_c = trz_to_rhop(ne_c.time, rmap_c, zmap_c, shot=shotnr, eq='FPP', squeeze=True)

    eq.Close()

    #Organize by increasing rho
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
                    