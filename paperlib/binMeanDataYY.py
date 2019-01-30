import numpy as np

def binMeanDataYY(XX, YY, binpts=100, right=True):
    """For XX and YY arrays, bins the YY data and averages XX data in the bins with actual data"""
    binYY = np.linspace(min(YY),max(YY), binpts)
    bins = np.digitize(YY, binYY, right=right)
    outbinXX = []
    outbinYY = []
    for i in range(binpts):
        dum = bins==i
        dum2 = YY[dum]
        if len(dum2)!=0:
            outbinXX.append(np.mean(XX[dum]))
            outbinYY.append(binYY[i])
    
    return np.array(outbinXX), np.array(outbinYY)