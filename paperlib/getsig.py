import dd
#import dd_20180216 as dd

def getsig(shotnr, diag, sig, tBegin=0.0, tEnd=10.0, exper='AUGD', edition=0):
    """Returns an object with the shotfile data
    
    Parameters
    -----------
    shotnr: int
        Number of the discharge
    diag: str
        Three letter code for the diagnostic
    sig: str
        Name of signal or signal group.
    tBegin: float
        Beginning of analysis window
    tEnd: float
        End of analysis window
    exper: str
        Experiment where to read the data. For private shotfiles it is the username in question.
    edition: int
        Version of the shotfile to be opened
        
    Returns
    -----------
    An object with the data respective to the signal or signal group. Look for more info in the DD library wrapper.
    """
    shtfl = dd.shotfile(diag, shotnr, experiment=exper, edition=edition)
    dta = shtfl(sig, tBegin=tBegin, tEnd=tEnd)
    shtfl.close()
    return dta
