import matplotlib.pyplot as plt

def plot4DiagPedSht(pedStruct,
                   rho=True,
                   zorder=[1,2,3,4],
                   labels=[r'$\mathrm{TS_{edge}}$', 'LIN', 'Ref. LFS', 'Ref. HFS']):
    """Takes a structure returned by pedshtf and plots it"""
    
    dotsize = 5
    dotsizeref = 4

    fig, ax = plt.subplots(figsize=(6,3), ncols=1, nrows=1, dpi=200)
    
    #pedStruct = pedshtf(shotnr=shotnr, exper=exper, edition=edition, nr_diags=4, elm_exper=elm_exper)
    ax.scatter(pedStruct.rad[pedStruct.indi[0]:pedStruct.indf[0]],pedStruct.dens[pedStruct.indi[0]:pedStruct.indf[0]]*1e-19,
               s=dotsize, zorder=2, label=labels[0],facecolors="C1", marker='o')
    #LIN
    ax.scatter(pedStruct.rad[pedStruct.indi[1]:pedStruct.indf[1]], pedStruct.dens[pedStruct.indi[1]:pedStruct.indf[1]]*1e-19,
               s=dotsize, zorder=3, label=labels[1],facecolors='C2', marker='D')
    #LFS
    lrad = pedStruct.rad[pedStruct.indi[2]:pedStruct.indf[2]]
    ldens = pedStruct.dens[pedStruct.indi[2]:pedStruct.indf[2]]*1e-19
    ax.scatter(lrad, ldens, s=dotsizeref, zorder=1, label=labels[2], facecolors='none', edgecolors="C0", marker='o')
    #Clean up HFS ref
    hrad = pedStruct.rad[pedStruct.indi[3]:pedStruct.indf[3]]
    hdens = pedStruct.dens[pedStruct.indi[3]:pedStruct.indf[3]]*1e-19
    ax.scatter(hrad, hdens, s=dotsizeref, zorder=0, label=labels[3],facecolors='none', edgecolors="C3", marker='o')

    #Separatrixes
    ax.axvline(pedStruct.fpgavg, color='k', lw=0.7)
    ax.axvspan(2.1, pedStruct.fpgavg, color='k', alpha=0.2)

    ax.set_xlim(left=2.1)
    ax.set_xlabel(r'$\mathrm{Major\,Radius\,[m]}$')
    
    
    ax.set_ylim(0,6)
    

    ax.set_title('t=[%0.1f,%0.1f]s'%(pedStruct.t1,pedStruct.t2), loc='left', fontsize=9)
    #ax.text(2.105,5.6, '', color='k')
    
    plt.show()
    return 0