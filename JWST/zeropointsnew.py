import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii

def zeropts(data1):

    F090W = Table.read('MATCHUP.XYMEEE.F090W', format = 'ascii')
    F150W = Table.read('MATCHUP.XYMEEE.F150W', format = 'ascii')
    A1 = (5.1125*F090W['xbar'] - 35680.375 > F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 < F090W['ybar'])
    A2 = (5.1125*F090W['xbar'] - 35680.375 > F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 > F090W['ybar'])
    A3 = (5.1125*F090W['xbar'] - 35680.375 < F090W['ybar']) & (5.7225352112676*F090W['xbar'] - 24000 > F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 < F090W['ybar'])
    A4 = (5.1125*F090W['xbar'] - 35680.375 < F090W['ybar']) & (5.7225352112676*F090W['xbar'] - 24000 > F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 > F090W['ybar'])
    B1 = (5.7225352112676*F090W['xbar'] - 8796 < F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 > F090W['ybar'])
    B2 = (5.7225352112676*F090W['xbar'] - 8796 < F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 < F090W['ybar'])
    B3 = (5.7225352112676*F090W['xbar'] - 8796 > F090W['ybar']) & (5.7225352112676*F090W['xbar'] - 24000 < F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 > F090W['ybar'])
    B4 = (5.7225352112676*F090W['xbar'] - 8796 > F090W['ybar']) & (5.7225352112676*F090W['xbar'] - 24000 < F090W['ybar']) & (-0.16752843846949*F090W['xbar'] + 5793.6608066184 < F090W['ybar'])
    mask = (F090W['xbar'] > 0) & (F150W['xbar'] > 0) & (F090W['msig'] < 9) & (F150W['msig'] < 9)
    CMD = F090W['mbar'][mask]-F150W['mbar'][mask]
    CMDA1 = F090W['mbar'][A1 & mask]-F150W['mbar'][A1 & mask]
    CMDA2 = F090W['mbar'][A2 & mask]-F150W['mbar'][A2 & mask]
    CMDA3 = F090W['mbar'][A3 & mask]-F150W['mbar'][A3 & mask]
    CMDA4 = F090W['mbar'][A4 & mask]-F150W['mbar'][A4 & mask]
    CMDB1 = F090W['mbar'][B1 & mask]-F150W['mbar'][B1 & mask]
    CMDB2 = F090W['mbar'][B2 & mask]-F150W['mbar'][B2 & mask]
    CMDB3 = F090W['mbar'][B3 & mask]-F150W['mbar'][B3 & mask]
    CMDB4 = F090W['mbar'][B4 & mask]-F150W['mbar'][B4 & mask]
    slicemaskA1 = (F150W['mbar'][A1 & mask] > -9) & (F150W['mbar'][A1 & mask] < -8)
    slicemaskA2 = (F150W['mbar'][A2 & mask] > -9) & (F150W['mbar'][A2 & mask] < -8)
    slicemaskA3 = (F150W['mbar'][A3 & mask] > -9) & (F150W['mbar'][A3 & mask] < -8)
    slicemaskA4 = (F150W['mbar'][A4 & mask] > -9) & (F150W['mbar'][A4 & mask] < -8)
    slicemaskB1 = (F150W['mbar'][B1 & mask] > -9) & (F150W['mbar'][B1 & mask] < -8)
    slicemaskB2 = (F150W['mbar'][B2 & mask] > -9) & (F150W['mbar'][B2 & mask] < -8)
    slicemaskB3 = (F150W['mbar'][B3 & mask] > -9) & (F150W['mbar'][B3 & mask] < -8)
    slicemaskB4 = (F150W['mbar'][B4 & mask] > -9) & (F150W['mbar'][B4 & mask] < -8)
    medianA = (np.median(CMDA1[slicemaskA1].flatten()), np.median(CMDA2[slicemaskA2].flatten()), np.median(CMDA3[slicemaskA3].flatten()), np.median(CMDA4[slicemaskA4].flatten()))
    medianB = (np.median(CMDB1[slicemaskB1].flatten()), np.median(CMDB2[slicemaskB2].flatten()), np.median(CMDB3[slicemaskB3].flatten()), np.median(CMDB4[slicemaskB4].flatten()))
    def diff(histogram):
        medianB1 = np.median(CMDB1[slicemaskB1].flatten())
        newhistogram = medianB1 - histogram
        return newhistogram
    difA = []
    difB = []
    for i in medianA:
        difA.append(diff(i))  
    print('The color shifts for chips A1-A4 are:', difA)
    for i in medianB:
        difB.append(diff(i))
    print('The color shifts for chips B1-B4 are:', difB)
    CMDA1s = CMDA1 + 0.11494999999999944
    CMDA2s = CMDA2 + 0.07399999999999984
    CMDA3s = CMDA3 + 0.12340000000000018
    CMDA4s = CMDA4 + 0.06705000000000005
    CMDB2s = CMDB2 + 0.005550000000001276
    CMDB3s = CMDB3 + 0.005350000000000854
    CMDB4s = CMDB4 - 0.041649999999998855
    
    fig, axs = plt.subplots(1, 2, sharex = True, sharey = True, figsize=(20, 10))
    axs[0].scatter(CMD, F150W['mbar'][mask], s = .01, color = 'black')
    axs[1].scatter(CMDA1s, F150W['mbar'][A1 & mask], s = .01, color = 'crimson', label = ("A1"))
    axs[1].scatter(CMDA2s, F150W['mbar'][A2 & mask], s = .01, color = 'yellow', label = ("A2"))
    axs[1].scatter(CMDA3s, F150W['mbar'][A3 & mask], s = .01, color = 'paleturquoise', label = ("A3"))
    axs[1].scatter(CMDA4s, F150W['mbar'][A4 & mask], s = .01, color = 'greenyellow', label = ("A4"))
    axs[1].scatter(CMDB1, F150W['mbar'][B1 & mask], s = .01, color = 'deepskyblue', label = ("B1"))
    axs[1].scatter(CMDB2s, F150W['mbar'][B2 & mask], s = .01, color = 'mediumslateblue', label = ("B2"))
    axs[1].scatter(CMDB3s, F150W['mbar'][B3 & mask], s = .01, color = 'mediumspringgreen', label = ("B3"))
    axs[1].scatter(CMDB4s, F150W['mbar'][B4 & mask], s = .01, color = 'violet', label = ("B4"))
    axs[1].set_xlim([-.5,1])
    axs[1].set_ylim([-12,-1])
    axs[1].invert_yaxis()
    axs[1].set_xlabel('Color', size = 20)
    axs[0].set_xlabel('Color', size = 20)
    axs[0].set_ylabel('Magnitude', size = 20)

