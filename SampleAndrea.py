import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column



def main():
    file="/Users/andreamaccarinelli/desktop/SDSS:DR7/gal_info_dr7_v5_2.fits"
    obj=fits.open(file)
    data=obj[1].data
    z=data["Z"]
    vdisp=data["V_DISP"]
    plt.scatter(z,vdisp,color="green",s=40,linestyle='None', marker="o",edgecolor='black', linewidth=2)
    vdisp2 = vdisp[np.where((vdisp != np.min(vdisp)) & (vdisp != np.max(vdisp)))]
    z2=z[np.where((vdisp != np.min(vdisp)) & (vdisp != np.max(vdisp)))]
    plt.scatter(z2,vdisp2,color="green",s=40,linestyle='None', marker="o",edgecolor='black', linewidth=2)
    plt.show()

if __name__ == "__main__":

 main ()
