import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table


from libreriaUNO import leggi_nomi_colonne, leggi_file_tsv
from libreriaDUE import plot_and_save_scatter, plot_scatter
from libreriaDUE import Aggiungi_Colonna, Aggiungi_Colonna6, AddColumn, Lettura_Colonne_RawDATA


def main():

    #apro il file con le info sugli SFR
    FILE_PATH_SFRs = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/Star_Formating_Ratios_Intersezione.txt'
    nomi_colonne_SFRs, RAWTabella_SFRs = leggi_file_tsv(FILE_PATH_SFRs)
    TableRAWTabella_SFRs = Table(RAWTabella_SFRs, names = nomi_colonne_SFRs)

    ARRNOMI_SFRs = np.array(nomi_colonne_SFRs)

   #apro il file con le info sugli sSFR
    FILE_PATH_sSFRs = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/Specific_Star_Formating_Ratios_Intersezione.txt'
    nomi_colonne_sSFRs, RAWTabella_sSFRs = leggi_file_tsv(FILE_PATH_sSFRs)
    TableRAWTabella_sSFRs = Table(RAWTabella_sSFRs, names = nomi_colonne_sSFRs)

    ARRNOMI_sSFRs = np.array(nomi_colonne_sSFRs)

    #apro il file con le info del file totLGM.txt
    FILE_PATH_TOTLGM = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/TOTLGM_Intersezione.txt'
    nomi_colonne_TOTLGM, RAWTABELLA_TOTLGM = leggi_file_tsv(FILE_PATH_TOTLGM)
    TableRAWTabella_TOTLGM = Table(RAWTABELLA_TOTLGM, names = nomi_colonne_TOTLGM)

    ARRNOMI_TOTLGM = np.array(nomi_colonne_TOTLGM)

    #apro il file infoLine.txt
    FILE_INFO = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfoLine.txt'
    nomi_colonne_INFO, RAWTabella_INFO = leggi_file_tsv(FILE_INFO)
    TableRAWTabella_INFO = Table(RAWTabella_INFO, names = nomi_colonne_INFO)

    ARRNOMI_INFO = np.array(nomi_colonne_INFO)


    #Realizzo i plot
    sigmab = TableRAWTabella_INFO['SIGMA_BALMER']
    sigmaf = TableRAWTabella_INFO['SIGMA_FORBIDDEN']

    FILE_PATH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfoLine.txt'
    nomi_colonne, RAWTabella = leggi_file_tsv(FILE_PATH)
    TableRAWTabella = Table(RAWTabella, names = nomi_colonne)

    ARRNOMI = np.array(nomi_colonne)



    #Indici di posizione delle colonne contenenti le righe che usiamo 
    Index_OIII_Flux = int(np.where(ARRNOMI == 'OIII_5007_FLUX')[0][0])
    Index_OIII_Flux_Err = int(np.where(ARRNOMI == 'OIII_5007_FLUX_ERR')[0][0])
    Index_OI_Flux = int(np.where(ARRNOMI == 'OI_6300_FLUX')[0][0])
    Index_OI_Flux_Err = int(np.where(ARRNOMI == 'OI_6300_FLUX_ERR')[0][0])
    Index_NII_Flux = int(np.where(ARRNOMI == 'NII_6584_FLUX')[0][0])
    Index_NII_Flux_Err = int(np.where(ARRNOMI == 'NII_6584_FLUX_ERR')[0][0])
    Index_SII6717_Flux = int(np.where(ARRNOMI == 'SII_6717_FLUX')[0][0])
    Index_SII6717_Flux_Err = int(np.where(ARRNOMI == 'SII_6717_FLUX_ERR')[0][0])
    Index_SII6731_Flux = int(np.where(ARRNOMI == 'SII_6731_FLUX')[0][0])
    Index_SII6731_Flux_Err = int(np.where(ARRNOMI == 'SII_6731_FLUX_ERR')[0][0])

    Index_HAlpha_Flux = int(np.where(ARRNOMI == 'H_ALPHA_FLUX')[0][0])
    Index_HAlpha_Flux_Err = int(np.where(ARRNOMI == 'H_ALPHA_FLUX_ERR')[0][0])

    Index_HBeta_Flux = int(np.where(ARRNOMI == 'H_BETA_FLUX')[0][0])
    Index_HBeta_Flux_Err = int(np.where(ARRNOMI == 'H_BETA_FLUX_ERR')[0][0])



    ha = np.ones(452)
    nii = np.ones(452)
    hb = np.ones(452)


    oiii = np.ones(452)
    for i in range(452):
        oiii[i] = TableRAWTabella[nomi_colonne[Index_OIII_Flux]][i]
        ha[i] = TableRAWTabella[nomi_colonne[Index_HAlpha_Flux]][i]
        nii[i] = TableRAWTabella[nomi_colonne[Index_NII_Flux]][i]
        hb[i] = TableRAWTabella[nomi_colonne[Index_HBeta_Flux]][i]
    ii=np.where((oiii > 0) & (hb > 0) & (ha > 0) & (nii > 0))[0]
    oiii=oiii[ii]
    hb=hb[ii]
    ha=ha[ii]
    nii=nii[ii]
    sigmab=sigmab[ii]
    sigmaf=sigmaf[ii]
    logoiiihb=np.log10(oiii/hb)
    logniiha=np.log10(nii/ha)
    i1 = np.where((logniiha < 0.3) & (logoiiihb > 0.8))[0]
    i2 = np.where((logniiha > 0.3) & (logoiiihb > 0.8))[0]
    sigmab_ion=sigmab[i1]
    sigmab_shock=sigmab[i2]
    sigmaf_ion=sigmaf[i1]
    sigmaf_shock=sigmaf[i2]
    """
    plt.hist(logniiha,20)
    plt.hist(logoiiihb,20)
    plt.show() 
    """

    fig,ax = plt.subplots()
    ax.hist(sigmab_ion, 10, color = 'red', label = 'Balmer:ion')
    ax.hist(sigmab_shock, 10, color = 'green', label = 'Balmer:shk')
    ax.legend()
    
    fig,ax = plt.subplots()
    ax.hist(sigmaf_ion, 10, color = 'blue', label = 'Forb:ion')
    ax.hist(sigmaf_shock, 10, color = 'black', label = 'Forb:shk')
    ax.legend()
    plt.show()
    
    
    print(len(sigmab_shock),len(sigmab_shock),len(sigmaf_ion),len(sigmaf_shock))

    

    

    







    

     




if __name__ == "__main__":

 main ()