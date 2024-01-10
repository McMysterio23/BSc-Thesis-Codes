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

from FUNZIONI import Kauffmann_03, Kewley_01

def main():

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

    #BPT Digagram del primo tipo : OIII/Hb versus NII/Ha
    """
    log([OIII]/Hb) = 0.61 / (log([NII]/Ha) - 0.05) + 1.3     (Kauffmann+03 line)
    log([OIII]/Hb) = 0.61 / (log([NII]/Ha) - 0.47) + 1.19    (Kewley+01 line)
    """

    #creo e popolo gli array con gli elementi da plottare
    Ordinate = np.ones(452)
    Ascisse = np.ones(452)
    for i in range(452):
        denom_HBeta = TableRAWTabella[ARRNOMI[Index_HBeta_Flux]][i]
        denom_HAlpha = TableRAWTabella[ARRNOMI[Index_HAlpha_Flux]][i]
        
        if denom_HBeta != 0 and denom_HAlpha != 0:
            Ordinate[i] = TableRAWTabella[ARRNOMI[Index_OIII_Flux]][i] / denom_HBeta
            Ascisse[i] = TableRAWTabella[ARRNOMI[Index_NII_Flux]][i] / denom_HAlpha
        else:
            Ordinate[i] = np.nan  # o un altro valore di default
            Ascisse[i] = np.nan   # o un altro valore di default


    ScalaBilogaritmica = True
    Overplotting = True
    plot_scatter(Ascisse, Ordinate, 'BPT Diagram BCG', 'NII/Halpha', 'OIII/Hbeta', ScalaBilogaritmica, Overplotting, Kauffmann_03, Kewley_01)






if __name__ == "__main__":

 main ()