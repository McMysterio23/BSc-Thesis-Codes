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

from FUNZIONI_BPTs import Kauffmann_03, Kewley_01, main_AGN_line, LINER_Sy2_line, BPT_type1, BPT_type2, BPT_type2_colored, BPTtype1_colored


def main():

    FILE_PATH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/INFOLINE_BCG.fits'
    

    #apro il file Fits
    hdul = fits.open(FILE_PATH)  
    #Matrice tutti dati
    data = hdul[1].data
    
    indici_colonne_interessate = [37, 51, 59, 67, 75, 83, 91]

    #print('I nomi delle colonne attualmente selezionate sono : \n', data.columns.names[37], data.columns.names[51], data.columns.names[59], data.columns.names[67], '\n', 
    #     data.columns.names[75], data.columns.names[83], data.columns.names[91])
    #print('Dovrebbero essere in questo ordine : hb, Oiii, oi, ha, nii, sii6717, sii6731')
    # Inizializza un array vuoto per affiancare le colonne
    
    TabellaSelezionata = np.empty((len(data), len(indici_colonne_interessate)))

    # Estrai e affianca le colonne una per una
    for i, indice_colonna in enumerate(indici_colonne_interessate):
        TabellaSelezionata[:, i] = data.field(indice_colonna)


    print(TabellaSelezionata.shape)
    #questi sono gli array con le righe !!
    oiii = TabellaSelezionata[:, 1]
    oi = TabellaSelezionata[:, 2]
    nii = TabellaSelezionata[:, 4]
    sii6731 = TabellaSelezionata[:, 6]
    sii6717 = TabellaSelezionata[:, 5]
    ha = TabellaSelezionata[:, 3]
    hb = TabellaSelezionata[:, 0]


    #PRIMA TIPOLOGIA DI BPT
    
    #Seleziono gli indici che hanno caratteristiche idonee per la realizzazione del BPT
    ii=np.where((oiii > 0) & (hb > 0) & (ha > 0) & (nii > 0))[0]

    oiii = oiii[ii]
    hb = hb[ii]
    ha = ha[ii]
    nii = nii[ii]


    
    log_oiiihb = np.log10(oiii/hb)
    log_niiha = np.log10(nii/ha)

    #plot_scatter(log_niiha, log_oiiihb, '', '')
    BPTtype1_colored(log_niiha, log_oiiihb)


    #SECONDA TIPOLOGIA DI BPT

    


    # S6731
    #Seleziono gli indici che hanno caratteristiche idonee per la realizzazione del BPT
    sii6731 = sii6731[:396]
    mask = (oiii > 0) & (hb > 0) & (ha > 0) & (sii6731 > 0)
    ii = np.where(mask)[0]
    #ii=np.where((oiii > 0) & (hb > 0) & (ha > 0) & (sii6731 > 0))[0]

    oiii = oiii[ii]
    hb = hb[ii]
    ha = ha[ii]
    sii6731 = sii6731[ii]


    
    log_oiiihb = np.log10(oiii/hb)
    log_sii6731ha = np.log10(sii6731/ha)

    #plot_scatter(log_sii6731ha, log_oiiihb, '', '')
    BPT_type2(log_sii6731ha, log_oiiihb)
 






    
     
            


if __name__ == "__main__":

 main ()