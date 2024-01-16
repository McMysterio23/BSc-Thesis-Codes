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
    print('I nomi delle colonne disponibili sono : \n', data.columns.names)  #Dai un'occhiata alle keywords
    #Per esempio, se vogliamo associate a sigma_bal la colonna chiamara "SIGMA_BALMER"
    #sigma_bal = data["SIGMA_BALMER"]

    Ascisse = np.ones(len(data['SIGMA_BALMER']))
    Ordinate = np.ones(len(data['SIGMA_BALMER']))

    IndiciSelezionati =[38, 52, 60, 68, 76, 84, 92] 
    TabellaSelezionata = data[:,IndiciSelezionati]
    print(TabellaSelezionata.columns.names)

    

    





    #BPT Diagram del primo tipo : OIII/Hb versus NII/Ha
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

    Comando = int(sys.argv[1])
    if Comando == 1:
       ScalaBilogaritmica = True
    else:
       ScalaBilogaritmica = False
    
    Overplotting = True
    #Faccio il logaritmo dei due array Ascisse e Ordinate 
    LOGAscisse = np.ones(452)
    LOGOrdinate = np.ones(452)

    for i in range(452):
        if Ascisse[i] > 0:
            LOGAscisse[i] = np.log10(Ascisse[i])
        else:
            LOGAscisse[i] = np.nan  # o un altro valore di default

        if Ordinate[i] > 0:
            LOGOrdinate[i] = np.log10(Ordinate[i])
        else:
            LOGOrdinate[i] = np.nan  # o un altro valore di default
 
    size = float(sys.argv[2])
    #plot_scatter(LOGAscisse, LOGOrdinate,Kauffmann_03, Kewley_01, 'BPT Diagram Type1 of the BCGs', size, 'log(NII/Halpha)', 'log(OIII/Hbeta)', ScalaBilogaritmica, Overplotting, 'Kauffman 03', 'Kewley 01')
    
    BPT_type1(LOGAscisse, LOGOrdinate, True)
    BPT_type1(LOGAscisse, LOGOrdinate)
    BPTtype1_colored(LOGAscisse, LOGOrdinate)
    




    #BPT Diagram del secondo tipo 
    """
    log([OIII]/Hb) = 0.72 / (log([SII]/Ha) - 0.32) + 1.30    (main AGN line)
    log([OIII]/Hb) = 1.89 log([SII]/Ha) + 0.76   (LINER/Sy2 line)
    """

    #creo e popolo gli array con gli elementi da plottare
   
    for i in range(452):
        denom_HBeta = TableRAWTabella[ARRNOMI[Index_HBeta_Flux]][i]
        denom_HAlpha = TableRAWTabella[ARRNOMI[Index_HAlpha_Flux]][i]
        
        if denom_HBeta != 0 and denom_HAlpha != 0:
            Ordinate[i] = TableRAWTabella[ARRNOMI[Index_OIII_Flux]][i] / denom_HBeta
            Ascisse[i] = TableRAWTabella[ARRNOMI[Index_SII6717_Flux]][i] / denom_HAlpha
        else:
            Ordinate[i] = np.nan  # o un altro valore di default
            Ascisse[i] = np.nan   # o un altro valore di default

    Comando = int(sys.argv[1])
    if Comando == 1:
       ScalaBilogaritmica = True
    else:
       ScalaBilogaritmica = False
    
    Overplotting = True
    #Faccio il logaritmo dei due array Ascisse e Ordinate 
    LOGAscisse = np.ones(452)
    LOGOrdinate = np.ones(452)

    for i in range(452):
        if Ascisse[i] > 0:
            LOGAscisse[i] = np.log10(Ascisse[i])
        else:
            LOGAscisse[i] = np.nan  # o un altro valore di default

        if Ordinate[i] > 0:
            LOGOrdinate[i] = np.log10(Ordinate[i])
        else:
            LOGOrdinate[i] = np.nan  # o un altro valore di default
 
    size = float(sys.argv[2])
    plot_scatter(LOGAscisse, LOGOrdinate,main_AGN_line, LINER_Sy2_line, 'BPT Diagram Type2 of the BCGs', size, 'log(SII-6717A/Halpha)', 'log(OIII/Hbeta)', ScalaBilogaritmica, Overplotting, 'Main AGN Line', 'LINER Sy2 LINE')
    BPT_type2(LOGAscisse, LOGOrdinate, True)
    BPT_type2(LOGAscisse, LOGOrdinate)
    BPT_type2_colored(LOGAscisse, LOGOrdinate)




    
#Secondo Type2 
    for i in range(452):
        denom_HBeta = TableRAWTabella[ARRNOMI[Index_HBeta_Flux]][i]
        denom_HAlpha = TableRAWTabella[ARRNOMI[Index_HAlpha_Flux]][i]
        
        if denom_HBeta != 0 and denom_HAlpha != 0:
            Ordinate[i] = TableRAWTabella[ARRNOMI[Index_OIII_Flux]][i] / denom_HBeta
            Ascisse[i] = TableRAWTabella[ARRNOMI[Index_SII6731_Flux]][i] / denom_HAlpha
        else:
            Ordinate[i] = np.nan  # o un altro valore di default
            Ascisse[i] = np.nan   # o un altro valore di default

    Comando = int(sys.argv[1])
    if Comando == 1:
       ScalaBilogaritmica = True
    else:
       ScalaBilogaritmica = False
    
    Overplotting =True
    #Faccio il logaritmo dei due array Ascisse e Ordinate 
    LOGAscisse = np.ones(452)
    LOGOrdinate = np.ones(452)

    for i in range(452):
        if Ascisse[i] > 0:
            LOGAscisse[i] = np.log10(Ascisse[i])
        else:
            LOGAscisse[i] = np.nan  # o un altro valore di default

        if Ordinate[i] > 0:
            LOGOrdinate[i] = np.log10(Ordinate[i])
        else:
            LOGOrdinate[i] = np.nan  # o un altro valore di default
 
    size = float(sys.argv[2])
    #plot_scatter(LOGAscisse, LOGOrdinate, main_AGN_line, LINER_Sy2_line, 'BPT Diagram Type2 of the BCGs', size, 'log(SII-6731A/Halpha)', 'log(OIII/Hbeta)', ScalaBilogaritmica, Overplotting, 'Main AGN Line', 'LINER Sy2 LINE') 
    BPT_type2(LOGAscisse, LOGOrdinate, True)
    BPT_type2(LOGAscisse, LOGOrdinate)
    BPT_type2_colored(LOGAscisse, LOGOrdinate)






    #BPT Diagram del Terzo tipo
    """
    log([OIII]/Hb) = 0.73 / (log([OI]/Ha) + 0.59) + 1.33    (main AGN line)
    log([OIII]/Hb) = 1.18 log([OI]/Ha) + 1.30  (LINER/Sy2 line)
    """









    
    
     
    






if __name__ == "__main__":

 main ()