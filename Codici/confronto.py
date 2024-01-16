"""
Questo codice va a trovare gli indici dei file con coordinate entro i 2 arcosecondi dall'archivio delle BCG
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column

from libreriaUNO import salva_array_senza_parentesi






def main():

    #Inizializzo i percorsi ai file che utilizzerò
    PATH="/Users/andreamaccarinelli/Desktop/SDSS/"
    All=PATH+"/gal_info_dr7_v5_2.fits"
    Bcg = PATH+"sampleC4.fits"


    #Raccolgo RA e DEC dal Catalogo SDSS/DR7
    obj=fits.open(All)
    data=obj[1].data
    ra,dec=data["RA"],data["DEC"]

    # Raccolgo RA E DEC dal Catalogo C4
    obj2=fits.open(Bcg)
    data2 = obj2[1].data
    rab,decb=data2["RAdeg"],data2["DEdeg"]

    #Inizializzo le variabili di zoom
    k1_initial = (1.5/3600)
    incr = (0.5/3600)

    
    #inizializzo una lista vuota per contenere gli indici in comune 
    indici_ottenuti = []

    #cerco dunque di implementare il tutto per l'intero catalogo C4 memorizzando gli indici in un array per il momento
    for i in range(len(rab)):
        k1 = k1_initial

        ival = np.where((ra >= rab[i]-k1) & (ra <= rab[i]+k1) & (dec >= decb[i]-k1) & (dec <= decb[i]+k1))[0]
        if(len(ival)==0 ):
         k1=k1+incr
         ival = np.where((ra >= rab[i]-k1) & (ra <= rab[i]+k1) & (dec >= decb[i]-k1) & (dec <= decb[i]+k1))[0]
         if (len(ival)==1):
            indici_ottenuti.append(ival)

         continue
        elif (len(ival)>1):
            k1=k1-incr
            ival = np.where((ra >= rab[i]-k1) & (ra <= rab[i]+k1) & (dec >= decb[i]-k1) & (dec <= decb[i]+k1))[0]
            if (len(ival)==1):
               indici_ottenuti.append(ival)

            continue
        elif (len(ival)==1):
           indici_ottenuti.append(ival)
        
    indici_ottenuti = np.array(indici_ottenuti)
    
    #Conclusione delle operazioni
    print('La dimensione della intersezione, entro 2 arcsec è di :', len(indici_ottenuti), 'elementi')
    salva_array_senza_parentesi('ElementiComuniC4.txt', indici_ottenuti )

    lengthC4 = len(rab)
    rapp = (len(indici_ottenuti))/(lengthC4)
    print('Considerando che il catalogo C4 comprendeva ',lengthC4,' elementi, esattamente il ', rapp*100, '% degli oggetti era presente in SDSS/DR7')



      
if __name__ == "__main__":

 main ()