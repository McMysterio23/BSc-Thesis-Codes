import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column

from libreriaUNO import save_array_to_txt





def main():

    PATH="/Users/andreamaccarinelli/Desktop/SDSS/"
    All=PATH+"/gal_info_dr7_v5_2.fits"
    #Bcg=PATH+"/DATA_BCG.txt" #versione con C4 in formato txt
    Bcg = PATH+"sampleC4.fits"

    obj=fits.open(All)
    data=obj[1].data
    # RA E DEC ALL SAMPLE
    ra,dec=data["RA"],data["DEC"]

    # RA E DEC IN BCG
    #rab,decb= np.loadtxt(Bcg, usecols=[1,2], unpack=True,dtype=float)
    obj2=fits.open(Bcg)
    data2 = obj2[1].data
    rab,decb=data2["RAdeg"],data2["DEdeg"]

    #ESEMPIO UN CASO i=0
    k1,k2=0.01,0.001
    incr = 0.0015
    i=0
    #TROVA INDICE (ival) CORRISPONDENTE A BCG IN FILE TOTALE SAMPLE
    ival = np.where((ra >= rab[i]-k2) & (ra <= rab[i]+k2) & (dec >= decb[i]-k2) & (dec <= decb[i]+k2))[0][0]
    #ival è un array bidimensionale con nella prima colonna l'elenco delle posizioni in cui la condizione dentro le parentesi tonde
    #della funzione np.where si verificano e nella seconda dove non si verificano.
    #così facendo la parte che compare sopra con le due parentesi quadre seleziona il primo elemento della tabella che si crea, 
    #nella fattispecie il primo elemento in cui la condizione si avvera 

    #CONFRONTA VALORI
    print(ra[ival],dec[ival]," == ",rab[i],decb[i])
    print("L'indice cui corrisponde il primo oggetto di C4, entro un intorno di 3.6 arcosecondi è ",ival)

    #inizializzo un vettore di zeri di dimensione pari a quella del catalogo C4 in cui cambiare i valori poi successivamente 
    indici_ottenuti = np.zeros(len(rab))

    #cerco dunque di implementare il tutto per l'intero catalogo C4 memorizzando gli indici in un array per il momento
    for i in range(len(rab)):
        ival = np.where((ra >= rab[i]-k1) & (ra <= rab[i]+k1) & (dec >= decb[i]-k1) & (dec <= decb[i]+k1))
        if(len(ival)==0 ):
         k2=k2+incr
        elif (len(ival)>1):
            k2=k2-incr
        #indici_ottenuti = ival



      #RESTA DA IMPLEMENTARE IL SISTEMA DI SALVATAGGIO DEGLI INDICI IN UN FILE .TXT
    #save_array_to_txt(indici_ottenuti,'ElementiComuniC4.txt' )

    #print(len(ival))
    #indici_ottenuti = ival[:,0]
    #save_array_to_txt(indici_ottenuti,'ElementiComuniC4.txt' )


      
if __name__ == "__main__":

 main ()