"""
Lo Scopo di questo programma è quello di andare ad aggiungere informazioni al file IntersezioneInfo.txt 
con tutte le informazioni relative all'intersezione provenienti dal file LINE, sono arrivato all'indice 16 passato da linea di comando
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table


from libreriaUNO import salva_array_senza_parentesi, leggi_nomi_colonne
from libreriaUNO import leggi_file_txt, leggi_file_tsv
from libreriaDUE import plot_and_save_scatter
from libreriaDUE import plot_scatter
from libreriaDUE import Intersection_Array_Type1, Aggiungi_Colonna, Aggiungi_Colonna6, AddColumn, Lettura_Colonne_RawDATA, Intersection_Array_Type2
from libreriaDUE import Aggiungi_Colonna6plus, CrossingColonna

def main():
    #Indice del file di dati grezzi da leggere 
    Indice_RAW = 1
    nomi_colonne_SDSS = Lettura_Colonne_RawDATA(Indice_RAW)
    
    lunghezza_lista_SDSS = len(nomi_colonne_SDSS)

    IndexColumnRequested = int(sys.argv[1])


    #STAMPE DI INFORMAZIONI UTILI !!!
    print('Ci sono un totale di ',lunghezza_lista_SDSS,' colonne in questo file !')
    print('Al momento stai cercando di analizzare la colonna con indice', IndexColumnRequested, 'ossia ', nomi_colonne_SDSS[IndexColumnRequested])


    
    Deposit= CrossingColonna(1, nomi_colonne_SDSS[IndexColumnRequested])
    ColonnaTrovata = Deposit[0]
    ARRDeposit = np.array(ColonnaTrovata)
    
    
    #print(np.shape(ARRDeposit))
    #TableDeposit = Table(Column(ARRDeposit, nomi_colonne_SDSS[IndexColumnRequested]))

    #print('L* Array creato ha dimensioni(Righe, Colonne) pari a :',Deposit.shape)
    
    
    

    #VADO DUNQUE A LEGGERE IL FILE CONTENENTE LE INFORMAZIONI BASALI : IntersecInfo.txt
    FILE_PATH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfo.txt'

    if(IndexColumnRequested > 2):
        FILE_PATH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfoLine.txt'
    nomi_colonne, RAWTabella = leggi_file_tsv(FILE_PATH)
    TableRAWTabella = Table(RAWTabella, names = nomi_colonne)

    TableRAWTabella = AddColumn(TableRAWTabella, ARRDeposit, nomi_colonne_SDSS[IndexColumnRequested])


        





    #Stampa dunque il tutto alla fine !
    salva_array_senza_parentesi('/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfoLine.txt', TableRAWTabella)














if __name__ == "__main__":

 main ()