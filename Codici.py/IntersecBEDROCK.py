"""
Lo Scopo di questo programma è quello di andare a creare un file .txt con tutte le informazioni relative all'intersezione
provenienti dal fiile Basale con il nome INFO contenente le informazioni basilari
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table

from libreriaUNO import salva_array_senza_parentesi, leggi_nomi_colonne
from libreriaUNO import leggi_file_txt
from libreriaDUE import plot_and_save_scatter
from libreriaDUE import plot_scatter
from libreriaDUE import Intersection_Array_Type1, Aggiungi_Colonna, Aggiungi_Colonna6, AddColumn, Lettura_Colonne_RawDATA, Intersection_Array_Type2
from libreriaDUE import Aggiungi_Colonna6plus, Crossing

def main():
    #Indice del file di dati grezzi da leggere 
    Indice_RAW = 0
    nomi_colonne_SDSS = Lettura_Colonne_RawDATA(Indice_RAW)
    
    lunghezza_lista_SDSS = len(nomi_colonne_SDSS)


    #STAMPE DI INFORMAZIONI UTILI !!!
    #print('Ci sono un totale di ',lunghezza_lista_SDSS,' colonne in questo file !')
    #print('Colonne disponibili nel file selezionato : \n', nomi_colonne_SDSS)

    
    Deposit, Columns_selected = Crossing(0, True)

    print('L* Array creato ha dimensioni(Righe, Colonne) pari a :',Deposit.shape)
    print('Al momento sto provvedendo un numero di ', Columns_selected, 'colonne all oggetto pandas')
    
    df = pd.DataFrame(Deposit, columns=Columns_selected)
    #print(df.columns)
    TableDeposit = Table.from_pandas(df)
    print(TableDeposit[1])


    #Stampa le colonne
    salva_array_senza_parentesi('/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneINFO.txt', TableDeposit)














if __name__ == "__main__":

 main ()