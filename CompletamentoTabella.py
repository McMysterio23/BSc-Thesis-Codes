"""
Lo scopo di questo programma è quello di andare a leggere il file .txt gia creato ed andare ad aggiungervi le informazioni relative
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table


from libreriaUNO import salva_array_senza_parentesi, leggi_nomi_colonne
from libreriaUNO import leggi_file_txt, leggi_file_completo_txt, leggi_file_csv, leggi_file_tsv
from libreriaDUE import plot_and_save_scatter
from libreriaDUE import plot_scatter
from libreriaDUE import Intersection_Array_Type1, Aggiungi_Colonna, Aggiungi_Colonna6, AddColumn, Lettura_Colonne_RawDATA, Intersection_Array_Type2
from libreriaDUE import Aggiungi_Colonna6plus, Crossing


def main():
    #RAWTabella = leggi_file_completo_txt('/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfo.txt')
    #Tabella = np.array(RAWTabella)
    FILE_PATH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfo.txt'
    nomi_colonne, RAWTabella = leggi_file_tsv(FILE_PATH)
    #Conviene Tenere una struttura di tipo Array di liste essendo che le colonne vengono identificate con le stringhe del relativo titolo
    #Questo è un esempio di come si risale al campo PLATEID del primo oggetto 
    #print(RAWTabella[0][nomi_colonne[1]])


    #creo un oggetto Table con quanto gia creato :
    mia_tabella_astropy = Table(RAWTabella, names=nomi_colonne)
    #NB quando vuoi accedere ad un elemento preciso dell'oggetto Table indica nella prima breket il nome della colonna e nel secondo l'indice di riga
    #print(mia_tabella_astropy[nomi_colonne[1]][0])
    

    #Indice del file di dati grezzi da leggere 
    Indice_RAW = 3
    nomi_colonne_SDSS = Lettura_Colonne_RawDATA(Indice_RAW)
    
    lunghezza_lista_SDSS = len(nomi_colonne_SDSS)
    #print(len(Intersection_Array_Type2(3, nomi_colonne_SDSS[1])))
    #print(nomi_colonne_SDSS)

    
    colonne = [Column(Crossing(Indice_RAW, nome_colonna), name=nome_colonna) for nome_colonna in nomi_colonne_SDSS]
    TabellaCreata = Table(colonne)
    
    
    # Creazione di colonne usando astropy.table.Column


    # Creazione della tabella usando astropy.table.Table
    TabellaCreata = Table(colonne)

    colonna1 = Crossing(3, nomi_colonne_SDSS[0])
    colonna2 = Crossing(3, nomi_colonne_SDSS[1])
    #print(colonna1[0], colonna2[0])

    # Stampa le colonne
    #print("Nomi delle colonne:", TabellaCreata.colnames)

    # Stampa le prime 5 righe della tabella
    #for i in range(5):
    #    print(f"Riga {i}: {TabellaCreata[i]}")

    
        

    salva_array_senza_parentesi('/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/Star_Formating_Ratios_Intersezione.txt', TabellaCreata)

    





if __name__ == "__main__":

 main ()