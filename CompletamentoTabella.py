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
from libreriaDUE import Aggiungi_Colonna6plus


def main():
    #RAWTabella = leggi_file_completo_txt('/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfo.txt')
    #Tabella = np.array(RAWTabella)
    FILE_PATH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfo.txt'
    nomi_colonne, RAWTabella = leggi_file_tsv(FILE_PATH)
    #Conviene Tenere una struttura di tipo Array di liste essendo che le colonne vengono identificate con le stringhe del relativo titolo
    #Questo è un esempio di come si risale al campo PLATEID del primo oggetto 
    print(RAWTabella[0][nomi_colonne[1]])

    #creo un oggetto Table con quanto gia creato :
    mia_tabella_astropy = Table(RAWTabella, names=nomi_colonne)
    #NB quando vuoi accedere ad un elemento preciso dell'oggetto Table indica nella prima breket il nome della colonna e nel secondo l'indice di riga
    print(mia_tabella_astropy[nomi_colonne[1]][0])


    





if __name__ == "__main__":

 main ()