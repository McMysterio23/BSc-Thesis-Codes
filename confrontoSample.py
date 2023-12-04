import sys
import numpy as np
from astropy.io import fits
from astropy.table import QTable
import matplotlib.pyplot as plt

from libConfronti import leggi_file_fits, confronta_files_fits, confronta_colonne, stampa_nomi_colonne


def main():

    #creo la path della cartella con dentro tutti i dati
    dataPATH = '/Users/andreamaccarinelli/Desktop/SDSS/'

    file1 = dataPATH + sys.argv[1]
    file2 = dataPATH + sys.argv[2]

    #memorizzo le due tabelle presenti nei prizm fits in due oggetti plt
    #data1 = file1[1].data
    #data2 = file2[1].data


    #colonna1_file1 = 'RA'
    #colonna2_file1 = 'DEC'
    #colonna1_file2 = 'RAdeg'
    #colonna2_file2 = 'DEdeg'

    #RA_1 = data1[colonna1_file1]






    #stampa_nomi_colonne(file1)
    #stampa_nomi_colonne(file2)

    #utilizzo dunque le funzioni delle cartella libConfronti per confrontare gli oggetti nei due sample

    #oggetti_comuni = confronta_files_fits(file1, file2)

    # Sostituisci 'NOME_COLONNA_1_FILE1', 'NOME_COLONNA_2_FILE1', 'NOME_COLONNA_1_FILE2', e 'NOME_COLONNA_2_FILE2' con i nomi reali delle colonne che vuoi confrontare nei tuoi file FITS
    colonna1_file1 = 'RA'
    colonna2_file1 = 'DEC'
    colonna1_file2 = 'RAdeg'
    colonna2_file2 = 'DEdeg'

    soglia = 4

    elementi_comuni_colonna1, elementi_comuni_colonna2 = confronta_colonne(file1, colonna1_file1, colonna2_file1, file2, colonna1_file2, colonna2_file2, soglia)

    print(f"Gli elementi comuni nella colonna '{colonna1_file1}' tra i due file sono:", elementi_comuni_colonna1)
    print(f"Gli elementi comuni nella colonna '{colonna2_file2}' tra i due file sono:", elementi_comuni_colonna2)








if __name__ == "__main__":
 if (len(sys.argv)<2):
        print ('usage: ', sys.argv[0], 'catalogoUNO.txt', sys.argv[1], 'catalogoDUE', sys.argv[2])
        sys.exit()

 main ()



