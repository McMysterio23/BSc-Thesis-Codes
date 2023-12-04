import sys
import numpy as np
from astropy.io import fits
from astropy.table import QTable
import matplotlib.pyplot as plt

from libConfronti import leggi_file_fits, confronta_files_fits, confronta_colonne, stampa_nomi_colonne, salva_su_file


def main():

    #creo la path della cartella con dentro tutti i dati
    dataPATH = '/Users/andreamaccarinelli/Desktop/SDSS/'

    file1 = dataPATH + sys.argv[1]
    file2 = dataPATH + sys.argv[2]
    soglia = float(sys.argv[3])

    colonna1_file1 = 'RA'
    colonna2_file1 = 'DEC'
    colonna1_file2 = 'RAdeg'
    colonna2_file2 = 'DEdeg'

    

    elementi_comuni = confronta_colonne(file1, colonna1_file1, colonna2_file1, file2, colonna1_file2, colonna2_file2, soglia)

    salva_su_file(elementi_comuni, 'output_colonne_unite.txt')

    if elementi_comuni:
        # Estrai RA e DEC dalle tuple
        ra_comuni, dec_comuni = zip(*elementi_comuni)

        # Plot dei dati
        plt.scatter(ra_comuni, dec_comuni, marker='o', label='Elementi Comuni')
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.title('Plot degli Elementi Comuni')
        plt.legend()
        plt.show()

    #stampa il messaggio di esito positivo dell'operazione
    print(f"Gli elementi comuni nelle colonne '{colonna1_file1}' e '{colonna2_file2}' entro la soglia sono stati salvati in 'output_colonne_unite.txt'")







if __name__ == "__main__":
 if (len(sys.argv)<3):
        print ('usage: ', sys.argv[0], 'catalogoUNO.txt', sys.argv[1], 'catalogoDUE', sys.argv[2], 'soglia di accettazione', sys.argv[3])
        sys.exit()

 main ()



