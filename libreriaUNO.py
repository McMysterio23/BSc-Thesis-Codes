from astropy.io import fits
import numpy as np


def leggi_file_fits(file_path):
    #restituisce un'array con il contenuto del file fits che gli viene fornito
    with fits.open(file_path) as hdul:
        # Assumiamo che l'informazione di interesse si trovi nella prima tabella del file FITS
        tabella = hdul[1].data
    return tabella


def stampa_nomi_colonne(file):
    with fits.open(file) as hdul:
        # Assumiamo che l'informazione di interesse si trovi nella prima tabella del file FITS
        tabella = hdul[1].data
        nomi_colonne = tabella.names

    print(f"Nomi delle colonne nel file {file}: {nomi_colonne}")

def salva_su_file(elementi_comuni, nome_file_output):
    np.savetxt(nome_file_output, elementi_comuni, fmt='%.18e', delimiter='\t')

def save_array_to_txt(array, filename):
    """
    Salva il contenuto di un array unidimensionale in un file di testo (.txt).

    Parameters:
    - array: array unidimensionale da salvare.
    - filename: nome del file di testo (.txt) in cui salvare l'array.
    """
    with open(filename, 'w') as file:
        for element in array:
            file.write(str(element) + '\n')
