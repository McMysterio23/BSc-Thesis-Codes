from astropy.io import fits
import numpy as np


def leggi_file_fits(file_path):
    with fits.open(file_path) as hdul:
        # Assumiamo che l'informazione di interesse si trovi nella prima tabella del file FITS
        tabella = hdul[1].data
    return tabella


def confronta_files_fits(file1, file2):
    # Leggi le tabelle dai file FITS
    tabella1 = leggi_file_fits(file1)
    tabella2 = leggi_file_fits(file2)

    # Estrai gli oggetti unici presenti in entrambi i file
    oggetti_comuni = set(tabella1['NOME_COLONNA_INTERESSATA']).intersection(set(tabella2['NOME_COLONNA_INTERESSATA']))

    return oggetti_comuni


def confronta_colonne(file1, colonna1_file1, colonna2_file1, file2, colonna1_file2, colonna2_file2, soglia):
    # Leggi i dati dalle colonne nei due file FITS
    dati_colonna1_file1 = leggi_colonna(file1, colonna1_file1)
    dati_colonna2_file1 = leggi_colonna(file1, colonna2_file1)
    dati_colonna1_file2 = leggi_colonna(file2, colonna1_file2)
    dati_colonna2_file2 = leggi_colonna(file2, colonna2_file2)

    # Trova gli elementi comuni tra le coppie di valori entro la soglia
    elementi_comuni = []
    for valore1_file1, valore2_file1, valore1_file2, valore2_file2 in zip(dati_colonna1_file1, dati_colonna2_file1, dati_colonna1_file2, dati_colonna2_file2):
        if abs(valore1_file1 - valore1_file2) < soglia and abs(valore2_file1 - valore2_file2) < soglia:
            elementi_comuni.append((valore1_file1, valore2_file1))

    return elementi_comuni

def leggi_colonna(file, colonna):
    with fits.open(file) as hdul:
        # Assumiamo che l'informazione di interesse si trovi nella prima tabella del file FITS
        tabella = hdul[1].data
        dati_colonna = tabella[colonna]

    return dati_colonna

def stampa_nomi_colonne(file):
    with fits.open(file) as hdul:
        # Assumiamo che l'informazione di interesse si trovi nella prima tabella del file FITS
        tabella = hdul[1].data
        nomi_colonne = tabella.names

    print(f"Nomi delle colonne nel file {file}: {nomi_colonne}")

def salva_su_file(elementi_comuni, nome_file_output):
    np.savetxt(nome_file_output, elementi_comuni, fmt='%.18e', delimiter='\t')

