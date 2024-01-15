from astropy.io import fits
import numpy as np
import pandas as pd
import csv


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



def leggi_file_txt(nome_file, numero_colonne):
    # Inizializza un array vuoto
    dati = []

    # Apri il file in modalità lettura
    with open(nome_file, 'r') as file:
        # Itera sulle righe del file
        for riga in file:
            # Suddividi la riga in colonne utilizzando lo spazio come delimitatore
            colonne = riga.split()

            # Verifica se il numero di colonne è coerente con quello specificato
            if len(colonne) == numero_colonne:
                # Converti le colonne in numeri float e aggiungi alla lista dati
                dati.append([float(colonna) for colonna in colonne])
            else:
                print(f"Attenzione: La riga '{riga}' non ha il numero corretto di colonne.")

    # Restituisci l'array con i dati
    return dati



def salva_array_senza_parentesi(file_path, array):
    # Salva l'array nel file di testo senza parentesi quadre
    #np.savetxt(file_path, array, delimiter=' ', fmt='%g')
    #fmt_list = ['%s', '%g']
    np.savetxt(file_path, array, delimiter='\t',header='\t'.join(list(array.dtype.names)), comments='', fmt='%g')


def leggi_nomi_colonne(file_fits):
    # Apri il file FITS
    hdul = fits.open(file_fits)

    # Se il file FITS contiene una tabella
    if isinstance(hdul[1], fits.BinTableHDU):
        # Estrai i nomi delle colonne
        nomi_colonne = hdul[1].columns.names
    else:
        # Se il file FITS contiene solo un'immagine, usa i nomi delle keywords dell'header
        nomi_colonne = hdul[0].header.keys()

    # Chiudi il file FITS
    hdul.close()

    return nomi_colonne

def leggi_file_completo_txt(nome_file):
    # Inizializza un array vuoto
    dati = []

    # Apri il file in modalità lettura
    with open(nome_file, 'r') as file:
        # Itera sulle righe del file
        for riga in file:
            # Suddividi la riga in colonne utilizzando lo spazio come delimitatore
            colonne = riga.split()

            # Converti le colonne in numeri float e aggiungi alla lista dati
            dati.append([float(colonna) for colonna in colonne])

    # Restituisci l'array con i dati
    return dati

def leggi_file_csv(file_path):
    """
    Legge File in formato csv, in cui la prima riga è costituita dai nomi delle colonne
    Viene Utilizzato come separatore la classica Virgola !!!
    """
    with open(file_path, 'r') as file:
        # Leggi la prima riga del file, NB il carattere separatore impostato per la prima riga è lo spazio \t
        colonne = file.readline().strip().split(',')
        
        # Utilizza il modulo csv per leggere il resto del file
        reader = csv.DictReader(file, fieldnames=colonne)
        
        # Inizializza un array per salvare i dati
        dati = []
        
        # Leggi le righe rimanenti e aggiungi i dati all'array
        for riga in reader:
            dati.append(riga)
    
    return colonne, dati







def leggi_file_tsv(file_path):
    """
    Legge un file in formato tsv, in cui la prima riga è costituita dai nomi delle colonne.
    Utilizza il carattere di tabulazione "\t" come separatore.
    Restituisce una lista eventualmente da trasformare 
    in un array numpy successivamente con i nomi delle colonne ed un np.array con i dati presenti nel file.
    """
    # Usa pandas per leggere il file in un DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Estrai nomi delle colonne e dati dal DataFrame
    colonne = df.columns.tolist()
    dati = df.to_numpy()

    return colonne, dati