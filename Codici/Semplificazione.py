import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table


from libreriaUNO import leggi_nomi_colonne, leggi_file_tsv
from libreriaDUE import plot_and_save_scatter, plot_scatter
from libreriaDUE import Aggiungi_Colonna, Aggiungi_Colonna6, AddColumn, Lettura_Colonne_RawDATA


def main():

    #apro il file con le info sugli SFR
    FILE_PATH_SFRs = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/Star_Formating_Ratios_Intersezione.txt'
    nomi_colonne_SFRs, RAWTabella_SFRs = leggi_file_tsv(FILE_PATH_SFRs)
    TableRAWTabella_SFRs = Table(RAWTabella_SFRs, names = nomi_colonne_SFRs)

    # Specifica il percorso del file .fits
    percorso_file_fits = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/SFR_BCG.fits'

    # Salva l'oggetto Table nel file .fits
    TableRAWTabella_SFRs.write(percorso_file_fits, format='fits', overwrite=True)

    print(f"Il file .fits è stato salvato con successo in {percorso_file_fits}")


    #apro il file con le info sugli sSFR
    FILE_PATH_sSFRs = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/Specific_Star_Formating_Ratios_Intersezione.txt'
    nomi_colonne_sSFRs, RAWTabella_sSFRs = leggi_file_tsv(FILE_PATH_sSFRs)
    TableRAWTabella_sSFRs = Table(RAWTabella_sSFRs, names = nomi_colonne_sSFRs)

    # Specifica il percorso del file .fits
    percorso_file_fits = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/sSFR_BCG.fits'

    # Salva l'oggetto Table nel file .fits
    TableRAWTabella_sSFRs.write(percorso_file_fits, format='fits', overwrite=True)

    print(f"Il file .fits è stato salvato con successo in {percorso_file_fits}")


    #apro il file con le info del file totLGM.txt
    FILE_PATH_TOTLGM = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/TOTLGM_Intersezione.txt'
    nomi_colonne_TOTLGM, RAWTABELLA_TOTLGM = leggi_file_tsv(FILE_PATH_TOTLGM)
    TableRAWTabella_TOTLGM = Table(RAWTABELLA_TOTLGM, names = nomi_colonne_TOTLGM)

    # Specifica il percorso del file .fits
    percorso_file_fits = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/TOTLGM_BCG.fits'

    # Salva l'oggetto Table nel file .fits
    TableRAWTabella_TOTLGM.write(percorso_file_fits, format='fits', overwrite=True)

    print(f"Il file .fits è stato salvato con successo in {percorso_file_fits}")

    #apro il file infoLine.txt
    FILE_INFO = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfoLine.txt'
    nomi_colonne_INFO, RAWTabella_INFO = leggi_file_tsv(FILE_INFO)
    TableRAWTabella_INFO = Table(RAWTabella_INFO, names = nomi_colonne_INFO)

    # Specifica il percorso del file .fits
    percorso_file_fits = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/INFOLINE_BCG.fits'

    # Salva l'oggetto Table nel file .fits
    TableRAWTabella_INFO.write(percorso_file_fits, format='fits', overwrite=True)

    print(f"Il file .fits è stato salvato con successo in {percorso_file_fits}")

    #apro il file FIBOH.txt
    FILE_PATH_FIBOH = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/FIBOH_Intersezione.txt'
    nomi_colonne_FIBOH, RAWTabella_FIBOH = leggi_file_tsv(FILE_PATH_FIBOH)
    TableRAWTabella_FIBOH = Table(RAWTabella_FIBOH, names = nomi_colonne_FIBOH)

    # Specifica il percorso del file .fits
    percorso_file_fits = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/FIBOH_BCG.fits'

    # Salva l'oggetto Table nel file .fits
    TableRAWTabella_FIBOH.write(percorso_file_fits, format='fits', overwrite=True)

    print(f"Il file .fits è stato salvato con successo in {percorso_file_fits}")

    #apro il file FIBLGM.txt
    FILE_PATH_FIBLGM = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/FIBLGM_Intersezione.txt'
    nomi_colonne_FIBLGM, RAWTabella_FIBLGM = leggi_file_tsv(FILE_PATH_FIBLGM)
    TableRAWTabella_FIBLGM = Table(RAWTabella_FIBLGM, names = nomi_colonne_FIBLGM)

    # Specifica il percorso del file .fits
    percorso_file_fits = '/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/DATA/FIBLGM_BCG.fits'

    # Salva l'oggetto Table nel file .fits
    TableRAWTabella_FIBLGM.write(percorso_file_fits, format='fits', overwrite=True)

    print(f"Il file .fits è stato salvato con successo in {percorso_file_fits}")



if __name__ == "__main__":

 main ()