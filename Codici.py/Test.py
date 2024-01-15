"""
Questo Codice deve essere eseguito dando da linea di comando la posizione nell'array degli indici dell'indice che si 
vuole stampare
"""

from libreriaUNO import leggi_file_txt, leggi_nomi_colonne
import numpy as np
import sys
from libreriaDUE import Aggiungi_Colonna6plus, Crossing, Lettura_Colonne_RawDATA


#Inizializzo i percorsi ai file che utilizzerò
PATH = "/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/"
PATHDATA="/Users/andreamaccarinelli/Desktop/SDSS/"
PATH_TESI_IMMAGINI ='/Users/andreamaccarinelli/Desktop/BSc-Thesis/Immagini'

INFO=PATHDATA+"/gal_info_dr7_v5_2.fits"
LINE = PATHDATA+'/gal_line_dr7_v5_2.fits'
INDX = PATHDATA+'/gal_indx_dr7_v5_2.fits'
SFRS = PATHDATA+'/gal_totsfr_dr7_v5_2.fits'
SPECSFR = PATHDATA+'/gal_totspecsfr_dr7_v5_2.fits'
TOTLGM = PATHDATA+'/totlgm_dr7_v5_2.fits'
SPECMAG = PATHDATA+'/gal_specmag_dr7_v5_2.fits'
FIBLGM = PATHDATA+'/fiblgm_dr7_v5_2.fits'
FIBOH = PATHDATA+'/gal_fiboh_dr7_v5_2.fits'
Lista_Di_Percorsi_Ai_Dati = [INFO, LINE, INDX, SFRS, SPECSFR, TOTLGM, SPECMAG, FIBLGM, FIBOH]
Intersezione = PATH+"Matricione.txt"
indices = PATH+"/Int_2_Asec.txt"



indiciint = leggi_file_txt(indices, 1)
indiciint2 = np.array(indiciint)
POV = int(sys.argv[1])
print(indiciint[POV])
#print(indiciint2[0])
#print(indiciint2[341])

"""
nomi_colonne_SDSS = Lettura_Colonne_RawDATA(3)
print(nomi_colonne_SDSS[0], nomi_colonne_SDSS[1])
colonna1 = Crossing(3, nomi_colonne_SDSS[0])
print(colonna1)
"""





"""

colonna2 = Crossing(3, nomi_colonne_SDSS[1])
if (colonna1[0] != colonna2[0]):
    print('Valori diversi, problema risolto')
    print(colonna1[0][0], colonna2[0][0])
else:
    print('Valori uguali, Problema NON risolto !')
    print(colonna1[0][0], colonna2[0][0])

"""

