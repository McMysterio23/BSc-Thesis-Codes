import numpy as np
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
from libreriaDUE import Aggiungi_Colonna6plus

def main():
    #Questo CODICE NON FUNZIONA !!!!!
    
   """
   Redshift = Intersection_Array_Type1('z')
    Errore_Redshift = Intersection_Array_Type1('z_err')
    Sigma = Intersection_Array_Type1('v_disp')
    Errore_Sigma = Intersection_Array_Type1('v_disp_err')
    SN_Medio = Intersection_Array_Type1('sn_median')
    Reddening = Intersection_Array_Type1('e_bv_sfd')




    
    #vado a creare due array distinti potenzialmente : uno privo dei dati elaborati ed un altro con essi inclusi
    #PROBLEMA DA RISOLVERE : CONTROLLARE CHE NON ESISTANO DEI NOMI DI COLONNE IN FILE DIVERSI CHE SONO DEGENERI !! 
    #CONTROLLA IN PARTICOLARE I FILE FIBOH, TOTLGM E COMPAGNIA BELLA 
    Structured = Aggiungi_Colonna6(Redshift, Errore_Redshift, Sigma, Errore_Sigma, SN_Medio, 
                                  Reddening, 'z', 'z_err','v_disp', 'v_disp_err',
                                  'sn_median', 'e_bv_sfd') 
    
   """ 
    #print(Structured[0])
    
    

    
    

    
    #salva_array_senza_parentesi('/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/IntersezioneInfo.txt', Structured)


    

    
    





    











if __name__ == "__main__":

 main ()