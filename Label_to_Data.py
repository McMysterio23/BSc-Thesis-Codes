import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column

from libreriaUNO import salva_array_senza_parentesi
from libreriaUNO import leggi_file_txt




def main():
    #Inizializzo i percorsi ai file che utilizzerò
    PATH = "/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/"
    PATHDATA="/Users/andreamaccarinelli/Desktop/SDSS/"
    All=PATHDATA+"/gal_info_dr7_v5_2.fits"
    Indices = PATH+"Int_2_Asec.txt"

    #Apro i file e ne memorizzo il contenuto in opportuni array

    obj=fits.open(All)
    data=obj[1].data
    plateID, fiberID = data["plateid"], data["fiberid"]

    Indici = leggi_file_txt(Indices, 1)
    
    #Popolo l'array di dimensioni len(indici) x 2 con i dati di plateID e fiberID relativi agli indici che mi interessano
    dimensioni = len(Indici)
    Matricione = np.ones((dimensioni,2))
    for i in range(dimensioni) :
       for j in range(2):
          if(j==0):
             Matricione[i,j] = plateID[i]
          elif(j==1):
             Matricione[i,j] = fiberID[i]
    
    
    #vado infine a creare un file .txt che contenga i dati che ho ottenuto 
   
    salva_array_senza_parentesi('Matricione.txt', Matricione)
   
             

if __name__ == "__main__":

 main ()