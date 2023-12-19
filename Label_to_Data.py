import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column

from libreriaUNO import save_array_to_txt
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
    





     
    






if __name__ == "__main__":

 main ()