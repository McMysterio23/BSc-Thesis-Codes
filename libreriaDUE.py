from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table


from libreriaUNO import leggi_file_txt


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
Intersezione = PATH+"Matricione.txt"
indices = PATH+"/Int_2_Asec.txt"

#conversione degli array di liste in array Numpy
RAW_intersection = leggi_file_txt(Intersezione, 2 )
coordinate = np.array(RAW_intersection)

indiciint = leggi_file_txt(indices, 1)
indiciint2 = np.array(indiciint)

dimensione = len(coordinate[:,0])

def plot_and_save_scatter(array1, array2, save_path, title="Scatter Plot", x_label="X", y_label="Y", nomefile = "Scatter_plot.png"):
    """
    Questa funzione realizza un plot scatter sulla base di due array passati da linea di comando, usala solo se non ti interessa
    visualizzare il plot sul momento ma solo salvarlo per consultarlo in un secondo momento !!!
    """

    # Crea uno scatter plot
    plt.scatter(array1, array2)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Salva l'immagine PNG nella cartella specificata dopo aver chiuso il plot
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    save_file_path = os.path.join(save_path, nomefile)
    plt.savefig(save_file_path)

    print(f"Scatter plot salvato in: {save_file_path}")



def plot_scatter(array1, array2, title="Scatter Plot", x_label="X", y_label="Y"):
    # Crea uno scatter plot
    plt.scatter(array1, array2)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Mostra il plot a video
    plt.show()


def Intersection_Array_Type1(Colonna="Nomecolonna"):
    """
    Crea un array basato sull'intersezione, ricevendo solo il nome della colonna che si vuole rendere array
    """
    #Apro i file fits che mi interessano 
    obj=fits.open(INFO)
    data=obj[1].data
    arr=data[Colonna]

    #creo l'array
    intersection = np.ones(dimensione)

    #Popolo l'array
    for i in range(dimensione):
       punto = int (indiciint2[i,0])
       intersection[i] = arr[punto]
    
    return intersection

def Aggiungi_Colonna(arr1, Arr2, arr3, arr4, arr5, arr6, arr7=np.ones(dimensione), arr8=np.zeros(dimensione), arr9=np.zeros(dimensione),
                     nome1='Nome Colonna1', nome2='Nome Colonna2', nome3='Nome Colonna3',
                     nome4='Nome Colonna4', nome5='Nome Colonna5', nome6='Nome Colonna6',
                     nome7='Nome Colonna7', nome8='Nome Colonna8', nome9='Nome Colonna9'):

    # Assicurati che gli array siano di tipo sequenza (es. liste o array NumPy)
    for arr in [arr1, Arr2, arr3, arr4, arr5, arr6]:
        if not isinstance(arr, (list, np.ndarray)):
            raise ValueError("Gli array devono essere di tipo sequenza (es. liste o array NumPy)")    


 #Sample_derived = Table([Column(indiciint2, name='Indici_INFO'), Column(coordinate[:,0], name='PlateID'), Column(coordinate[:,1], name='FiberID'),
                            #Column(arr1, name=nome1), Column(Arr2, name= nome2), Column(arr3, name = nome3), Column(arr4, name = nome4),
                            #Column (arr5, name = nome5), Column(arr6, name = nome6), Column(arr7, name= nome7), Column(arr8, name = nome8),
                            #Column(arr9, name = nome9)])
    #return Sample_derived

    
