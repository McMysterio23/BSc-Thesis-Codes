from astropy.io import fits
from os.path import exists
import numpy as np

#CARICA fits file
file="/Users/andreatravascio/Desktop/Tesisti/TesiAndrea/gal_line_dr7_v5_2.fits"
hdul = fits.open(file)  
#Matrice tutti dati
data = hdul[1].data
print(data)  #Dai un'occhiata alle keywords
#Per esempio, se vogliamo associate a sigma_bal la colonna chiamara "SIGMA_BALMER"
sigma_bal = data["SIGMA_BALMER"]

#Per salvare i file ti consiglio questa procedura:
FileToSave= "Path/File.txt"
fmt = "%.5e" #Se hai tutti numeri in notazione scientifica
fmt = "%f" #se hai tutti numeri floats #Vedi google per altri formati
#Se hai colonne di diversi type fai per esempio così:
fmt=["%i","%f","%.5e"] #
#a, b, c e d non devono essere list ma array (colonne)
#Se per esempio a è una lista, che può essere comodo quando non sai a priori la grandezza dell'array e vuoi implementralo con append()
#allora puoi farlo diventare un array nel seguente modo: a= np.asarray(a)
data = np.column_stack((a,b,c,d))
if os.path.exists(FileToSave) == False:#Se il file già esiste non runnare questo comando
    np.savetxt(FileToSave, data, fmt=fmt)