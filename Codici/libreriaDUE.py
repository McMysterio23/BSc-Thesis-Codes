from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import QTable
from astropy.table import Column
from astropy.table import Table


from libreriaUNO import leggi_file_txt, leggi_nomi_colonne


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
Intersezione = PATH+"/INDICI_BCG/Matricione.txt"
indices = PATH+"/INDICI_BCG/Int_2_Asec.txt"

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



def plot_scatter(array1, array2, funzione1, funzione2, title="Scatter Plot", s=15, x_label="X", y_label="Y", Bilogaritmico = False, Overplotting = False,
                  NomeFunzione1 = 'Inserire il nome della prima funzione', 
                  NomeFunzione2 = 'Inserire il nome della seconda funzione'):
    


    Xmin, Xmax         = -1.2, 1.2   # Define the maximum and minimum limit in X-axis
    Ymin, Ymax         = -1.5, 1.0   # Define the maximum and minimum limit in Y-axis
    # Crea uno scatter plot
    plt.scatter(array1, array2)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if Bilogaritmico:
        # Impostazione della scala logaritmica sugli assi
        plt.xscale('log')
        plt.yscale('log')

    if Overplotting:
        # Dati di esempio per le funzioni
        x_func = np.linspace(0.001, 1, 1000)  
        x_fun2 = np.linspace(0.95, 1.35, 1000)  
        y_func1 = funzione1(x_func)
        y_func2 = funzione2(x_func)
        # Sovrapposizione del grafico delle funzioni
        plt.plot(x_func, y_func1, label=NomeFunzione1)
        plt.plot(x_func, y_func2, label=NomeFunzione2)
        plt.xlim(Xmin, Xmax)
        plt.ylim(Ymin, Ymax)
        plt.legend()
    


    # Mostra il plot a video
    plt.show()




def Aggiungi_Colonna(arr1, Arr2, arr3, arr4, arr5, arr6, arr7=np.ones(dimensione), arr8=np.zeros(dimensione), arr9=np.zeros(dimensione),
                     nome1='Nome Colonna1', nome2='Nome Colonna2', nome3='Nome Colonna3',
                     nome4='Nome Colonna4', nome5='Nome Colonna5', nome6='Nome Colonna6',
                     nome7='Nome Colonna7', nome8='Nome Colonna8', nome9='Nome Colonna9'):

    # Assicurati che gli array siano di tipo sequenza (es. liste o array NumPy)
    for arr in [arr1, Arr2, arr3, arr4, arr5, arr6]:
        if not isinstance(arr, (list, np.ndarray)):
            raise ValueError("Gli array devono essere di tipo sequenza (es. liste o array NumPy)")    


    Sample_derived = Table([Column(indiciint2, name='Indici_INFO'), Column(coordinate[:,0], name='PlateID'), Column(coordinate[:,1], name='FiberID'),
                            Column(arr1, name=nome1), Column(Arr2, name= nome2), Column(arr3, name = nome3), Column(arr4, name = nome4),
                            Column (arr5, name = nome5), Column(arr6, name = nome6), Column(arr7, name= nome7), Column(arr8, name = nome8),
                            Column(arr9, name = nome9)])
    return Sample_derived

def Aggiungi_Colonna6(arr1, arr2, arr3, arr4, arr5, arr6, nome1='Nome Colonna1', nome2='Nome Colonna2', nome3='Nome Colonna3',
                     nome4='Nome Colonna4', nome5='Nome Colonna5', nome6='Nome Colonna6'):

    # Assicurati che gli array siano di tipo sequenza (es. liste o array NumPy)
    for arr in [arr1, arr2, arr3, arr4, arr5, arr6]:
        if not isinstance(arr, (list, np.ndarray)):
            raise ValueError("Gli array devono essere di tipo sequenza (es. liste o array NumPy)")    

    Sample_derived = Table([Column(indiciint2, name='Indici_INFO'), Column(coordinate[:, 0], name='PlateID'),
                            Column(coordinate[:, 1], name='FiberID'), Column(arr1, name=nome1),
                            Column(arr2, name=nome2), Column(arr3, name=nome3), Column(arr4, name=nome4),
                            Column(arr5, name=nome5), Column(arr6, name=nome6)])
    
    return Sample_derived
   
def AddColumn(table_esistente, nuova_colonna, nome_colonna):
    """
    Aggiunge una nuova colonna a un oggetto Table esistente.

    Parametri:
    - table_esistente: Oggetto Table a cui aggiungere la colonna.
    - nuova_colonna: Array numpy o lista rappresentante la nuova colonna da aggiungere.
    - nome_colonna: Nome della nuova colonna.

    Restituisce:
    - Un nuovo oggetto Table con la colonna aggiunta.
    """

    # Verifica se table_esistente è un oggetto Table
    if not isinstance(table_esistente, Table):
        raise ValueError("L'argomento 'table_esistente' deve essere un oggetto Table.")

    # Verifica se nuova_colonna è una sequenza (array o lista)
    if not isinstance(nuova_colonna, (np.ndarray, list)):
        raise ValueError("La nuova colonna deve essere una sequenza (array o lista).")

    # Verifica se la lunghezza di nuova_colonna è la stessa delle righe in table_esistente
    if len(nuova_colonna) != len(table_esistente):
        raise ValueError("La lunghezza della nuova colonna deve essere la stessa delle righe in table_esistente.")

    # Aggiunge la nuova colonna a table_esistente
    table_esistente[nome_colonna] = Column(nuova_colonna, name=nome_colonna)

    return table_esistente




    


def Lettura_Colonne_RawDATA(indice):
    """
    Questa funzione restituisce un array contenente i nomi di tutte le colonne presenti nel file fits di dati Grezzi
    presente all'indice i della lista di percorsi Lista_Di_Percorsi_Ai_Dati
    """

    nomi_colonne = leggi_nomi_colonne(Lista_Di_Percorsi_Ai_Dati[indice])
    #print('I nomi delle colonne nel file :', indice, 'sono :')
    #àprint(nomi_colonne)
    return nomi_colonne

def Aggiungi_Colonna6plus(arr1, arr2, arr3, arr4, arr5, arr6, indiciint2, coordinate, nome1='Nome Colonna1', nome2='Nome Colonna2', nome3='Nome Colonna3',
                     nome4='Nome Colonna4', nome5='Nome Colonna5', nome6='Nome Colonna6'):

    # Assicurati che gli array siano di tipo sequenza (es. liste o array NumPy)
    for arr in [arr1, arr2, arr3, arr4, arr5, arr6]:
        if not isinstance(arr, (list, np.ndarray)):
            raise ValueError("Gli array devono essere di tipo sequenza (es. liste o array NumPy)")

    # Crea una tabella Astropy
    coo1 = np.array(coordinate[:,0])
    coo2 = np.array(coordinate[:,1])

    Sample_derived = Table([Column(indiciint2, name='Indici_INFO'), Column(coo1, name='PlateID'),
                            Column(coo2, name='FiberID'), Column(arr1, name=nome1),
                            Column(arr2, name=nome2), Column(arr3, name=nome3), Column(arr4, name=nome4),
                            Column(arr5, name=nome5), Column(arr6, name=nome6)])

    # Converte la tabella in un DataFrame di pandas
    dataframe = Sample_derived.to_pandas()

    # Restituisci un array NumPy
    return dataframe.values



def Crossing(Indice_file, Base = False, Intersezione=True):
    """
    Questa funzione restituisce l'array contenente i dati di tutte le colonne del file analizzato secondo l'indice dell'array
    Lista_di_percorQuesta funzione restituisce l'array contenente i dati di tutte le colonne del file analizzato secondo l'indice di ell'array
    Lista_di_percorsi. Al momento l'implementazione della selezione Intersezione si/no non è ancora implementata !
    Ultima modifica, restituisce in seconda istanza anche un array con il nome delle colonne effettivamente selezionate !!!
    
    Ricordare che :
    Indice = 1 ) File Fits LINE
    Indice = 2 ) File Fits INDX
    Indice = 3 ) File Fits SFRS
    Indice = 4 ) File Fits SPECSFR
    Indice = 5 ) File Fits TOTLGM
    Indice = 6 ) File Fits SPECMAG
    Indice = 7 ) File Fits FIBLGM
    Indice = 8 ) File Fits FIBOH
    """

    # Apro il file FITS che ti interessa
    obj = fits.open(Lista_Di_Percorsi_Ai_Dati[Indice_file])
    data = obj[1].data


    #print(data.columns.names)
    
    
    # Filtra solo le colonne di tipo float e lunghezza 1
    if (Base == True):
        colonne_selezionate = [colonna for colonna in data.columns.names if (data[colonna].dtype.kind == 'f' or np.issubdtype(data[colonna].dtype, np.integer)) and data[colonna].shape == (len(data),)
                           and colonna not in ['RELEASE', 'MJD']]
    else:
       #Filtra solo le colonne di tipo float e lunghezza 1 con nome diverso da PLATEID e FIBERID
        colonne_selezionate = [colonna for colonna in data.columns.names if
                      (data[colonna].dtype.kind == 'f' or np.issubdtype(data[colonna].dtype, np.integer)) and
                       data[colonna].shape == (len(data),) and
                       colonna not in ['PLATEID', 'FIBERID']] 
    
   
    
    # Ottieni la dimensione dai dati
    dimensione = len(indiciint2)
    print('Sto Selezionando le seguenti colonne con i criteri attuali :\n',colonne_selezionate)

    # Creazione di un array 2D per contenere i valori estratti da tutte le colonne per ciascun indice

    intersection = np.ones((dimensione, len(colonne_selezionate)))
    

    # Popola l'array con i valori delle colonne filtrate
    for i in range(dimensione):
        for j, colonna in enumerate(colonne_selezionate):
         # Sostituisci il valore eventualmente mancante con la media della colonna
         data[colonna] = np.where(np.isnan(data[colonna]), np.nanmean(data[colonna]), data[colonna])
         punto = int(indiciint2[i, 0])
         rra = data[colonna]
         intersection[i, j] = rra[punto]
    return intersection, colonne_selezionate




def CrossingColonna(Indice_file, Colonna_da_analizzare, Base=False, Intersezione=True):
    """
    Questa funzione restituisce una colonna contenente i dati della colonna specificata
    del file FITS analizzato secondo l'indice dell'array Lista_di_percorsi.

    Parametri:
    - Indice_file: Indice del file FITS da analizzare.
    - Colonna_da_analizzare: Nome della colonna da estrarre.
    - Base: Se True, usa una base diversa di colonne selezionate.
    - Intersezione: Non ancora implementato.

    Restituisce:
    - Una colonna con i dati della colonna specificata.
    - Un array con i nomi delle colonne effettivamente selezionate.
    """

    # Apro il file FITS che ti interessa
    obj = fits.open(Lista_Di_Percorsi_Ai_Dati[Indice_file])
    data = obj[1].data

    # Filtra solo le colonne di tipo float e lunghezza 1
    if Base:
        colonne_selezionate = [colonna for colonna in data.columns.names if
                               (data[colonna].dtype.kind == 'f' or np.issubdtype(data[colonna].dtype, np.integer)) and
                               data[colonna].shape == (len(data),) and
                               colonna not in ['RELEASE', 'MJD']]
    else:
        # Filtra solo le colonne di tipo float e lunghezza 1 con nome diverso da PLATEID e FIBERID
        colonne_selezionate = [colonna for colonna in data.columns.names if
                               (data[colonna].dtype.kind == 'f' or np.issubdtype(data[colonna].dtype, np.integer)) and
                               data[colonna].shape == (len(data),) and
                               colonna not in ['PLATEID', 'FIBERID']]

    # Verifica se la colonna specificata è presente nelle colonne selezionate
    if Colonna_da_analizzare not in colonne_selezionate:
        raise ValueError(f"La colonna {Colonna_da_analizzare} non è presente nelle colonne selezionate.")

    # Ottieni la dimensione dai dati
    dimensione = len(indiciint2)
    print('Sto Selezionando la colonna con il nome:', Colonna_da_analizzare)

    # Creazione di un array per contenere i valori estratti dalla colonna per ciascun indice
    colonna_selezionata = np.ones(dimensione)

    # Popola l'array con i valori della colonna filtrata
    for i in range(dimensione):
        punto = int(indiciint2[i, 0])
        rra = data[Colonna_da_analizzare]

        # Sostituisci il valore eventualmente mancante con la media della colonna
        data[Colonna_da_analizzare] = np.where(np.isnan(data[Colonna_da_analizzare]),
                                               np.nanmean(data[Colonna_da_analizzare]),
                                               data[Colonna_da_analizzare])

        colonna_selezionata[i] = rra[punto]

    return colonna_selezionata, colonne_selezionate