"""
Questo file contiene tutte le componenti necessarie per la creazione dei BPT diagrams. tutto però è fermo al momento prima dell'accelerata di metà gennaio!!!
"""

import numpy as np
import scipy as sc
import pylab as plt
from astropy.io import fits


def Kauffmann_03(x):
    # Prima funzione per i BPT del primo tipo
    return 0.61 / (np.log10(x) - 0.05) + 1.3


def Kauffmann_03_specific(x):
    """
    Questa funzione è da intendersi solamente per la realizzazione della colorazione del BPT type1 in quanto, riceve in ingresso
    numeri di cui è già stato fatto il logaritmo !!!
    """
    return 0.61 / (x - 0.05) + 1.3



def Kewley_01(x):
    # Seconda Funzione per i BPT del primo tipo
    return 0.61 / (np.log(x) - 0.47) + 1.19

def main_AGN_line(x):
    """
    Calcola log(y) per la main AGN line.

    Parametri:
    - x: Valore dell'ascissa.

    Restituisce:
    - log(y) per la main AGN line.
    """
    return 0.72 / (np.log(x) - 0.32) + 1.30

def LINER_Sy2_line(x):
    """
    Calcola log(y) per la LINER/Sy2 line.

    Parametri:
    - x: Valore dell'ascissa.

    Restituisce:
    - log(y) per la LINER/Sy2 line.
    """
    return 1.89 * np.log(x) + 0.76

def BPT_type1(LOGarr_x, LOGarr_y, Histogram = False):
   if Histogram: 
        f      = plt.figure(1)

        

        # Specifying integer values for bins
        bins_X, bins_Y = 40, 40




        ax                 = f.add_subplot(1,1,1)
        #bins_X, bins_Y     =  60., 60.   # Define the number of bins in X- and Y- axis
        Xmin, Xmax         = -1.2, 1.2   # Define the maximum and minimum limit in X-axis
        Ymin, Ymax         = -1.5, 1.0   # Define the maximum and minimum limit in Y-axis
        Nlevels            = 4           # Define the number of levels of isocontour


        hist,xedges,yedges = np.histogram2d(LOGarr_x,LOGarr_y,bins=(bins_X, bins_Y),range=[[Xmin,Xmax],[Ymin,Ymax]])
        masked             = np.ma.masked_where(hist==0, hist)
        plotting           = ax.imshow(masked.T,extent=[Xmin, Xmax, Ymin, Ymax],interpolation='nearest',origin='lower',cmap=plt.cm.gray_r)
        levels             = np.linspace(0., np.log10(masked.max()), Nlevels)[1:]
        CS                 = ax.contour(np.log10(masked.T), levels, colors='k',linewidths=1,extent=[Xmin,Xmax,Ymin,Ymax])


            # Kewley+01 ------------------------------------------
        X_kewley = np.linspace(0.001, 1)
        Y_kewley = (0.61 / (np.log10(X_kewley) - 0.47)) + 1.19

        # Schawinski+07 --------------------------------------
        X_schawinski = np.linspace(0.65, 2)
        Y_schawinski = 1.05 * np.log10(X_schawinski) + 0.45

            # Kauffmann+03 ---------------------------------------
        X_kauffmann = np.linspace(0.001, 1)
        Y_kauffmann = 0.61 / (np.log10(X_kauffmann) - 0.05) + 1.3
        # Regions --------------------------------------------
        ax.plot(X_kewley,   Y_kewley, '-' , color='blue', lw=3, label='Kewley+01'    ) # Kewley+01
        ax.plot(X_schawinski, Y_schawinski, '-', color='red', lw=5, label='Schawinski+07') # Schawinski+07
        ax.plot(X_kauffmann, Y_kauffmann, '--', color='blue', lw=5, label='Kauffmann+03' ) # Kauffmann+03

        # Axi name here ...
        Nsize = 25
        ax.set_xlabel(r'log([NII] $\lambda$ 6583/H$\alpha$)',fontsize=Nsize)
        ax.set_ylabel(r'log([OIII] $\lambda$ 5007/H$\beta$)',fontsize=Nsize)
        ax.tick_params(labelsize = Nsize)
        ax.set_ylim(Ymin, Ymax)
        ax.set_xlim(Xmin, Xmax)
        plt.legend()

        plt.show()
   else:
     f = plt.figure(1)

     ax = f.add_subplot(1, 1, 1)

        # Definizione del numero di bin negli assi X e Y
     Xmin, Xmax = -0.75, 1.2  # Limiti massimi e minimi sull'asse X
     Ymin, Ymax = -0.75, 1.0  # Limiti massimi e minimi sull'asse Y

        # Kewley+01 ------------------------------------------
     X_kewley = np.linspace(0.001, 1, 100)
     Y_kewley = 0.61 / (np.log10(X_kewley) - 0.47) + 1.19

     # Schawinski+07 --------------------------------------
     X_schawinski = np.linspace(0.3823, 2, 100)
     Y_schawinski = 1.05 * np.log10(X_schawinski) + 0.45

        # Kauffmann+03 ---------------------------------------
     X_kauffmann = np.linspace(0.001, 1)
     Y_kauffmann = 0.61 / (np.log10(X_kauffmann) - 0.05) + 1.3
     ax.scatter(LOGarr_x, LOGarr_y, color='black', marker='.', label='Data Points')

        # Plot delle linee delle regioni
     ax.plot(X_kewley, Y_kewley, '-', color='blue', lw=0.5, label='Kewley+01')  # Kewley+01
     ax.plot(X_schawinski, Y_schawinski, '-', color='red', lw=1, label='Schawinski+07')  # Schawinski+07
     ax.plot(X_kauffmann, Y_kauffmann, '--', color='yellow', lw=3, label='Kauffmann+03')  # Kauffmann+03

        # Etichette degli assi
     Nsize = 20
     ax.set_xlabel(r'log([NII] $\lambda$ 6583/H$\alpha$)', fontsize=Nsize)
     ax.set_ylabel(r'log([OIII] $\lambda$ 5007/H$\beta$)', fontsize=Nsize)
     ax.tick_params(labelsize=Nsize)
     ax.set_ylim(Ymin, Ymax)
     ax.set_xlim(Xmin, Xmax)
     plt.legend()

     plt.show()
     





def BPT_type2(LOGarr_sii_ha, LOGarr_oiii_hb, Histogram = False):
    

    if Histogram :
        f = plt.figure(1)

        # Specifica dei valori interi per i bin
        bins_X, bins_Y = 60, 60

        ax = f.add_subplot(1, 1, 1)

        # Definizione del numero di bin negli assi X e Y
        Xmin, Xmax = -1.5, 1.5  # Limiti massimi e minimi sull'asse X
        Ymin, Ymax = -1.5, 1.5  # Limiti massimi e minimi sull'asse Y
        Nlevels = 6  # Numero di livelli degli isocontorni

        # Creazione dell'istogramma bidimensionale
        hist, xedges, yedges = np.histogram2d(LOGarr_sii_ha, LOGarr_oiii_hb, bins=(bins_X, bins_Y),
                                            range=[[Xmin, Xmax], [Ymin, Ymax]])
        # Mascheramento delle regioni con conteggio zero
        masked = np.ma.masked_where(hist == 0, hist)

        # Plot dell'istogramma come immagine
        plotting = ax.imshow(masked.T, extent=[Xmin, Xmax, Ymin, Ymax], interpolation='nearest', origin='lower',
                            cmap=plt.cm.gray_r)

        # Creazione delle linee di isocontorno
        levels = np.linspace(0., np.log10(masked.max()), Nlevels)[1:]
        CS = ax.contour(np.log10(masked.T), levels, colors='k', linewidths=1, extent=[Xmin, Xmax, Ymin, Ymax])

        # Nuove linee
        X_main_agn = np.linspace(0.08, 2, 100)  
        Y_main_agn = 0.72 / (np.log10(X_main_agn) - 0.32) + 1.30

        X_liner_sy2 = np.linspace(0.49, 2, 100)
        Y_liner_sy2 = 1.89 * np.log10(X_liner_sy2) + 0.76


        # Plot delle nuove linee
        ax.plot(X_main_agn, Y_main_agn, '-', color='blue', lw=3, label='Main AGN line') 
        ax.plot(X_liner_sy2, Y_liner_sy2, '-', color='red', lw=5, label='LINER/Sy2 line') 

        # Etichette degli assi
        Nsize = 25
        ax.set_xlabel(r'log([SII] $\lambda\lambda$ 6716, 6731/H$\alpha$)', fontsize=Nsize)
        ax.set_ylabel(r'log([OIII] $\lambda$ 5007/H$\beta$)', fontsize=Nsize)
        ax.tick_params(labelsize=Nsize)
        ax.set_ylim(Ymin, Ymax)
        ax.set_xlim(Xmin, Xmax)
        plt.legend()

        plt.show()
    else :
        
        f = plt.figure(1)

        ax = f.add_subplot(1, 1, 1)

        # Definizione del numero di bin negli assi X e Y
        Xmin, Xmax = -1.5, 1.5  # Limiti massimi e minimi sull'asse X
        Ymin, Ymax = -1.5, 1.5  # Limiti massimi e minimi sull'asse Y
        Nlevels = 6  # Numero di livelli degli isocontorni

        # Nuove linee
        X_main_agn = np.linspace(0.08, 2, 100)  
        Y_main_agn = 0.72 / (np.log10(X_main_agn) - 0.32) + 1.30

        X_liner_sy2 = np.linspace(0.49, 2, 100)
        Y_liner_sy2 = 1.89 * np.log10(X_liner_sy2) + 0.76

        # Plot delle nuove linee
        ax.plot(X_main_agn, Y_main_agn, '-', color='blue', lw=1.5, label='Main AGN line') 
        ax.plot(X_liner_sy2, Y_liner_sy2, '-', color='red', lw=1, label='LINER/Sy2 line') 

        # Plot dei punti
        ax.scatter(LOGarr_sii_ha, LOGarr_oiii_hb, color='black', marker='.', label='Data Points')

        # Etichette degli assi
        Nsize = 25
        ax.set_xlabel(r'log([SII] $\lambda\lambda$ 6716, 6731/H$\alpha$)', fontsize=Nsize)
        ax.set_ylabel(r'log([OIII] $\lambda$ 5007/H$\beta$)', fontsize=Nsize)
        ax.tick_params(labelsize=Nsize)
        ax.set_ylim(Ymin, Ymax)
        ax.set_xlim(Xmin, Xmax)
        plt.legend()

        plt.show()



def BPT_type2_colored(LOGarr_sii_ha, LOGarr_oiii_hb):
    f = plt.figure(1)

    ax = f.add_subplot(1, 1, 1)

    # Definizione del numero di bin negli assi X e Y
    Xmin, Xmax = -1.5, 1.5  # Limiti massimi e minimi sull'asse X
    Ymin, Ymax = -1.5, 1.5  # Limiti massimi e minimi sull'asse Y
    Nlevels = 6  # Numero di livelli degli isocontorni

    # Nuove linee
    X_main_agn = np.linspace(-1.5, 0.5)
    Y_main_agn = 0.72 / (np.log10(X_main_agn) - 0.32) + 1.30

    X_liner_sy2 = np.linspace(-1.5, 1.5)
    Y_liner_sy2 = 1.89 * np.log10(X_liner_sy2) + 0.76

    # Plot delle nuove linee
    ax.plot(X_main_agn, Y_main_agn, '-', color='blue', lw=3, label='Main AGN line') 
    ax.plot(X_liner_sy2, Y_liner_sy2, '-', color='red', lw=5, label='LINER/Sy2 line') 

    # Definizione delle regioni
    region_main_agn = (LOGarr_oiii_hb > 0.72 / (np.log10(LOGarr_sii_ha) - 0.32) + 1.30)
    region_liner_sy2 = (LOGarr_oiii_hb > 1.89 * np.log10(LOGarr_sii_ha) + 0.76)

    # Plot dei punti colorati
    ax.scatter(LOGarr_sii_ha[~region_main_agn & ~region_liner_sy2], 
               LOGarr_oiii_hb[~region_main_agn & ~region_liner_sy2], 
               color='black', marker='.', label='Other Points')

    ax.scatter(LOGarr_sii_ha[region_main_agn], 
               LOGarr_oiii_hb[region_main_agn], 
               color='blue', marker='.', label='Main AGN Points')

    ax.scatter(LOGarr_sii_ha[region_liner_sy2], 
               LOGarr_oiii_hb[region_liner_sy2], 
               color='red', marker='.', label='LINER/Sy2 Points')

    # Etichette degli assi
    Nsize = 25
    ax.set_xlabel(r'log([SII] $\lambda\lambda$ 6716, 6731/H$\alpha$)', fontsize=Nsize)
    ax.set_ylabel(r'log([OIII] $\lambda$ 5007/H$\beta$)', fontsize=Nsize)
    ax.tick_params(labelsize=Nsize)
    ax.set_ylim(Ymin, Ymax)
    ax.set_xlim(Xmin, Xmax)
    plt.legend()

    plt.show()

def BPTtype1_colored(LOGarr_x, LOGarr_y):
    f = plt.figure(1)

    ax = f.add_subplot(1, 1, 1)

    # Definizione del numero di bin negli assi X e Y
    Xmin, Xmax = -1.2, 1.2  # Limiti massimi e minimi sull'asse X
    Ymin, Ymax = -1.5, 1.0  # Limiti massimi e minimi sull'asse Y
    #Nlevels = 4  # Numero di livelli degli isocontorni

    # Kewley+01 ------------------------------------------
    X_kewley = np.linspace(0.001, 1)
    Y_kewley = (0.61 / (np.log10(X_kewley) - 0.47)) + 1.19

    # Schawinski+07 --------------------------------------
    X_schawinski = np.linspace(-0.180, 1.5)
    Y_schawinski = 1.05 * X_schawinski + 0.45

    # Kauffmann+03 ---------------------------------------
    X_kauffmann = np.linspace(0.001, 1)
    Y_kauffmann = 0.61 / (np.log10(X_kauffmann) - 0.05) + 1.3
    # Plot delle linee delle regioni
    ax.plot(X_kewley, Y_kewley, '-', color='blue', lw=1, label='Kewley+01')  # Kewley+01
    ax.plot(X_schawinski, Y_schawinski, '-', color='red', lw=2, label='Schawinski+07')  # Schawinski+07
    ax.plot(X_kauffmann, Y_kauffmann, '--', color='orange', lw=2, label='Kauffmann+03')  # Kauffmann+03

    # Definizione delle regioni

    # Inizializza un array per la popolazione
    popolazione = np.full_like(LOGarr_x, "Others", dtype=np.object_)

    # Imposta la popolazione 1 (AGNs)
    popolazione[(LOGarr_x > 0) & (LOGarr_y > Kauffmann_03_specific(LOGarr_x))] = "AGNs"

    # Imposta la popolazione 2 (SF)
    popolazione[(LOGarr_x < 0) & (LOGarr_y < Kauffmann_03_specific(LOGarr_x))] = "SF"

    # Crea un plot delle funzioni (sostituisci con le tue funzioni reali)
    # Ad esempio, plot delle funzioni y = tua_funzione(LOGarr_x) e y = -tua_funzione(LOGarr_x)
    #plt.plot(LOGarr_x, Kauffmann_03_specific(LOGarr_x), label="Tua Funzione")
    #plt.plot(LOGarr_x, Kauffmann_03_specific(LOGarr_x), label="-Tua Funzione")

    # Crea un scatter plot per i punti delle popolazioni
    plt.scatter(LOGarr_x[popolazione == "AGNs"], LOGarr_y[popolazione == "AGNs"], color='red', label='AGNs')
    plt.scatter(LOGarr_x[popolazione == "SF"], LOGarr_y[popolazione == "SF"], color='blue', label='SF')
    plt.scatter(LOGarr_x[popolazione == "Others"], LOGarr_y[popolazione == "Others"], color='green', label='Others')

    

    # Plot dei punti colorati

    """
    ax.scatter(LOGarr_x[~region_kewley & ~region_schawinski & ~region_kauffmann], 
               LOGarr_y[~region_kewley & ~region_schawinski & ~region_kauffmann], 
               color='black', marker='.', label='Other Points')
               """
    """
    ax.scatter(LOGarr_x[region_kewley], 
               LOGarr_y[region_kewley], 
               color='blue', marker='.', label='Kewley+01 Points')
               """

    """
    ax.scatter(LOGarr_x[region_schawinski], 
               LOGarr_y[region_schawinski], 
               color='red', marker='.', label='LINERs')
               """
    """
    ax.scatter(LOGarr_x[region_kauffmann], 
               LOGarr_y[region_kauffmann], 
               color='blue', marker='.', label='AGNs')
               """

    # Etichette degli assi
    Nsize = 25
    ax.set_xlabel(r'log([NII] $\lambda$ 6583/H$\alpha$)', fontsize=Nsize)
    ax.set_ylabel(r'log([OIII] $\lambda$ 5007/H$\beta$)', fontsize=Nsize)
    ax.tick_params(labelsize=Nsize)
    ax.set_ylim(Ymin, Ymax)
    ax.set_xlim(Xmin, Xmax)
    plt.legend()

    plt.show()
