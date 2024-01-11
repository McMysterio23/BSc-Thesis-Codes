import numpy as np

def Kauffmann_03(x):
    # Prima funzione per i BPT del primo tipo
    return 0.61 / (np.log(x) - 0.05) + 1.3

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