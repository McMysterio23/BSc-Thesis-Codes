import numpy as np

def Kauffmann_03(x):
    # Prima funzione per i BPT del primo tipo
    return 0.61 / (np.log(x) - 0.05) + 1.3

def Kewley_01(x):
    # Seconda Funzione per i BPT del primo tipo
    return 0.61 / (np.log(x) - 0.47) + 1.19