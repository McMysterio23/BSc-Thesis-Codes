import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Column

from libreriaUNO import salva_array_senza_parentesi
from libreriaUNO import leggi_file_txt
from libreriaDUE import plot_and_save_scatter
from libreriaDUE import plot_scatter
from libreriaDUE import Intersection_Array_Type1

def main():

    Redshift = Intersection_Array_Type1('z')
    Errore_Redshift = Intersection_Array_Type1('z_err')
    Sigma = Intersection_Array_Type1('v_disp')
    Errore_Sigma = Intersection_Array_Type1('v_disp_err')
    SN_Medio = Intersection_Array_Type1('sn_median')


    

    
    





    











if __name__ == "__main__":

 main ()