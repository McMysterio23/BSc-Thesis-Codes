from astropy.io import fits
from os.path import exists
import numpy as np
import pandas as pd


#Inizializzo i percorsi ai file che utilizzerò
PATH = "/Users/andreamaccarinelli/Desktop/BSc-Thesis-Codes/"
PATHDATA="/Users/andreamaccarinelli/Desktop/SDSS/"
PATH_TESI_IMMAGINI ='/Users/andreamaccarinelli/Desktop/BSc-Thesis/Immagini'

INFO=PATHDATA+"/gal_info_dr7_v5_2.fits"
LINE = PATHDATA+'/gal_line_dr7_v5_2.fits'
INDX = PATHDATA+'/gal_indx_dr7_v5_2.fits'
SFRS = PATHDATA+'gal_totsfr_dr7_v5_2.fits'
SPECSFR = PATHDATA+'/gal_totspecsfr_dr7_v5_2.fits'
TOTLGM = PATHDATA+'/totlgm_dr7_v5_2.fits'
SPECMAG = PATHDATA+'/gal_specmag_dr7_v5_2.fits'
FIBLGM = PATHDATA+'/fiblgm_dr7_v5_2.fits'
FIBOH = PATHDATA+'/gal_fiboh_dr7_v5_2.fits'

file = SFRS

#CARICA fits file
hdul = fits.open(file)  
#Matrice tutti dati
data = hdul[1].data

sigma_bal = data["SIGMA_BALMER"]

