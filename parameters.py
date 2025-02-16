import pandas as pd
import numpy as np
#user inputs: 
calc_type = 'SNR' #limmag|spec_per_hour|SNR

#SNR
spec_type = 'bb' #bb|flat|stellar|WD
stellar_type = None
Teff = None 
logg = None

bb_temp = 20000
AB_mag_renorm = 18
wl=5000 #wavelength for limiting magnitude
T_norm = 1200 #exposure time for SNR calc
n_tel_group = 4 #number of telescopes observing the same target for SNR calculation

#limmag
sigma_limit = 10 #SNR limit for limiting magnitude
n_tel_arr=[1,4,20] #numer of telescopes observing the same target for limiting magnitudes
best_sky_brightness_surface_den=20.9 #mag/arcsec**2
delta_sky = 0 # difference from Dark
Sky_brightness_surface_den = best_sky_brightness_surface_den + delta_sky 
T_exp=2000 #max exposure time for plots
Type = 'per pixel'

#spec_per_hour
mag_analyze = np.linspace(13,20,15)       
overhead_sec = 300
SNR_arr = [10,20,50]
n_tel=20 #numer of telescopes per full array. 


## static parameters
# spectrograph parameters
#DeepSpec
magnification = 42/25 #camera to fiber magnification
fiber=25*magnification # projected fiber size in microns
DC=0.0005 # Dark current 
read_noise=2.9 # Read noise
mu=13 #microns per detector pixels
binning=[3,1] #pixel binning
bias = 4000 #bias count in electrons
saturation = 95000 #saturation count in electrons
#telescope parameters 
r_cm=61/2 #cm
eff_area = np.pi*r_cm**2*0.95
F_num=3 #Telescope focal ratio
f_mm=2*F_num*r_cm*10 #focal length in mm
Plate_Scale=206265*mu/1000/f_mm #arcsec/pixel
N_pack=1 #number of MAST arrays

#astronomical conditions
#Plate_Scale=0.44
bkg2src=1

#As built throughputs
res_wl = pd.read_csv("resolution_as_built.csv",sep = ',')
throughput= pd.read_csv("slit_to_e_thorughput.csv",sep = ',')
res_wl['lambda'] = res_wl['lambda']*10
throughput['lambda'] = throughput['lambda']*10
#tel_ref = 0.85 #minimum acceptence critetia in the 400-700 nm range
#fiber_eff =0.86
#Eff_wl = np.array(throughput['T']) * fiber_eff *  tel_ref
Eff_wl = np.array(throughput['T_sky'])
throughput['Eff_tot'] = Eff_wl

#calculate resolution, efficiency and size of resolution element in ang
R =  np.interp(wl,res_wl['lambda'],res_wl['Res_25um'])
R = int(R)
Eff = np.interp(wl,throughput['lambda'],throughput['Eff_tot'])
d_lam=wl/R



#plotting parameters
T_exp_vec=10**np.linspace(0,np.log10(T_exp),100)
lam = np.linspace(3700,9000,1800)
vega_path = r"alpha_lyr_004.dat"
V_band = np.loadtxt('./johnson_v.txt',delimiter = ',')
V_band[:,1] = V_band[:,1]/np.trapz(V_band[:,1],V_band[:,0])

#paths for stellar spectra which can be selected in the GUI for SNR calculations
stellar_path_dic = {'O5 V':'./PicklesStellarSpec/o5v.mat',
            'O9 V':'./PicklesStellarSpec/o9v.mat',
            'B0 V':'./PicklesStellarSpec/b0v.mat',
            'B1 V':'./PicklesStellarSpec/b1v.mat',
            'B3 V':'./PicklesStellarSpec/b3v.mat',
            'B8 V':'./PicklesStellarSpec/b8v.mat',
            'B9 V':'./PicklesStellarSpec/b9v.mat',
            'A0 V':'./PicklesStellarSpec/a0v.mat',
            'A2 V':'./PicklesStellarSpec/a2v.mat',
            'A3 V':'./PicklesStellarSpec/a3v.mat',
            'A5 V':'./PicklesStellarSpec/a5v.mat',
            'A7 V':'./PicklesStellarSpec/a7v.mat',
            'F0 V':'./PicklesStellarSpec/f0v.mat',
            'F2 V':'./PicklesStellarSpec/f2v.mat',
            'F5 V':'./PicklesStellarSpec/f5v.mat',
            'F6 V':'./PicklesStellarSpec/f6v.mat',
            'F8 V':'./PicklesStellarSpec/f8v.mat',
            'G0 V':'./PicklesStellarSpec/g0v.mat',
            'G2 V':'./PicklesStellarSpec/g2v.mat',
            'G5 V':'./PicklesStellarSpec/g5v.mat',
            'G8 V':'./PicklesStellarSpec/g8v.mat',
            'K0 V':'./PicklesStellarSpec/k0v.mat',
            'K2 V':'./PicklesStellarSpec/k2v.mat',
            'K3 V':'./PicklesStellarSpec/k3v.mat',
            'K4 V':'./PicklesStellarSpec/k4v.mat',
            'K5 V':'./PicklesStellarSpec/k5v.mat',
            'K7 V':'./PicklesStellarSpec/k7v.mat',
            'M0 V':'./PicklesStellarSpec/m0v.mat',
            'M1 V':'./PicklesStellarSpec/m1v.mat',
            'M2 V':'./PicklesStellarSpec/m2v.mat',
            'M3 V':'./PicklesStellarSpec/m3v.mat',
            'M4 V':'./PicklesStellarSpec/m4v.mat',
            'M5 V':'./PicklesStellarSpec/m5v.mat',
            'M6 V':'./PicklesStellarSpec/m6v.mat'}

import glob
import os
WD_names = glob.glob('./KoesterWDmodels/*.txt')
WD_dic = {os.path.basename(x).split('.')[0]:x for x in WD_names}

Teff = [int(key.split('_')[0].split('da')[1]) for key in WD_dic.keys()]
logg = [float(key.split('_')[1].split('.')[0])/100 for key in WD_dic.keys()]

base_names = list(WD_dic.keys())
WD_path_dic = {(Teff[i],logg[i] ): WD_dic[base_names[i]] for i in range(len(base_names)) }


#label details for the stellar spectra which can be selected in the GUI. 
#Value corresponds to the key for stellar_path_dic, and label is what will be displayed in the GUI dropdown menu
counts = None                                