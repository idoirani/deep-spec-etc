
from astropy import constants
import numpy as np
import os
import sys
#from PhotoUtils import * 
import matplotlib.pyplot as plt
from scipy.optimize  import broyden1,broyden2,newton_krylov
import math
from matplotlib.widgets import TextBox
import pandas as pd
import plotly.express as px
from SNR_func import *

## spectrum functions
def bb_F(lam_AA,T_K,r_cm = 1e14,d_cm = 3.08e19):
    '''
    Calculates blackbody flux for a given wavelength, temperature, radius and distance
    IN:         lam_AA - wavelength in Angstrom,
                T_K - temperature in Kelvin,
                r_cm - radius in cm,
                d_cm - distance in cm
    OUT:       flux - flux in erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$
    '''
    flux=1.191/((lam_AA*1e-4)**(5))/(np.exp((1.9864/1.38064852)/((lam_AA*1e-4)*(T_K*1e-4)))-1)*(np.pi*((r_cm/d_cm)*1e10)**2)*1e-13 
    return flux  

def flux_vec_to_mag(lam,flux):
    '''
    Converts flux vector to AB magnitude vector
    IN:        lam - wavelength vector in Angstrom,
               flux - flux array in erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$
    OUT:       m_array - magnitude vector
    '''
    cAA = 3e18 # in AA/s
    m_array = np.zeros_like(lam)
    ZPflux=3631 #Jy
    ZPflux=ZPflux*1e-23
    ZPflux_vec=ZPflux*cAA*(lam**(-2))
    m_array = -2.5*np.log10(flux/ZPflux_vec)
    return m_array 
    
def mag_to_flux(mag,lam):
    '''
    Converts AB magnitude to flux
    IN:         mag - magnitude array,
                lam - wavelength array in Angstrom (of same length as mag)
    OUT:        flux - flux in erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$
    '''
    cAA = 3e18 # in AA/s
    ZPflux=3631 #Jy
    ZPflux=ZPflux*1e-23
    ZPflux_vec=ZPflux*cAA*(lam**(-2))
    flux = ZPflux_vec*10**(-0.4*mag)
    return flux


def generate_bb_mag_spec(lam,T,r = 1e14, d = 3.08e19):
    '''
    Generates blackbody spectrum in magnitudes for a given wavelength, temperature, radius and distance
    IN:         lam - wavelength vector in Angstrom,
                T - temperature in Kelvin,
                r - radius in cm,
                d - distance in cm
    OUT:        m_array - magnitude vector
    '''
    flux = bb_F(lam,T,r,d)                                                        
    m_array = flux_vec_to_mag(lam,flux)
    return m_array 


def generate_flat_spec(lam, norm = 3631*1e-9):
    '''
    Generates flat spectrum in f_nu for a given wavelength and normalization
    IN:         lam - wavelength vector in Angstrom,
                norm - normalization factor (defualt is the AB magnitude zero point flux of 3631 Jy))
    OUT:        flux_vec - flux vector in erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$

    '''
    cAA = 3e18
    flux=norm*1e-23
    flux_vec=flux*cAA*(lam**(-2))
    return flux_vec 

def generate_flat_spec_mag(lam, norm = 3631*1e-9):
    '''
    Generates flat spectrum in magnitudes for a given wavelength and normalization
    IN:         lam - wavelength vector in Angstrom,
                norm - normalization factor (defualt is the AB magnitude zero point flux of 3631 Jy))
    OUT:        m_array - magnitude vector
    '''
    flux  = generate_flat_spec(lam, norm = norm)                                                
    m_array =flux_vec_to_mag(lam,flux)
    return m_array 


def generate_stellar_spec(lam,spec_type):
    import scipy.io
    mat = scipy.io.loadmat(parameters.path_dic[spec_type])
    k = spec_type.replace(' ','').lower()
    tab = mat[k]
    ltab = tab[:,0]
    stab = tab[:,1]
    spec = np.interp(lam,ltab,stab)
    return spec 


def spec_per_hour(SNR_target = 10, mag_analyze = np.linspace(14,20,25), overhead_sec = 300,Xin = [14,17]): 
    '''
    Calculates number of spectra per hour for a given SNR target and source brightness
    IN:         SNR_target - SNR target,
                mag_analyze - source brightness vector,
                overhead_sec - overhead time in seconds
    OUT:        spec_per_hour - number of spectra per hour, 
                T_exp_optimal - optimal exposure time per target, 
                n_tel_optimal - optimal number of telescopes per target
    '''
    T_exp_check = np.array([10,50,150,200,300,700,1200,2400])
    T_exp_check_high = np.logspace(1,np.log10(2400),50)
    T_exp_optimal = np.zeros_like(mag_analyze)
    n_tel_optimal = np.zeros_like(mag_analyze)
    spec_per_hour = np.zeros_like(mag_analyze)
    for i,mag in enumerate(mag_analyze): 
        limmag = []
        for n_tel_try in range(1,21):
            limag = find_limiting_mag(T_exp_check,SNR_func,SNR_lim=SNR_target,n_tel=n_tel_try,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,Xin=Xin,N_pack=1,binning=parameters.binning)[0]
            limmag.append(limag)
        limmag = np.array(limmag)     
        cond = limmag>mag
        cond_telescope = np.sum(cond,axis = 1)>0 
        if np.sum(cond_telescope) == 0: 
            spec_per_hour[i]  = np.nan
        else: 
            n_tel_best = np.min(np.arange(1,21)[cond_telescope])
            if n_tel_best>10: 
                n_tel_best = 20
            limag = find_limiting_mag(T_exp_check_high,SNR_func,SNR_lim=SNR_target,n_tel=n_tel_best,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,Xin=Xin,N_pack=1,binning=parameters.binning)[0]
            best_T_Exp = T_exp_check_high[limag>mag][0]
            spec_per_hour[i] = (20//n_tel_best)*3600/(best_T_Exp+overhead_sec)
            n_tel_optimal[i] = n_tel_best
            T_exp_optimal[i] = best_T_Exp
    return spec_per_hour, T_exp_optimal,   n_tel_optimal  


def SynPhot_fast(Lambda,Flux,Filter):
    cAA = 3e18
    ZPflux=3631
    ZPflux=ZPflux*1e-23
    ZPflux_vec=ZPflux*cAA*(Lambda**(-2))
    Filter = Filter.astype(float)
    lam_filt = Filter[:,0]
    T_filt = Filter[:,1]
    T_filt_interp=np.interp(Lambda,lam_filt,T_filt,left=0,right=0)
    trans_flux=np.trapz(Flux*Lambda*T_filt_interp,Lambda)
    norm_flux=np.trapz(ZPflux_vec*T_filt_interp*Lambda,Lambda)
    mag=-2.5*np.log10(trans_flux/norm_flux)
    return mag  

def renorm_spec(lambda_arr, spec,mV_norm):
    '''
    Renormalizes spectrum to a given AB magnitude
    IN:         lambda_arr - wavelength array in Angstrom,
                spec - flux spectrum,
                mV_norm - new V band AB magnitude
    OUT:        spec_norm - renormalized spectrum
    '''

    mag_cur = SynPhot_fast(lambda_arr,spec,parameters.V_band)
    delta_mag = mV_norm - mag_cur
    newspec = spec*10**(-0.4*delta_mag)
    return newspec





def generate_spec(spec_type, stellar_type = None):
    if spec_type == 'bb':
        spec = bb_F(parameters.lam, parameters.bb_temp,r_cm = 1e14,d_cm = 1e24)
    elif spec_type == 'flat':
        spec = generate_flat_spec(parameters.lam)
    elif spec_type == 'stellar':
        spec = generate_stellar_spec(parameters.lam,stellar_type)
    return spec
    


def SNR_sequence(lam, 
                 spec, 
                 AB_mag_renorm,
                 T_exp = parameters.T_norm,
                 binning=parameters.binning, 
                 n_tel=parameters.n_tel_group, 
                 Sky_brightness= parameters.Sky_brightness_surface_den):
    '''
    Generates the projected SNR array and a simulated spectrum for a set of exposure and target parameters
    IN:         lam - wavelength vector in Angstrom,
                spec - flux spectrum of source
                AB_mag_renorm - AB magnitude for spectrum renormalization
                T_exp - exposure time in seconds,
                binning - pixel binning,
                n_tel - number of telescopes,
                Sky_brightness - sky brightness in mag/arcsec^2

    OUT:        spec_sim - simulated spectrum,
                SNR_proj - projected SNR array
    '''
    
    renormed_spec = renorm_spec(lam, spec,AB_mag_renorm)
    renormed_magspec = flux_vec_to_mag(lam,renormed_spec)
    SNR_proj = SNR_in_wl_array(lam, renormed_magspec,T_exp = T_exp,binning=binning, n_tel=n_tel, Sky_brightness= Sky_brightness)
    noise = np.random.normal(size =np.size(renormed_spec) )/SNR_proj
    noise = np.random.normal(size =np.size(renormed_spec) )/SNR_proj
    spec_err = renormed_spec * noise
    spec_sim = renormed_spec+spec_err
    return spec_sim, SNR_proj