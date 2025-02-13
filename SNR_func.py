
from astropy import constants
import numpy as np
import os
import sys
#from PhotoUtils import * 
from scipy.optimize  import broyden1,broyden2,newton_krylov
import math
import pandas as pd
import parameters




def Sky_bkg_density_per_pixel(mu,Sky_brightness_surface_den,f_mm=1000):
    '''
    Calculates integrated sky background in AB magnitude over the fiber area
    IN: mu - pixel size in microns, Sky_brightness_surface_den - sky brightness in AB mag/arcsec**2, f_mm - telescope focal length in mm
    OUT: Sky_brightness_mag - sky brightness in AB mag
    '''
    Plate_Scale=206265*mu/1000/f_mm #arcsec/pixela
    Sky_brightness_flux_den=10**(-0.4*(Sky_brightness_surface_den))
    Sky_brightness_flux=Sky_brightness_flux_den*Plate_Scale**2
    Sky_brightness_mag=-2.5*np.log10(Sky_brightness_flux)
    return Sky_brightness_mag

def N_phot_in_element(M_ab,T_exp=500,d_lam=parameters.d_lam,r_cm=parameters.r_cm,wl_AA=parameters.wl):
    '''
    Calculates number of photons in a given resolution element given the brightness of the source and telescope parameters
    IN: M_ab - source brightness in AB mag, T_exp - exposure time in seconds, d_lam - resolution element in Angstrom, r_cm - telescope radius in cm, wl_AA - wavelength in Angstrom
    OUT: N_phot - number of photons in the resolution element
    '''
    h=6.62607015e-27
    c=3e10
    wl_cm=wl_AA*1e-8
    nu_cgs=c/wl_cm
    photon_energy=h*nu_cgs
    f_lam=mag2flux(M_ab,wl_AA)
    luminosity=f_lam*d_lam*np.pi*r_cm**2
    photon_flux=luminosity/photon_energy
    N_phot=photon_flux*T_exp
    return N_phot

def SNR_func(M_src,T_exp=500,
             Sky_brightness= parameters.Sky_brightness_surface_den,
             DC=parameters.DC,
             Read=parameters.read_noise,
             G=2,
             Sig_ADU=0.289,
             eff=parameters.Eff,
             n_tel=parameters.n_tel,
             mu=parameters.mu,
             fiber=parameters.fiber,
             d_lam=parameters.d_lam,
             r_cm=parameters.r_cm,
             bkg2src=parameters.bkg2src,
             N_pack=1,binning=[1,1],
             wl_AA=parameters.wl):
    '''
    Calculates SNR for a given source brightness and telescope/spectrograph parameters
    IN:         M_src - source brightness in AB mag, 
                T_exp - exposure time in seconds, 
                Sky_brightness - sky brightness in AB mag/arcsec**2, 
                DC - dark current in electrons/pixel/s, 
                Read - read noise in electrons, 
                G - gain, 
                Sig_ADU - digitization noise, 
                eff - throughput, 
                n_tel - number of telescopes, 
                mu - pixel size in microns, 
                fiber - projected fiber size in microns, 
                d_lam - resolution element in Angstrom, 
                r_cm - telescope radius in cm, 
                bkg2src - background to source ratio, 
                N_pack - number of MAST arrays, 
                binning - pixel binning factor (2 element list, [binX,binY]), 
                wl_AA - wavelength in Angstrom
    OUT:        SNR - signal to noise ratio
                noise_comp - array of noise components (source,backround,dark current, read noise)
    '''

    w=math.ceil(fiber/mu)
    if w<3:
        w=3
    fill_factor=1 #ratio of used to unused pixels which are read out

    N_pix_tot=n_tel*w**2*(1/fill_factor)

    bin_factor=binning[0]*binning[1]
    N_pix_read = N_pix_tot/bin_factor
    N2_readout_Gain=(Read**2+G**2*Sig_ADU**2)*N_pix_read
    N2_Dark=N_pix_tot*DC*T_exp

    pixel_sky_mag=Sky_bkg_density_per_pixel(mu,Sky_brightness,f_mm=parameters.f_mm)
    N_bkg=eff*n_tel*N_phot_in_element(pixel_sky_mag,T_exp=T_exp,d_lam=d_lam,r_cm=r_cm,wl_AA=wl_AA)
    N_src=eff*n_tel*N_phot_in_element(M_src,T_exp=T_exp,d_lam=d_lam,r_cm=r_cm,wl_AA=wl_AA)

    noise=np.sqrt(N_src+(1+1/bkg2src)*(N_bkg+N2_Dark+N2_readout_Gain))
    SNR=N_src/noise
    SNR=SNR*np.sqrt(N_pack)

    noise_comp=[N_src,(1+1/bkg2src)*N_bkg, (1+1/bkg2src)*N2_Dark, (1+1/bkg2src)*N2_readout_Gain]
    return SNR, noise_comp



def find_limiting_mag(T_vec,
                      f_SNR,
                      SNR_lim,
                      n_tel=parameters.n_tel,
                      mu = parameters.mu,
                      Sky_brightness= parameters.Sky_brightness_surface_den,
                      Xin=[15,19], 
                      N_pack=1,binning=[2,2],
                      wl_AA=parameters.wl, 
                      read_noise = parameters.read_noise, 
                      fiber = parameters.fiber ):
    '''
    Finds limiting magnitude for a given SNR limit
    IN:         T_vec - exposure time vector,
                f_SNR - SNR function,
                SNR_lim - SNR limit,
                n_tel - number of telescopes in array,
                mu - pixel size in microns,
                Sky_brightness - sky brightness in AB mag/arcsec**2,
                X_in - initial guess for the rnage in which the limiting magnitude resides,
                N_pack - number of MAST arrays,
                binning - pixel binning factor (2 element list, [binX,binY]),
                wl_AA - wavelength in Angstrom,
                read_noise - read noise in electrons,
                fiber - projected fiber size in microns
    '''
    R =  np.interp(parameters.wl,parameters.res_wl['lambda'],parameters.res_wl['Res_25um'])
    R = parameters.R
    Eff = np.interp(parameters.wl,parameters.throughput['lambda'],parameters.throughput['Eff_tot'])
    d_lam=parameters.d_lam
    def func(M,T):
        res=f_SNR(M_src=M,T_exp=T,n_tel=n_tel,mu=mu,fiber = fiber,Sky_brightness=Sky_brightness,N_pack=N_pack,binning=binning,wl_AA=wl_AA, d_lam=d_lam, eff = Eff,  Read = read_noise)[0]-SNR_lim
        return res
    
    if (isinstance(T_vec,float)|(isinstance(T_vec,int))):
        f=lambda m: func(m,T=T_vec)
        if T_vec>50:
            xin=Xin[1]

        elif T_vec<50:
            xin=Xin[0]
        try:
            try:
                sol=newton_krylov(f,xin=xin,f_tol=1e-2)
            except:
                sol=newton_krylov(f,xin=xin-3,f_tol=1e-2)
        except:
            sol=newton_krylov(f,xin=xin+3,f_tol=1e-2)
        _ , noise_comp =f_SNR(M_src=sol,T_exp=T_vec,n_tel=n_tel,mu=mu,fiber = fiber,Sky_brightness=Sky_brightness,N_pack=N_pack,binning=binning,wl_AA=wl_AA, d_lam=d_lam, eff = Eff,  Read = read_noise)
    else:
        sol=[]
        noise_comp=np.empty_like(np.array([0,0,0,0]))
        for T in T_vec:
            f=lambda m: func(m,T=T)
            if T>50:
                xin=Xin[1]
            elif T<50:
                xin=Xin[0]
            try:
                try:
                    res=newton_krylov(f,xin=xin,f_tol=1e-2)
                except:
                    res=newton_krylov(f,xin=xin-3,f_tol=1e-2)
            except:
                res=newton_krylov(f,xin=xin+3,f_tol=1e-2)
            sol.append(res)
            last=np.array(f_SNR(M_src=sol[0],T_exp=T,n_tel=n_tel,mu=mu,fiber = fiber,Sky_brightness=Sky_brightness,N_pack=N_pack,binning=binning,wl_AA=wl_AA, d_lam=d_lam, eff = Eff,  Read = read_noise)[1])
            noise_comp=np.vstack([noise_comp,last])
        sol=np.array(sol)

    noise_comp=np.delete(noise_comp,0,0)
    return sol, noise_comp


def zp_piv_wl_ab(piv_wl):
    '''
    Calculates zero point magnitude for a given pivot wavelength
    IN:         piv_wl - pivot wavelength in Angstrom
    OUT:        zp_band_ab - zero point magnitude in AB system
    '''
    # in angstrom
    c=constants.c.value*1e10
    zp_band_ab=-2.5*np.log10((piv_wl**2)/c)-48.6
    return  zp_band_ab
def mag2flux(mag,piv_wl):
    zp=zp_piv_wl_ab(piv_wl)
    flux=10**(-0.4*(mag-zp)) 
    return flux
def flux2mag(flux,piv_wl):
    zp=zp_piv_wl_ab(piv_wl)
    mag=-2.5*np.log10(flux)+zp 
    return mag
def magerr2fluxerr(magerr,flux):
    fluxerr=np.abs(-2.303/2.5*magerr*flux)
    return fluxerr
def fluxerr2magerr(fluxerr,flux):
    magerr=np.abs(-2.5/2.303*fluxerr/flux)
    return magerr

def SNR_in_wl_array(lam_array, 
                    M_src_array,
                    T_exp,binning=[1,1], 
                    n_tel = 20, 
                    Sky_brightness= parameters.Sky_brightness_surface_den):
    '''
    Calculates SNR for pairs of source brightness and wavelengths provided as arrays. Uses the as built throughput and resolution at the specified wavelengths to calculate the SNR
    IN:         lam_array - wavelength array,
                M_src_array - source brightness array of the same length as lam_array,
                T_exp - exposure time in seconds,
                binning - pixel binning factor (2 element list, [binX,binY]),
                n_tel - number of telescopes in array
                Sky_brightness - sky brightness surface density in AB mag/sq arcsec
    OUT:        SNR_array - array of SNR values of the same length as lam_array and M_src_array
    '''
    SNR_array = np.zeros_like(lam_array)
    for i,lam in enumerate(lam_array):
        M_src = M_src_array[i]
        R_lam =  np.interp(lam,parameters.res_wl['lambda'],parameters.res_wl['Res_25um'])
        Eff_lam = np.interp(lam,parameters.throughput['lambda'],parameters.throughput['Eff_tot'])
        d_l=lam/R_lam
        SNR_array[i] = SNR_func(M_src,
                                T_exp=T_exp, 
                                Sky_brightness=Sky_brightness,
                                eff=Eff_lam,
                                binning=binning,
                                wl_AA=lam, 
                                d_lam = d_l, 
                                n_tel = n_tel)[0]
    return SNR_array

def create_limmag_sequence(n_tel_arr = [1,4,20], type = 'per pixel',wl_AA=5500):
    '''
    Creates limiting magnitude sequence for all telescopes, a single telecesope, four telescope and all telescopes per pixel readout for a given SNR limit
    IN:         n_tel_arr - array like. Number of telescopes in array, for each value of n_tel, a limiting magnitude sequence is calculated
                type - 'per pixel' or 'per element'
                wl_AA - wavelength in Angstrom. The limiting magnitude is calculated at this wavelength, taking the as built resolution and throughput at this wavelength into account. 
    OUT:        limmag - dictionary with limiting magnitudes as a function of exposure time for every value of n_tel

    '''

    pixel_factor = parameters.fiber//parameters.mu/parameters.binning[1]
    Xin = [14,17]
    if type == 'per pixel':
       sigma_factor= np.sqrt(pixel_factor)
    elif type == 'per element':
       sigma_factor= 1

    if parameters.R>5000:
        Xin = [12,15]
    if (parameters.R>5000)&parameters.sigma_limit>10:
        Xin = [10,15]
    limmag = {}
    for n in n_tel_arr:
        try:
            l, nc    =find_limiting_mag(parameters.T_exp_vec,SNR_func,SNR_lim=parameters.sigma_limit*sigma_factor,n_tel=n,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,Xin=Xin,N_pack=1,binning=parameters.binning, wl_AA = wl_AA)
        except:
               try:
                   Xin = [Xin[0]-3,Xin[1]-3]
                   l, nc    =find_limiting_mag(parameters.T_exp_vec,SNR_func,SNR_lim=parameters.sigma_limit*sigma_factor,n_tel=n,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,Xin=Xin,N_pack=1,binning=parameters.binning, wl_AA = wl_AA)
               except:
                   Xin = [Xin[0]-3,Xin[1]-3]
                   l, nc    =find_limiting_mag(parameters.T_exp_vec,SNR_func,SNR_lim=parameters.sigma_limit*sigma_factor,n_tel=n,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,Xin=Xin,N_pack=1,binning=parameters.binning, wl_AA = wl_AA)
        limmag[n] = l
    return limmag


