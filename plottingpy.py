
from astropy import constants
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import parameters
from SNR_func import *
from spec_func import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots




def plot_limmag_vs_exp_time(T_exp_vec, 
                             f_SNR=SNR_func, 
                             binning=parameters.binning, 
                             n_tel_arr=[1,4,20], 
                             type = 'per pixel',
                             wl_AA=parameters.wl, 
                             sigma_limit=10):
    """
    Plot42/s limiting magnitude as a function of exposure time for 20, 4, and 1 telescopes using Plotly.
    
    Parameters:
      T_exp_vec   : array-like, exposure time vector.
      f_SNR       : SNR function.
      binning     : pixel binning factor [binX, binY].
      n_tel       : array like: number of telescopes in the array. For every value of n, a curve will be created
      type        : 'per pixel' or 'per element' (default 'per pixel').
      wl_AA       : wavelength in Angstrom.
      sigma_limit : SNR limit (default 10).
    """
    # Generate the data using your existing routine.
    limmag = create_limmag_sequence(n_tel_arr = n_tel_arr, type = type, wl_AA=wl_AA)
    R =  np.interp(wl_AA,parameters.res_wl['lambda'],parameters.res_wl['Res_25um'])
    R = int(R)
    l = [T_exp_vec]
    header = {'calc_type':'Limiting magnitude',
              'n_tel':n_tel_arr,
              'binning':binning,
              'type':type,
              'wl_AA':wl_AA,
              'sigma_limit':sigma_limit,
              'Avg. resolution':R,
              'SNR type': type,
              'Sky brightness': parameters.Sky_brightness_surface_den,
              'col1':'T_exp'}
    for i,n in enumerate(n_tel_arr):
        l.append(limmag[n_tel_arr[i]])
        key = 'col'+str(i+1)
        header[key] = f'limmag for n_tel = {n}'

    output = np.array(l).T


    # Create a Plotly figure.
    fig = go.Figure()
    for n in n_tel_arr:

        fig.add_trace(go.Scatter(
            x=T_exp_vec,
            y=limmag[n],
            mode='lines',
            name=f'{n} telescopes',
            line=dict(dash='dash', width=3)
        ))

    # Define axis ranges.
    # Note: For log axes, Plotly expects the range values in log10 space.
    x_range = [np.log10(50), np.log10(1.1*parameters.T_exp)]
    mmin = np.min([np.min(val) for val in limmag.values()])+4
    mmax = np.max([np.max(val) for val in limmag.values()])+1
    y_range = [mmax, mmin]  # Inverted (high value at top, low value at bottom)
    # Update layout properties.
    if type == 'per pixel':
        title=dict(text=f'Limiting magnitude at {wl_AA} Ang (R = {R}), per pixel', font=dict(size=20))
    elif type == 'per element':
        title=dict(text=f'Limiting magnitude at {wl_AA} Ang (R = {R}), per resolution element', font=dict(size=20))
    else:
        raise ValueError(f'Invalid type: {type}. type must be either "per pixel" or "per element".')
    Tmax = np.max(T_exp_vec)
    logmax = np.floor(np.log10(Tmax))
    tickmax = Tmax/10**logmax
    tickvals=[]
    ticktext=[]    

    for i in range(int(np.log10(Tmax)//1)):
        tickvals.append(10**i)        
        tickvals.append(2*10**i)
        tickvals.append(5*10**i)
        ticktext.append(str(10**i))
        ticktext.append(str(2*10**i))
        ticktext.append(str(5*10**i))

    if tickmax>=0.1:
        tickvals.append(10**logmax)        
        ticktext.append(str(int(10**logmax)))
        tickvals.append(2*10**logmax)        
        ticktext.append(str(int(2*10**logmax)))
    if tickmax>=0.2:
        tickvals.append(5*10**logmax)        
        ticktext.append(str(int(5*10**logmax)))

    fig.update_layout(
        title=title,
        xaxis=dict(
                title=dict(text='Exposure time (s)', font=dict(size=20)),
                type='log',
                range=x_range,
                tickfont=dict(size=20),
                gridcolor='rgba(77,77,77,0.3)',  # approximately (0.3, 0.3, 0.3, 0.3)
                gridwidth=0.5,
                zeroline=False,
                tickmode="array",
                tickvals = tickvals, 
                ticktext = ticktext
        ),
        yaxis=dict(
            title=dict(text=f'{sigma_limit:.0f}-sigma limiting magnitude (AB)', font=dict(size=20)),
            autorange='reversed',  # reverses the y-axis (largest value at bottom)
            range=y_range,
            tickfont=dict(size=20),
            gridcolor='rgba(77,77,77,0.3)',
            gridwidth=0.5,
            zeroline=False
        ),
        # Add a second (twin) y-axis with the same range on the right-hand side.
        yaxis2=dict(
            overlaying='y',
            side='right',
            range=y_range,
            tickfont=dict(size=20),
            showgrid=False
        ),
        legend=dict(font=dict(size=20)),
        width=1000,
        height=800,
        margin=dict(l=80, r=80, t=80, b=80)
    )

    return fig, output, header


def plot_SNR_simspec(lam, simspec, SNR_proj):
    """
    Plots the simulated spectrum and its SNR on a single panel with dual y-axes.
    
    Parameters:
      lam      : array-like
                 Wavelength vector in Angstrom.
      simspec  : array-like
                 Simulated spectrum flux.
      SNR_proj : array-like
                 SNR vector (per pixel).

                 
    The left (primary) y-axis corresponds to the simulated flux (red),
    and the right (secondary) y-axis corresponds to the SNR (blue).
    """
    out = np.array([lam,simspec, SNR_proj]).T
    header = {'calc_type':'SNR/simulated spectrum',
              'n_tel':parameters.n_tel_group,
              'binning':parameters.binning,
              'T_exp': parameters.T_norm,
              'SNR type': 'per pixel',
              'Sky brightness': parameters.Sky_brightness_surface_den, 
              'V_mag_source': parameters.AB_mag_renorm
              ,'spec_type': parameters.spec_type}
    if parameters.spec_type == 'bb':
        header['T_bb'] = parameters.bb_temp
    elif parameters.spec_type == 'WD':
        header['Teff'] = parameters.Teff
        header['logg'] = parameters.logg
    elif parameters.spec_type == 'stellar':
        header['stellar_type'] = parameters.stellar_type

    # Compute y-limits for the simulated spectrum using points where SNR > 3.
    mask = SNR_proj > 3
    if np.any(mask):
        ymax = np.max(simspec[mask]) * 1.2
        ymin = np.min(simspec[mask]) * 0.8
    else:
        ymax = np.max(simspec) * 1.2
        ymin = np.min(simspec) * 0.8

    # Create the figure.
    fig = go.Figure()

    # Add the simulated spectrum trace (primary y-axis) in red.
    fig.add_trace(go.Scatter(
        x=lam,
        y=simspec,
        mode='lines',
        name='Simulated Spectrum',
        line=dict(color='rgba(255,0,0,0.75)')  # red with 75% opacity
    ))

    # Add the SNR trace (secondary y-axis) in blue.
    fig.add_trace(go.Scatter(
        x=lam,
        y=SNR_proj,
        mode='lines',
        name='SNR',
        line=dict(color='blue'),
        yaxis='y2'
    ))

    # Update layout to use dual y-axes with matching colors.
    fig.update_layout(
        title=dict(text='Simulated Spectrum and SNR', font=dict(size=20)),
        xaxis=dict(
            title=dict(text='Wavelength [Å]', font=dict(size=18)),
            tickfont=dict(size=18)
        ),
        yaxis=dict(
            title=dict(text='Simulated Flux [erg cm⁻² s⁻¹ Å⁻¹]', font=dict(size=18, color='red')),
            range=[ymin, ymax],
            tickfont=dict(size=18, color='red')
        ),
        yaxis2=dict(
            title=dict(text='SNR ratio (per pixel)', font=dict(size=18, color='blue')),
            overlaying='y',
            side='right',
            tickfont=dict(size=18, color='blue')
        ),
        margin=dict(l=80, r=80, t=80, b=80)
    )

    # Display the figure if requested.
    return fig, out, header


def plot_spec_per_hour(mag_analyze, SNR = [10,20,50], overhead_sec = 300):
    """
    Plots the number of spectra per hour versus V-band brightness for various SNR targets.
    
    Parameters:
      mag_analyze : array-like, V-band brightness values (in AB mag).
      SNR         : array-like, a curve will be created for each SNR target. 
      overhead_sec: float, overhead time in seconds.
    """
    # Compute spectra per hour for various SNR targets.
    SNR_per_hour = {}
    fig = go.Figure()
    out = [mag_analyze]

    for snr in SNR:
        s_per_hour, _, _ = spec_per_hour(SNR_target=snr, mag_analyze=mag_analyze, overhead_sec=overhead_sec)
        SNR_per_hour[snr] = s_per_hour
        out.append(s_per_hour)
    out = np.array(out).T
    header = {'calc_type':'Spectra per hour',
              'n_tel':parameters.n_tel,
              'binning':parameters.binning,
              'overhead': overhead_sec,
              'SNR type': 'per element',
              'Sky brightness': parameters.Sky_brightness_surface_den, 
              'sigma limit': SNR}

    # Create a Plotly figure.
    for snr in SNR:
    # Add traces for the different SNR targets.
        fig.add_trace(go.Scatter(
            x=mag_analyze,
            y=SNR_per_hour[snr],
            mode='lines',
            name=str(snr)+'σ'
        ))
    

    # Determine the x-axis limits from the provided magnitude data.
    x_min = np.min(mag_analyze)
    x_max = np.max(mag_analyze)

    # Calculate the overhead limit value.
    overhead_limit = parameters.n_tel * 3600 / overhead_sec

    Nmax = np.max(overhead_limit*2)
    logmax = np.floor(np.log10(Nmax))
    tickmax = Nmax/10**logmax
    tickvals=[]
    ticktext=[]    

    for i in range(int(np.log10(Nmax)//1)):
        tickvals.append(10**i)        
        tickvals.append(2*10**i)
        tickvals.append(5*10**i)
        ticktext.append(str(10**i))
        ticktext.append(str(2*10**i))
        ticktext.append(str(5*10**i))

    if tickmax>=0.1:
        tickvals.append(10**logmax)        
        ticktext.append(str(int(10**logmax)))
        tickvals.append(2*10**logmax)        
        ticktext.append(str(int(2*10**logmax)))
    if tickmax>=0.2:
        tickvals.append(5*10**logmax)        
        ticktext.append(str(int(5*10**logmax)))


    # Add a horizontal line (as a scatter trace) for the overhead limit.
    fig.add_trace(go.Scatter(
        x=[x_min, x_max],
        y=[overhead_limit, overhead_limit],
        mode='lines',
        line=dict(color='black', dash='dash'),
        name='overhead limit'
    ))

    # Update the layout with logarithmic y-axis.
    fig.update_layout(
        width=1000,
        height=800,
        title=dict(text='Spectra per Hour', font=dict(size=20)),
        xaxis=dict(
            title=dict(text='V band brightness [AB mag]', font=dict(size=20)),
            tickfont=dict(size=20)
        ),
        yaxis=dict(
            title=dict(text='Spectra/hour', font=dict(size=20)),
            type='log',
            # Specify the range in log10 space: [log10(min), log10(max)]
            range=[np.log10(0.9), np.log10(overhead_limit*2)],
            tickfont=dict(size=20),
            tickmode="array",
            tickvals = tickvals, 
            ticktext = ticktext
        ),
        legend=dict(font=dict(size=20)),
        margin=dict(l=80, r=80, t=80, b=80)
    )

    # Display the figure.
    return fig, out, header


### not implemented in the ETC

def plot_noise_comp(T_exp_vec,mag,f_SNR = SNR_func,binning=parameters.binning,n_tel = parameters.n_tel, R = parameters.R, wl_AA = parameters.wl,plot_only=-1, **kwargs):
    '''
    Plots noise components as a function of exposure time
    IN:         T_exp_vec - exposure time vector,
                mag - source brightness vector,
                f_SNR - SNR function,
                binning - pixel binning factor (2 element list, [binX,binY]),
                n_tel - number of telescopes in array,
                R - spectral resolution,
                wl_AA - wavelength in Angstrom,
                plot_only - noise component to plot (0-source, 1-background, 2-dark current, 3-read noise
                **kwargs - keyword arguments accepted by the matplotlib plot function
    '''
    d_lam=wl/R
    noise_comp=np.empty_like(np.array([0,0,0,0]))
    
    for T,M in zip(T_exp_vec,mag):
        f_SNR(M_src=M,T_exp=T,n_tel=n_tel,mu=mu,fiber = fiber,Sky_brightness=Sky_brightness_surface_den,N_pack=1,binning=binning)
        last=np.array(f_SNR(M_src=M,T_exp=T,n_tel=n_tel,mu=mu,fiber = fiber,Sky_brightness=Sky_brightness_surface_den,N_pack=1,binning=binning, d_lam = d_lam ,wl_AA=wl)[1])
        noise_comp=np.vstack([noise_comp,last])
    noise_comp=np.delete(noise_comp,0,0)
    if plot_only==-1:
        plt.plot(T_exp_vec,noise_comp[:,0],'b',label='source')
        plt.plot(T_exp_vec,noise_comp[:,1],'k',label='background')
        plt.plot(T_exp_vec,noise_comp[:,2],'m',label='dark current')
        plt.plot(T_exp_vec,noise_comp[:,3],'r',label='read noise')
    else:
        plt.plot(T_exp_vec,noise_comp[:,plot_only],**kwargs)
    pass


def plot_noise_comp_resolution(T_exp_vec,mag,f_SNR,binning=[1,1], label = None,n_tel = parameters.n_tel, R = parameters.R, wl_AA = parameters.wl):
    '''
    Plots noise components as a function of spectral resolution
    IN:         T_exp_vec - exposure time vector,
                mag - source brightness,
                f_SNR - SNR function,
                binning - pixel binning factor (2 element list, [binX,binY]),
                label - plot label,
                n_tel - number of telescopes in array,
                R - spectral resolution vector,
                wl_AA - wavelength in Angstrom
     '''
    d_lam=wl/R
    noise_comp=np.empty_like(np.array([0,0,0,0]))
    
    for i in range(len(R)):
        f_SNR(M_src=mag[i],T_exp=T_exp_vec[i],n_tel=n_tel,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,N_pack=1,binning=binning, d_lam = d_lam[i] ,wl_AA=wl_AA[i] )
        last=np.array(f_SNR(M_src=mag[i],T_exp=T_exp_vec[i],n_tel=n_tel,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,N_pack=1,binning=binning, d_lam = d_lam[i] ,wl_AA=wl_AA[i])[1])
        noise_comp=np.vstack([noise_comp,last])
    noise_comp=np.delete(noise_comp,0,0)

    plt.plot(R,noise_comp[:,0],'b',label='source')
    plt.plot(R,noise_comp[:,1],'k',label='background')
    plt.plot(R,noise_comp[:,2],'m',label='dark current')
    plt.plot(R,noise_comp[:,3],'r',label='read noise')

    pass

def plot_noise_comp_wl(T_exp_vec,mag,f_SNR,binning=[1,1], label = None,n_tel = parameters.n_tel, R = parameters.R, wl_AA = parameters.wl):
    '''
    Plots noise components as a function of wavelength
    IN:         T_exp_vec - exposure time vector,
                mag - source brightness,
                f_SNR - SNR function,
                binning - pixel binning factor (2 element list, [binX,binY]),
                label - plot label,
                n_tel - number of telescopes in array,
                R - spectral resolution,
                wl_AA - wavelength vector in Angstrom
    '''
    d_lam=wl/R
    noise_comp=np.empty_like(np.array([0,0,0,0]))
    
    for i in range(len(parameters.res_wl)):
        f_SNR(M_src=mag[i],T_exp=T_exp_vec[i],n_tel=n_tel,mu=mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,N_pack=1,binning=binning, d_lam = d_lam[i] ,wl_AA=wl_AA[i] )
        last=np.array(f_SNR(M_src=mag[i],T_exp=T_exp_vec[i],n_tel=n_tel,mu=parameters.mu,fiber = parameters.fiber,Sky_brightness=parameters.Sky_brightness_surface_den,N_pack=1,binning=binning, d_lam = d_lam[i] ,wl_AA=wl_AA[i])[1])
        noise_comp=np.vstack([noise_comp,last])
    noise_comp=np.delete(noise_comp,0,0)

    plt.plot(wl_AA,noise_comp[:,0],'b',label='source')
    plt.plot(wl_AA,noise_comp[:,1],'k',label='background')
    plt.plot(wl_AA,noise_comp[:,2],'m',label='dark current')
    plt.plot(wl_AA,noise_comp[:,3],'r',label='read noise')
    pass    
