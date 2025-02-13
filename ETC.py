import parameters
from SNR_func import *
from spec_func import *
from plottingpy import *

def run_ETC():
    """
    Runs the exposure time calculator logic based on the current parameters
    and returns a Plotly figure.
    """
    if parameters.calc_type == 'limmag':
        # Assume plot_limmag_vs_exp_time returns a Plotly figure when show=False.
        fig = plot_limmag_vs_exp_time(parameters.T_exp_vec, wl_AA=parameters.wl, n_tel_arr=parameters.n_tel_arr, type=parameters.Type, sigma_limit = parameters.sigma_limit)
    elif parameters.calc_type == 'SNR':
        spec = generate_spec(parameters.spec_type, stellar_type=parameters.stellar_type, Teff = parameters.Teff, log_g = parameters.logg)
        spec_sim, SNR_proj = SNR_sequence(parameters.lam,
                                         spec, parameters.AB_mag_renorm,
                                         T_exp = parameters.T_norm,
                                         binning=parameters.binning, 
                                         n_tel=parameters.n_tel_group, 
                                         Sky_brightness= parameters.Sky_brightness_surface_den)

        fig = plot_SNR_simspec(parameters.lam, spec_sim, SNR_proj)
    elif parameters.calc_type == 'spec_per_hour':
        fig = plot_spec_per_hour(parameters.mag_analyze, SNR=parameters.SNR_arr)
    else:
        raise ValueError("calc_type not recognized")
    return fig

if __name__ == '__main__':
    # For debugging outside Dash.
    fig = run_ETC()
    fig.show()