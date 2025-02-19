import parameters
from SNR_func import *
from spec_func import *
from plottingpy import *

def run_ETC(spec = None):
    """
    Runs the exposure time calculator logic based on the current parameters
    and returns a Plotly figure.
    """
    if parameters.calc_type == 'limmag':
        # Assume plot_limmag_vs_exp_time returns a Plotly figure when show=False.
        fig, out, header = plot_limmag_vs_exp_time(parameters.T_exp_vec, wl_AA=parameters.wl, n_tel_arr=parameters.n_tel_arr, type=parameters.Type, sigma_limit = parameters.sigma_limit)
    elif parameters.calc_type == 'SNR':
        if spec is None:
            spec = generate_spec(parameters.spec_type, stellar_type=parameters.stellar_type, Teff = parameters.Teff, log_g = parameters.logg)
        else: 
            spec = prep_user_spec(spec)


        spec_sim, SNR_proj, total_counts = SNR_sequence(parameters.lam,
                                         spec, parameters.AB_mag_renorm,
                                         T_exp = parameters.T_norm,
                                         binning=parameters.binning, 
                                         n_tel=parameters.n_tel_group, 
                                         Sky_brightness= parameters.Sky_brightness_surface_den
                                         , Type = 'per pixel')
        parameters.counts = total_counts
        fig, out, header = plot_SNR_simspec(parameters.lam, spec_sim, SNR_proj)
    elif parameters.calc_type == 'spec_per_hour':
        fig, out, header = plot_spec_per_hour(parameters.mag_analyze, SNR=parameters.SNR_arr)
    else:
        raise ValueError("calc_type not recognized")
    parameters.output = out
    parameters.header = header

    return fig, out, header

def save_output(output, header, path):
    """
    Saves the output of the ETC to a text file. writes the header line by line and then writes the array line by line
    """
    with open(path, "w") as f:
        f.write('#DeepSpec/MAST ETC output generator \n')
        f.write('#By Ido Irani \n')
        for key in header.keys():
            f.write('#'+ key + ':' + str(header[key])+'\n')
        for i in range(np.shape(output)[0]):
            line = output[i,:]
            string = f"{line[0]:.4g}"
            for i in range(1,len(line)):
                string += f",{line[i]:.4g}"
            string += "\n"
            f.write(string)
    pass   



if __name__ == '__main__':
    # For debugging outside Dash.
    fig, output,header = run_ETC()
    fig.show()