


import camb, numpy as np
import pdb

def websky_cosmology():
#stolen from ~https://mocks.cita.utoronto.ca/data/websky/v0.0/cosmology.py
    output = {}
    output['omega_b'] = 0.049
    output['omega_c'] = 0.261
    output['omega_m'] = output['omega_b'] + output['omega_c']
    output['h']      = 0.68
    output['n_s']     = 0.965
    # sigma8 = 0.81

    output['A_s'] = 2.022e-9 #note this gets me s8 = 0.81027349, pretty close to the specified 0.81

    return output

def websky_cmb_spectra():


    websky_params = websky_cosmology()

    #Set up a new set of parameters for CAMB
    pars = camb.CAMBparams()
    #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
    pars.set_cosmology(
        H0 = websky_params['h'] * 100,
        ombh2 = websky_params['omega_b'] * websky_params['h']**2,
        omch2 = websky_params['omega_c'] * websky_params['h']**2,
        mnu = 0.,
        omk = 0,
        tau = 0.055)

    pars.InitPower.set_params(
        websky_params['A_s'],
        ns = websky_params['n_s'],
        r=0)

    pars.set_for_lmax(10000, lens_potential_accuracy=0)


    results = camb.get_results(pars)

    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
    for name in powers: print(name)

    to_use = powers['unlensed_scalar']

    input_power = np.zeros((3, 3, to_use.shape[0]))

    ells = np.arange(to_use.shape[0])
    camb_factor = np.append([0], 2. * np.pi / (ells[1:] * (ells[1:] + 1) ))


    input_power[0,0,:] = to_use[:,0] * camb_factor #TT
    input_power[1,1,:] = to_use[:,1] * camb_factor #EE
    input_power[1,0,:] = to_use[:,3] * camb_factor #TE
    input_power[0,1,:] = to_use[:,3] * camb_factor #EE
    
    return input_power

