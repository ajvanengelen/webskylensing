import matplotlib
matplotlib.use('Agg')
import get_cmb_powerspectra
import numpy as np, matplotlib.pyplot as plt





specs = get_cmb_powerspectra.websky_cmb_spectra(return_tensors = True)
# import pickle
# pickle.save(specs['tensor'], open("websky_BB_from_tensors_r_eq_1.py"))

np.save("/global/cscratch1/sd/engelen/websky_BB_from_tensors_r_eq_1.npz", specs['tensor'])

