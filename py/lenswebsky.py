import healpy

import flipper.flipperDict as flipperDict

import matplotlib.pyplot as plt

from pixell import  curvedsky, enmap, utils, lensing, sharp
import numpy as np
import get_cmb_powerspectra


p = flipperDict.flipperDict()
p.read_from_file('lenswebsky.dict')

do_all = True
if do_all:


    kappa_map_hp = healpy.read_map(p['websky_dir'] + p['kappa_map_name'])




    kappa_alms = healpy.sphtfunc.map2alm(kappa_map_hp)

    shape, wcs = enmap.fullsky_geometry(p['PIX_SIZE_ARCMIN']*utils.arcmin)

    # phi_alms =   healpy.sphtfunc.almxfl(kappa_alms, 1. / (ells * (ells + 1) / 2.))

    # phi_map = 

    kappa_map_car = enmap.enmap(np.zeros(shape), wcs)
    curvedsky.alm2map(kappa_alms, kappa_map_car)

    lmax = healpy.Alm.getlmax(len(kappa_alms))
    ainfo = sharp.alm_info(lmax)

    ellvals = np.arange(lmax)

    #factor to change kappa into phi
    func = 2. / (ellvals * (ellvals + 1.))
    func[0] = 0

    phi_alm = healpy.sphtfunc.almxfl(kappa_alms, func)
    
    cmb_powers = get_cmb_powerspectra.websky_cmb_spectra()

    cmb_alm, cmb_ainfo = curvedsky.rand_alm(cmb_powers, lmax = lmax, seed = 0, return_ainfo = True)

    unlensed_map, lensed_map = lensing.lens_map_curved((3,) + shape, wcs, phi_alm, cmb_alm, ainfo, output = 'ul')


    cmb_alm_lensed = curvedsky.map2alm(lensed_map, ainfo = cmb_ainfo)
    

    healpy.fitsfunc.write_alm(p['outdir'] + 'lensed_alm.fits',
                              np.complex64(cmb_alm_lensed), overwrite = True)

    healpy.fitsfunc.write_alm(p['outdir'] + 'unlensed_alm.fits',
                              np.complex64(cmb_alm), overwrite = True)

    

degs = 3
crange = [-300, 300]

unl_cutout = unlensed_map[0, shape[0] / 2 - int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN'])) :
                          shape[0] / 2 + int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN'])),
                          shape[1] / 2 - int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN'])):
                          shape[1] / 2 + int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN']))]

len_cutout = lensed_map[0, shape[0] / 2 - int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN'])) :
                         shape[0] / 2 + int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN'])),
                         shape[1] / 2 - int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN'])):
                         shape[1] / 2 + int(np.rint(degs * 60 / p['PIX_SIZE_ARCMIN']))]




plt.figure('unl_map', figsize = (10,10))
plt.clf()
plt.imshow(np.flipud(unl_cutout), clim = crange)
plt.savefig('../plot/unl_cutout.png', dpi = 500)
plt.colorbar()

plt.figure('len_map', figsize = (10,10))
plt.clf()
plt.imshow(np.flipud(len_cutout), clim = crange)
plt.savefig('../plot/len_cutout.png', dpi = 500)
plt.colorbar()

