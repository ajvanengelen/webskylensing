import healpy

import flipper.flipperDict as flipperDict

import matplotlib.pyplot as plt

from pixell import  curvedsky, enmap, utils, lensing, sharp
import numpy as np
import get_cmb_powerspectra

plotloc = '/project/projectdirs/sobs/www/users/engelen/'

p = flipperDict.flipperDict()
p.read_from_file('lenswebsky.dict')

do_all = True
if do_all:

    print 'grabbing kappa map from', \
        healpy.read_map(p['websky_dir'] + p['kappa_map_name'])
    
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
    
    cmb_powers = get_cmb_powerspectra.websky_cmb_spectra()['unlensed_scalar']

    cmb_alm, cmb_ainfo = curvedsky.rand_alm(cmb_powers, lmax = lmax, seed = 1, return_ainfo = True)

    unlensed_map, lensed_map = lensing.lens_map_curved((3,) + shape, wcs, phi_alm, cmb_alm, ainfo, output = 'ul')


    cmb_alm_lensed = curvedsky.map2alm(lensed_map, ainfo = cmb_ainfo)
    

    healpy.fitsfunc.write_alm(p['outdir'] + 'lensed_alm.fits',
                              np.complex128(cmb_alm_lensed), overwrite = True)

    healpy.fitsfunc.write_alm(p['outdir'] + 'unlensed_alm.fits',
                              np.complex128(cmb_alm), overwrite = True)

    

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
plt.savefig(plotloc + 'unl_cutout.png', dpi = 500)
plt.colorbar()

plt.figure('len_map', figsize = (10,10))
plt.clf()
plt.imshow(np.flipud(len_cutout), clim = crange)
plt.savefig(plotloc + 'len_cutout.png', dpi = 500)
plt.colorbar()

cmb_cl_unlensed = healpy.alm2cl(cmb_alm)
cmb_cl_lensed = healpy.alm2cl(cmb_alm_lensed)

plt.figure('powerspec', figsize = (10,10))
plt.clf()
plt.subplot(2, 1, 1)

theory_powers = get_cmb_powerspectra.websky_cmb_spectra()

# for si, spec in enumerate([cmb_cl_lensed, cmb_cl_unlensed]):
for ii in range(2):
    line, = plt.plot(cmb_cl_lensed[ii] / cmb_cl_unlensed[ii], 
                        marker = 'o', markersize = .2, linestyle = 'None')

    plt.plot(theory_powers['lensed_scalar'][ii, ii, :] / theory_powers['unlensed_scalar'][ii, ii, :],
                 color = line.get_color(), label = ['TT', 'EE'][ii])

    plt.axhline(1., linestyle = 'dashed', color = '.5')


plt.ylim([.8, 1.4])
plt.xlim([0,5000])
plt.legend()
plt.ylabel('lensed power / unlensed power')
plt.xlabel('l')

plt.subplot(2, 1, 2)
ii = 2
line, = plt.plot(cmb_cl_lensed[ii] , 
                        marker = 'o', markersize = .2, linestyle = 'None', color = 'red')

plt.plot(theory_powers['lensed_scalar'][ii, ii, :] ,
         color = line.get_color(), label = 'BB')


plt.xlim([0,5000])
plt.legend()
plt.ylabel('$C_l^{BB}$ (uK$^2$)')
plt.xlabel('l')

plt.ylim([0, (6. * (np.pi / 180. / 60.))**2])

plt.savefig(plotloc + 'lens_power_check.png', dpi = 300)

