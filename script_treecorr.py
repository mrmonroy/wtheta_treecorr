import numpy as np
import os
import healpy as hp
import healpix_util as hu
import sys
import time

return_var = False
return_corr = False

nside = 512

nthetabins = 14
thetamin = 9.95267926/60.
thetamax = 250./60.
bin_slop = 0.01
bin_type = None #None or Linear

sysweights = True
fracweights = True

jointmask = #Vector of the observed pixels. True/False means masked/observed pixel

N1 = # map1 data (evaluated where there is data, i.e. with no hp.UNSEEN)
fracdet1 = # fractional coverage of the observed pixels in map1
scale1 = # rescaling factor applied to map1 observed pixels. This should be 1./fracdet1, such that the weighted mean of N1, N1_bar, is N1_bar = sum((N1_i/fracdet1_i) * fracdet1_i)/ sum(fracdet1_i)
sw1 = # systematic weight map applied to map1 observed pixels 

N2 = # map2 data (evaluated where there is data, i.e. with no hp.UNSEEN)
fracdet2 = # fractional coverage of the observed pixels in map2
scale2 = # rescaling factor applied to map2 observed pixels. This should be 1./fracdet2, such that the weighted mean of N2, N2_bar, is N2_bar = sum((N2_i/fracdet2_i) * fracdet2_i)/ sum(fracdet2_i)
sw2 = # systematic weight map applied to map2 observed pixels

kk = treecorr.KKCorrelation(min_sep=thetamin, max_sep=thetamax, nbins=nthetabins, sep_units='degrees', bin_slop=bin_slop, bin_type=bin_type)

print ('corr2pt using bin_slop = {0}'.format(kk.bin_slop))
pixindex = np.arange(hp.nside2npix(nside))
pixels = pixindex[jointmask == False]
hpix = hu.HealPix('ring',nside)
ra,dec = hpix.pix2eq(pixels)

assert len(N1) == len(N2)

if systweights:
	N1 = N1*sw1

if scale1 is not None:
	N1 = N1*scale1

if delta_input == False:
	N1_bar = np.average(N1, weights = fracdet1)
	delta1 = (N1-N1_bar)/N1_bar
else:
	delta1 = N1

if fracweights:
	fracweights1 = fracdet1
else: 
	fracweights1 = np.ones(len(delta1))

start = time.time()

cat1 = treecorr.Catalog(ra=ra,dec=dec,ra_units='degrees',dec_units='degrees', k=delta1, w=fracweights1)


if (N1 == N2).all() == True:
	print ('corr2pt: auto')
	kk.process(cat1,  num_threads=num_threads)
else:
	if systweeights:
		N2 = N2*sw2
	
	if scale2 is not None:
		N2 = N2*scale2
	
	if delta_input == False:
		N2_bar = np.average(N2, weights = fracdet2)
		delta2 = (N2-N2_bar)/N2_bar
	else:
		delta2 = N2


	if fracweights:
		fracweights2 = fracdet2
	else:
		fracweights2 = np.ones(len(delta2))

	cat2 = treecorr.Catalog(ra=ra,dec=dec,ra_units='degrees',dec_units='degrees', k=delta2, w=fracweights2)
	print ('corr2pt: cross')
	kk.process(cat1,cat2, num_threads=num_threads)
	
finish = time.time()
print ("corr2pt took {0}s".format(finish - start))

if return_var == True:
	return np.exp(kk.logr), kk.xi, kk.varxi
if return_corr == True:
	return np.exp(kk.logr), kk.xi, kk
else:
	return np.exp(kk.logr), kk.xi







