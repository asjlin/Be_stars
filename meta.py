import numpy as np
import matplotlib.pyplot as plt
import oifits
from astropy.io import fits

import extract_data as ext


# location of input files
# IMAGE/SPEC from BeRay
# oifits file
root = r'D:\Andrea\Documents\research\59Cyg_combined\\'

# file names
photfile = 'phot_59Cyg.txt'
intffile = 'oiprep_59_Cyg.2014_2015_tcoh75ms_GHS.XCHAN.cal_fix.AVG.oifits'

# stellar parameters
name = '59Cyg'
temp = '14880'
dist = 434.8 * 3.086e18 # pc * cm/pc

# parameters
powlist = np.array(['2d0', '2d25', '2d5', '2d75', '3d0','3d25'])
denlist = np.array(['3d0em11','4d0em11','5d0em11','6d0em11','7d0em11','8d0em11'])
inclist = np.array(['56', '58', '60', '62', '64', '66', '68', '70',\
'72', '74', '76', '78', '80', '82', '84', '86', '88', '90'])

# minimum error (default 5%)
minerr = 0.05

# -----------------------------------------------------------------------------

# dict of magnitude zero points (erg s-1 cm-2 A-1)
# U B V R I J H K - Vega zeropoints, Bessell et al. 1998
# L - Bessell 1990
# u v b y - R.O. Gray 1998
# W1 W2 W3 W4 - Jarrett et al. 2011 (F_lambda(iso))
zeropt = {\
'U':417.5e-11, 'B':632e-11, 'V':363.1e-11, 'R':217.7e-11, 'I':112.6e-11,\
'J':31.47e-11, 'H':11.38e-11, 'K':3.961e-11, 'L':0.708e-11,\
'u':11.72e-9, 'v':8.66e-9, 'b':5.89e-9, 'y':3.73e-9,\
'W1':0.81787e-11, 'W2':0.24150e-11, 'W3':0.0065151e-11, 'W4':0.0005090e-11}

# dict of effective wavelengths (angstroms)
lam_eff = {\
'U':3600, 'B':4380, 'V':5450, 'R':6410, 'I':7980,\
'J':12200, 'H':16300, 'K':21900, 'L':34500,\
'u':3491, 'v':4111, 'b':4662, 'y':5456,\
'W1':33530, 'W2':46030, 'W3':115610, 'W4':220880,\
'S9W':90000, 'L18W':180000}

# -----------------------------------------------------------------------------
# photometry
# -----------------------------------------------------------------------------

phot = np.genfromtxt(root+photfile, dtype=None, usecols=(0,1,2), skip_footer=8, comments='#')

mags_list = []
emag_list = []
wave_list = []
lfl_list = []
err_list = []
zpt_list = []

# calculate lfl & error
for x in np.arange(0,len(phot)):
	if lam_eff.has_key(phot[x][0]):
		if phot[x][0] in ['S9W','L18W']:
			zpt = 0
			leff = lam_eff[phot[x][0]]
			
			lfl = phot[x][1] * 1e-23 * 2.998e18 / leff
			
			# use larger of min error or photometry error
			err = lfl * minerr
			err_phot = phot[x][2] * 1e-23 * 2.998e18 / leff
			if err_phot > err:
				err = err_phot

		else:
			zpt = zeropt[phot[x][0]]
			leff = lam_eff[phot[x][0]]
			
			lfl = leff * zpt * 10**(-phot[x][1]/2.5)
			
			# use larger of min error or photometry error
			err = lfl * minerr
			err_phot = (leff * zpt * 10**(-(phot[x][1] - phot[x][2])/2.5)) - lfl
			if err_phot > err:
				err = err_phot
			
		# add to lists
		mags_list.append(phot[x][1])
		emag_list.append(phot[x][2])
		wave_list.append(leff)
		lfl_list.append(lfl)
		err_list.append(err)
		zpt_list.append(zpt)
				
	else:
		print phot[x][0], "is not a recognized waveband"

# change types to np.arrays
mags_array = np.array(mags_list)
emag_array = np.array(emag_list)
wave_array = np.array(wave_list)
lfl_array = np.array(lfl_list)
err_array = np.array(err_list)
zpt_array = np.array(zpt_list)


def chisquarephot(root, name, temp, dist, powlist, denlist, inclist, verbose=0):
	
	photmags, photemag, photwave, photlfl, photerr, photzpt =\
	mags_array, emag_array, wave_array, lfl_array, err_array, zpt_array
	
	cspdata = np.zeros((len(powlist),len(denlist)))
	bestinc = np.zeros((len(powlist),len(denlist)))
	
	# loop through densities & power law indexes
	for r in np.arange(0,len(powlist)):
		for c in np.arange(0,len(denlist)):
			
			incchi = []
			
			currbestchisq = 1000000 # arbitrary, just has to be very large
			currbestpow = 100 # arbitrary
			currbestden = 100 # arbitrary
			currbestinc = 100 # arbitrary
			
			for k in np.arange(0,len(inclist)):
				specfile = open(root + 'SPEC.' + name + '_t' + temp + '_w0d0_j0h0_' + powlist[r] + '_' + denlist[c] + '_rd12p0_cont_alp0p0_i' + inclist[k], 'r')
				
				# ignore header
				specfile.readline()

				specwave = []
				specflux = []
				speclfl = []
				
				for line in specfile:
					specwave.append(float(line[6:18]))
					specflux.append(float(line[20:32]))
	
				for x in np.arange(0,len(specflux)):
					speclfl.append(specwave[x]*1e-8 * specflux[x] * (2.998e10 / (specwave[x]*1e-8)**2) / dist**2)
	
				specfile.close()
				
				# lfl at the wavelengths of the photometry (linear interpolation)
				modellfl = np.interp(photwave, specwave, speclfl)
				
				# chi square statistic (magnitude, log space)
				
				modelmags = np.zeros(len(photwave))
				
				for x in np.arange(0,len(modellfl)):
					if photwave[x] in [lam_eff['S9W'], lam_eff['L18W']]:
						modelmags[x] = modellfl[x] * photwave[x] / (1e-23 * 2.998e18)
					else:
						modelmags[x] = -2.5 * np.log10(modellfl[x]/(photwave[x]*photzpt[x]))
				
				# print photmags - modelmags
				# print photemag
				chisq = 1/float(len(photmags)-3) * sum((photmags - modelmags)**2 / (photemag)**2)
								
				# chi square statistic (flux)
				# chisq = 1/float(len(photlfl)-3) * sum(((modellfl - photlfl)**2/photerr**2))			
				
				incchi.append(chisq)
				
				# pick best inclination
				if np.abs(chisq - 1) < np.abs(currbestchisq - 1):
					currbestchisq = chisq
					currbestpow = powlist[r]
					currbestden = denlist[c]
					currbestinc = inclist[k]
			
			if verbose==1:
				print powlist[r], denlist[c], '\n', incchi, '\n'
			
			# enter chi square and best inclination
			cspdata[r,c] = np.min(incchi)
			bestinc[r,c] = np.float(currbestinc)
	
	print "\nchi square grid"
	print np.flipud(cspdata)
	
	# table of densities vs power law indexes
	fig, csp = plt.subplots()
	heatmap = csp.pcolor(np.abs(1-1./cspdata), cmap=plt.cm.hot_r)
	# gray_r : white is low chisq (good), black is high chisq (bad)
	# hot_r : light is low, dark is high
	# rainbow_r : red is low, blue is high
		
	csp.set_xticks(np.arange(cspdata.shape[1])+0.5, minor=False)
	csp.set_yticks(np.arange(cspdata.shape[0])+0.5, minor=False)
	
	csp.set_xticklabels(denlist, fontweight='bold', minor=False)
	csp.set_yticklabels(powlist, fontweight='bold', minor=False)
	
	# print chi-squared statistic and best inclination on each block
	for r in range(cspdata.shape[0]):
		for c in range(cspdata.shape[1]):
			csp.text(c+0.5, r+0.5, '{0:.6f}'.format(cspdata[r,c]), horizontalalignment='center', verticalalignment='bottom')
			csp.text(c+0.5, r+0.5, 'Inc = ' + str(bestinc[r,c]), horizontalalignment='center', verticalalignment='top')
	
	plt.title('Chi-Square Plot - Visible/IR Photometry',fontsize=16,fontweight='bold')
	plt.xlabel('Initial Density (rho_0)', fontsize=12,fontweight='bold')
	plt.ylabel('Power Law Index', fontsize=12, fontweight='bold')
	plt.show()
	
	return photwave
	
def plotSED(root, name, temp, pow, den, inc):
	spec = open(root + 'SPEC.' + name + '_t' + temp + '_w0d0_j0h0_' + pow + '_' + den + '_rd12p0_cont_alp0p0_i' + inc, 'r')
		
	# model spectra
	wavelength = []
	intensity = []

	spec.readline()	# read and toss header

	for line in spec:
		#print line
		wavelength.append(float(line[6:18]))
		intensity.append(float(line[20:32]))

	lfl_spec = np.zeros(len(intensity))

	for i in np.arange(0, len(intensity)):
		lfl_spec[i] = wavelength[i]*1e-8 * intensity[i] * (2.998e10 / (wavelength[i]*1e-8)**2) / dist**2

	spec.close()
	
	# photometry w/errors
	plt.scatter(wave_array, lfl_array)
	plt.errorbar(wave_array, lfl_array, yerr=err_array, fmt='o')
	
	# model
	plt.plot(wavelength, lfl_spec, color = "red")
	

	specfile = open(root + 'SPEC.' + name + '_t' + temp + '_w0d0_j0h0_' + pow + '_' + den + '_rd12p0_cont_alp0p0_i' + inc, 'r')
	
	# ignore header
	specfile.readline()

	specwave = []
	specflux = []
	speclfl = []
	
	for line in specfile:
		specwave.append(float(line[6:18]))
		specflux.append(float(line[20:32]))

	for x in np.arange(0,len(specflux)):
		speclfl.append(specwave[x]*1e-8 * specflux[x] * (2.998e10 / (specwave[x]*1e-8)**2) / dist**2)

	specfile.close()
	
	# lfl at the wavelengths of the photometry (linear interpolation)
	modellfl = np.interp(pwave, specwave, speclfl)
		
	# plt.scatter(pwave, modellfl, color = 'red')
	
	plt.title(name + ' ' + pow + ' ' + den + ' ' + inc)
	plt.xlabel("wavelength (Angstroms)", fontsize=12, fontweight='bold')
	plt.ylabel("lambda F_lambda (ergs/s/cm^2)", fontsize=12, fontweight='bold')
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(1e-12, 1e-2)
	plt.xlim(5e2, 5e5)
	plt.show()

# -----------------------------------------------------------------------------
# interferometry
# -----------------------------------------------------------------------------

intf = oifits.open(root + intffile)

vis2 = intf.vis2
wave = intf.wavelength

# extracting data into arrays
vis2data = np.zeros((len(vis2),1))
vis2sfu = np.zeros(np.shape(vis2data))

u = []
v = []

for i in np.arange(0,len(vis2)):
	# ONLY H-BAND FOR NOW
	vis2data[i] = vis2[i].vis2data[2]
	key = wave.keys()[wave.values().index(vis2[i].wavelength)] # gets key from values
	wavelist = wave[key].eff_wave
	
	u.append(vis2[i].ucoord / wavelist[2]) # 1/rad
	v.append(vis2[i].vcoord / wavelist[2]) # 1/rad
	
	sfu = np.sqrt(u[i]**2 + v[i]**2)
	vis2sfu[i] = sfu

u = np.array(u)
v = np.array(v)
vis2data = np.array(vis2data)


print np.max(lam_eff), np.min(lam_eff)
print np.shape(u), np.shape(v), np.shape(vis2data)

plt.scatter(vis2sfu[:,:],vis2data[:,:])
plt.xlabel('SFU')
plt.ylabel('vis2')
plt.show()


def chisquareintf(root, name, temp, dist, powlist, denlist, inclist, verbose=0):
	
	cspdata = np.zeros((len(powlist),len(denlist)))
	bestinc = np.zeros((len(powlist),len(denlist)))
	
	# loop through densities & power law indexes
	for r in np.arange(0,len(powlist)):
		for c in np.arange(0,len(denlist)):
			
			incchi = []
			
			currbestchisq = 1000000 # arbitrary, just has to be very large
			currbestpow = 100 # arbitrary
			currbestden = 100 # arbitrary
			currbestinc = 100 # arbitrary
			
			for k in np.arange(0,len(inclist)):
				imagefile = open(root + 'IMAGE.' + name + '_t' + temp + '_w0d0_j0h0_' + powlist[r] + '_' + denlist[c] + '_rd12p0_cont_alp0p0_i' + inclist[k], 'r')
				
				templine = imagefile.readline()	
				# get dimensions
				width = int(templine[18:23])
				height = int(templine[26:])

				imagefile.readline()
				imagefile.readline()
				
				image = np.zeros((height,width))
				
				for line in imagefile:
					w = int(line[0:5])
					h = int(line[6:10])
					
					intensity = line[65:75]
					
					image[h-1][w-1] = intensity
				
				imagefile.close()
				
				# norm to 1
				image = image / np.sum(image)
				# chi square statistic
				
				realdata, imagdata = ext.extract_data(image, u, v, im_scale=0.007, cubic=1, nohan=1)
				modelvis2 = realdata**2 + imagdata**2
				chisq = 1/float(np.shape(u)[0]-3) * sum((vis2data - modelvis2)**2 / (0.05*vis2data)**2)
				incchi.append(chisq)
				
				# pick best inclination
				if np.abs(chisq - 1) < np.abs(currbestchisq - 1):
					currbestchisq = chisq
					currbestpow = powlist[r]
					currbestden = denlist[c]
					currbestinc = inclist[k]
			
			if verbose==1:
				print powlist[r], denlist[c], '\n', incchi, '\n'
			
			# enter chi square and best inclination
			cspdata[r,c] = np.min(incchi)
			bestinc[r,c] = np.float(currbestinc)
	
	print "\nchi square grid"
	print np.flipud(cspdata)
	
	# table of densities vs power law indexes
	fig, csp = plt.subplots()
	heatmap = csp.pcolor(np.abs(1-1./cspdata), cmap=plt.cm.hot_r)
	# gray_r : white is low chisq (good), black is high chisq (bad)
	# hot_r : light is low, dark is high
	# rainbow_r : red is low, blue is high
		
	csp.set_xticks(np.arange(cspdata.shape[1])+0.5, minor=False)
	csp.set_yticks(np.arange(cspdata.shape[0])+0.5, minor=False)
	
	csp.set_xticklabels(denlist, fontweight='bold', minor=False)
	csp.set_yticklabels(powlist, fontweight='bold', minor=False)
	
	# print chi-squared statistic and best inclination on each block
	for r in range(cspdata.shape[0]):
		for c in range(cspdata.shape[1]):
			csp.text(c+0.5, r+0.5, '{0:.6f}'.format(cspdata[r,c]), horizontalalignment='center', verticalalignment='bottom')
			csp.text(c+0.5, r+0.5, 'Inc = ' + str(bestinc[r,c]), horizontalalignment='center', verticalalignment='top')
	
	plt.title('Chi-Square Plot - H-band Interferometry', fontsize=16, fontweight='bold')
	plt.xlabel('Initial Density (rho_0)', fontsize=12, fontweight='bold')
	plt.ylabel('Power Law Index', fontsize=12, fontweight='bold')
	plt.show()

# -----------------------------------------------------------------------------

# pwave = chisquarephot(root, name, temp, dist, powlist, denlist, inclist, verbose=0)

# plotSED(root, name, temp, '2d25', '4d0em11', '78')
# plotSED(root, name, temp, '2d75', '8d0em11', '76')

chisquareintf(root,name,temp,dist,powlist,denlist,inclist,verbose=1)
'''
imagefile = open(root + 'IMAGE.' + name + '_t' + temp + '_w0d0_j0h0_2d25_3d0em11_rd12p0_cont_alp0p0_i74', 'r')

templine = imagefile.readline()	
# get dimensions
width = int(templine[18:23])
height = int(templine[26:])

imagefile.readline()
imagefile.readline()

image = np.zeros((height,width))

for line in imagefile:
	w = int(line[0:5])
	h = int(line[6:10])
	
	intensity = line[65:75]
	
	image[h-1][w-1] = intensity

imagefile.close()

plt.imshow(image, cmap='hot')
'''
plt.show()