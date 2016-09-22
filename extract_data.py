import numpy as np
import scipy.special as sp
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits

# This routine extracts the fourier amps and phases of a given image.
# East is to the LEFT, North UP
'''
testimage = np.zeros((100,100))
for i in np.arange(0,np.shape(testimage)[0]):
	for j in np.arange(0,np.shape(testimage)[1]):
		# cos (radial)
		# testimage[i,j] = 100 * np.cos(2*np.pi*np.sqrt((i-np.shape(testimage)[0]/2.)**2 + (j-np.shape(testimage)[1]/2.)**2)/scale)
		# bessel fxn j0
		# testimage[i,j] = 100 * sp.j0(2*np.pi*np.sqrt((i-np.shape(testimage)[0]/2.)**2 + (j-np.shape(testimage)[1]/2.)**2)/scale)
		# step fxn
		if np.sqrt((i-np.shape(testimage)[0]/2.+0.5)**2 + (j-np.shape(testimage)[1]/2.+0.5)**2) < 5:
			testimage[i,j] = 10.
		# if np.sqrt((i-np.shape(testimage)[0]/2.+0.5)**2 + (j-np.shape(testimage)[1]/2.+0.5)**2) < 1:
			# testimage[i,j] = 1000

testimage = testimage / np.sum(testimage)
plt.imshow(testimage)
plt.show()

testu = np.array(np.arange(-2e8,2e8,1e7))
testv = np.array(np.arange(-2e8,2e8,1e7))

testu = np.random.randint(-2e8,2e8+1,20)
testv = np.random.randint(-2e8,2e8+1,20)

testu = np.array([float(i) for i in testu])
testv = np.array([float(j) for j in testv])
'''
# extract H band, norm to 1
# calculate vis2 data at the us, vs in the oifits
# chisq against model

def extract_data(image1, u, v, im_scale=1, cubic=0, help=0, nohan=0, nopad=0, get=0):
	#vis, phases,, u, v, \
	#fft_scale=scale, amp=a, theta=t, noprint=0,\
	#cvis_real=vreal, cvis_imag=vimag):
	if cubic == 0:
		print "cubic=0 is not an available option."
	
	if help == 1:
		print 'image, vis, phases (degs), scale=im_scale (mas/pixel), u=u, v=v, /cubic, /nohan, /nopad'
		
	ydim, xdim = np.shape(image1)

	if nohan == 0:
		print "nohan=0 is not an available option."
		return

	im_scale_rad = im_scale / 206265000. # mas/pix * rad/mas

	if nopad == 1:
		factor = 1
	else:
		factor = 16
	
	# scale of u,v
	scale = 1.0 / (xdim*factor*im_scale_rad)
	
	if get == 1:
		print "get=1 is not an available option."
		return
		
	# calculate u/v from ucoord/vcoord
	'''
	sfus = np.sqrt(np.square(u) + np.square(v))
	print sfus
	maxsfu = np.amax(sfus/1e5)
	minsfu = np.amin([x/1e5 for x in sfus if x > 0])
	'''
	image = np.zeros((ydim*factor, xdim*factor))
	image[0:ydim,0:xdim] = image1 # would multiply by Hanning fxn here if nohan=0

	
	# inverse 2D FFT
	# when shifting, use python roll array
	image_fft = float((xdim*factor)*(ydim*factor)) * np.fft.ifft2(np.roll(np.roll(image,-xdim/2, axis=1),-ydim/2, axis=0))
	# image_fft = np.roll(np.roll(image_fft,xdim*factor/2+1, axis=1), ydim*factor/2+1, axis=0)
	
	# separate into real and imag
	realimage = image_fft.real
	imagimage = image_fft.imag
	
	# print np.fft.fftfreq(xdim*factor)

	'''
	plt.figure(1)
	plt.imshow(image1)
	plt.figure(2)
	plt.imshow(np.roll(np.roll(image,-xdim/2+xdim*factor/2, axis=1),-ydim/2+ydim*factor/2, axis=0), cmap=cm.gray) 
	plt.figure(3)
	# plt.imshow(np.sqrt(realimage**2 + imagimage**2), cmap=cm.gray)
	plt.imshow(np.roll(np.roll(realimage, xdim*factor/2, axis=1),ydim*factor/2, axis=0))
	plt.figure(4)
	plt.imshow(np.roll(np.roll(imagimage, xdim*factor/2, axis=1),ydim*factor/2, axis=0))
	plt.show()
	'''
	
	# print len(testu), len(testv)
	# print np.size(realimage)
	# interpolate
	
	realroll = np.roll(np.roll(realimage, xdim*factor/2, axis=1),ydim*factor/2, axis=0)
	imagroll = np.roll(np.roll(imagimage, xdim*factor/2, axis=1),ydim*factor/2, axis=0)
	
	realfxn = interp.interp2d(np.arange(-xdim*factor/2, xdim*factor/2), np.arange(-ydim*factor/2, ydim*factor/2), realroll, kind='cubic')
	imagfxn = interp.interp2d(np.arange(-xdim*factor/2, xdim*factor/2), np.arange(-ydim*factor/2, ydim*factor/2), imagroll, kind='cubic')
	
	realdata = []
	imagdata = []
	
	# this loop is SLOW
	for x in np.arange(0,np.shape(u)[0]):
		realval = realfxn((-u[x]/scale), (v[x]/scale))
		imagval = imagfxn((-u[x]/scale), (v[x]/scale))
		
		realdata.append(realval)
		imagdata.append(imagval)
		
	realdata = np.array(realdata)
	imagdata = np.array(imagdata)
	'''
	realdata = realfxn((u/scale), (v/scale))
	imagdata = imagfxn((u/scale), (v/scale))
	'''	
	'''
	# plot vis2 against sqrt(u^2+v^2) as a check
	sfus = []
	vis2 = []
	
	# for i in np.arange(0,len(u)):
		# for j in np.arange(0,len(v)):
			# sfus.append(np.sqrt(u[i]**2 + v[j]**2))
			# vis2.append(np.sqrt(realdata[j][i]**2 + imagdata[j][i]**2))
	
	sfus.append(np.sqrt(u**2 + v**2))
	vis2.append(realdata**2 + imagdata**2)
	
	print np.shape(sfus), np.shape(vis2)
	plt.scatter(sfus, vis2)
	plt.show()
	
	
	# convert to polar (np.angle takes quadrants into account)
	amp = np.abs(np.conj(image_fft))
	phases = np.angle(np.conj(image_fft))
	
	print np.min(amp), np.max(amp)
	print np.min(phases), np.max(phases)
	'''
	return realdata, imagdata
	
	# may have to convert from ucoord/vcoord to u/v, divide coord by wavelength
'''
extract_data(testimage, testu, testv, im_scale=0.01, cubic=1, nohan=1)
'''
# delta u = 1/theta (total size of image, radians)
# plot sqrt (u^2 + v^2) vs abs(FFT)