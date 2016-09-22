import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import oifits

file = oifits.open(r'D:\Andrea\Documents\research\phiper_allepochs_SPLIT_middle.oifits')

vis2 = file.vis2
t3 = file.t3
wave = file.wavelength

# extracting data into arrays
vis2data = np.zeros((len(vis2),len(vis2[0].vis2data)))
t3amp = np.zeros((len(t3),len(t3[0].t3amp)))
t3phi = np.zeros((len(t3),len(t3[0].t3phi)))

vis2sfu = np.zeros(np.shape(vis2data))

for i in np.arange(0,len(vis2)):
	vis2data[i] = vis2[i].vis2data
	
	key = wave.keys()[wave.values().index(vis2[i].wavelength)] # gets key from values
	wavelist = wave[key].eff_wave
	
	u = vis2[i].ucoord / wavelist # 1/rad
	v = vis2[i].vcoord / wavelist # 1/rad
	
	sfu = np.sqrt(np.square(u) + np.square(v))
	vis2sfu[i] = sfu

print np.shape(vis2sfu)

# haven't done closure phases yet
'''
for i in np.arange(0,len(t3)):
	t3amp[i] = t3[i].t3amp
	t3phi[i] = t3[i].t3phi
'''
# quick plot to confirm vis2's are there
plt.scatter(vis2sfu[:,:],vis2data[:,:])
plt.xlabel('SFU')
plt.ylabel('vis2')
plt.show()

'''
print vis2[0].vis2data
print vis2[0].vis2err
'''
target = file.target
print target[0]

wave = file.wavelength
print "effective wavelengths\n", wave["H_PRISM"].eff_wave