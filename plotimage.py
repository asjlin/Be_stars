import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plotimage(filename):
	file = open(filename, 'r')
	
	templine = file.readline()	# get dimensions
	width = int(templine[18:23])
	height = int(templine[26:])
	
	templine = file.readline()	# get num of wavelengths
	numwaves = int(templine[26:])
	
	# print "dimensions", width, "by", height
	# print numwaves, "wavelengths"
		
	file.readline()	# ignore units line
	
	# cols 3&4 look like position data
	# what does column 5 do?
	
	templine = file.readline()	# get column headings (for labels later)
	headings = []
		
	for i in np.arange(0,numwaves):
		begin = (10*i)+45
		end = (10*i)+55
		headings.append(float(templine[begin:end]))
	
	# print headings
	
	# create [numwaves] arrays of zeros
	meta = []
	
	for i in np.arange(0,numwaves+1):
		meta.append(np.zeros((height,width)))
	
	# input all the table values into meta array
	for line in file:
		w = int(line[0:5])
		h = int(line[6:10])
		
		intensity = []
		
		for i in np.arange(0,numwaves):
			begin = (10*i)+45
			end = (10*i)+55
			intensity.append(float(line[begin:end]))
		
		for i in np.arange(0,numwaves):
			meta[i][h-1][w-1] = intensity[i]
	
	file.close()
	
	# plot all the wavelengths
	
	prows = int(np.floor(np.sqrt(numwaves)))
	pcols = int(np.ceil(numwaves/prows))
	
	f, fig = plt.subplots(prows, pcols)
	
	plotnum = 0
	
	for i in np.arange(0,prows):
		for j in np.arange(0,pcols):
			fig[i,j].imshow(meta[plotnum], extent=[0,width,0,height], cmap = cm.hot)
			fig[i,j].set_title(str(headings[plotnum]) + " microns")
			plotnum += 1
			fig[i,j].xaxis.set_visible(False)
			fig[i,j].yaxis.set_visible(False)
	plt.show()
	
# need to input full path as raw string
plotimage(r'D:\research\59Cyg_combined\IMAGE.59Cyg_t14880_w0d0_j0h0_2d25_3d0em11_rd12p0_cont_alp0p0_i74')