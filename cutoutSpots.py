#!/usr/bin/env python3
import argparse, sys, json, datetime, os, re
import matplotlib.pyplot
import matplotlib
import numpy, math
from astropy.io import fits
import astropy, photutils, scipy
from QuickLook import QL
from photutils.segmentation import deblend_sources
from photutils.segmentation import SourceCatalog
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources

def boostImageData(imageData):
    """ Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255 """
    hi = 99
    lo = 20
    data = imageData
    max = data.max()
    dataArray = data.flatten()
    pHi = numpy.percentile(dataArray, hi)
    pLo = numpy.percentile(dataArray, lo)
    range = pHi - pLo
    scale = range/255
    data = numpy.clip(data, pLo, pHi)
    data-= pLo
    data/=scale
    return data


class exposure:
	def __init__(self, filename, doplot = False):
		self.filename = filename
		self.headers = {}
		self.sources = []
		# Load the left data
		hdulist = fits.open(filename)  # open the FITS file
		hdr = hdulist[0].header  # the primary HDU header
		card = hdulist[0]
		for key in card.header.keys():
			self.headers[key] = card.header[key]
		hdulist.close(output_verify='ignore')

		QLook = QL(biascorr=True, plot=doplot, verbose=False)
		QLook.combineCCDs(filename)
		self.data = QLook.newarray
		self.arm = QLook.arm
		self.shutter = self.headers['HART']
		self.VPH = self.headers['VPH']
		self.focusA = self.headers['FOCUSMTA']
		self.focusB = self.headers['FOCUSMTB']
		self.RUN = self.headers['IRAFNAME']
		self.fpmode = self.headers['FPMODE']
		self.dimensions = numpy.shape(self.data)
		self.cutoutActive = False
		self.fullData = self.data
		self.thresh = -1
		#print("Loaded image %s.  %s arm"%(self.filename, self.arm))

	def summary(self):
		return "[%s] %s Arm: %s  Shutter: %s, Grating: %s Focus A: %dum Focus B: %dum"%(self.RUN, self.filename, self.arm, self.shutter, self.VPH, self.focusA, self.focusB)

	def subtractBackground(self):
		### subtract background
		bkg_estimator = MedianBackground()
		if self.arm == 'red':
			self.bkg = Background2D(self.data, (500, 20), filter_size=(3, 3),
							bkg_estimator=bkg_estimator)
		else:
			self.bkg = Background2D(self.data, (10,10), filter_size=(3, 3),
						bkg_estimator=bkg_estimator)
		self.data -= self.bkg.background 

	def reset(self):
		self.data = self.fullData
		self.dimensions = numpy.shape(self.data)
		self.sources = []

	def cutout(self, cutoutDimensions):
		halfwidth = cutoutDimensions[2] >> 1
		halfheight = cutoutDimensions[3] >> 1
		x1 = cutoutDimensions[0] - halfwidth
		x2 = cutoutDimensions[0] + halfwidth
		y1 = cutoutDimensions[1] - halfheight
		y2 = cutoutDimensions[1] + halfheight
		self.data = self.data[y1:y2, x1:x2]
		self.dimensions = numpy.shape(self.data)
		
	def segmentSources(self):
		### segment
		thresh  = self.thresh
		if thresh < 0:
			if self.arm == 'red':
				thresh = 5
			else:
				thresh = 20
		print("\tusing threshold value of:", thresh)
		data = self.data
		threshold = thresh * self.bkg.background_rms
		segment_map = detect_sources(data, threshold, npixels=10)
		cat = SourceCatalog(data, segment_map)
		tbl = cat.to_table()

		for t in tbl:
			source = { "x" :  round(t['xcentroid'], 3), "y" : round(t['ycentroid'], 3), "peak" : round(t['max_value'], 3), 'flux' : t['segment_flux']}
			self.sources.append(source)
		print("\t%d sources found."%len(self.sources))

	def plotSources(self):
		self.figure = matplotlib.pyplot.figure()
		xValues = [ s['x'] for s in self.sources ]
		yValues = [ s['y'] for s in self.sources ]
		cValues = [ s['peak'] for s in self.sources ]
		matplotlib.pyplot.scatter(xValues, yValues, c=cValues, cmap="hot")
		matplotlib.pyplot.draw()
		matplotlib.pyplot.pause(0.1)

		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads the left and right Hartmann images and looks for spots. This version cuts out regions of the detector.')
	parser.add_argument('-l', '--left', type=str, help='left shutter in Hartmann exposure. (right shutter open)')
	parser.add_argument('-r', '--right', type=str, help='right shutter in Hartmann exposure. (left shutter open)')
	parser.add_argument("--plot", action="store_true", help='Hold the plots on the display.')
	parser.add_argument("--plotsources", action="store_true", help='Plots the found sources as a schematic.')
	parser.add_argument("--version", action="store_true", help='Show the versions of the libraries used. Then exits.')
	parser.add_argument("--cutouts", help="A JSON file defining the cutouts (aka sub-regions) of the image to use.")
	parser.add_argument("--thresh", default=-1, type=int,  help="Threshhold for segmentation (default 5 red and 20 blue).")
	args = parser.parse_args()
	(execPath, exec) = os.path.split(__file__)
	
	if args.version: 
		print("Python version:", sys.version)
		print("Astropy version:", astropy.__version__)
		print("Photutils version:", photutils.__version__)
		print("Scipy version:", scipy.__version__)
		print("Matplotlib version:", matplotlib.__version__)
		sys.exit()

	if (args.left is None) or (args.right is None):
		print("You must specify a left and right input image.")
		sys.exit()

	if (args.cutouts is None):
		cutoutsFilename = os.path.join(execPath, "cutouts.json")
	else:
		cutoutsFilename = args.cutouts

	# Load the cutouts
	cutoutsDefFile = open(cutoutsFilename, "rt")
	cutoutsDef = json.load(cutoutsDefFile)
	cutoutsDefFile.close()
	
	left_exposure = exposure(args.left, doplot = args.plot)
	if args.thresh!=-1: left_exposure.thresh = args.thresh
	print(left_exposure.summary())
	right_exposure = exposure(args.right, doplot = args.plot)
	if args.thresh!=-1: right_exposure.thresh = args.thresh
	print(right_exposure.summary())

	notFound = True
	for c in cutoutsDef:
		if c['arm'] == left_exposure.arm and c['VPH'] == left_exposure.VPH: 
			print("Found cutout def")
			notFound = False
			cutouts = c['cutouts']
	if notFound:
		print("Could not find cutout info for arm: %s and VPH: %s\nPlease check the 'cutouts.json' file."%(left_exposure.arm, left_exposure.VPH))
		sys.exit()
	
	i=0
	j=0
	dx_matrix = [
		[ 0, 0, 0], 
		[ 0, 0, 0], 
		[ 0, 0, 0]
	]
	dy_matrix = [
		[ 0, 0, 0], 
		[ 0, 0, 0], 
		[ 0, 0, 0]
	]  
	quiverLength = 1
		

	for cutout in cutouts:
		print(cutout)	
		print("Processing left...")
		print("\tSubtracting background...")
		left_exposure.cutout([cutout['x'], cutout['y'], cutout['width'], cutout['height']])	
		left_exposure.subtractBackground()
		matplotlib.pyplot.figure(figsize=(4,4))
		matplotlib.pyplot.imshow(1 - boostImageData(left_exposure.data), cmap='gray', origin='lower', aspect='equal')
		matplotlib.pyplot.title("centre x: %d y: %d"%(cutout['x'], cutout['y']))
		#matplotlib.pyplot.pause(0.1)
		filename = "%02d_left.png"%cutout['id']
		matplotlib.pyplot.savefig(filename)
		print("\tFinding bright spots...")
		left_exposure.segmentSources()


		print("Processing right...")
		print("\tSubtracting background...")
		right_exposure.cutout([cutout['x'], cutout['y'], cutout['width'], cutout['height']])
		right_exposure.subtractBackground()
		matplotlib.pyplot.figure(figsize=(4,4))
		matplotlib.pyplot.imshow(1 - boostImageData(right_exposure.data), cmap='gray', origin='lower', aspect='equal')
		matplotlib.pyplot.title("centre x: %d y: %d"%(cutout['x'], cutout['y']))
		#matplotlib.pyplot.pause(0.1)
		filename = "%02d_right.png"%cutout['id']
		matplotlib.pyplot.savefig(filename)
		print("\tFinding bright spots...")
		right_exposure.segmentSources()

		if args.plotsources:
			left_exposure.plotSources()
			right_exposure.plotSources()
			matplotlib.pyplot.show(block=False)

		# match sources
		distanceThreshold = 10
		matches = []
		for l in left_exposure.sources:
			#print(l)
			closestMatch = 1E6
			matched = False	
			match = {}
			for r in right_exposure.sources:
				distance = math.sqrt( (l['x'] - r['x'])**2 + (l['y'] - r['y'])**2)
				peakRatio = l['peak'] / r['peak']
				fluxRatio = l['flux'] / r['flux']
				if l['peak']>60000 or r['peak']>60000: 
					print("Saturated spot!")
					continue
				if distance<distanceThreshold and distance<closestMatch:
					closestMatch = distance
					print("peak ratio: peak: %.2f flux ratio: %f %.2f  %.2f  (%.0f, %.0f)"%(peakRatio, fluxRatio, l['peak'], r['peak'], l['x'], l['y']))	
					if (fluxRatio>0.8 and fluxRatio<1.2):
						match = { "leftSource": l, "rightSource": r, "distance": distance, "dx" : l['x'] - r['x'], "dy" : l['y'] - r['y'], 'flux' : numpy.mean([l['flux'], r['flux']])}
						matched = True
						print("...matched")
			if matched: matches.append(match)

		
		print("%d matches found."%len(matches))
		
		# Write them to file
		outputFile = open("%02d_amatches.json"%cutout['id'], "wt")
		json.dump(matches, outputFile, indent = 4)
		outputFile.close()

		xValues = [ m['leftSource']['x'] for m in matches]
		yValues = [ m['leftSource']['y'] for m in matches]
		vValues = [ m['dy'] for m in matches]
		uValues = [ m['dx'] for m in matches]
		fluxValues = [ m['flux'] for m in matches]
		
		median_dx = numpy.median(uValues)
		median_dy = numpy.median(vValues)
		uValues.sort()
		lowerP = 0.166666 * 100
		upperP = 0.833333 * 100
		lowerPercentile = numpy.percentile(uValues, lowerP)
		upperPercentile = numpy.percentile(uValues, upperP)
		#print("Percentile: %.0f%% %.2f"%(lowerP, lowerPercentile))
		#print("Percentile: %.0f%% %.2f"%(upperP, upperPercentile))
		#print("total points:", len(uValues))
		countInside = 0
		for u in uValues:
			if u<=upperPercentile and u>=lowerPercentile: countInside+=1
		print("%d points lie between these values or %.2f%%"%(countInside, countInside/len(uValues)*100))
		cutout['samples'] = len(uValues)
		cutout['stats'] = { "percentage" : round(countInside/len(uValues)*100,2),
							"samples_in_range" : countInside }
		cutout['stats']['range_lower'] = round(lowerPercentile, 2)
		cutout['stats']['range_upper'] = round(upperPercentile, 2)

		dx_matrix[i][j] = round(median_dx,2)
		dy_matrix[i][j] = round(median_dy,2)
		j=j+1
		if j==3:
			j=0
			i=i+1	

		cutout['median_dx'] = round(median_dx,2)
		cutout['median_dy'] = round(median_dy,2)
		print("median values {:.2f} dx and {:.2f} dy.".format(median_dx, median_dy))

		uValues = [ m['dx'] for m in matches]
		
		fig, ax = matplotlib.pyplot.subplots(figsize=(4,4))
		q = ax.quiver(xValues, yValues, uValues, vValues, scale=10)

		# Draw circles at the base of the arrows
		for index in range(len(xValues)):
			circle1 = matplotlib.pyplot.Circle((xValues[index], yValues[index]), numpy.log10(fluxValues[index]), color='b', fill=False)
			ax.add_patch(circle1)
			print(index, fluxValues[index], numpy.log10(fluxValues[index]))

		# quiverLength = round(numpy.mean(uValues), 0)
		ax.quiverkey(q, X=0.3, Y=-.10, U=quiverLength, label='Quiver key, length = %.0f pixels'%quiverLength, labelpos='E')
		# Add median quiver
		q = ax.quiver([250], [250], [median_dx], [median_dy], color='r', scale=10)
		ax.set_xlim(left=0, right=left_exposure.dimensions[1])
		ax.set_ylim(bottom=0, top=left_exposure.dimensions[0])
		matplotlib.pyplot.title("median dx: %.2f dy: %.2f (red arrow)"%(median_dx, median_dy))
		matplotlib.pyplot.draw()
		matplotlib.pyplot.savefig("%02d_zquiver.png"%cutout['id'])
		if args.plot: 
			matplotlib.pyplot.show(block=False)
			matplotlib.pyplot.pause(0.1)
		else: matplotlib.pyplot.close()
		left_exposure.reset()
		right_exposure.reset()
		# end of the loop through each of the 9 sections

	results = { "left_image" : { "filename": left_exposure.filename, "arm": left_exposure.arm, "VPH" : left_exposure.VPH, "shutter" : left_exposure.shutter, 
				"fpmode" : left_exposure.fpmode,
				"focusMTA": left_exposure.focusA,
				"focusMTB": left_exposure.focusB
				},
				"right_image" : { "filename": right_exposure.filename, "arm": right_exposure.arm, "VPH" : right_exposure.VPH, "shutter" : right_exposure.shutter,
				"fpmode" : right_exposure.fpmode,
				"focusMTA": right_exposure.focusA,
				"focusMTB": right_exposure.focusB
				 }
	}
	
	results['dx_matrix'] = dx_matrix
	results['dy_matrix'] = dy_matrix

	print("median of 9 dx values:", numpy.median(dx_matrix))
	tilt = [ 
		dx_matrix[0][2] - dx_matrix[0][0], 
		dx_matrix[1][2] - dx_matrix[1][0], 
		dx_matrix[2][2] - dx_matrix[2][0] 
	 ]
	print("tilt:", tilt)
	tip = [ 
		dx_matrix[0][0] - dx_matrix[2][0], 
		dx_matrix[1][1] - dx_matrix[2][1], 
		dx_matrix[0][2] - dx_matrix[2][2] 
	 ]
	print("tip:", tip)
	
	print("median tilt:", round(numpy.median(tilt), 3))
	print("median tip:", round(numpy.median(tip), 3))

	results["median_focus"] = round(numpy.median(dx_matrix), 3)
	results["median_tilt"] = round(numpy.median(tilt), 3)
	results["median_tip"] = round(numpy.median(tip), 3)

	results['cutouts'] = cutouts

	# Produce an overall quiver plot
	fig, ax = matplotlib.pyplot.subplots()
	fig.set_figheight(6)
	fig.set_figwidth(12)
	matplotlib.pyplot.imshow(1-boostImageData(left_exposure.data), cmap='gray', origin='lower', aspect='equal')
	xValues = [ c['x'] for c in cutouts ]
	yValues = [ c['y'] for c in cutouts ]
	uValues = [ c['median_dx'] for c in cutouts]
	vValues = [ c['median_dy'] for c in cutouts]
	q = ax.quiver(xValues, yValues, uValues, vValues, color='w', scale=10)
	ax.quiverkey(q, X=0.3, Y=-.10, U=quiverLength, label='Quiver key, length = %.0f pixels'%quiverLength, labelpos='E')
		
	ax.set_xlim(left=0, right=left_exposure.dimensions[1])
	ax.set_ylim(bottom=0, top=left_exposure.dimensions[0])
	median_dx = numpy.median(uValues)
	median_dy = numpy.median(vValues)
	matplotlib.pyplot.title("median dx: %.2f dy: %.2f"%(median_dx, median_dy))
	matplotlib.pyplot.draw()
	matplotlib.pyplot.savefig("full_quiver.png")
	if args.plot: 
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(0.1)
	else: matplotlib.pyplot.close()
	
	outfile = open("results.json", "wt")
	json.dump(results, outfile, indent=4)
	outfile.close()

	# Write the json to a <script> javascript file
	outfile = open("results.js", "wt")
	outfile.write("var results = ")
	outfile.write(json.dumps(results, indent=4))
	outfile.write(";")
	outfile.close()

	# Copy the .html template to the local directory

	from shutil import copyfile
	sourceFile = os.path.join(execPath, "results.html")
	destinationFile = "results.html"
	copyfile(sourceFile, destinationFile)
	
	if args.plot: input("press a key to quit")

	print("Results page is ready. Open it by typing 'firefox results.html'")


	sys.exit()
	
	
