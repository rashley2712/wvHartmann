#!/usr/bin/env python3
import argparse, sys, json, datetime, os, re
import matplotlib.pyplot
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
		print(self.thresh)
		thresh  = self.thresh
		if thresh < 0:
			if self.arm == 'red':
				thresh = 5
			else:
				thresh = 20
		print("threshold is", thresh)
		data = self.data
		threshold = thresh * self.bkg.background_rms
		segment_map = detect_sources(data, threshold, npixels=10)
		cat = SourceCatalog(data, segment_map)
		tbl = cat.to_table()

		for t in tbl:
			source = { "x" :  round(t['xcentroid'], 3), "y" : round(t['ycentroid'], 3), "peak" : round(t['max_value'], 3)}
			self.sources.append(source)
		print("\t%d sources found."%len(self.sources))

	def plotSources(self):
		self.figure = matplotlib.pyplot.figure()
		xValues = [ s['x'] for s in self.sources]
		yValues = [ s['y'] for s in self.sources]
		cValues = [ s['peak'] for s in self.sources]
		matplotlib.pyplot.scatter(xValues, yValues, c=cValues, cmap="hot")
		matplotlib.pyplot.draw()
		matplotlib.pyplot.pause(0.1)

		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads the left and right Hartmann images and looks for spots. This version cuts out regions of the detector.')
	parser.add_argument('--left', default="r3000987.fit", type=str, help='left hartmann exposure.')
	parser.add_argument('--right', default="r3000989.fit", type=str, help='right hartmann exposure.')
	parser.add_argument("--plot", action="store_true", help='Hold the plots on the display.')
	parser.add_argument("--plotsources", action="store_true", help='Plots the found sources as a schematic.')
	parser.add_argument("--version", action="store_true", help='Show the versions of the libraries used.')
	parser.add_argument("--cutouts", default="../cutouts.json", help="A JSON file defining the cutouts (aka sub-regions) of the image to use.")
	parser.add_argument("--thresh", default=-1, type=int,  help="Threshhold for segmentation (default 5 red and 20 blue).")
	
	args = parser.parse_args()

	if args.version: 
		print("Python version:", sys.version)
		print("Astropy version:", astropy.__version__)
		print("Photutils version:", photutils.__version__)
		print("Scipy version:", scipy.__version__)

	# Load the cutouts
	cutoutsDefFile = open(args.cutouts, "rt")
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
		print("Could not find cutout info for arm: %s and VPH: %s"%(left_exposure.arm, left_exposure.VPH))
	
	i=0
	j=0
	dx_matrix = [
		[ 0, 0, 0], 
		[ 0, 0, 0], 
		[ 0, 0, 0]
	] 

	for cutout in cutouts:
		print(cutout)	
		print("Processing left...")
		print("\tSubtracting background...")
		left_exposure.cutout([cutout['x'], cutout['y'], cutout['width'], cutout['height']])	
		left_exposure.subtractBackground()
		matplotlib.pyplot.figure()
		matplotlib.pyplot.imshow(boostImageData(left_exposure.data), cmap='gray', origin='lower', aspect='equal')
		#matplotlib.pyplot.pause(0.1)
		filename = "%02d_left.png"%cutout['id']
		matplotlib.pyplot.savefig(filename)
		print("\tFinding bright spots...")
		left_exposure.segmentSources()


		print("Processing right...")
		print("\tSubtracting background...")
		right_exposure.cutout([cutout['x'], cutout['y'], cutout['width'], cutout['height']])
		right_exposure.subtractBackground()
		matplotlib.pyplot.figure()
		matplotlib.pyplot.imshow(boostImageData(right_exposure.data), cmap='gray', origin='lower', aspect='equal')
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

				if distance<distanceThreshold and (peakRatio>0.8 and peakRatio<1.2) and distance<closestMatch:
					match = { "leftSource": l, "rightSource": r, "distance": distance, "dx" : l['x'] - r['x'], "dy" : l['y'] - r['y']}
					closestMatch = distance
					matched = True
			if matched: matches.append(match)

		
		print("%d matches found."%len(matches))
		
		# Write them to file
		outputFile = open("%02d_amatches.json"%cutout['id'], "wt")
		json.dump(matches, outputFile, indent = 4)
		outputFile.close()

		xValues = [ m['leftSource']['x'] for m in matches]
		yValues = [ m['leftSource']['y'] for m in matches]
		uValues = [ m['dx'] for m in matches]
		vValues = [ m['dy'] for m in matches]

		median_dx = numpy.median(uValues)
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

		dx_matrix[i][j] = round(median_dx,3)
		j=j+1
		if j==3:
			j=0
			i=i+1	

		median_dy = numpy.median(vValues)
		cutout['median_dx'] = round(median_dx,2)
		cutout['median_dy'] = round(median_dy,2)
		print("median values {:.2f} dx and {:.2f} dy.".format(median_dx, median_dy))

		fig, ax = matplotlib.pyplot.subplots()
		q = ax.quiver(xValues, yValues, uValues, vValues)
		ax.quiverkey(q, X=0.3, Y=1.05, U=2.5,
				label='Quiver key, length = 2.5 pixel', labelpos='E')
		ax.set_xlim(left=0, right=left_exposure.dimensions[1])
		ax.set_ylim(bottom=0, top=left_exposure.dimensions[0])
		

		#matplotlib.pypl/--ot.quiver(xValues, yValues, uValues, vValues)
		#ax.quiverkey(q, X=0.3, Y=1.1, U=10,label='Quiver key, length = 10', labelpos='E')
		matplotlib.pyplot.draw()
		matplotlib.pyplot.savefig("%02d_zquiver.png"%cutout['id'])
		if args.plot: 
			matplotlib.pyplot.show(block=False)
			matplotlib.pyplot.pause(0.1)
		else: matplotlib.pyplot.close()
		left_exposure.reset()
		right_exposure.reset()
	
	print(json.dumps(cutouts, indent=4))

	results = { "left_image" : { "filename": left_exposure.filename, "arm": left_exposure.arm, "VPH" : left_exposure.VPH, "shutter" : left_exposure.shutter },
				"right_image" : { "filename": right_exposure.filename, "arm": right_exposure.arm, "VPH" : right_exposure.VPH, "shutter" : right_exposure.shutter }
	}
	
	results['dx_matrix'] = dx_matrix

	print(dx_matrix)
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

	outfile = open("results.json", "wt")
	json.dump(results, outfile, indent=4)
	outfile.close()

	outfile = open("results.js", "wt")
	outfile.write("var results = ")
	outfile.write(json.dumps(results, indent=4))
	outfile.write(";")
	outfile.close()
	
	if args.plot: input("press a key to quit")



	sys.exit()
	
	