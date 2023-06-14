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

	def segmentSources(self):
		### segment
		if self.arm == 'red':
			thresh = 15
		else:
			thresh = 20
		threshold = thresh * self.bkg.background_rms
		segment_map = detect_sources(self.data, threshold, npixels=10)
		cat = SourceCatalog(self.data, segment_map)
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
	parser = argparse.ArgumentParser(description='Loads the left and right Hartmann images and looks for spots.')
	parser.add_argument('--left', default="r3000987.fit", type=str, help='left hartmann exposure.')
	parser.add_argument('--right', default="r3000989.fit", type=str, help='right hartmann exposure.')
	parser.add_argument("--plot", action="store_true", help='Show the image')
	parser.add_argument("--version", action="store_true", help='Show the versions of the libraries used.')
	args = parser.parse_args()

	if args.version: 
		print("Python version:", sys.version)
		print("Astropy version:", astropy.__version__)
		print("Photutils version:", photutils.__version__)
		print("Scipy version:", scipy.__version__)


	left_exposure = exposure(args.left, doplot = args.plot)
	print(left_exposure.summary())
	right_exposure = exposure(args.right, doplot = args.plot)
	print(right_exposure.summary())

	metadata = {
		"left" : left_exposure.filename,
		"right" : right_exposure.filename,
		"width" : left_exposure.dimensions[1],
		"height" : left_exposure.dimensions[0],
		"arm" : left_exposure.arm, 
		"VPH" : left_exposure.VPH
	}
	outputMeta = open("metadata.json", "wt")
	json.dump(metadata, outputMeta)
	outputMeta.close()
	
	print("Processing left...")
	print("\tSubtracting background...")
	left_exposure.subtractBackground()
	print("\tFinding bright spots...")
	left_exposure.segmentSources()
	
	print("Processing right...")
	print("\tSubtracting background...")
	right_exposure.subtractBackground()
	print("\tFinding bright spots...")
	right_exposure.segmentSources()

	print("Stitched detector dimensions:", left_exposure.dimensions)
	
	metadata = {
		"left" : left_exposure.filename,
		"right" : right_exposure.filename,
		"width" : left_exposure.dimensions[1],
		"height" : left_exposure.dimensions[0],
		"arm" : left_exposure.arm, 
		"VPH" : left_exposure.VPH
	}
	outputMeta = open("metadata.json", "wt")
	json.dump(metadata, outputMeta)
	outputMeta.close()

	if args.plot:
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
	
	# Write them to file.
	outputFile = open("matches.json", "wt")
	json.dump(matches, outputFile, indent = 4)
	outputFile.close()

	#figure = matplotlib.pyplot.figure()
	xValues = [ m['leftSource']['x'] for m in matches]
	yValues = [ m['leftSource']['y'] for m in matches]
	uValues = [ m['dx'] for m in matches]
	vValues = [ m['dy'] for m in matches]

	fig, ax = matplotlib.pyplot.subplots()
	q = ax.quiver(xValues, yValues, uValues, vValues)
	ax.quiverkey(q, X=0.3, Y=1.05, U=2.5,
             label='Quiver key, length = 2.5 pixel', labelpos='E')
	ax.set_xlim(left=0, right=left_exposure.dimensions[1])
	ax.set_ylim(bottom=0, top=left_exposure.dimensions[0])
	

	#matplotlib.pyplot.quiver(xValues, yValues, uValues, vValues)
	#ax.quiverkey(q, X=0.3, Y=1.1, U=10,label='Quiver key, length = 10', labelpos='E')
	matplotlib.pyplot.draw()
	matplotlib.pyplot.savefig("tempquiver.png")
	matplotlib.pyplot.show()
	#input("Press Enter to continue...")
	
	sys.exit()
	
	