#!/usr/bin/env python3
import argparse, sys, json, datetime, os, re
import matplotlib.pyplot
import numpy
from astropy.io import fits
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Searches an obsdata folder for Hartmann tests (looking in the FITS headers).')
	parser.add_argument('-d', '--directory', default="/obsdata/whta/20230418", type=str, help='Root folder for the obsdata. Default is ''/obsdata/whta/[yesterday]''.')
	args = parser.parse_args()
	
	folder = os.listdir(args.directory)
	
	folderList = []
	for f in folder:
		if re.search("r[0-9]{7}.fit", f): folderList.append(f)

	folderList.sort(reverse=True)

	limit = 4E6
	counter = 1
	saved = []
	for f in folderList:
		print(counter, f)
		allHeaders = {}
		# Get the FITS headers
		hdulist = fits.open(os.path.join(args.directory, f))  # open a FITS file
		hdr = hdulist[0].header  # the primary HDU header
		#for h in hdulist:
			#print(h)
		#	if type(h.data) is numpy.ndarray:
		#		imageObject = {}
				# imageObject['data'] = h.data
		#		imageObject['size'] = numpy.shape(h.data)
			
				# print(imageObject)

		#for card in hdulist:
		card = hdulist[0]
		for key in card.header.keys():
			allHeaders[key] = card.header[key]
		
		hdulist.close(output_verify='ignore')
		try:
			print("HART shutters: ", allHeaders['HART'])
			if allHeaders['HART']!="BOTH_OPEN":
				save = {
					"filename": f,
					"HART" : allHeaders['HART']
				}
				saved.append(save)
		except KeyError:
			print("No HART header info... skipping")
		if counter==limit: break
		counter+=1

	print("Found the following HARTMANN tests:")
	
	copyScript = "#!/bin/bash\n"

	for s in saved:
		print(s['filename'])
		copyScript+="cp %s .\n"%(os.path.join(args.directory, s['filename']))

	print(copyScript)

	writer = open("copyscript.bash", "wt")
	writer.write(copyScript)
	writer.close()

	