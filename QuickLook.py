#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

class QL:
    def __init__(self, biascorr=False, plot=False, verbose=False):
        self.biascorr = biascorr
        self.plot = plot
        self.verbose = verbose

    def quadrant(self, hdulist, chip, trimsec, biassec=None):
        # return array containing just bias-corrected quadrant
        # input FITS-standard x and y in trimsec and biassec
        x0 = trimsec[2]-1
        x1 = trimsec[3]
        y0 = trimsec[0]-1
        y1 = trimsec[1]
        if self.verbose: print('in quadrant: chip={} x0={} x1={} y0={} y1={}'.format(chip,x0,x1,y0,y1))
        quad = hdulist[chip].data[x0:x1,y0:y1]
        if self.biascorr:
            if biassec != None:
                biaslevel = np.median((hdulist[chip].data[biassec[2]-1:biassec[3],biassec[0]-1:biassec[1]]).flatten())
                quad = quad - biaslevel
                if self.verbose:
                    print('Chip {} quadrant bias level: {}'.format(chip,biaslevel))
        return quad

    def combineCCDs(self, fitsfile):
        with fits.open(fitsfile) as hdulist:
            self.header = hdulist[0].header
            if hdulist[0].header['DETECTOR'][:3] == 'WVB':
                self.arm = 'blue'
                if self.verbose: print('File {} Blue arm'.format(fitsfile))
                # blue 1
                newarray1 = np.zeros((6170,6164))

                # quadrant 1: blue1, TRIMSEC2= '[51:3122,1:3080] ', BIASSEC2= '[10:48,100:3079] '
                trimsec1 = [51, 3122, 1, 3080]
                biassec1 = [10, 48, 100, 3079]
                newarray1[0:trimsec1[3]-trimsec1[2]+1, 20:20+trimsec1[1]-trimsec1[0]+1] = self.quadrant(hdulist, 1, trimsec1, biassec1)

                # quadrant 2: blue1, TRIMSEC1= '[3155:6226,1:3080] ', '[6229:6267,100:3079] ' 
                trimsec2 = [3155, 6226, 1, 3080]
                biassec2 = [6229, 6267, 100, 3079]
                newarray1[0:trimsec2[3]-trimsec2[2]+1, 20+trimsec2[1]-trimsec2[0]+1:20+(trimsec2[1]-trimsec2[0]+1)+trimsec2[1]-trimsec2[0]+1] \
                    = self.quadrant(hdulist, 1, trimsec2, biassec2)

                # quadrant 3: blue1, TRIMSEC3= '[51:3122,3113:6192] ', BIASSEC3= '[10:48,3114:6093] '
                trimsec3 = [51, 3122, 3113, 6192]
                biassec3 = [10, 48, 3114, 6093]
                newarray1[trimsec2[3]-trimsec2[2]+1:(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1, 20:20+trimsec3[1]-trimsec3[0]+1] \
                    = self.quadrant(hdulist, 1, trimsec3, biassec3)

                # quadrant 4: blue1, TRIMSEC4= '[3155:6226,3113:6192] ', BIASSEC4= '[6229:6267,3114:6093] ' 
                trimsec4 = [3155, 6226, 3113, 6192]
                biassec4 = [6229, 6267, 3114, 6093]
                newarray1[trimsec2[3]-trimsec2[2]+1:(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1,
                          20+trimsec2[1]-trimsec2[0]+1:20+(trimsec2[1]-trimsec2[0]+1)+trimsec4[1]-trimsec4[0]+1] \
                          = self.quadrant(hdulist, 1, trimsec4, biassec4)

                newarray1 = newarray1[:,::-1]

                # blue 2
                newarray2 = np.zeros((6170,6144))

                # quadrant 1: blue2, TRIMSEC2= '[51:3122,1:3080] ', BIASSEC2= '[10:48,100:3079] '
                trimsec1 = [51, 3122, 1, 3080]
                biassec1 = [10, 48, 100, 3079]
                newarray2[9:trimsec1[3]-trimsec1[2]+1+9, 0:trimsec1[1]-trimsec1[0]+1] = self.quadrant(hdulist, 2, trimsec1, biassec1)

                # quadrant 2: blue2, TRIMSEC1= '[3155:6226,1:3080] ', '[6229:6267,100:3079] ' 
                trimsec2 = [3155, 6226, 1, 3080]
                biassec2 = [6229, 6267, 100, 3079]
                newarray2[9:trimsec2[3]-trimsec2[2]+1+9, trimsec2[1]-trimsec2[0]+1:(trimsec2[1]-trimsec2[0]+1)+trimsec2[1]-trimsec2[0]+1] \
                    = self.quadrant(hdulist, 2, trimsec2, biassec2)

                # quadrant 3: blue2, TRIMSEC3= '[51:3122,3113:6192] ', BIASSEC3= '[10:48,3114:6093] '
                trimsec3 = [51, 3122, 3113, 6192]
                biassec3 = [10, 48, 3114, 6093]
                newarray2[9+trimsec2[3]-trimsec2[2]+1:(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1+9, 0:trimsec3[1]-trimsec3[0]+1] \
                    = self.quadrant(hdulist, 2, trimsec3, biassec3)

                # quadrant 4: blue2, TRIMSEC4= '[3155:6226,3113:6192] ', BIASSEC4= '[6229:6267,3114:6093] ' 
                trimsec4 = [3155, 6226, 3113, 6192]
                biassec4 = [6229, 6267, 3114, 6093]
                newarray2[9+trimsec2[3]-trimsec2[2]+1:(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1+9,
                          trimsec2[1]-trimsec2[0]+1:(trimsec2[1]-trimsec2[0]+1)+trimsec4[1]-trimsec4[0]+1] \
                          = self.quadrant(hdulist, 2, trimsec4, biassec4)            

                newarray2 = newarray2[:,::-1]

                self.newarray = np.hstack((newarray1,newarray2))

            elif hdulist[0].header['DETECTOR'][:3] == 'WVR':
                self.arm = 'red'
                if self.verbose: print('File {} Red arm'.format(fitsfile))
                # red 1
                newarray1 = np.zeros((6170,6144))

                # quadrant 1: red1, TRIMSEC2= '[51:3122,1:3080] ', BIASSEC2= '[10:48,100:3079] '
                trimsec1 = [51, 3122, 1, 3080]
                biassec1 = [10, 48, 100, 3079]
                newarray1[2:2+trimsec1[3]-trimsec1[2]+1, 0:0+trimsec1[1]-trimsec1[0]+1] = self.quadrant(hdulist, 1, trimsec1, biassec1)

                # quadrant 2: red1, TRIMSEC1= '[3155:6226,1:3080] ', '[6229:6267,100:3079] ' 
                trimsec2 = [3155, 6226, 1, 3080]
                biassec2 = [6229, 6267, 100, 3079]
                newarray1[2:2+trimsec2[3]-trimsec2[2]+1, 0+trimsec2[1]-trimsec2[0]+1:0+(trimsec2[1]-trimsec2[0]+1)+trimsec2[1]-trimsec2[0]+1] \
                    = self.quadrant(hdulist, 1, trimsec2, biassec2)

                # quadrant 3: red1, TRIMSEC3= '[51:3122,3113:6192] ', BIASSEC3= '[10:48,3114:6093] '
                trimsec3 = [51, 3122, 3113, 6192]
                biassec3 = [10, 48, 3114, 6093]
                newarray1[2+trimsec2[3]-trimsec2[2]+1:2+(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1, 0:0+trimsec3[1]-trimsec3[0]+1] \
                    = self.quadrant(hdulist, 1, trimsec3, biassec3)

                # quadrant 4: red1, TRIMSEC4= '[3155:6226,3113:6192] ', BIASSEC4= '[6229:6267,3114:6093] ' 
                trimsec4 = [3155, 6226, 3113, 6192]
                biassec4 = [6229, 6267, 3114, 6093]
                newarray1[2+trimsec2[3]-trimsec2[2]+1:2+(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1,
                          0+trimsec2[1]-trimsec2[0]+1:0+(trimsec2[1]-trimsec2[0]+1)+trimsec4[1]-trimsec4[0]+1] \
                          = self.quadrant(hdulist, 1, trimsec4, biassec4)            

                # red 2
                newarray2 = np.zeros((6170,6164))

                # quadrant 1: red2, TRIMSEC2= '[51:3122,1:3080] ', BIASSEC2= '[10:48,100:3079] '
                trimsec1 = [51, 3122, 1, 3080]
                biassec1 = [10, 48, 100, 3079]
                newarray2[0:trimsec1[3]-trimsec1[2]+1, 0:trimsec1[1]-trimsec1[0]+1] = self.quadrant(hdulist, 2, trimsec1, biassec1)

                # quadrant 2: red2, TRIMSEC1= '[3155:6226,1:3080] ', '[6229:6267,100:3079] ' 
                trimsec2 = [3155, 6226, 1, 3080]
                biassec2 = [6229, 6267, 100, 3079]
                newarray2[0:trimsec2[3]-trimsec2[2]+1, trimsec2[1]-trimsec2[0]+1:(trimsec2[1]-trimsec2[0]+1)+trimsec2[1]-trimsec2[0]+1] \
                    = self.quadrant(hdulist, 2, trimsec2, biassec2)

                # quadrant 3: red2, TRIMSEC3= '[51:3122,3113:6192] ', BIASSEC3= '[10:48,3114:6093] '
                trimsec3 = [51, 3122, 3113, 6192]
                biassec3 = [10, 48, 3114, 6093]
                newarray2[0+trimsec2[3]-trimsec2[2]+1:(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1,
                          0:trimsec3[1]-trimsec3[0]+1] = self.quadrant(hdulist, 2, trimsec3, biassec3)

                # quadrant 4: red2, TRIMSEC4= '[3155:6226,3113:6192] ', BIASSEC4= '[6229:6267,3114:6093] ' 
                trimsec4 = [3155, 6226, 3113, 6192]
                biassec4 = [6229, 6267, 3114, 6093]
                newarray2[0+trimsec2[3]-trimsec2[2]+1:(trimsec2[3]-trimsec2[2]+1)+trimsec3[3]-trimsec3[2]+1,
                          trimsec2[1]-trimsec2[0]+1:(trimsec2[1]-trimsec2[0]+1)+trimsec4[1]-trimsec4[0]+1] \
                          = self.quadrant(hdulist, 2, trimsec4, biassec4)

                self.newarray = np.hstack((newarray2,newarray1))

        if self.plot:
            plt.figure()
            plt.imshow(boostImageData(self.newarray), cmap='gray', origin='lower', aspect='equal')
            #plt.colorbar()
            plt.show(block=True)

    def writeArray(self, outfile):
        fits.writeto(outfile, self.newarray, self.header)

def boostImageData(imageData):
    """ Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255 """
    hi = 99
    lo = 20
    data = imageData
    max = data.max()
    dataArray = data.flatten()
    pHi = np.percentile(dataArray, hi)
    pLo = np.percentile(dataArray, lo)
    range = pHi - pLo
    scale = range/255
    data = np.clip(data, pLo, pHi)
    data-= pLo
    data/=scale
    return data

if __name__ == "__main__":
    QLook = QL(biascorr=True, plot=False, verbose=False)
    
    fitslist = sorted(glob.glob('r2*.fit'))
    rcounter=0
    bcounter=0
    for fitsfile in fitslist:
        QLook.combineCCDs(fitsfile)
        if QLook.arm == 'blue':
            QLook.writeArray('blue{}.fits'.format(bcounter))
            print("{} --> blue{}.fits".format(fitsfile, bcounter))
            bcounter = bcounter+1
        else:
            QLook.writeArray('red{}.fits'.format(rcounter))
            print("{} --> red{}.fits".format(fitsfile, rcounter))
            rcounter = rcounter+1
