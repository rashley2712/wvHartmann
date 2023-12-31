# wvHartmann
Python scripts for performing shuttered Hartmann tests on the WEAVE spectrograph

## Installing this code

Perform a git clone into a local directory on your machine.
```
git clone https://github.com/rashley2712/wvHartmann
```

Make sure that you have a relatively up to date version of the scipy, astropy and photutils Python libraries.

You can check your versions of these libraries by running 
```
./cutoutSpots.py --version  
```
The versions of the libraries that work for me are shown below
```
Python version: 3.10.6 (main, May 29 2023, 11:10:38) [GCC 11.3.0]
Astropy version: 5.2.2
Photutils version: 1.7.0
Scipy version: 1.10.1
```

Upgrade your instance of these libraries using the following command. 
```
pip3 install scipy astropy photutils --upgrade 
```
## Running the script

Get help for the command line parameters by typing
```
./cutoutSpots.py -h
```
1. Create a folder for your work. Note that the script will create a bunch of json, html and png files when it runs. This will mean clutter in your working directory. 
2. Copy the two .FITS files, your left and right Hartmann exposures, into this folder. 
3. Run the script by running
   ```
      [pathto]/cutoutSpots.py -l r4000332.fit -r r4000331.fit
   ```
4. Once the script has run, view your results by typing
      ```
         firefox results.html
      ```
      
