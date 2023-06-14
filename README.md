# wvHartmann
Python scripts for performing shuttered Hartmann tests on the WEAVE spectrograph

## Installing this code

Perform a git clone into a local directory on your machine.
```
git clone https://github.com/rashley2712/wvHartmann
```

Make sure that you have a relatively up to date version of the scip

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
