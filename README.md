# download-nldas-forcings
This function will help users download [NLDAS](http://ldas.gsfc.nasa.gov/nldas/) data forcings to their local computer using a Matlab script

### Requirements

Matlab R2008a+. You can verify the version of Matlab by typing:

      version

You must first install [nctoolbox](https://github.com/nctoolbox/nctoolbox). Follow the instructions on that page to install.

You must have java version 7 or higher. You can verify the version of Javay used by Matlab by typing:

      version('-java')

The version returned should start with 'Java 1.7.' If it doesn't, you can try updating the Matlab JVM: http://www.mathworks.com/support/solutions/en/data/1-1812J/

### Demo

Use the file, 'callGetNldasForcing.m' to test out the downloading of two days of hourly NLDAS data. Locations to be downloaded are in the input file, 'inFile.txt.'

This file will create an output directory, 'forcingFromNasa,' with text files holding each location's NLDAS forcing data.
