# download-nldas-forcings

## Overview
The Matlab function, 'getNldasForcing.m' downloads [NLDAS](http://ldas.gsfc.nasa.gov/nldas/) primary forcings (A) at any number of point locations and saves them as individual text files that contain the corresponding meteorological timeseries. If a requested location is not in the NLDAS domain, the script will issue a warning and download the closest NLDAS pixel to the requested location.

## Requirements

You must first install [nctoolbox](https://github.com/nctoolbox/nctoolbox). Follow the instructions on that page to install.

Matlab R2008a+. You can verify the version of Matlab by typing:

      version

Java version 7 or higher. You can verify the version of Java used by Matlab by typing:

      version('-java')

The version returned should start with 'Java 1.7.' If it doesn't, you can try updating the Matlab JVM: http://www.mathworks.com/support/solutions/en/data/1-1812J/

## Demo

The script, 'callGetNldasForcing.m' will
* Read the site name, latitude, and longitude of two arbitrary locations found in the input file, 'inFile.txt.'
* Define the range of dates over which to download. This is initially set to only 2 days for testing, which should take about a minute to run.
* Place output in a directory it creates called 'forcingFromNasa.' Each location's forcing data will be in a tab-delimited .txt file named after the site name. Details about the site's location are provided in a header in each file.

The above options can be changed in the 'callGetNldasForcing.m' script.

## Contact
If you have any questions or edits, please [contact me](mailto:peter.shellito@colorado.edu) or create a new pull request.
