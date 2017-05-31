# download-nldas-forcings

## Overview
The Matlab script, 'callGetNldasForcing.m' will call the function, 'getNldasForcing.m,' which downloads [NLDAS](http://ldas.gsfc.nasa.gov/nldas/) primary forcings (A) at any number of point locations and saves them as individual text files that contain the corresponding meteorological timeseries. If a requested location is not in the NLDAS domain, the script will issue a warning and download the closest NLDAS pixel to the requested location.

Use 'inFile_test.txt' to select site locations to save.

## Requirements

You must first install [nctoolbox](https://github.com/nctoolbox/nctoolbox). Follow the instructions on that page to install.

Matlab R2008a+. You can verify the version of Matlab by typing:

      version

Java version 7 or higher. You can verify the version of Java used by Matlab by typing:

      version('-java')

The version returned should start with 'Java 1.7.' If it doesn't, you can try updating the Matlab JVM: http://www.mathworks.com/support/solutions/en/data/1-1812J/

You must register with Earthdata and authorize NASA GESDICS DATA ARCHIVE Data Access in Earthdata Login. You must then set up .netrc and create a cookie file. To do both these steps, follow #2 in the "Procedure" section [here](https://disc.sci.gsfc.nasa.gov/recipes/?q=recipes/How-to-Download-Data-Files-from-HTTP-Service-with-wget)
## Demo

The script, 'callGetNldasForcing.m' will
* Read the site name, latitude, and longitude of two arbitrary locations found in the input file, 'inFile.txt.'
* Define the range of dates over which to download. This is initially set to only 2 days for testing, which should take about a minute to run.
* Place output in a directory it creates called 'forcingFromNasa.' Each location's forcing data will be in a tab-delimited .txt file named after the site name. Details about the site's location are provided in a header in each file.

The above options can be changed in the 'callGetNldasForcing.m' script.

Using a CU internet connection and a 2015 Macbook Pro, it takes about 4 hours to download one year of forcing data.

## Contact
If you have any questions or edits, please [contact me](mailto:peter.shellito@colorado.edu) or create a new pull request.
