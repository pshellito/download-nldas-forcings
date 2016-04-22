# download-nldas-forcings
The Matlab function, 'getNldasForcing.m' will help users download [NLDAS](http://ldas.gsfc.nasa.gov/nldas/) data forcings to their local computer.

### Requirements

Matlab R2008a+. You can verify the version of Matlab by typing:

      version

You must first install [nctoolbox](https://github.com/nctoolbox/nctoolbox). Follow the instructions on that page to install.

You must have java version 7 or higher. You can verify the version of Java used by Matlab by typing:

      version('-java')

The version returned should start with 'Java 1.7.' If it doesn't, you can try updating the Matlab JVM: http://www.mathworks.com/support/solutions/en/data/1-1812J/

### Demo

The script, 'callGetNldasForcing.m' will
* Read the site name, latitude, and longitude of two arbitrary locations found in the input file, 'inFile.txt.'
* Define the range of dates over which to download. This is initially set to only 2 days for testing, which should take about a minute.
* Place output in a directory it creates called 'forcingFromNasa.' Each location's forcing data will be in a tab-delimited .txt file named after the site name. Details about the site's location are provided in a header in each file.

The above options can be changed in the 'callGetNldasForcing.m' script.

### Contact
If you have any questions or edits, please [contact me](mailto:peter.shellito@colorado.edu) or create a new pull request.
