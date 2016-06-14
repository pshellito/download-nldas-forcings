function [ outDir ] = getNldasForcing(qNames, qLat, qLon, qStart, qEnd, outDir)
% GETNLDASFORCING This script will download NLDAS primary forcing data from
%       Nasa's servers.
% Created by Peter J. Shellito 2/16/16
% 
% qNames: cell array of site names
% qLat: latitude corresponding to each site
% qLon: longitude corresponding to each site
% qStart: a vector [yyyy, mm, dd] specifying the day to start downloading.
%       Start date must be [1979, 1, 2] or later. If set to 'apnd,'
%       script will look for an existing file with the
%       same site name and start at the last date in that file and
%       append new data to it.
% qEnd: a vector [yyyy, mm, dd] specifying the day to stop downloading. End
%       date must be 4 days before today or earlier. Neither qStart nor 
%       qEnd support a starting hour and minute.
% outDir: Directory to place the output files. If nothing is provided,
%       default is to create a directory in the present directory titled
%       './outFiles/'
% 
% For your own reference, or if something goes wrong, the directory where
% you can find all the NLDAS data:
% http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings
% And the location of a sample file:
% fileDir = 'ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/2015/001/NLDAS_FORA0125_H.A20150101.0600.002.grb';
% The file containing elevation data is from here:
% http://ldas.gsfc.nasa.gov/nldas/NLDASelevation.php
% The file condtaining soil data is from here:
% http://ldas.gsfc.nasa.gov/nldas/NLDASsoils.php
% 
% Initial testing revealed that processing one month took 20 minutes using
% a CU internet connection.
% 
% ========================================================================
% NB: This script requires the user to have downloaded and installed the
% nctoolbox: http://nctoolbox.github.io/nctoolbox/
% Extract nctoolbox-1.1.3, and then in the matlab command line, navigate 
% to the nctoolbox: cd('Path/to/toolbox/nctoolbox-1.1.3')
% Then type 'setup_nctoolbox' and you should be able to use the toolbox.
% It may be necessary to run setup_nctoolbox every time you start a new 
% matlab session. To automate this, add the '/Path/to/toolbox/' to your 
% matlab path. Then edit (or create) a startup.m file that is also in your
% matlab path. Add the following to startup.m: setup_nctoolbox;
% ========================================================================
% 
% Copyright (c) 2016, Peter J. Shellito, University of Colorado Boulder
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% -------------------------------------------------------------------------
% If no output directory was provided
if nargin<6
    outDir = './outFiles';
end

% -------------------------------------------------------------------------
% Some initial checks

% Names must be a cell array of strings
if ~iscellstr(qNames)
    error('Input names must be a cell array of strings.')
end

% If a start date is provided
if isfloat(qStart)
    % Do not append files
    apnd = false;
    % The start datenum is provided in qStart
    dnStart = floor(datenum(qStart));
    if dnStart < datenum(1979,1,2)
        dnStart = datenum(1979,1,2);
    end
% If the string 'apnd' is provided, use the date of the last record
elseif strcmp(qStart,'apnd')
    % Append these files, starting with the datenum after the last line,
    % which should be at hour 2300.
    apnd = true;
    % Open the last file
    fid = fopen([outDir '/' qNames{end} '.txt'],'r');
    if fid == -1
        error(['There is no file by name ' outDir '/' qNames{end} '.txt to which to append data'])
    end
    % Count the lines in that file
    nLines = 0;
    tline = fgetl(fid);
    while ischar(tline)
        tline = fgetl(fid);
        nLines = nLines+1;
    end
    % Read the last line to get the last date
    frewind(fid)
    lastLine = textscan(fid, '%f', 'headerlines', nLines-1);
    % Assign the last date in the file as the date to start with
    qStart = [lastLine{1}(1:5)' 0];
    dnStart = ceil(datenum(qStart));
    
else
    error('qStart must be either a start date or the string ''apnd'' if you wish to append to an existing record')
end

% Start date is not after end date
dnEnd = floor(datenum(qEnd));
if dnEnd > floor(datenum(now))-4
    dnEnd = floor(datenum(now))-4;
end
if dnEnd < dnStart
    error('Start date cannot be after end date')
end


% -------------------------------------------------------------------------
% Set up some variables
% Number of sites requested
nSites = length(qNames);
% Initialize vectors for lat/lon idcs
latIdcs = nan(nSites,1);
lonIdcs = nan(nSites,1);
% The years and days of years of the query
qDatenums = dnStart:dnEnd;
% Initialize a day of year vector
qDoy = nan(1,length(qDatenums));
% The years, months, and days of the query
[qYears, qMonths, qDays] = datevec(qDatenums);
% The unique years in the query
qUniqueYears = unique(qYears);
% For each unique year
for yy = 1:length(qUniqueYears)
    % Find all the indices where qYears matches this unique year
    sameYrIdcs = find(qYears == qUniqueYears(yy));
    % Subtract the datenum of day zero of the year from the datenums in the
    % year to get day of year
    qDoy(sameYrIdcs) = datenum(qDatenums(sameYrIdcs)) - datenum(qUniqueYears(yy),0,0);
end % Loop through each unique year
% Convert qYears to string
qYearStr = num2str(qYears');
% Convert qMonths to string
qMonthStr = num2str(qMonths', '%02d');
% Convert qDays to string
qDayStr = num2str(qDays', '%02d');
% Convert qDoy to string
qDoyStr = num2str(qDoy', '%03d');
% Create hour strings
qHourStr = num2str((0:100:2300)','%04d');
% The directory where nldas forcings are held
ftpBaseDir = '/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/';
% The local directory where nldas forcings will be placed
localBaseDir = [pwd '/data'];
% If there is already a directory named data in the working direcotry, do not continue because
% it will be deleted at the end of this script and I don't want to delete
% anything that this script itself did not create.
if exist(localBaseDir, 'dir') == 7
    error('%s\n%s\n%s','This function requires use of the local directory named ''./data.'',', ...
        'and you already have a directory with that name. We are stopping here',...
        'because the end of this function will DELETE ./data and all files within it.')
end
% The beginning of the file name of forcings
ftpBaseFn = 'NLDAS_FORA0125_H.A';
% The end of the file name of forcings
ftpEndFn = '.002.grb';
% Strings of variables to read
latStr = 'lat'; % Center of 1/8 degree pixel
lonStr = 'lon'; % Center of 1/8 degree pixel
pressureStr = 'Pressure_surface'; % Surface pressure [Pa]
windStrU = 'u-component_of_wind_height_above_ground'; % 10-m above ground Zonal wind speed [m/s]
windStrV = 'v-component_of_wind_height_above_ground'; % 10-m above ground Meridional wind speed [m/s]
temperatureStr = 'Temperature_height_above_ground'; % 2-m above ground Temperature [K]
specificHumidityStr = 'Specific_humidity_height_above_ground'; % 2-m above ground Specific humidity [kg/kg]
shortwaveStr = 'Downward_shortwave_radiation_flux_surface'; % Shortwave radiation flux downwards (surface) [W/m^2]
longwaveStr = 'Downward_longwave_radiation_flux_surface'; % Longwave radiation flux downwards (surface) [W/m^2]
precipStr = 'Total_precipitation_surface_1_Hour_Accumulation'; % Precipitation hourly accumulated total [kg/m^2]
potentialEvapStr = 'Potential_evaporation_surface_1_Hour_Accumulation'; % Potential evaporation hourly total [kg/m^2]
convFractStr = 'VAR7-12-130-153_surface_1_Hour_Accumulation'; % Fraction of total precipitation that is convective [unitless]
capeStr = 'Convective_available_potential_energy_CAPE_layer_between_two_pressure_difference_from_ground_layer'; % 180-0 mb above ground Convective Available Potential Energy [J/kg] 
% Put all the strings into one cell
allStrings = {pressureStr; windStrU; windStrV; temperatureStr; ...
    specificHumidityStr; shortwaveStr; longwaveStr; precipStr; ...
    potentialEvapStr; convFractStr; capeStr};
% The number of variables I'm getting
nVars = length(allStrings);
% Names and units for the variables (for headers)
namesStrAll = 'Pressure  WindU       WindV       Temperature   SpecificHumid   Shortwave  Longwave  Precip      PE          ConvFrac    CAPE';
unitsStrAll = '[Pa]      (10m)[m/s]  (10m)[m/s]  (2m)[K]       (2m)[kg/kg]     [W/m2]     [W/m2]    [kg/m2]     [kg/m2]     [unitless]  [J/kg]';
% Variable formats
varFmt = '%08.1f %+10f  %+10f   %5.2f %17.4e %11.2f %10.2f %13.4e %11.4e %11.4e %11.4e\n';
% Date formats
dateFmt = '%04d   %02d    %02d  %02d   %02d     ';
% -------------------------------------------------------------------------
% Create a directory to hold output data
if exist(outDir, 'dir') ~= 7
    disp(['Making an output directory here: ' outDir])
    mkdir(outDir);
end

% -------------------------------------------------------------------------
% Obtain soil texture data
% Name of the file holding dominant layer 1 soil textures
soilFileName = './STEX_TAB.01.txt';
% Descriptions of each texture number
soilTxtDescription = {'sand', 'loamy sand', 'sandy loam', 'silt loam', ...
    'silt', 'loam', 'sandy clay loam', 'silty clay loam', 'clay loam', ...
    'sandy clay', 'silty clay', 'clay', 'organic materials', 'water', ...
    'bedrock', 'other'};
% Read soil texture file
soilDf = load(soilFileName);
% Assign lats and lons to a vector each
allSoilLat = soilDf(:,4);
allSoilLon = soilDf(:,3);
% Assign unique lats and lons to a vector each
soilLat = unique(allSoilLat);
soilLon = unique(allSoilLon);
% Predominant soil textures
[~, allSoil] = max(soilDf(:,5:20), [], 2);
% Remove soili data frame from memory
clear soilDf
% Initialize a vector for soil texture at each site
soilTxt = nan(nSites,1);
% Loop through each site to get the predom soil texture (indices of soil
% texture are different from indices from elevation and forcing data, so I
% have to do this separately)
for ss = 1:nSites
    % lat/lon idcs of the site
    [~, soilLatIdx] = min(abs(soilLat-qLat(ss)));
    [~, soilLonIdx] = min(abs(soilLon-qLon(ss)));
    % The NLDAS lat/lon of this site
    nearestSoilLat = soilLat(soilLatIdx);
    nearestSoilLon = soilLon(soilLonIdx);
    % The overall index of this lat/lon
    soilTxtIdx = find(all([allSoilLat==nearestSoilLat, allSoilLon==nearestSoilLon],2));
    % The dominant soil type of this site
    soilTxt(ss) = allSoil(soilTxtIdx);
end % ss loop through each site

% -------------------------------------------------------------------------
% Obtain elevation data
% Name of the file holding NLDAS mean elevation of the
% topography from the GTOP30 dataset.
elevFileName = './gtopomean15k.asc';
% Read elevation file
elevDf = load(elevFileName);
% Assign unique lats and lons to a vector each
lat = unique(elevDf(:,3));
lon = unique(elevDf(:,4));
% Size of these vectors
nLat = length(lat);
nLon = length(lon);
% Assign mean elevation to a vector
allElev = elevDf(:,5);

% -------------------------------------------------------------------------
% Set up output files. Write headers with lat/lon data in them. Get lat/lon
% idcs.
% Initialize a vector of file IDs
fid = nan(1,nSites);
% Loop through each site
for ss = 1:nSites
    disp(qNames{ss})
    % Get lat/lon idcs
    [latDiff latIdcs(ss)] = min(abs(lat-qLat(ss)));
    [lonDiff lonIdcs(ss)] = min(abs(lon-qLon(ss)));
    disp(latIdcs(ss))
    disp(lonIdcs(ss))
    % Get the rounded lat and lon. This is the center point of the NLDAS
    % cell.
    nearestLat = lat(latIdcs(ss));
    nearestLon = lon(lonIdcs(ss));
    % Get the elevation
    meanElevation = allElev((latIdcs(ss)-1)*nLon+lonIdcs(ss));
    % If the difference between the queried and actual lat or lon is too
    % big, display a warning
    if latDiff>0.125 || lonDiff>0.125
        warning('%s site does not have a nearby nldas cell. \nQueried lat/lon: %f %f\nClosest NLDAS lat/lon: %f %f',...
            qNames{ss}, qLat(ss), qLon(ss), nearestLat, nearestLon)
    end
    % The name for this output file
    outFile = [outDir '/' qNames{ss} '.txt'];

    % If we are not appending data, write a header
    if ~apnd
        % Open the output file for the first time
        fid(ss) = fopen(outFile,'w');
        % Print header lines to the file
        fprintf(fid(ss),['%% Site: ' qNames{ss} '\n']);
        fprintf(fid(ss),['%% Site lat/lon:\n']);
        fprintf(fid(ss),['%% ' num2str([qLat(ss) qLon(ss)], '%12.8f\t') '\n']);
        fprintf(fid(ss),['%% Closest NLDAS pixel (1/8 degree) center:\n']);
        fprintf(fid(ss),['%% ' num2str([nearestLat nearestLon], '%8.4f\t') '\n']);
        fprintf(fid(ss),['%% Avg NLDAS pixel elevation (m):\n']);
        fprintf(fid(ss),['%% ' num2str(meanElevation, '%9.4f') '\n']);
        fprintf(fid(ss),['%% Dominant surface soil texture number:\n']);
        fprintf(fid(ss),['%% ' num2str(soilTxt(ss), '%d') '\n']);
        fprintf(fid(ss),['%% Dominant surface soil texture description:\n']);
        fprintf(fid(ss),['%% ' soilTxtDescription{soilTxt(ss)} '\n']);
        fprintf(fid(ss),['%% File created on ' date '.\n']);
        fprintf(fid(ss),['%% Date                       ' namesStrAll '\n']);
        fprintf(fid(ss),['%% year month day hour minute ' unitsStrAll '\n']);
    else
        % If we are appending data, open the output file with append
        % permission only.
        fid(ss) = fopen(outFile,'a');
    end
end % Loop through each site

% -------------------------------------------------------------------------
% Set up ftp connection
% The Nasa host that holds nldas forcings
nasaHost = 'hydro1.sci.gsfc.nasa.gov';
% Create an ftp object (open the connection)
ftpObj = ftp(nasaHost);

% -------------------------------------------------------------------------
% Loop through each day in the record
for dd = 1:length(qDatenums)
    % Location of this day's data on the local machine
    localDir = [pwd ftpBaseDir qYearStr(dd,:) '/' qDoyStr(dd,:)];
    % Loop through each hour of the day
    for hh = 1:24
        % Create strings and such to define where the files are located on the host
        qFileName = [ftpBaseDir qYearStr(dd,:) '/' qDoyStr(dd,:) '/' ftpBaseFn qYearStr(dd,:) qMonthStr(dd,:) qDayStr(dd,:) '.' qHourStr(hh,:) ftpEndFn];
        % Get the file from Nasa's server
        disp(['Getting ' qFileName '...'])
        localFileName = mget(ftpObj,qFileName);
        % Create ncgeodataset object
        geo = ncdataset(localFileName{1});
        % Initialize a vector to hold the timestep's data (all variables) for entire domain
        domainData = nan(nLat, nLon, nVars);
        % Loop through each variable
        for vv = 1:nVars
            % Get the data for each variable
            domainData(:,:,vv) = squeeze(geo.data(allStrings{vv}));
        end
        % Loop through each site 
        for ss = 1:nSites
            % Extract the point data at each site from the domain data
            siteMetData = squeeze(domainData(latIdcs(ss), lonIdcs(ss), :));
            % Print the data. Hour is hh-1 because hours are listed from
            % 00:00 to 23:00.
            fprintf(fid(ss), [dateFmt varFmt], [qYears(dd) qMonths(dd) qDays(dd) (hh-1) 0 siteMetData']);
        end % loop through each site
    end % Loop through each hour of the day
    % Delete this day's directory and all data within
    rmdir(localDir, 's')
end % Loop through each day in the record
% -------------------------------------------------------------------------
% Clean up and close files
% Delete all data stored in this session
disp(['Cleaning up...'])
rmdir(localBaseDir, 's')
% Loop through each site and close the files
for ss = 1:nSites
    outFile = qNames{ss};
    fclose(fid(ss));
end
% Close the ftp connection
close(ftpObj)

end