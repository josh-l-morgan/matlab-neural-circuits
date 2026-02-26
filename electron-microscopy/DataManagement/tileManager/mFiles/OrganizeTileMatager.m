

%{



checkSectionOrdering = ?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Automatic quality checks

checkFileQual = basic quality check of single file

findBestQualityFromLog: use logfile to identify pull out quality ratings
for tiles
    input: logBook
    output: logQulas

checkSelectQuals: use checkFileQual to check quality of images identified
by findBestQualityFromLog
    input:  [logQuals] = findBestQualityFromLog 
    output: chkQual

checkUnLoggedQuals: check quality of images with no log file
    input: [logQuals] = findBestQualityFromLog
    output: [noLogQuals]


reviewQuals:  Decide on best tile and retakes based on chkQual and
noLogQuals


%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% Manual Quality Checks




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% presence of data

tileListFromUTSL: produce list of expected tiles based on UTSL
    input: UTSL
    output: utslTiles?

listDataLocations:  manually input list of data folders
    output: montageDirectories

tifMapMontageDirectories: describe location of all tiffs at select server
locations
    input: montageDirectories
    output: tifMap

checkTileProgress:  check if tiles from UTSL exist on server

lineCheckTiles:  Read bottom line of ever tile to see if data is present
    input: [tifMap] = tifMapMontageDirectories

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% format Data

ID2tile
sub2tile
tif2tile

makeIDs: Make single ID list from utsl map

