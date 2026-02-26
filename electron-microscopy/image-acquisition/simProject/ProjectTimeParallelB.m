clear all

%% Data Set variables
sectionThickness = .025;
xyResolution = .004;
zLengths = 100: 10: 1000; %length in microns
xLengths = zLengths;%[ 1 1 1 1 1] * 300; %xy length in microns
yLengths = zLengths;%[ 1 1 1 1 1 ] * 300;

%% Tissue Prep Numbers
BlockPrepTime = 7 * 24 * 60 * 60;
SectionsCutPerSec = 0.02; %% Cutting speed
SectionsStainedPerSec = 100000000; %% Currently set large for en block staining. 
WaferConstructionTime = 30 * 60;
SectionsPerWafer = 200;

%% Imaging Prep Time
SectionMapTime = 60; %% Seconds required to map each section
WaferChangeTime = 15 * 60 ; % Seconds for changing Wafers

%% Imaging Time
PixelDwell = 500 * 10 ^ -9; %% seconds
PixNumX = 2000; %% X Dimension of subTile
PixNumY = 3000;
BeamNumber = 61;
BackLashTime = 0.001; % ?

Section2SectionTime = 1; % seconds between Planes
SettlingTime = 1;
MovingTime = 0.01;
AutoFocusTime = 6;
AutoFocusFrequency = 1;
AutoStigTime = 6;
AutoStigFrequency = 0.01;
Tile2TileSoftwareTime = 0;

RetakeFrequency = 0.01;

for i = 1:length(zLengths);
%% Tissue Size
xLength = xLengths(i);
yLength = yLengths(i);
zLength = zLengths(i);
xyPixels = (xLength / xyResolution) * (yLength / xyResolution);  % number of pixels per plane
numPlane = zLength / sectionThickness; % number of planes

%% SEM How long
SEMRawScanMinutesPerPlane= (xyPixels/SEMpixelsPerSecond)/60;
SEMminutesPerPlane =  SEMRawScanMinutesPerPlane +  ((xyPixels/SEMtileSize)*SEMsecondsBetweenTiles)/60;
SEMScanHours(i) = (SEMminutesPerPlane * numPlane)/60;

SEMSectioningHours = numPlane/SEMcutsPerHour + numPlane/SEMplanesStainedPerHour;
SEMhoursBetweenSections = (numPlane *SEMminBetweenPlanes)/60;
SEMSectionHours(i) = SEMSectioningHours + SEMhoursBetweenSections;

%% TEM How long
TEMRawScanMinutesPerPlane= (xyPixels/TEMpixelsPerSecond)/60;
TEMminutesPerPlane =  TEMRawScanMinutesPerPlane +  ((xyPixels/TEMtileSize)*TEMsecondsBetweenTiles)/60;
TEMScanHours(i) = (TEMminutesPerPlane * numPlane)/60;

TEMSectioningHours = numPlane/TEMcutsPerHour + numPlane/TEMplanesStainedPerHour;
TEMhoursBetweenSections = (numPlane *TEMminBetweenPlanes)/60;
TEMSectionHours(i) = TEMSectioningHours + TEMhoursBetweenSections;

end

%% Show data

plot(zLengths,(SEMScanHours/24 + SEMSectionHours/24),'b')
hold on
plot(zLengths,SEMSectionHours/24,'c')
%plot(zLengths,(TEMScanHours + TEMSectionHours)/24/30,'r')
%plot(zLengths,TEMSectionHours/24/30,'m')
hold off

% 
% scatter(zLengths,(SEMScanHours + SEMSectionHours)/24/30,'b')
% hold on
% scatter(zLengths,SEMSectionHours/24/30,'c')
% scatter(zLengths,(TEMScanHours + TEMSectionHours)/24/30,'r')
% scatter(zLengths,TEMSectionHours/24/30,'m')
% hold off


%ylim([ 0 36])

