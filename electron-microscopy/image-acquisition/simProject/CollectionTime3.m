clear all

%% Tissue variables
sectionThickness = .025;
xyResolution = .004;
zLengths = [ 1 1 1 1 1 ] * 200; %length in microns
xLengths = [ 1:5 ] * 100; %xy length in microns
yLengths = [ 1:5 ] * 100;

%% SEM numbers

SEMcutsPerHour = 240; %% Cutting speed
SEMplanesStainedPerHour = 100000; %% number of planes stained per day

SEMpixelsPerSecond = 8 * 10^6; %% Imaging speed

SEMminBetweenPlanes = 5; %% minutes between Planes
SEMtileSize = 16000 * 16000; %% size of SEM tile
SEMsecondsBetweenTiles = 30; %% seconds between tiles within a plane

%% TEM numbers

TEMcutsPerHour = 5; %% Cutting speed
TEMplanesStainedPerHour = 100; %% number of planes stained per day

TEMpixelsPerSecond = 12 * 10^6; %% Imaging speed
TEMminBetweenPlanes = 20; %% minutes between Planes
TEMtileSize = 5000 * 3000; %% size of SEM tile
TEMsecondsBetweenTiles = 0; %% seconds between tiles within a plane


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

plot(zLengths,(SEMScanHours + SEMSectionHours)/24/30,'b')
hold on
plot(zLengths,SEMSectionHours/24/30,'c')
plot(zLengths,(TEMScanHours + TEMSectionHours)/24/30,'r')
plot(zLengths,TEMSectionHours/24/30,'m')
hold off

% 
% scatter(zLengths,(SEMScanHours + SEMSectionHours)/24/30,'b')
% hold on
% scatter(zLengths,SEMSectionHours/24/30,'c')
% scatter(zLengths,(TEMScanHours + TEMSectionHours)/24/30,'r')
% scatter(zLengths,TEMSectionHours/24/30,'m')
% hold off


%ylim([ 0 36])

