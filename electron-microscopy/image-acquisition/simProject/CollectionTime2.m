clear all

%% Tissue variables
sectionThickness = .040;
xyResolution = .004;
zLengths = 100: 100:500; %length in microns
xyLengths = [500 300]; %xy length in microns

%% SEM numbers

SEMcutsPerHour = 60; %% Cutting speed
SEMplanesStainedPerHour = 100; %% number of planes stained per day

SEMpixelsPerSecond = 2 * 10^6; %% Imaging speed

SEMminBetweenPlanes = 5; %% minutes between Planes
SEMtileSize = 5000 * 3000; %% size of SEM tile
SEMsecondsBetweenTiles = 1; %% seconds between tiles within a plane

%% TEM numbers

TEMcutsPerHour = 5; %% Cutting speed
TEMplanesStainedPerHour = 100; %% number of planes stained per day

TEMpixelsPerSecond = 12 * 10^6; %% Imaging speed
TEMminBetweenPlanes = 20; %% minutes between Planes
TEMtileSize = 5000 * 3000; %% size of SEM tile
TEMsecondsBetweenTiles = 1; %% seconds between tiles within a plane


for i = 1:length(zLengths);
%% Tissue Size
xyLength = [xyLengths(1) xyLengths(2)];
zLength = zLengths(i);
xyPixels = prod(xyLength ./ xyResolution);  % number of pixels per plane
numPlane = zLength / sectionThickness; % number of planes

%% SEM How long
SEMSectioningHours(i) = numPlane/SEMcutsPerHour + numPlane/SEMplanesStainedPerHour;
SEMRawScanMinutesPerPlane= (xyPixels/SEMpixelsPerSecond)/60;
SEMminutesPerPlane =  SEMRawScanMinutesPerPlane +  ((xyPixels/SEMtileSize)*SEMsecondsBetweenTiles)/60;
SEMhoursPerBlock(i) = (SEMminutesPerPlane * numPlane + numPlane * SEMminBetweenPlanes)/60;


%% TEM How long

TEMSectioningHours(i) = numPlane/TEMcutsPerHour + numPlane/TEMplanesStainedPerHour;
TEMRawScanMinutesPerPlane = (xyPixels/TEMpixelsPerSecond)/60;
TEMminutesPerPlane = TEMRawScanMinutesPerPlane +  ((xyPixels/TEMtileSize)*TEMsecondsBetweenTiles)/60;
TEMhoursPerBlock(i) = (TEMminutesPerPlane/60) * numPlane + numPlane * TEMminBetweenPlanes/60;


end

%% Show data

plot(zLengths,(SEMhoursPerBlock + SEMSectioningHours)/24/30,'b')
hold on
plot(zLengths,SEMSectioningHours/24/30,'c')
plot(zLengths,(TEMhoursPerBlock + TEMSectioningHours)/24/30,'r')
plot(zLengths,TEMSectioningHours/24/30,'m')
hold off



%ylim([ 0 36])

