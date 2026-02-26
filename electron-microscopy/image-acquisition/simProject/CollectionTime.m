clear all

%% Tissue variables
sectionThickness = .040;
xyResolution = .004;
zLengths = 300; %100: 100:500; %length in microns
xyLengths = 300; %xy length in microns

%% SEM numbers

SEMplanesPerHour = 600; %% Cutting speed
SEMplanesStainedPerHour = 10000000; %% number of planes stained per day

SEMpixelsPerSecond = 1 * 10^6; %% Imaging speed
SEMminBetweenPlanes = .1; %% minutes between Planes
SEMtileSize = 5000 * 5000; %% size of SEM tile
SEMsecondsBetweenTiles = 1; %% seconds between tiles within a plane

%% TEM numbers

TEMplanesPerHour = 60; %% Cutting speed
TEMplanesStainedPerHour = 100; %% number of planes stained per day

TEMpixelsPerSecond = 100 * 10^6; %% Imaging speed
TEMminBetweenPlanes = .1; %% minutes between Planes
TEMtileSize = 5000 * 5000; %% size of SEM tile
TEMsecondsBetweenTiles = 1; %% seconds between tiles within a plane


for i = 1:length(zLengths);
%% Tissue Size
xyLength = [xyLengths(1) xyLengths(1)];
zLength = zLengths(i);
xyPixels = prod(xyLength ./ xyResolution);  % number of pixels per plane
numPlane = zLength / sectionThickness; % number of planes

%% SEM How long
SEMSectioningHours = numPlane/SEMplanesPerHour + numPlane/SEMplanesStainedPerHour;
SEMminutesPerPlane = ((xyPixels/SEMpixelsPerSecond) +  (xyPixels/SEMtileSize)*SEMsecondsBetweenTiles)/60;
SEMhoursPerBlock = (SEMminutesPerPlane/60) * numPlane + numPlane * SEMminBetweenPlanes/60;

SEMRawDaysPerPlane = xyPixels/SEMpixelsPerSecond/60/60/24;
SEMTotalHours = SEMSectioningHours + SEMhoursPerBlock;
SEMTotalMonths(i) = SEMTotalHours/12/30;

%% TEM How long

TEMSectioningHours = numPlane/TEMplanesPerHour + numPlane/TEMplanesStainedPerHour;
TEMminutesPerPlane = ((xyPixels/TEMpixelsPerSecond) +  (xyPixels/TEMtileSize)*TEMsecondsBetweenTiles)/60;
TEMhoursPerBlock = (TEMminutesPerPlane/60) * numPlane + numPlane * TEMminBetweenPlanes/60;

TEMRawDaysPerPlane = xyPixels/TEMpixelsPerSecond/60/60/24;
TEMTotalHours = TEMSectioningHours + TEMhoursPerBlock;
TEMTotalMonths(i) = TEMTotalHours/12/30;



end

SEMTotalMonths(length(SEMTotalMonths))
TEMTotalMonths(length(TEMTotalMonths))

plot(zLengths,SEMTotalMonths,'b')
hold on
plot(zLengths,TEMTotalMonths,'r')
hold off
