clear all

%% Data Set Dimensions
sectionThickness = .025;
xyResolution = .004;
zLengths = 100: 10: 1000; %length in microns
xLengths = zLengths;%[ 1 1 1 1 1] * 300; %xy length in microns
yLengths = zLengths;%[ 1 1 1 1 1 ] * 300;

%% Tissue Prep Time Variables
BlockPrepTime = 7 * 24 * 60 * 60;
SectionsCutPerSec = 0.02; %% Cutting speed
SectionsStainedPerSec = 100000000; %% Currently set large for en block staining. 
WaferConstructionTime = 60 * 60; %half an hour per wafer
SectionsPerWafer = 200;

%% Imaging Prep Time Varibles
SectionMapTime = 60; %% Seconds required to map each section
WaferChangeTime = 20 * 60 ; % Seconds for changing Wafers

%% Tile Time Variables
PixelDwell = 100 * 10 ^ -9; %% seconds
PixNumX = 2592; %% X Dimension of subTile
PixNumY = 2131;
BeamNumber = 61;
BeamReturnTime = 1 * 10^ -10; % ?
TilePix = BeamNumber * PixNumX * PixNumY;

%% Tile to Tile Time Variables
Section2SectionTime = 1; % seconds between Planes
SettlingTime = 0.1;
MovingTime = 0.01;
AutoFocusTime = 3; % delay for tile scan due to autofocus
AutoFocusFrequency = .1; % number of times to autofocus / num tiles
AutoStigTime = 3;
AutoStigFrequency = 0.1;
Tile2TileSoftwareTime = 0;
Overhead = 0; % Mystery/Software time cost for each tile

RetakeFrequency = 0.01;
DownTime = 0.2; 

%% Calculate Project Times

%%Tissue Size
xyPixels = (xLengths / xyResolution) .* (yLengths / xyResolution);  % number of pixels per plane
Sections = zLengths / sectionThickness; % number of planes

%%TissuePrepTime
TissuePrepTime = BlockPrepTime ...
    + Sections / SectionsCutPerSec ...
    + Sections / SectionsStainedPerSec ...
    + WaferConstructionTime * Sections/SectionsPerWafer;

%%ImagePrepTime
ImagingPrepTime = SectionMapTime * Sections ...
    + WaferChangeTime;

%%ImagingTime
NumTiles = (fix(xyPixels / TilePix)+1) .* Sections; % Assuming Maximum efficiency ROI
%should be improved

TileTime = PixelDwell * PixNumX * PixNumY ...
    + BeamReturnTime * PixNumY ...
    + Overhead;

Tile2TileTime = SettlingTime + MovingTime ...
    + Tile2TileSoftwareTime ...
    + AutoFocusTime * AutoFocusFrequency ...
    + AutoStigTime * AutoStigFrequency;

RetakeTime = NumTiles * RetakeFrequency ...  
    * (TileTime + Tile2TileTime + Section2SectionTime ...
    + AutoFocusTime + AutoStigTime); %Assumes automated failure detection

ImagingTime = TileTime * NumTiles ...
    + (Section2SectionTime ) * Sections ...
    + Tile2TileTime * NumTiles...
    + RetakeTime;

ScopeTime = (ImagingTime + ImagingPrepTime) * ( 1 + DownTime );



%% Show Times
subplot(1,5,2:5)
sec2day = 1/60/60/24;
plot(zLengths,TissuePrepTime * sec2day,'r')
hold on
plot(zLengths,ImagingPrepTime * (1 + DownTime) * sec2day,'g')
plot(zLengths,(ScopeTime)* sec2day,'b')
hold off

subplot(1,5,1)
bar(TileTime + Tile2TileTime,'r')
hold on
bar(TileTime,'b')
hold off

%% Show data
