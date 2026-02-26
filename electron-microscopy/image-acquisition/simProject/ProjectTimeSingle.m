clear all

%% Cost
costPerTer = 100; %dollars per terrabyte


%% Data Set Dimensions
sectionThickness = .035;
xyResolution = .004;
zLengths = 300;%100: 10: 500; %length in microns
xLengths = 300;%zLengths;%./zLengths * 500; %xy length in microns
yLengths = 500;%zLengths;%./zLengths * 500;

%% Tissue Prep Time Variables
BlockPrepTime = 7 * 24 * 60 * 60;
SectionsCutPerSec = 0.01; %% Cutting speed
SectionsStainedPerSec = 100000000; %% Currently set large for en block staining. 
WaferConstructionTime = 60 * 60; %half an hour per wafer
SectionsPerWafer = 200;

%% Imaging Prep Time Varibles
SectionMapTime = 60; %% Seconds required to map each section
WaferChangeTime = 20 * 60 ; % Seconds for changing Wafers

%% Tile Time Variables
CycleTime = 0;  % If cycle time is known, Use Cycle time instead of Calculating TileTime
PixelDwell = 100 * 10 ^ -9; %% seconds
PixNumX = 16000; %% X Dimension of subTile
PixNumY = 16000;
BeamNumber = 1;
BeamReturnTime = 1 * 10^ -10; % ?
TilePix = BeamNumber * PixNumX * PixNumY;

%% Tile to Tile Time Variables
Section2SectionTime = 3; % seconds between Planes
SettlingTime = 2;
MovingTime = 0.1;
AutoFocusTime = 6; % delay for tile scan due to autofocus
AutoFocusFrequency =1; % number of times to autofocus / num tiles
AutoStigTime = 3;
AutoStigFrequency = 1;
Tile2TileSoftwareTime = 5;
Overhead = 0; % Mystery/Software time cost for each tile

RetakeFrequency = 0.01;
DownTime = 0.1; 

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

%%Cost
TileCost = (PixNumX * PixNumY)/(10^12)* costPerTer;
storeCost = NumTiles * TileCost;

TileTime = PixelDwell * PixNumX * PixNumY ...
    + BeamReturnTime * PixNumY ...
    + Overhead;

if CycleTime >0
    TileTime = CycleTime;
end

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
hold off
subplot(6,6,[2:6 8:12 14:18 20:24])
sec2day = 1/60/60/24/30;
plot(zLengths,TissuePrepTime * sec2day,'r')
hold on
plot(zLengths,ImagingPrepTime * (1 + DownTime) * sec2day,'g')
plot(zLengths,(ScopeTime)* sec2day,'b')
subplot(6,6,[32:36])
plot(zLengths,storeCost)
hold off

subplot(1,6,1)
bar(TileTime + Tile2TileTime,'r')
hold on
bar(TileTime,'b')
hold off

%% Show data
