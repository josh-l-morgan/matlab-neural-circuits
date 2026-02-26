%% simulate the response properties of the VG3, 
% given a simplified version of its connectivity

%% basic params
fieldDim=1600;

%% create the empty field
blank=zeros(fieldDim,fieldDim);
bpcDensity=1;

%% create the bipolar cells locations

%start with a simple grid
gridLocs=fieldDim/(100/bpcDensity)/2:fieldDim/(100/bpcDensity):fieldDim;
BPClocs=nchoosek(gridLocs,2);
