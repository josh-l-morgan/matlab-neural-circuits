function[RGBs] = returnMaxFromOIF(CPN,stkChan)

if ~exist('stkChan','var')
stkChan = {[1 2 3] [4 5 6]};
end

useSlices = [];%[1:34];


%% Parse names
dCPN = dir([CPN '*.tif']);
inams = {dCPN.name};
chans = [];
slice = [];
for i  = 1:length(inams)
    
    nam = inams{i};
    
    c = regexp(nam,'C');
    d = regexp(nam,'.tif');
    z = regexp(nam,'Z');
    
    if isempty(z)
            chans(i) = str2num(nam(c+1:d-1));
            slice(i) = 1;
    else
    chans(i) = str2num(nam(c+1:z-1));
    slice(i) = str2num(nam(z+1:d-1));
    end
    
    
end

%% Choose slices and colors
numCol = max(chans);

if isempty(useSlices)
    slices = unique(slice);
else
    slices = useSlices;
end

%% Get max from each channel
testI = imread([CPN inams{1}]);
[ys xs zs] = size(testI);
IcCols = zeros(ys,xs,numCol,'uint8');
for c = 1:numCol

    Ic = zeros(ys,xs,length(slices));
    for i = 1:length(slices)

        targ = find((slice == slices(i)) & ( chans == c));
        if ~isempty(targ)
            Ic(:,:,i) = imread([CPN inams{targ}]);
        end

    end
    IcMax = double(max(Ic,[],3));
    IcCols(:,:,c) = IcMax*255/max(IcMax(:));
end


for i = 1:length(stkChan)
    RGBs{i} = zeros(ys,xs,3,'uint8');
    for ci = 1:length(stkChan{i})
        c = stkChan{i}(ci);
        if c<=numCol
            RGBs{i}(:,:,ci) = uint8(IcCols(:,:,c));
        end
    end
end












