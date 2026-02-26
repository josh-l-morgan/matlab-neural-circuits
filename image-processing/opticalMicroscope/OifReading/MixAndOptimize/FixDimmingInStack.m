
%Parietal eye
SPN = 'H:\FV1000\AxU\Anole6\parietalEye40x3D_full.oif.files\'

%IPL dendrites
%SPN = 'H:\FV1000\AxV\Anole5\fovea_13(Direct File Access Mode)3.oif.files\'

%%Strange fovea
SPN = 'H:\FV1000\AxV\Anole5\fovea(Direct File Access Mode)2.oif.files\'

%% spiny dendrite
SPN = 'H:\FV1000\AxV\Anole5_reslice\image_09.oif.files\'

%% Manual decisions
getZ = [];
getChan = [1 3];
putChan = [3 3 1];



%% Parse file names
dSPN = dir([SPN '*.tif']);
iNams = {dSPN.name}';
cs = zeros(length(iNams),1);
zs = zeros(length(iNams),1);
for i = 1:length(iNams)
    nam = iNams{i};
    c = str2num(nam(4:6));
    z = str2num(nam(8:10));
    cs(i) = c;
    zs(i) = z;
end

if isempty(getZ), getZ = [1:max(zs)];end

%% Read in raw data
ifo = imfinfo([SPN iNams{1}]);
maxC = length(getChan);
maxZ = length(getZ);
maxX = ifo.Width;
maxY = ifo.Height;
for c = 1:length(getChan)
    I{c} = zeros(maxY,maxX,length(getZ));
end


for c = 1:maxC
    isChan = find(cs == getChan(c));
    checkZs = zs(isChan);
    [a idx] = intersect(checkZs,getZ);
    getI = isChan(idx);
    for i = 1:length(getI)
        disp(sprintf('for chan %d of %d, reading image %d of %d',c,maxC,i,length(getI)))
        I{c}(:,:,zs(i)) = imread([SPN iNams{isChan(i)}]);
    end
end

%% Median filter
if 0
for c = 1:maxC
    for z = 1:maxZ
        disp(sprintf('median filtering channel %d of %d, section %d of %d',c,maxC,z,maxZ))
        I{c}(:,:,z) = medfilt2(I{c}(:,:,z),[2 2]);
    end
end
end

%% Find histograms
numVals = maxX * maxY;
low = round(numVals * .01);
high = round(numVals * .99);
atLow = zeros(maxZ,maxC);
atHigh = atLow;
for c = 1:maxC
    for z = 1: maxZ
        disp(sprintf('finding vals of c%d z%d',c,z))
        plane = I{c}(:,:,z);
        vals = sort(plane(:),'ascend');
        atLow(z,c) = vals(low);
        atHigh(z,c) = vals(high);
    end
end

plotCols = [0 0 1; 0 1 0; 1 0 0];
clf
sp1 = subplot(2,1,1);
sp1.NextPlot = "add";
sp2 = subplot(2,1,2);
sp2.NextPlot = "add";
for c = 1:maxC
    plot(sp1,atLow(:,c),'color',plotCols(c,:))
    plot(sp2,atHigh(:,c),'color',plotCols(c,:))
end

allLow = min(atLow,[],1)
allHigh = max(atHigh,[],1)


%% Manual decisions
refMax = [10 200; 760 1689]; % Channel 1 reference points


%% Adjust Brightness for depth
refRange = refMax;
refRange(2,:) = refRange(2,:) - allLow(1);
rise = (refRange(2,2)-refRange(2,1))/refRange(2,1); %fraction of max missing from lower sections
run = refRange(1,2)-refRange(1,1);
fixDepth = rise/run;


for z = 1:maxZ
    disp(sprintf('fixing depth brightness for plane %d of %d',z,maxZ))
    for c = 1:maxC
        newPlane = I{c}(:,:,z) - allLow(c); %set floor to 0
        newPlane = allLow(c) + newPlane * (1 + (maxZ-z) * fixDepth); % increase brightness of low sections
        I{c}(:,:,z) = newPlane;
    end
end


%% Adjust Channels
targetFloor = 10;
targetRange = 200;
for c = 1:maxC
    disp(sprintf('adjusting range of channel %d of %d',c,maxC))
    vals = sort(I{c}(:),'ascend');
    numVal = length(vals);
    lowVal = vals(round(numVal*.001));
    highVal = vals(round(numVal*.999));
    currentRange = highVal-lowVal;
    I{c} = ((I{c}-lowVal) * targetRange/currentRange) + targetFloor;
end


%% gamma
if 0
chanGamma = [1,1,1];
for c = 1:maxC
    I{c}(I{c}<0) = 0;
    I{c} = ((I{c}/255).^chanGamma(c)) * 255;
end


%%Adjust Channels again
targetFloor = 10;
targetRange = 200;
for c = 1:maxC
    disp(sprintf('again adjusting range of channel %d of %d',c,maxC))
    vals = sort(I{c}(:),'ascend');
    numVal = length(vals);
    lowVal = vals(round(numVal*.001));
    highVal = vals(round(numVal*.999));
    currentRange = highVal-lowVal;
    I{c} = ((I{c}-lowVal) * targetRange/currentRange) + targetFloor;
end
end

%% Create 8 bit images
Ic = zeros(maxY,maxX,maxC,maxZ,'uint8');
Ic(:,:,1,:) = I{2};
Ic(:,:,2,:) = I{2}*.5;
Ic(:,:,3,:) = I{1};


%% Show

clf
for i = 1:maxZ
    disp(sprintf('showing plane %d of %d',i,maxZ))
    image(uint8(Ic(:,:,:,i)));
    drawnow
    %pause
end

%% write

TPN = [SPN(1:end-1) '_Tweaked\'];
if ~exist(TPN,'dir'),mkdir(TPN);end
for i = 1:maxZ
    disp(sprintf('writing plane %d of %d',i,maxZ))
    fileName = sprintf('s_z%03f.tif');
    imwrite(Ic(:,:,:,i),[TPN fileName],'tif')
end
