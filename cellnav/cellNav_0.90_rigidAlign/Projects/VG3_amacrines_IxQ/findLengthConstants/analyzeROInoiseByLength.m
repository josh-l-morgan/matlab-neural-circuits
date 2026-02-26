
%%Check for unassigned bipolar cells
%%Determine if polarity of Ca responses are consistent with other cells. Is
%%it SNR dependent.
%%See if changing masks changes bimodality or polarity spread.
%%Consolidate redundant ROIs


if 0
    global glob
    datFold = [glob.datDir 'Analysis\Data\preproc\'];
    SPN =datFold;
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'GOI.mat']);
end

clf

filterBySNR = 0;%0.96; %remove rois with low signal
filterByEdge = 1;
weightErrors = 0;
standardize = 0;



%% Get ROIs
roiCids = ptDat(:,3);
numRoi = length(roiCids);
goodRoiCid = zeros(numRoi,1);


%% Get calcium data
realCal = ROI.Polarity1;
SNR = ROI.qual;

roiNearEdge = [66 67 68 69 71 120 144 70 198 ]; % manually identified by jm 1/6/2022 (GOI?)

targetType = 'real'; %Keep original calcium data or change to test model
%%real = keep observed values
%%allMean = replace values with the mean
%%meanPlusNorm = use mean + normal noise * caNoise
%%evenRand = use random distribution evenly spread between -1 and 1
%%meanPlusNorm = use mean + normal noise * caNoise

caPol = realCal;
caStd = std(caPol);
caNoise = caStd;

meanCaPol = mean(caPol(:));
if strcmp(targetType,'real')
    caPol = realCal;
elseif strcmp(targetType,'allMean')
    caPol = randn(size(caPol))*0 + meanCaPol; % all values equal mean
elseif strcmp(targetType,'noisyMean')
    caPol = randn(size(caPol)) * caNoise + meanCaPol;
elseif strcmp(targetType,'evenRand')
    caPol = rand(size(caPol))*2-1; % even random distribution between -1 and 1
elseif strcmp(targetType,'scramble')
    caPol = randsample(caPol,length(caPol)); % even random distribution between -1 and 1
end

%%Correct for weird errors that create outliers
caPol = caPol * -1;
caPol(caPol>1) = 1;
caPol(caPol<-1) = -1;


SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');
minBip = 5;


%% Remove bad frames
frames = ptDat(:,1);
allFrames = unique(frames); %all
useFrames = allFrames;
%useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 2006]; %all
%useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 ]; %all but worst
%useFrames = [1005 1006 1007 1008 1009 1010 ]; % only first stack
%useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 ]; %all
removeFrames = setdiff(allFrames,useFrames);

%% Choose cells
vCids = SOI.cids;
vCids = [2 3 4 5 13 14];
goodCid = zeros(length(vCids),1);



%% Get weight

numOn = zeros(length(vCids),1);
numOff = zeros(length(vCids),1);
clear d caDif
d = zeros(numRoi,numRoi) + inf;
caDif = abs(caPol-caPol');
for v = 1:length(vCids)
    vCid = vCids(v); %Get cid of VG3 to run
    vTarg = find(SOI.cid == vCid);
    rTarg = find(runCids == vCid);
    isCid = find(roiCids==vCid);  % Find all grouped rois for the cid being run
    useNodes = SOI.closeNode(vTarg); % regtrieve list of skeleton nodes for rois on cid
    D = sms(rTarg).sm.skel2skel.linDist;
    %d = SOI.cell(vTarg).d(:,useNodes); % Get matrix of distiances from all synapse nodes (Y) to all functional roi nodes (X)
    d(isCid,isCid) = D(useNodes,useNodes);
    goodRoiCid(isCid) = 1;
    
end


%% Filter Rois
%useRoi = find(~isnan(sum(allPred,2)));

goodRoi = goodRoiCid;

if filterBySNR
    goodRoi = goodRoi & (SNR>=filterBySNR);
end
meanGoodRoi = mean(goodRoi)

if filterByEdge
goodRoi(roiNearEdge) = 0;
meanGoodRoi = mean(goodRoi)
end


if exist('useFrames')
    for i = 1:length(removeFrames)
        for g = 1:length(goodRoi)
            gFrames = frames(g);
            if sum(gFrames==removeFrames(i))
                goodRoi(g) = 0;
            end
        end
    end
end
meanGoodRoi = mean(goodRoi)

useRoi = find(goodRoi>0);
numUse = length(useRoi);



%% compare d and caDif

clf
uD = d(useRoi,useRoi);
uCaDif = caDif(useRoi,useRoi);
sameCell = uD<inf;

yRoi = repmat(useRoi,[1 length(useRoi)]);
xRoi = repmat(useRoi',[length(useRoi) 1]);
difRoi = yRoi>xRoi;

compRoi = sameCell & difRoi;
compX = xRoi(compRoi);
compY = yRoi(compRoi);

compRois = find(compRoi);
dListSame = uD(compRoi);
cListSame = uCaDif(compRoi);
%scatter(dListSame,cListSame)
%scatter(uD(uD<inf),uCaDif(uD<inf)))


%% Look at differences as specific distance
hRange = [0:1:150];
binWidth = 3;
meanDif = hRange * 0;
numDifs = hRange * 0;

for i = 1:length(hRange);

    isD = find((dListSame< (hRange(i)+binWidth/2) ) & ...
        (dListSame >= (hRange(i) - binWidth/2)));
    caHits = cListSame(isD);
    meanDif(i) = mean(caHits);
    numDifs(i) = length(caHits);

end

fewDifs = find(numDifs<50);
lastUse = min(fewDifs);

plot(hRange(1:fewDifs),meanDif(1:fewDifs));


%% Look at differences as range
subplot(2,1,1)
hRange = [0:1:150];
meanDifR = hRange * 0;
numDifsR = hRange * 0;

for i = 1:length(hRange);

    isD = find((dListSame<= (hRange(i)) )& (dListSame>0));
    caHits = cListSame(isD);
    meanDifR(i) = mean(caHits);
    numDifsR(i) = length(caHits);

end


plot(hRange,meanDifR);
ylim([0 .5])

%% find different close ROIs
subplot(2,1,2)
maxSamp = 3;
isClose = (dListSame < maxSamp);
[cSort idx] = sort(cListSame(isClose),'descend');
dSort = dListSame(isClose);
dSort = dSort(idx);

xSort = compX(isClose);
xSort = xSort(idx);
ySort = compY(isClose);
ySort = ySort(idx);
xySort = [xSort ySort]

binRad = 1;
sRange = [binRad:.1:maxSamp];
meanS = sRange * 0;
for i = 1:length(sRange);
    inBin = (dSort>=(sRange(i)-binRad)) & (dSort <(sRange(i)+binRad));
    meanS(i) = mean(cSort(inBin));
end


sFrames = frames(xySort);
uFrames = unique(frames);

subplot(2,1,2)
cla
hold on
fCol = hsv(length(uFrames));
for i = 1:length(uFrames);
    isFrame = sFrames == uFrames(i);
    scatter(dSort(isFrame(:,1)),cSort(isFrame(:,1)),150,...
        'markerfacecolor',fCol(i,:));
    scatter(dSort(isFrame(:,2)),cSort(isFrame(:,2)),...
        '^','markerfacecolor',fCol(i,:));
end

% scatter(dSort,cSort)
% hold on
plot(sRange,meanS,'linewidth',2);
hold off
ylim([0 1])

frames(xySort)

medCaDifClose = mean(cSort)
N = length(cSort)



