

%%Run getSynapseNumberWithLocalExclusion.m to get cellCounts
runBipOrRgc = 0; %1 = RGCs
minDist = 10;%.1;









%% load data
if 0

    global glob tis

    figure
    SPN = [glob.datDir 'Analysis\Data\preproc\'];
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'GOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load([SPN 'COI.mat']);



    %%Load all sms
    smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
    clear sms
    roiCid = ptDat(:,3);
    runCids = unique(roiCid);%MOI.cids;
    for i = 1:length(runCids);
        cid = runCids(i);
        %%Get distances between nodes
        disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
        fileName = sprintf('sm_cid%d.mat',cid);
        useSM(i) = 1;
        %sm = load([smDir fileName],'skel2skel','nep','syn2skel','syn');
        load([smDir fileName]);
        sm.nep.swcS = nep2swc(sm.nep);
        sms(i).sm = sm;
    end
    useVgc = runCids(useSM>0);

    s = makeSynapseClassifyer(COI); %make structure describing types of synapses


end


useVgc =[ 2 3 4 5 13 14];
for i = 1:length(sms)
    smCids(i) = sms(i).sm.cid;
end


if runBipOrRgc == 1
    cellCounts = getSynapseNumberWithLocalExclusion(sms,minDist,1);% RGCs
else
    cellCounts = getPreSynapseNumberWithLocalExclusion(sms,minDist,7);% Bipolars

end
NormalizeStratToOneEach = 1;



%% Import Eyewire cells
loadBool=1;
%set the source directories
jsonDirBPC='Z:\Active\MorganLab\karlsRetina\eyewire_data\bpc\';
jsonDirRGC='Z:\Active\MorganLab\karlsRetina\eyewire_data\rgc\';
%get the file lists
jsonDirStructBPC=dir(jsonDirBPC);
jsonDirStructRGC=dir(jsonDirRGC);
%trim to just json files
jsonFileListBPC=string([]);
jsonFileListRGC=string([]);

for k=1:length(jsonDirStructBPC)
    curFileName=string(jsonDirStructBPC(k).name);
    if endsWith(curFileName,'json')
        jsonFileListBPC=[jsonFileListBPC; curFileName];
    end
end


for l=1:length(jsonDirStructRGC)
    curFileName=string(jsonDirStructRGC(l).name);
    if endsWith(curFileName,'json')
        jsonFileListRGC=[jsonFileListRGC; curFileName];
    end
end

%%Put everything into a structure
%load in bpc by type

%make a struct to hold all the bpc data
global cellDat;
cellDat=struct;
global ewIDList;
ewIDList=[];
cellIt=1;
%load bpc id, type, and stratDat
if runBipOrRgc
    useSubTypes = jsonFileListRGC;
    useSubDir = jsonDirRGC;
else
    useSubTypes = jsonFileListBPC;
    useSubDir = jsonDirBPC;
end

for k=1:length(useSubTypes)
    curFileName=[useSubDir+useSubTypes(k)];
    rawTxt=fileread(curFileName);
    curJ=jsondecode(rawTxt);
    for curCell=1:length(curJ)
        curCellDat=curJ(curCell);
        cellDat(cellIt).id=curCellDat.id;
        ewIDList=[ewIDList curCellDat.id];
        cellDat(cellIt).type=curCellDat.type;
        cellDat(cellIt).strat=curCellDat.stratification;
        cellIt=cellIt+1;
    end
end

%% Correct cell names
for i = 1:length(cellDat)
    nam = cellDat(i).type;
    if strcmp(nam(1:2),'37')
        cellDat(i).type = '37';
    end
end

%% colect data on cells

allPos = [];
allVoxPos = [];
allLengths = [];
allDepths = [];
vArborLengths = [];
for v = 1:length(useVgc)

    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    vc = find(smCids==cid,1);
    sm = sms(vc).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    nodeLengths = sm.nep.props.nodeLength;


    allVoxPos = cat(1,allVoxPos,sm.arbor.subs);
    allPos = cat(1,allPos,sm.nep.pos);
    allLengths = cat(1,allLengths(:),nodeLengths(:));


    [zGCL,zINL,nDepth1] = getIPLdepth(sm.nep.pos(:,3),sm.nep.pos(:,1),sm.nep.pos(:,2),[],[]);
    notCB = nDepth1 > 0.1;
    vArborLengths(v) = sum(nodeLengths(notCB));
    allDepths = cat(1,allDepths,nDepth1);

end

%% Build RGC probability bins
zBins = cellDat(1).strat(:,1);
binWidth = mean(abs(zBins(1:end-1) - zBins(2:end)));
allEWTypes = {cellDat.type};
uEWTypes = unique(allEWTypes);
uEWTypes = {uEWTypes{:}}
typeNum = length(uEWTypes);
bipCol = hsv(typeNum) * .75;
bipCol = flipud(bipCol);

bcTypeZH = zeros(length(zBins),length(uEWTypes));

for i = 1:length(cellDat)
    targ = find(strcmp(uEWTypes,cellDat(i).type));
    bcTypeZH(:,targ) = bcTypeZH(:,targ) + cellDat(i).strat(:,2);
end

%% Count hits

realCount = zeros(1,length(uEWTypes));
for i = 1:size(cellCounts,1)
    targ = find(strcmp(uEWTypes,cellCounts{i,1}));
    if ~isempty(targ)
        realCount(targ) = cellCounts{i,2};
    end
end



%% Remove primaries
for i = 1:typeNum
    typeH = bcTypeZH(:,i);
    typeThresh = max(typeH) * .1;
    typeH = typeH - typeThresh;
    typeH(typeH<0) = 0;
    bcTypeZH(:,i) = typeH;
end
bcTypeZH = bcTypeZH / max(bcTypeZH(:));


if NormalizeStratToOneEach % each RGC has equal occupation (summed) of IPL
    for i = 1:typeNum
        typeH = bcTypeZH(:,i);
        bcTypeZH(:,i) = typeH/sum(typeH);
    end
end


clf
hold on
for i = 1:typeNum
    typeH = bcTypeZH(:,i);
    % plot(typeH,'color',bipCol(i,:))
    plot(zBins,typeH,'color',bipCol(i,:));
end
legend( uEWTypes{:})


%% build VGC probability bins
clf
hold on
vgcZH = zeros(size(zBins));
vgcBindWidth = binWidth * 1;
for i = 1:length(zBins)
    isDeep = (allDepths>= (zBins(i)-vgcBindWidth/2)) & ...
        (allDepths< (zBins(i)+vgcBindWidth/2));
    vgcZH(i) = sum(allDepths(isDeep));
end

vgcZH = vgcZH/ max(vgcZH);

plot(zBins,vgcZH,'k')
legend(uEWTypes)
xlim([0 1])
ylim([0 1.1])

%% Find overlaps

overLaps = bcTypeZH .* repmat(vgcZH,[1 length(uEWTypes)]);

for i = 1:typeNum
    %plot(zBins,bcTypeZH(:,i),'color',bipCol(i,:))
    plot(zBins,overLaps(:,i),'color',bipCol(i,:))
end

typeP = sum(overLaps,1)/sum(overLaps(:));
legend({'VGC' uEWTypes{:}})
xlim([0 1])
ylim([0 max(overLaps(:) * 1.1)])

%% Run Monte
numSyn = sum(realCount);
reps = 10000;
randHR = zeros(reps,typeNum);
randHS = zeros(reps,typeNum);
for r = 1:reps
    randSamp = randsample(typeNum,numSyn,true,typeP);
    randH = histcounts(randSamp,1:typeNum+1);
    randHR(r,:) = randH/numSyn;
    randHS(r,:) = sort(randH/numSyn,'descend');
end

rHBounds = zeros(typeNum,4);
rSBounds = zeros(typeNum,4);
rHBounds2 = zeros(typeNum,4);
rSBounds2 = zeros(typeNum,4);
for i = 1:typeNum
    rHBounds(i,:) = bounds95(randHR(:,i),0.999);
    rSBounds(i,:) = bounds95(randHS(:,i),0.999);
    rHBounds2(i,:) = bounds95(randHR(:,i),0.50);
    rSBounds2(i,:) = bounds95(randHS(:,i),0.50);
end

realSort = sort(realCount,'descend')/sum(realCount);

clf
subplot(1,1,1)
hold on
xRange = 1:length(realSort);
w = .3;
% 
% bar(xRange ,rHBounds(:,3),'facecolor',[.5 .5 .5],'barwidth',.8,'edgecolor','k');
% bar(xRange ,rHBounds(:,1),'facecolor',[1 1 1],'barwidth',.8,'edgecolor','k');

% 
% x = 1;
% y = [.2 .3 .4 .5];
% 
% plot([x-w *.6 x+w *.6],[y(1) y(1)],'k')
% plot([x-w x+w],[y(2) y(2)],'k')
% plot([x-w x+w],[y(3) y(3)],'k')
% plot([x-w *.6 x+w *.6],[y(4) y(4)],'k')
% plot([x-w x-w],[y(2) y(3)],'k')
% plot([x+w x+w],[y(2) y(3)],'k')
% plot([x x],[y(3) y(4)],'k')
% plot([x x],[y(1) y(2)],'k')

clf

bar(xRange,realCount/sum(realCount),'facecolor',[1 .4 .4],'barwidth',.8,'edgecolor','k');
xticks(1:typeNum)
xticklabels(uEWTypes)
for i = 1:typeNum
    hold on
    b1 = bounds95(randHR(:,i),0.999);
    b2 = bounds95(randHR(:,i),0.5);
    y = [b1(1) b2(1) b1(2) b2(3) b1(3)];
    x = i;
    xb = [x-w *.6 x+w *.6;x-w x+w; x-w x+w; x-w *.6 x+w *.6;x-w x-w;x+w x+w;x x;x x; x-w x+w]
    yb = [y(1) y(1); y(2) y(2);y(4) y(4); y(5) y(5); y(2) y(4); y(2) y(4); y(4) y(5); y(1) y(2); y(3) y(3)];
    plot(xb',yb','k','linewidth',2)
end




% 
% 
% subplot(2,1,2)
% hold on
% bar(rSBounds(:,3),'facecolor',[.5 .5 .5]);
% bar(rSBounds(:,1),'facecolor',[1 1 1]);
% bar(realSort/sum(realSort),'facecolor',[1 0 0],'barwidth',.5);
% xticks(1:typeNum)
% xticklabels(1:typeNum)




checkTypes = [12     ];

uEWTypes(checkTypes)
realCount(checkTypes)
N = sum(realCount)
rHBounds(checkTypes,:)*N

 bounds95(randHR(:,checkTypes),0.95)*N

rSBounds(1,:)*N










