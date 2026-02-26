

%%Run getSynapseNumberWithLocalExclusion.m to get cellCounts
runBipOrRgc = 1; %1 = RGCs
minDist = .01;
lc = 16; %Length constant default = 16


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


synV = getSynapsePropertyForType(sms);

if runBipOrRgc == 1
    cellCounts = getSynapseNumberWithLocalExclusion(sms,minDist,1);% RGCs
else
    cellCounts = getPreSynapseNumberWithLocalExclusion(sms,minDist,7);% Bipolars

end
NormalizeStratToOneEach = 0;



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

zBins = cellDat(1).strat(:,1);
binWidth = mean(abs(zBins(1:end-1) - zBins(2:end)));
onZ = zBins * 0; %record total on influence at depth
offZ = zBins * 0;

%% Correct cell names
for i = 1:length(cellDat)
    nam = cellDat(i).type;
    if strcmp(nam(1:2),'37')
        cellDat(i).type = '37';
    end
    if strcmp(nam(1:2),'7i')
        cellDat(i).type = '7i';
    end
end



%% colect data on cells

allPos = [];
allVoxPos = [];
allLengths = [];
allDepths = [];
vArborLengths = [];
clear onZv offZv
onZa = cell(length(zBins),1);
offZa = cell(length(zBins),1);
offWv = {};
onWv  ={};
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
    %Get average polarity at depth.
    vgcBindWidth = binWidth * 1;
    soiTarg = find(SOI.cids==cid);
     W = exp(-skelD/lc);
     offW = sum(W(SOI.cell(soiTarg).preSign==1,:),1);
     onW = sum(W(SOI.cell(soiTarg).preSign==2,:),1);
     offWv{v} = offW;
     onWv{v} = onW;
     onZv{v} = zBins *0;
     offZv{v} = zBins *0;
     for i = 1:length(zBins)
         isDeep = (nDepth1>= (zBins(i)-vgcBindWidth/2)) & ...
             (nDepth1< (zBins(i)+vgcBindWidth/2));
         %          onZ(i) = onZ(i) + sum(SOI.cell(soiTarg).sumOn(isDeep>0));
         %         offZ(i) = offZ(i) + sum(SOI.cell(soiTarg).sumOff(isDeep>0));
         onZv{v}(i) = sum(onW(isDeep>0));
         offZv{v}(i) = sum(offW(isDeep>0));

         onZa{i} = [onZa{i} onW(isDeep>0)];
         offZa{i} = [offZa{i} offW(isDeep>0)];

     end

end


%% Show pol
onZ = zBins *0;
offZ = zBins * 0;
for v = 1:length(onZv);
    clf, hold on
    plot(onZv{v},'r')
    plot(offZv{v},'b')
    onZ = onZ + onZv{v};
    offZ = offZ + offZv{v};
    drawnow

end
clf, hold on
plot(onZ,'r')
plot(offZ,'b')


%% Build RGC probability bins
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

%% Filter curvss, bcTypeZH, vgcZH


clf, hold on
maxInf = max([onZ(:); offZ(:)]);
plot(vgcZH/max(vgcZH),'k')
plot(onZ/maxInf,'r')
plot(offZ/maxInf,'b')


%%Filter
filtOrd = 11;
vgcZHf = medfilt1(vgcZH,filtOrd);
onZf = medfilt1(onZ,filtOrd);
offZf = medfilt1(offZ,filtOrd);
bcTypeZHf = medfilt1(bcTypeZH,filtOrd,0,1)
% 
% 
% w = gausswin(1,1);
% vgcZHf = conv(vgcZH,w,'same');
% onZf = conv(onZ,w,'same');
% offZf = conv(offZ,w,'same');
% bcTypeZHf = convn(bcTypeZH,w,'same')

maxInf = max([onZf(:); offZf(:)]);



onBiasZf = (onZf-offZf)./(onZf+offZf);
onBiasZf(isnan(onBiasZf)) = 0;


clf,hold on
plot(zBins,vgcZHf/max(vgcZHf),'k')
plot(zBins,onZf/maxInf,'g')
plot(zBins,offZf/maxInf,'r')
plot(zBins,onBiasZf,'b')





%% Calculate bias for RGC types


overLapsStruc = bcTypeZHf .* repmat(vgcZHf,[1 length(uEWTypes)]);
overLapsON = bcTypeZHf .* repmat(onZf,[1 length(uEWTypes)]);
overLapsOFF = bcTypeZHf .* repmat(offZf,[1 length(uEWTypes)]);
overLapsPol = (overLapsON-overLapsOFF)./(overLapsON + overLapsOFF);

overLapsSumPol = (sum(overLapsON,1)-sum(overLapsOFF,1))./(sum(overLapsON,1) + sum(overLapsOFF,1));


maxPol = max(abs(overLapsPol(:)));

%checkRGCLabels =  COI.checkRGCLabels;
checkRGCLabels = {'4on' '4i' '4ow' '37' '5ti' '63' '6sw'};


checkUEW = zeros(length(uEWTypes),1);
runUEW = [];
for i = 1:length(checkRGCLabels)
    targType =  find(strcmp(uEWTypes,checkRGCLabels{i}));
    if ~isempty(targType)
        runUEW = [runUEW targType];
    end
end

listPredPol = overLapsSumPol(runUEW);

clf,hold on
for i = 1:length(runUEW)
    plot(zBins,overLapsPol(:,runUEW(i))+i/10,'b')
end

%% run Monte

%%Calculate real type polarities by summing off and on influence for all
%%synapses
typePols = zeros(length(checkRGCLabels),1);
clear realCount
for i = 1:length(checkRGCLabels)
    isSubType  = strcmp(synV.subTypeNames, checkRGCLabels{i});
    isSubType = isSubType' & (synV.type == 1);
    targTypes = find(isSubType);
    offS = [];
    onS = [];
    c = 0;
    for f = 1:length(targTypes)
        vCid = synV.vCid(targTypes(f));
        v = find(useVgc==vCid);
        n = synV.synID(targTypes(f));
        if ~isempty(v)
            c = c+1;
            offS(c) = offWv{v}(n);
            onS(c) = onWv{v}(n);

        end
    end
    realCount(i) = c;
    typePols(i) = (sum(onS)-sum(offS))/((sum(onS)+sum(offS)));
    realSynOn{i} = onS;
    realSynOff{i} = offS;
    realSynPol{i} = (onS-offS)./(onS+offS);
end

%% hist real
clf
pickA = [1 2 3];
pickB = [7];
realPolA = cat(2,realSynPol{pickA});
realPolB = cat(2,realSynPol{pickB});

%wDallV{g1,2}

realPolAm = mean(realPolA,2);
realPolBm = mean(realPolB,2);
realPolDif = realPolBm - realPolAm;

meanA = mean(realPolAm);
seA = std(realPolAm)/sqrt(length(realPolAm));
[meanA-seA meanA+seA]


meanB = mean(realPolBm);
seB = std(realPolBm)/sqrt(length(realPolBm));
[meanB - seB meanB+seB]

hRange = [-1:.2:1];
realhA = hist(realPolA,hRange);
realhB = hist(realPolB,hRange);
bar(hRange,[realhA/sum(realhA);realhB/sum(realhB)],'barwidth',1.8)




%%Run Monte Carlo
numSyn = sum(realCount);
reps = 10000;

numStrat = size(overLapsStruc,1);

clear sPol
for i = 1:length(onZa)
    Ls(i) = length(onZa{i}); %get number of VG3 nodes in each stratification bin
end
for u = 1:length(runUEW);
    disp(sprintf('running %d of %d',u,length(runUEW)));
    pickNum = realCount(u);
    ovLap = overLapsStruc(:,runUEW(u));
    sPol{u} = zeros(reps,pickNum);
    for r = 1:reps
        randSamp = randsample(numStrat,pickNum,true,ovLap);
        for s = 1:length(randSamp)
            sOn{u}(r,s) = onZa{randSamp(s)}(ceil(rand*Ls(randSamp(s))));
            sOff{u}(r,s) = offZa{randSamp(s)}(ceil(rand*Ls(randSamp(s))));
        end
    end
    sPol{u} = (sum(sOn{u},2)-sum(sOff{u},2)) ./ (sum(sOn{u},2)+sum(sOff{u},2)) ;
end


%% Show results
checkNum = length(checkRGCLabels);

rHBounds = zeros(checkNum,4);
rSBounds = zeros(checkNum,4);
rHBounds2 = zeros(checkNum,4);
rSBounds2 = zeros(checkNum,4);
for u = 1:checkNum
    rHBounds(u,:) = bounds95(sPol{u},0.999);
    rSBounds(u,:) = bounds95(sPol{u},0.999);
    rHBounds2(u,:) = bounds95(sPol{u},0.50);
    rSBounds2(u,:) = bounds95(sPol{u},0.50);
end

% 
% clf
% subplot(1,1,1)
% hold on
% 
% bar(xRange ,rHBounds(:,3),'facecolor',[.5 .5 .5],'barwidth',.8,'edgecolor','k');
% bar(xRange ,rHBounds(:,1),'facecolor',[1 1 1],'barwidth',.8,'edgecolor','k');

% 
% x = 1;
% y = [.2 .3 .4 .5];
% w = .3;
% plot([x-w *.6 x+w *.6],[y(1) y(1)],'k')
% plot([x-w x+w],[y(2) y(2)],'k')
% plot([x-w x+w],[y(3) y(3)],'k')
% plot([x-w *.6 x+w *.6],[y(4) y(4)],'k')
% plot([x-w x-w],[y(2) y(3)],'k')
% plot([x+w x+w],[y(2) y(3)],'k')
% plot([x x],[y(3) y(4)],'k')
% plot([x x],[y(1) y(2)],'k')

clf
xRange = 1:length(realCount);

w = .25;
for i = 1:checkNum
    hold on
    b1 = bounds95(sPol{i},0.995);
    b2 = bounds95(sPol{i},0.5);
    b3 = bounds95(sPol{i},0.95);
    y = [b1(1) b2(1) b1(2) b2(3) b1(3)];
    x = i;
    xb = [x-w *.8 x+w *.8;x-w x+w; x-w x+w; x-w *.8 x+w *.8;x-w x-w;x+w x+w;x x;x x; x-w x+w]
    yb = [y(1) y(1); y(2) y(2);y(4) y(4); y(5) y(5); y(2) y(4); y(2) y(4); y(4) y(5); y(1) y(2); y(3) y(3)];
    plot(xb',yb','k','linewidth',2)
%     xc = [x-w *.7 x+w * .7; x-w *.7 x+w * .7]
%     yc = [b3(1) b3(1); b3(3) b3(3)];
%     plot(xc',yc','k','linewidth',2)
end

scatter(xRange,typePols,800,'markerfacecolor',[1 .2 .2],'markeredgecolor','k',...
    'markerfacealpha',.6);
xticks(1:checkNum)
xticklabels(checkRGCLabels)

if 0
    epsName = 'Z:\Active\morganLab\PUBLICATIONS\VG3\Figure\Draft4\Pics3\MontePolarityTypePrediction.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
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

%% Compare 4s to 6s

for u = 1:length(sOn)
   synPol{u} =  (sOn{u}-sOff{u})./(sOn{u}+sOff{u});

end


clf, hold on
aPol = synPol;
polA = cat(2,aPol{pickA});
polB = cat(2,aPol{pickB});

polAm = mean(polA,2);
polBm = mean(polB,2);
polDif = polBm - polAm;
sPolDif = sort(polDif);
mc95CI = [sPolDif(round(reps*0.025)) sPolDif(round(reps*0.975))];
median(sPolDif)



meanA = mean(polA);
seA = std(polA)/sqrt(length(polA));
[meanA-seA meanA+seA];

meanB = mean(polB);
seB = std(polB)/sqrt(length(polB));
[meanB - seB meanB+seB];


subplot(2,1,1)
hold on
realhA = hist(realPolA,hRange);
realhB = hist(realPolB,hRange);
bar(hRange,realhA/max(realhA),'b','barwidth',1,'facealpha',.4)
bar(hRange,realhB/max(realhB),'r','barwidth',1,'facealpha',.4)


subplot(2,1,2)
hold on
hA = hist(polA(:),hRange);
hB = hist(polB(:),hRange);
%bar(hRange,[hA/sum(hA);hB/sum(hB)],'barwidth',1.8)
% plot(hRange,hA/max(hA),'b')
% plot(hRange,hB/max(hB),'r')
bar(hRange,hA/max(hA),'b','barwidth',1,'facealpha',.4);%'edgecolor','none'
bar(hRange,hB/max(hB),'r','barwidth',1,'facealpha',.4);


bootDifCI(realPolA,realPolB,.95)
standardDifCI(realPolA,realPolB,.95)





if 0
    epsName = 'Z:\Active\morganLab\PUBLICATIONS\VG3\Figure\Draft4\Pics3\4v6r.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
end

%% Other
% 
% checkTypes = [12     17      25];
% realCount(checkTypes)
% 
% uEWTypes(checkTypes)
% realCount(checkTypes)
% N = sum(realCount)
% rHBounds(checkTypes,:)*N
% 
% rSBounds(1,:)*N










