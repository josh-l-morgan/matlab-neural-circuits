
clear all
nVcName = 'nVc.mat';
%nVcName = 'nVcT.mat';

doSmooth = 1;

f = figure;


%% Set options for what to compare
vCids = [2 3 4 5 13 14];
%vCids = [20];

if doSmooth
    nVcName = 'nVcSm.mat';
end


%%Define Must Be conditions for each variable.  Empty = no filter
%%Filtering by x labels so most values are in log10

showType = 3
if showType ==1 % no filter
    cMB = [];   eMB = [];   iMB = [];   osMB = [];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = []; osMin = [];
    dimOrder = [2 3 1 4  ]; % define or comment out

elseif showType ==2 % no inhibition
    cMB = [];   eMB = [];   iMB = [];   osMB = [];
    cMax = [];  eMax = [0];  iMax = [-5];  osMax = [];
    cMin = []; eMin= []; iMin = []; osMin = [];
    dimOrder = [1 2 3 4 ]; % define or comment out

elseif showType == 3; % Best Realistic
    cMB = [-4];   eMB = [-3];   iMB = [-3.65];   osMB = [.5];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = []; osMin = [];
    dimOrder = [2 3 1 4   ]; % define or comment out
elseif showType == 4; % Realistic conductance, no inhibition
    cMB = [];   eMB = [-3];   iMB = [];   osMB = [];
    cMax = [0];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = [-8]; osMin = [];
    dimOrder = [1 4 2 3 ]; % define or comment out
elseif showType == 5; % field of inhib
    cMB = [-6];   eMB = [-3];   iMB = [];   osMB = [];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = [-9]; osMin = [];
    dimOrder = [3 4 1 2]; % define or comment out
elseif showType == 6; % field of inhib
    cMB = [-4];   eMB = [];   iMB = [];   osMB = [];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = [-9]; osMin = [];
    dimOrder = [2 3 1 4]; % define or comment out
elseif showType ==7 % PUBLISH? nice map at c = -4
    cMB = [-4];   eMB = [-3];   iMB = [];   osMB = [];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = [-6]; osMin = [];
    dimOrder = [3 4 1 2 ]; % define or comment out
elseif showType ==8 % PUBLISH? vary synapse strength at c = -4, os fixed at 0.5
    cMB = [-4];   eMB = [];   iMB = [];   osMB = [.5];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = [-8]; osMin = [];
    dimOrder = [2 3 1 4 ]; % define or comment out
elseif showType ==9 % PUBLISH? vary synapse strength at c = -4, os fixed at 0.5
    cMB = [-6];   eMB = [];   iMB = [];   osMB = [1];
    cMax = [];  eMax = [];  iMax = [];  osMax = [];
    cMin = []; eMin= []; iMin = [-8]; osMin = [];
    dimOrder = [2 3 1 4 ]; % define or comment out
end


%% Set options for filtering functional ROIs
filterBySNR = 0.96; %remove rois with low signal
filterByEdge = 1;
weightErrors = 0;
standardize = 0;
weightRois = 1;


%% Load up cell data (non neuron)

SPN = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\'
load([SPN 'ptDat.mat'])
load([SPN 'GOI.mat']);
load([SPN nVcName])

dat = nVc.dat;
caPol = nVc.caPol;
nPol = nVc.nPol;
roiCids = nVc.roiCids;
goodRoiCid = nVc.goodRoiCid;
gChecked = nVc.gChecked;

SNR = GOI.qual;

useV = unique(roiCids(gChecked>0));
numCon = size(nPol,2);


allNn2Roi = cat(1,nVc.nn2Roi{:});
uNn2Roi = unique(allNn2Roi);
goodRoiCid(uNn2Roi) = 1;



%% Filter conditions
useC = ones(numCon,1,"logical");

mB = {cMB eMB iMB osMB}
for i = 1:length(mB)
    cMB = mB{i};
    if ~isempty(cMB)
        isOK = useC * 0;
        for c = 1:length(cMB)
            isOK = isOK | (dat(i).x == cMB(c));
        end
        useC = useC & isOK;
    end
end


mB = {cMax eMax iMax osMax}
for i = 1:length(mB)
    cMB = mB{i};
    if ~isempty(cMB)
        useC = useC & (dat(i).x <= cMB);
    end
end

mB = {cMin eMin iMin osMin}
for i = 1:length(mB)
    cMB = mB{i};
    if ~isempty(cMB)
        useC = useC & (dat(i).x >= cMB);
    end
end

useCfrac = mean(useC);
disp(sprintf('using %0.1f%% conditions',useCfrac * 100))

for c = 1:length(dat);
    dat(c).x = dat(c).x(useC);
    dat(c).r = dat(c).r(useC);
    dat(c).ux = unique(dat(c).x);
    dat(c).u = unique(dat(c).r);
    dat(c).num = length(dat(c).u);
end

nPol = nPol(:,useC);

numCon = size(nPol,2);


%% Remove bad frames
frames = ptDat(:,1);
allFrames = unique(frames); %all
%useFrames = allFrames;
%useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 2006]; %all
useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 ]; %all but worst
%useFrames = [1005 1006 1007 1008 1009 1010 ]; % only first stack
removeFrames = setdiff(allFrames,useFrames);

roiNearEdge = [66 67 68 69 71 120 144 70 198 ]; % manually identified by jm 1/6/2022 (GOI?)


%% Filter Rois
%useRoi = find(~isnan(sum(allPred,2)));
goodNNroi = goodRoiCid & gChecked;

goodRoi = goodRoiCid;
goodRoi = goodRoi & goodNNroi;

if filterBySNR
    goodRoi = goodRoi & (SNR>=filterBySNR);
end
meanGoodRoi = mean(goodRoi)


goodRoi = goodRoi & (abs(caPol)<=1);

if filterByEdge
    for r = 1:length(GOI.roiID)
        hit = intersect(GOI.roiID{r},roiNearEdge);
        if ~isempty(hit)
            goodRoi(r) = 0;
        end
    end
end
meanGoodRoi = mean(goodRoi)

if exist('useFrames')
    for i = 1:length(removeFrames)
        for g = 1:length(GOI.roiID)
            grIDs = GOI.roiID{g};
            gFrames = frames(grIDs);
            if sum(gFrames==removeFrames(i))
                goodRoi(g) = 0;
            end
        end
    end
end
meanGoodRoi = mean(goodRoi)

useRoi = find(goodRoi>0);
numUse = length(useRoi);

uC = caPol(useRoi);
uP = nPol(useRoi,:);
uV = roiCids(useRoi);


%% weight rois by polarity
if weightRois
    uC(uC<-1.1) = -1.1;
    uC(uC>1.1) = 1.1;
    caDif = abs(uC - uC');

    c = .05;
    x = caDif;
    caGau = exp( -1 * x.^2/ (2 * c^2));
    caSim = sum(caGau,2);
    gWeights = 1./caSim;

    caRange = [-1:.1:1];
    bin = .1;
    clear caWeight caCount;
    for i = 1:length(caRange)
        hit = (uC>=(caRange(i)-bin)) & (uC<=caRange(i)+bin);
        caCount(i) = mean(hit);
        caWeight(i) = sum(hit.*gWeights)/sum(gWeights);
    end
    clf, hold on
    plot(caRange,caCount,'r')
    plot(caRange,caWeight,'b')
else
    gWeights = uC*0+1;
end

%% Standardize polarities
uCN = uC;
if standardize
    uCN = uCN - mean(uCN);
    %uCN = uCN./std(uCN,1);
end
uPN = uP;
if standardize
    uPN = uPN - repmat(mean(uPN,2),[1 size(uPN,2)]);
    %uPN = uPN ./ repmat(std(uPN,1,2),[1 size(uPN,2)]);
end

uCNmat = repmat(uCN',[size(uPN,1) 1]);


%% find differencs
clear useError binError errorOfBins errorOfMeans errorOfVars rms predCor meanError meanErrorW

useWeights = gWeights;
for c = 1:numCon
    uPNs = uPN(:,c);
    meanError(c) = mean(abs(uPNs-uCN));
    meanErrorW(c) = sum((abs(uPNs-uCN).*useWeights))/sum(useWeights);
    binError(c) = mean(abs((uPNs>0)-(uCN>0)));
    errorOfBin(c) = abs(mean(uPNs>0)-mean(uCN>0));
    errorOfMeans(c) = abs(mean(uPNs) - mean(uCN));
    errorOfVars(c) = abs(var(uPNs) - var(uCN));
    rms(c) = sqrt(mean((uPNs-uCN).^2));
    cc = corrcoef(uPNs,uCN) ;
    predCor(c) = cc(1,2);
    %useError(c) = cc(1,2);
    %useError(c) = mean(sqrt((uP-uC).^2));
end
useError = meanErrorW;

%%Find error for each cell
clear meanErrorV rmsV predCorV errorOfMeansV errorOfBinV binErrorV errorOfVarsV
for v = 1:length(vCids)
    for c = 1:numCon
        isV = (uV == vCids(v));

        uPNs = uPN(:,c);
        meanErrorV{v}(c) = mean(abs(uPNs(isV)-uCN(isV)));
        binErrorV{v}(c) = mean(abs((uPNs(isV)>0)-(uCN(isV)>0)));
        errorOfBinV{v}(c) = abs(mean(uPNs(isV)>0)-mean(uCN(isV)>0));
        errorOfMeansV{v}(c) = abs(mean(uPNs(isV)) - mean(uCN(isV)));
        errorOfVarsV{v}(c) = abs(var(uPNs(isV)) - var(uCN(isV)));
        rmsV{v}(c) = sqrt(mean((uPNs(isV)-uCN(isV)).^2));
        cc = corrcoef(uPNs(isV),uCN(isV));
        predCorV{v}(c) = cc(1,2);
        %useError(c) = cc(1,2);
        %useError(c) = mean(sqrt((uP-uC).^2));
    end
end


errorMat = zeros(dat(1).num,dat(2).num,dat(3).num,dat(4).num);
idMat = cell(dat(1).num,dat(2).num,dat(3).num,dat(4).num);
for id1 = 1:dat(1).num
    for id2 = 1:dat(2).num
        for id3 = 1:dat(3).num
            for id4 = 1:dat(4).num
                targ = (dat(1).r == dat(1).u(id1)) & (dat(2).r == dat(2).u(id2)) & ...
                    (dat(3).r == dat(3).u(id3)) & (dat(4).r == dat(4).u(id4));

                idMat{id1,id2,id3,id4} = find(targ);
                sME = useError(targ);
                errorMat(id1,id2,id3,id4) = mean(sME);
                if length(sME)>1
                    disp('found multiple results for combination of conditions')
                end
            end
        end
    end
end



%% Show differences
if ~exist('dimOrder','var')
    sd = [1 2 3 4];
    if showType == 1
        sd = [ 1 2 3 4];
    elseif showType == 2
        sd = [3 4 1 2];
    elseif showType == 3
        sd = [1 2 3 4];
    end
else
    sd = dimOrder;
end

for d = 1:length(dat)

    dat(d).x(dat(d).x == -Inf) = -10;
    dat(d).ux(dat(d).ux == -Inf) = -10;
end
d1 = dat(sd(1));
d2 = dat(sd(2));
d3 = dat(sd(3));
d4 = dat(sd(4));
em = permute(errorMat,sd);


clf

prop = useError;
prop = prop-min(prop);
prop = round(prop * 99/max(prop)+1);
prop(isnan(prop)) = 1;

propMat = em;
propMat = propMat-min(propMat(:));
propMat = round(propMat .* 99/max(propMat(:))+1);

pCol = jet(100);
[mX mY] = meshgrid(d2.ux,d1.ux);


minError = min(useError);
bInd = find(useError==minError,1);
bX1 = d1.x(bInd);
bX2 = d2.x(bInd);
bX3 = d3.x(bInd);
bX4 = d4.x(bInd);

b1 = find(d1.ux == bX1);
b2 = find(d2.ux == bX2);
b3 = find(d3.ux == bX3);
b4 = find(d4.ux == bX4);
disp(sprintf('min error %0.4f',minError))
bestStr = sprintf('best %s = %0.2f, %s = %0.2f, %s = %0.2f, %s = %0.2f,',...
    d1.lab,bX1,d2.lab,bX2,d3.lab,bX3,d4.lab,bX4);
disp(bestStr)

bestPlane = mX * 0 + minError;

if 0 % show each plane

    for dC = 1:d4.num
        for dR = 1:d3.num
            sp = subplot(d3.num,d4.num,(dR-1)*d4.num+dC); cla(sp); hold on
            %sp = subplot(1,1,1); cla(sp); hold on

            isCond = find((d4.x == d4.ux(dC) ) & ( d3.x == d3.ux(dR)));
            isCondX1 = d1.x(isCond);
            isCondX2 = d2.x(isCond);
            isCondY = useError(isCond);


            %plot3(d1.x,d2.x,em(:,tC),'k')
            %m = mesh(mg,em(:,:,dR,dC));
            emS = em(:,:,dR,dC);
            pM = propMat(:,:,dR,dC);
            isNum = ~isnan(emS);
            suM = surf(mX,mY,bestPlane,'facecolor','flat','faceAlpha',.2);

            vq = griddata(mX(isNum),mY(isNum),emS(isNum),mX,mY);
            if ~isempty(vq)
                vProp = vq(~isnan(vq));
                vProp = vProp-min(vProp);
                vProp = round(vProp .* 99/max(vProp(:))+1);
                su = surf(mX,mY,vq,'facecolor','interp');
                su.CData(~isnan(vq)) =vProp;
            end

            % sc = scatter3(mX(isNum),mY(isNum),emS(isNum),'k','filled');
            % sc.CData =pCol(pM(isNum),:);

            sc = scatter3(isCondX2,isCondX1,isCondY,'k','filled');
            sc.CData =pCol(prop(isCond),:);
            %sc = scatter3(d2.x,d1.x,useError,'k','filled');
            %sc.CData =pCol(prop,:);
            zlim([0 max(useError)])
            ylabel(d1.lab)
            xlabel(d2.lab)
            zlabel('error')
            view(sp1,0,90)
            title(sprintf('%s = %0.2f, %s = %0.2f',d3.lab,d3.ux(dR),d4.lab,d4.ux(dC)));

            drawnow
            if (d4.ux(dC) == bX4) & (d3.ux(dR) == bX3)
                scatter3(d2.x(bInd),d1.x(bInd),minError,100,'o','r','linewidth',3)
            end
        end
    end


end

%% Show best plane
sp1 = subplot(1,2,1),cla(sp1),hold on

dC = find(d4.ux == bX4);
dR = find(d3.ux == bX3);
isCond = find((d4.x == bX4 ) & ( d3.x == bX3));
isCondX1 = d1.x(isCond);
isCondX2 = d2.x(isCond);
isCondY = useError(isCond);


%plot3(d1.x,d2.x,em(:,tC),'k')
%m = mesh(mg,em(:,:,dR,dC));
emS = em(:,:,dR,dC);
if (size(emS,1)>1) & (size(emS,2)>1)
    pM = propMat(:,:,dR,dC);
    isNum = ~isnan(emS);
    suM = surf(mX,mY,bestPlane,'facecolor','flat','faceAlpha',.2);

    vq = griddata(mX(isNum),mY(isNum),emS(isNum),mX,mY);
    if ~isempty(vq)
        vProp = vq(~isnan(vq));
        vProp = vProp-min(vProp);
        vProp = round((vProp .* 99 / max(vProp(:)))+1);
        su = surf(mX,mY,vq,'facecolor','interp');
        %su.CDataMapping = 'direct';
        su.CData(~isnan(vq)) =vProp;
        %colormap(sp1,jet(100))
        %colorbar(sp1)
    end
end
% sc = scatter3(mX(isNum),mY(isNum),emS(isNum),'k','filled');
% sc.CData =pCol(pM(isNum),:);
%scK = scatter3(isCondX2,isCondX1,isCondY,'k','filled','MarkerEdgeColor','k','linewidth',3);
sc = scatter3(isCondX2,isCondX1,isCondY,'k','filled','MarkerEdgeColor','w');
sc.CData =pCol(prop(isCond),:);
%sc = scatter3(d2.x,d1.x,useError,'k','filled');
%sc.CData =pCol(prop,:);
zlim([0 max(useError)])
ylabel(d1.lab)
xlabel(d2.lab)
zlabel('error')
view(-8,20)
title(sprintf('%s = %0.2f, %s = %0.2f',d3.lab,d3.ux(dR),d4.lab,d4.ux(dC)));

drawnow
scatter3(sp1,bX2,bX1,minError,100,'o','r','linewidth',3)



%% Match errors for each cell
%%error measures = useError errorOfMeans errorOfVars rms predCor errorOfBin binError

errorMeasureV = meanErrorV;

sp2 = subplot(1,2,2); cla(sp2),hold on

plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
plot([-1 1],[-1 1],'k')
axis 'equal'

vs = unique(uV);
vCol = jet(length(vs));
bCStr = sprintf('%s    %s    %s   %s',d1.lab(1:2), d2.lab(1:2), d3.lab(1:2), d4.lab(1:2));
for v = 1:length(useV)
    isV = (uV == useV(v));
    [sortError idxE] = sort(errorMeasureV{v},'ascend');
    bInd = idxE(1)
    bestScat = uPN(isV,bInd);
    scB = scatter(uCN(isV),bestScat,50,'markerfacecol',vCol(v,:),'markeredgealph',0,...
        'markerfacealpha',1);
    bCStr = cat(2,bCStr,...
        sprintf('\n%0.2f %0.2f %0.2f %0.2f',...
        d1.x(bInd),d2.x(bInd),d3.x(bInd),d4.x(bInd)));
end
title('Best fit each cell')

%% Match errors for all cells
sp2 = subplot(1,2,2); cla(sp2),hold on

errorMeasure = useError;
[sortError idxE] = sort(errorMeasure,'ascend');
bInd = idxE(1);

bestScat = uPN(:,idxE(1));

plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
plot([-1 1],[-1 1],'k')
axis 'equal'

vs = unique(uV);
vCol = jet(length(vs));
for v = 1:length(vs)
    isV = (uV == vs(v));
    scB = scatter(uCN(isV),bestScat(isV),50,'markerfacecol',vCol(v,:),'markeredgealph',0,...
        'markerfacealpha',1);
end


bestStr = sprintf('best %s = %0.2f, %s = %0.2f, %s = %0.2f, %s = %0.2f,',...
    d1.lab,bX1,d2.lab,bX2,d3.lab,bX3,d4.lab,bX4);
bestStrS = sprintf('best %s = %0.2f, %s = %0.2f, %s = %0.2f, %s = %0.2f,',...
    d1.lab(1:2),bX1,d2.lab(1:2),bX2,d3.lab(1:2),bX3,d4.lab(1:2),bX4);
disp(bestStr)
errorStr = sprintf('min error %0.4f, cc = %0.3f, bin = %0.3f',minError,predCor(bInd),(1-binError(bInd))*100);
disp(errorStr)
titleStr = sprintf('%s\n%s',errorStr,bestStrS);
title(titleStr)



%% get length constants

%subplot(2,2,3),cla,hold on
uConlog10 = dat(1).ux;
conNames = {d1.lab d2.lab d3.lab d4.lab};
isGm = find(strcmp(conNames,'conductances'));
bestCon = eval(sprintf('bX%d',isGm));

bestCon = -4
bestConInd = find(uConlog10==bestCon,1);

r = 0.25; % Median radius
%logGm = [-4 -2.75 -2.5 -2];
plotLogGm = uConlog10;
Gm = 10.^ plotLogGm;
%Gm = sortCon;
Cm = 1;
Ri = 100;

rCm = r / 10^4;
Rm = 1./Gm;
ra = Ri / (pi * rCm^2);
rm = Rm ./ (2 * pi * rCm);
lcCm = sqrt(rm./ra);


%lcCm = sqrt((2*rCm*Rm)/(4 * Ri));

lc = lcCm * 10^4;
logLc = log10(lc);
bestLogLc = logLc(bestConInd);
disp(sprintf('tested logLc %s',num2str(logLc')))
disp(sprintf('best logLc is %0.2f. %0.1f um',bestLogLc,10^bestLogLc))

if 0
    clf, hold on
    plot(uConlog10,logLc,'k')
    scatter(uConlog10,logLc,'k')

    scatter(bestCon,bestLogLc,'r')
end

%
% plot(plotLogGm,logLc,'k')
% scatter(plotLogGm,logLc,'r','filled')
% %yticklabels(10.^[0:5])
% ylim([0 5])
%


%% navigate
if 1

    f = gcf;
 fPos = [-1851         175        1605         688];
    f.Position = fPos;

    view(sp1,0,90)
    view(sp1,-63,47)
    sd2 = sd;

    minError = min(useError);
    bInd = find(useError==minError,1);
    bX1 = d1.x(bInd);
    bX2 = d2.x(bInd);
    bX3 = d3.x(bInd);
    bX4 = d4.x(bInd);

    b1 = find(d1.ux == bX1);
    b2 = find(d2.ux == bX2);
    b3 = find(d3.ux == bX3);
    b4 = find(d4.ux == bX4);

    [s1 s2 s3 s4] = size(em);
    disp(sprintf('min error %0.4f',minError))
    bestStr = sprintf('best %s = %0.2f, %s = %0.2f, %s = %0.2f, %s = %0.2f,',...
        d1.lab,bX1,d2.lab,bX2,d3.lab,bX3,d4.lab,bX4);
    disp(bestStr)

    sp2 = subplot(2,2,4); cla(sp2),hold on

    errorMeasure = useError;
    [sortError idxE] = sort(errorMeasure,'ascend');
    bInd = idxE(1);

    bestScat = uPN(:,idxE(1));

    plot([-1 1],[0 0],'k')
    plot([0 0],[-1 1],'k')
    plot([-1 1],[-1 1],'k')
    axis 'equal'




    while 1



    [mX mY] = meshgrid(d2.ux,d1.ux);
    bestPlane = mX * 0 + minError;

    if 0
        %%Show available variables
        h1 = histc(d1.x,d1.ux);
        sh1 = subplot(8,2,2),cla(sh1); hold on;
        bar(d1.ux,h1,'k');
        scatter(d1.ux(b1),0,'r','filled')
        title(d1.lab);
        h2 = histc(d2.x,d2.ux);
        sh2 = subplot(8,2,4), cla(sh2); hold on;
        bar(d2.ux,h2,'k');
        scatter(d2.ux(b2),0,'r','filled')
        title(d2.lab);
        h3 = histc(d3.x,d3.ux);
        sh3 = subplot(8,2,6), cla(sh3); hold on;
        bar(d3.ux,h3,'k');
        scatter(d3.ux(b3),0,'r','filled')

        title(d3.lab);
        h4 = histc(d4.x,d4.ux);
        sh4 = subplot(8,2,8), cla(sh4); hold on;
        bar(d4.ux,h4,.4,'k');
        scatter(d4.ux(b4),0,'r','filled')

        title(d4.lab);
    end



        %%Show best plane
        cla(sp1),

        isCond = find((d4.x == bX4 ) & ( d3.x == bX3));
        isCondX1 = d1.x(isCond);
        isCondX2 = d2.x(isCond);
        isCondY = useError(isCond);

        %plot3(d1.x,d2.x,em(:,tC),'k')
        %m = mesh(mg,em(:,:,dR,dC));
        emS = em(:,:,b3,b4);
        nums = 0;
        if (size(emS,1)>1) & (size(emS,2)>1)
            pM = propMat(:,:,b3,b4);
            isNum = ~isnan(emS);
            nums = sum(isNum(:));

            if nums
                suM = surf(sp1,mX,mY,bestPlane,'facecolor','flat','faceAlpha',.2);

                vq = griddata(mX(isNum),mY(isNum),emS(isNum),mX,mY);
                if ~isempty(vq)
                    vProp = vq(~isnan(vq));
                    vProp = vProp-min(vProp);
                    vProp = round((vProp .* 99 / max(vProp(:)))+1);
                    su = surf(sp1,mX,mY,vq,'facecolor','interp');
                    %su.CDataMapping = 'direct';
                    su.CData(~isnan(vq)) =vProp;
                    %colormap(sp1,jet(100))
                    %colorbar(sp1)
                end
            end
        end
        % sc = scatter3(mX(isNum),mY(isNum),emS(isNum),'k','filled');
        % sc.CData =pCol(pM(isNum),:);
        %scK = scatter3(isCondX2,isCondX1,isCondY,'k','filled','MarkerEdgeColor','k','linewidth',3);
        sc = scatter3(sp1,isCondX2,isCondX1,isCondY,'k','filled','MarkerEdgeColor','w');
        sc.CData =pCol(prop(isCond),:);
        %sc = scatter3(d2.x,d1.x,useError,'k','filled');
        %sc.CData =pCol(prop,:);
        zlim(sp1,[0 max(useError)]);
        ylabel(sp1,d1.lab);
        xlabel(sp1,d2.lab);
       
        zlabel(sp1,'error');
        titleStr1 = sprintf('%s = %0.2f, %s = %0.2f',d3.lab,d3.ux(b3),d4.lab,d4.ux(b4));
        title(sp1,titleStr1);

        %%Match errors for all cells
        sp2 = subplot(2,2,4); cla(sp2),hold on
        cn = find((d1.x==bX1) & (d2.x==bX2) &(d3.x==bX3) &(d4.x==bX4),1);

        if nums & ~isempty(vq)
            mE = vq(b1,b2);
            sn = scatter3(sp1,bX2,bX1,mE,100,'o','g','linewidth',3);
            sn2 = scatter3(sp1,[bX2 bX2],[bX1 bX1],mE+[.01 -.01],10,'o','g','linewidth',3);
        end


       sp2 = subplot(1,2,2);cla(sp2); hold on

        if ~isempty(cn)
            bE = useError(cn);

            bestScat = uPN(:,cn);

            plot([-1 1],[0 0],'k')
            plot([0 0],[-1 1],'k')
            plot([-1 1],[-1 1],'k')
            axis 'equal'

            vs = unique(uV);
            vCol = jet(length(vs));
            for v = 1:length(vs)
                isV = (uV == vs(v));
                scB = scatter(sp2,uCN(isV),bestScat(isV),50,'markerfacecol',vCol(v,:),'markeredgealph',0,...
                    'markerfacealpha',1);
            end


            bestStr = sprintf('best %s = %0.2f, %s = %0.2f, %s = %0.2f, %s = %0.2f,',...
                d1.lab,bX1,d2.lab,bX2,d3.lab,bX3,d4.lab,bX4);
            bestStrS = sprintf('best %s = %0.2f, %s = %0.2f, %s = %0.2f, %s = %0.2f,',...
                d1.lab(1:2),bX1,d2.lab(1:2),bX2,d3.lab(1:2),bX3,d4.lab(1:2),bX4);
            disp(bestStr)
            errorStr = sprintf('min error %0.4f, cc = %0.3f, bin = %0.3f',useError(cn),predCor(cn),(1-binError(cn))*100);
            disp(errorStr);
            titleStr = sprintf('%s\n%s',errorStr,bestStrS);
            title(sp2,titleStr);
        end


        drawnow

        %%get button press for navigation
        disp(sprintf('wating for button press %0.2f',rand));
        set(f,'CurrentCharacter',' ');
        w = waitforbuttonpress;
        k = get(f,'CurrentCharacter');

        switch(lower(k))
            case 'd'
                b2 = b2+1;
                b2 = min(b2,s2);
                bX2 = d2.ux(b2);
            case 'a'
                b2 = b2-1;
                b2 = max(b2,1);
                bX2 = d2.ux(b2);
            case 'w'
                b1 = b1+1;
                b1 = min(b1,s1);
                bX1 = d1.ux(b1);
            case 's'
                b1 = b1-1;
                b1 = max(b1,1);
                bX1 = d1.ux(b1);
            case 'e'
                b3 = b3+1;
                b3 = min(b3,s3);
                bX3 = d3.ux(b3);
            case 'q'
                b3 = b3-1;
                b3 = max(b3,1);
                bX3 = d3.ux(b3);
            case 'z'
                b4 = b4+1;
                b4 = min(b4,s4);
                bX4 = d4.ux(b4);
            case 'c'
                b4 = b4-1;
                b4 = max(b4,1);
                bX4 = d4.ux(b4);
            case 'x'
                bs(sd2) = [b1 b2 b3 b4];

                %sd2 = circshift(sd2,1,2);
                sd2 = randperm(4,4);
                d1 = dat(sd2(1));
                d2 = dat(sd2(2));
                d3 = dat(sd2(3));
                d4 = dat(sd2(4));
                em = permute(errorMat,sd2);

                propMat = em;
                propMat = propMat-min(propMat(:));
                propMat = round(propMat .* 99/max(propMat(:))+1);
                b1  = bs(sd2(1));
                b2  = bs(sd2(2));
                b3  = bs(sd2(3));
                b4  = bs(sd2(4));

                bX1 = d1.ux(b1);
                bX2 = d2.ux(b2);
                bX3 = d3.ux(b3);
                bX4 = d4.ux(b4);




        end

    end
end


%% print figure
if 0

    % sp1.YLim = [-5 -2];
    % sp1.XLim = [0 1.25];

    printPath = 'D:\WorkDocs\Publications\VG3\Revisions\Neuron\Pics\'
    printName = 'realisticInhibMap4'


   
    if ~exist(printPath,'dir'),mkdir(printPath); end
    epsName = sprintf('%s%s.eps',printPath,printName);
    print(gcf, epsName, '-depsc2','-painters','-r300')


end

