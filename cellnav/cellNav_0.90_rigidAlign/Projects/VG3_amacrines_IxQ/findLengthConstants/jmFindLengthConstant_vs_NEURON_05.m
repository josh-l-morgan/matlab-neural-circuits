
%%Check for unassigned bipolar cells
%%Determine if polarity of Ca responses are consistent with other cells. Is
%%it SNR dependent.
%%See if changing masks changes bimodality or polarity spread.
%%Consolidate redundant ROIs
clear all


%% Set options for what to compare
vCids = [2 3 4 5 13 14];
%vCids = [20];

compare = 'conduct'; % 'conduct' , 'inhib'


%% Set options for filtering functional ROIs
filterBySNR = 0.96; %remove rois with low signal
filterByEdge = 1;
weightErrors = 0;
standardize = 0;



%% Choose model outputs to display

expName = 'testOFFON_Cond07';
expName = 'testOFFON_changeInhib_02'; %Stronglly OFF biased
expName = 'exOFFON_CondInhib_03'; % Nice fitting. no minima.
expName = 'testOFFON_Cond08';
expName = 'exOFFON_CondVsInhib_01';
expName = 'exOFFON_CondInhib_02'; % Nice fitting. no minima.
expName = 'testOFFON_Cond09';
expName = 'exp_01_a'; %First to vary conductance, inhibition weight, on/off inhibition ratio, excitation weight
expName = 'exp_01_c'; %
expName = 'exp_01_d'; %
expName = 'exp_01_e'; %

%% Experiments with inhibition
expNames = {'test_exp_01_c'};
expNamesWithInhib = {'exp_01_c' 'exp_01_d' 'exp_01_e'  'exp_01_g' 'exp_01_h' };
expNamesNoInhib = {'exp_02_a' 'exp_02_b' 'exp_02_c' 'exp_02_d' 'exp_02_e' 'exp_02_g' 'exp_02_i' 'exp_02_j' 'exp_02_k'};
expNames = cat(2,expNamesNoInhib,expNamesWithInhib);

%% Load up cell data (non neuron)
if 1
    %global glob
    % datFold = [glob.datDir 'Analysis\Data\preproc\'];
    % SPN =datFold;
    SPN = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\'
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'GOI.mat']);
end



%% Get ROIs
roiCids = GOI.roiCids;%ptDat(:,3);
numRoi = length(roiCids);
goodRoiCid = zeros(numRoi,1);


%% Get calcium data

realCal = GOI.Polarity1;
SNR = GOI.qual;

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
% caPol(caPol>1) = 1;
% caPol(caPol<-1) = -1;


SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');
minBip = 5;


%% Remove bad frames
frames = ptDat(:,1);
allFrames = unique(frames); %all
%useFrames = allFrames;
%useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 2006]; %all
useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 ]; %all but worst
%useFrames = [1005 1006 1007 1008 1009 1010 ]; % only first stack
removeFrames = setdiff(allFrames,useFrames);


%% Get neuron predictions for ROIs
clear nn gnn nnOff nnOn nn2Roi trackC
conductances = [];
inhibWeight  = [];
exciteWeight  = [];
offInhibScale = [];

aC = 0; % Count ALL Conditions
for f = 1:length(expNames)
    expName = expNames{f}
    for v = 1:length(vCids)
        vCid = vCids(v); %Get cid of VG3 to run
        isCid = find(GOI.roiCids==vCid);  % Find all grouped rois for the cid being run
        useNodes = GOI.closeNode(isCid); % regtrieve list of skeleton nodes for rois on cid
        gPos = GOI.pos(isCid,:);

        runCid = vCid;

        nnDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\nn\';
        expDir = [nnDir expName '\'];
        cellDir = sprintf('%scid%d\\',expDir,runCid);

        gnn.nnON(size(GOI.pos,1),1) = 0;
        gnn.nnOFF(size(GOI.pos,1),1) = 0;
        if  exist([cellDir 'expInfo.mat'],'file')
            load([cellDir 'expInfo.mat'])
            % conductancesF = expInfo.conductances(:);
            % inhibWeightF  = expInfo.inhibWeight(:);
            % exciteWeightF = expInfo.exciteWeight(:);
            % offInhibScaleF = expInfo.offInhibScale(:);
            startV(f) = expInfo.nneuron.params.v_init;


            swcDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\swc\';
            swcFile = sprintf('%scid%d.swc',swcDir,runCid);
            mtrFile = sprintf('%scid%d.mtr',swcDir,runCid);

            tree = load_tree(mtrFile);

            %%Find ROIs on NEURON tree
            clf
            scatter3(tree{1}.X,tree{1}.Y,tree{1}.Z,'.','k')
            hold on
            scatter3(gPos(:,1),gPos(:,2),gPos(:,3),'o','r')
            hold off
            drawnow

            for i = 1:length(useNodes)
                d = sqrt((tree{1}.X-gPos(i,1)).^2 + (tree{1}.Y-gPos(i,2)).^2 + (tree{1}.Z-gPos(i,3)).^2);
                minD = min(d);
                nTarg = find(d == minD,1);
                gnn.vcid(isCid(i)) = vCid;
                gnn.nnNode(isCid(i)) = nTarg;
                goodRoiCid(isCid(i)) = 1;
            end

            dCellDir = dir([cellDir '*.mat']);
            eNames = {dCellDir.name};

            for con = 1:length(eNames);

                expFileNameOFF = sprintf('%sex%d_con%02.0f.mat',cellDir,1,con);
                expFileNameON = sprintf('%sex%d_con%02.0f.mat',cellDir,2,con);


                if exist(expFileNameON,'file')

                    aC = aC + 1;
                    filesRead{aC} = expFileNameOFF;

                    trackV(aC) = v;
                    a =  load(expFileNameOFF);
                    % nn(aC).exOff = a.nnRes;
                    b = load(expFileNameON);
                    r.nnON{aC} = b.nnRes.maxV(gnn.nnNode(isCid))-startV(f);
                    r.nnOFF{aC} = a.nnRes.maxV(gnn.nnNode(isCid))-startV(f);

                    % nn(aC).exOn = b.nnRes;

                    r.conductances(aC) = b.nnRes.nneuron.mech{1}.all.pas.g;
                    r.inhibWeight(aC)  = a.nnRes.syn2.weight;
                    r.exciteWeight(aC) = b.nnRes.syn1.weight;
                    r.nn2Roi{aC} = isCid;

                    iWa =  a.nnRes.nneuron.con(2).weight;
                    iWb =  b.nnRes.nneuron.con(2).weight;
                    oS = iWb/iWa;
                    r.offInhibScale(aC) = oS;

                end
            end
        end
    end
end

offInhibScale(isnan(offInhibScale)) = 1;

%% Start nVc

nVc.read.aC = aC;
nVc.read.files =  filesRead;
nVc.read.trackV = trackV;
nVc.read.nnON = nnON;
nVc.read.nnOF = nnOFF;
nVc.read.conductances = conductances;
nVc.




%% recombine results by condition


%%Collect unique values for comparison
uIW = unique(inhibWeight)';
uCon = unique(conductances)';
uEW = unique(exciteWeight)';
uOS = unique(offInhibScale)';

%uOS = 0.5;

cON = zeros(length(uCon),length(uEW),length(uIW),length(uOS),length(goodRoiCid));
cOFF = cON;
isRes = cON;
gChecked = zeros(length(goodRoiCid),1);
for d1 = 1:length(uCon);
    for d2 = 1:length(uEW)
        for d3 = 1:length(uIW)
            for d4 = 1:length(uOS)
                isCond = find((conductances==uCon(d1)) & (exciteWeight==uEW(d2)) ....
                    & (inhibWeight == uIW(d3)) & (offInhibScale == uOS(d4)));
                if ~isempty(isCond)
                    n2r = cat(1,nn2Roi{isCond});
                    nON = cat(2,nnON{isCond});
                    nOFF = cat(2,nnOFF{isCond});
                    uG = unique(n2r);
                    if length(uG) ~= length(n2r)
                        disp('multiple results for same condition')
                    end
                    for g = 1:length(uG)
                        gTo = uG(g);
                        gFrom = find(n2r==gTo);
                        cON(d1,d2,d3,d4,gTo) = mean(nON(gFrom));
                        cOFF(d1,d2,d3,d4,gTo) = mean(nOFF(gFrom));
                        isRes(d1,d2,d3,d4,gTo) = 1;
                        gChecked(gTo) = 1;
                    end
                end
            end
        end
    end
end

useV = unique(roiCids(gChecked>0));
numRes = sum(isRes,5);
[s1 s2 s3 s4 ] = ind2sub(size(numRes),find(numRes));
cPol = (cOFF-cON)./(cOFF+cON);

nPol = zeros(size(isRes,5),length(s1));
for i = 1:length(s1)
    nPol(:,i) = cPol(s1(i),s2(i),s3(i),s4(i),:);
end
nPol(isnan(nPol)) = 0;


numCon = size(nPol,2);


%% Pull out data to use
%%Sort variables into numbered data structure


dat(1).r = uCon(s1);
dat(1).x = log10(dat(1).r);
dat(1).ux = unique(dat(1).x);
dat(1).u = unique(dat(1).r);
dat(1).lab = 'conductances';
dat(1).num = length(dat(1).u);

dat(2).r = uEW(s2);
dat(2).x = log10(dat(2).r);
dat(2).ux = unique(dat(2).x);
dat(2).u = unique(dat(2).r);
dat(2).lab = 'exciteWeight';
dat(2).num = length(dat(2).u);


dat(3).r = uIW(s3);
dat(3).x = log10(dat(3).r);
dat(3).ux = unique(dat(3).x);
dat(3).u = unique(dat(3).r);
dat(3).lab = 'inhibWeight';
dat(3).num = length(dat(3).u );


dat(4).r = uOS(s4);
dat(4).x = dat(4).r;
dat(4).ux = unique(dat(4).x);
dat(4).u = unique(dat(4).r);
dat(4).lab = 'offInhibScale';
dat(4).num = length(dat(4).u);

%% save neuron vs calcium data
nVc.caPol = caPol;
nVc.dat = dat;
nVc.cPol = cPol;
nVc.nPol = nPol;
nVc.roiCids = roiCids
nVc.goodRoiCid = goodRoiCid;
nVc.gChecked = gChecked;
save([SPN 'nVc.mat'],'nVc');
disp('finished saving nVc.mat')

return

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
clear meanError binError errorOfBins errorOfMeans errorOfVars rms predCor
for c = 1:numCon
    uPNs = uPN(:,c);
    meanError(c) = mean(abs(uPNs-uCN));
    binError(c) = mean(abs((uPNs>0)-(uCN>0)));
    errorOfBin(c) = abs(mean(uPNs>0)-mean(uCN>0));
    errorOfMeans(c) = abs(mean(uPNs) - mean(uCN));
    errorOfVars(c) = abs(var(uPNs) - var(uCN));
    rms(c) = sqrt(mean((uPNs-uCN).^2));
    cc = corrcoef(uPNs,uCN);
    predCor(c) = cc(1,2);
    %meanError(c) = cc(1,2);
    %meanError(c) = mean(sqrt((uP-uC).^2));
end

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
        %meanError(c) = cc(1,2);
        %meanError(c) = mean(sqrt((uP-uC).^2));
    end
end





errorMat = zeros(dat(1).num,dat(2).num,dat(3).num,dat(4).num);
idMat = cell(dat(1).num,dat(2).num,dat(3).num,dat(4).num);
for d1 = 1:dat(1).num
    for d2 = 1:dat(2).num
        for d3 = 1:dat(3).num
            for d4 = 1:dat(4).num
                targ = (dat(1).r == dat(1).u(d1)) & (dat(2).r == dat(2).u(d2)) & ...
                    (dat(3).r == dat(3).u(d3)) & (dat(4).r == dat(4).u(d4));

                idMat{d1,d2,d3,d4} = find(targ);
                sME = meanError(targ);
                errorMat(d1,d2,d3,d4) = mean(sME);
                if length(sME)>1
                    disp('found multiple results for combination of conditions')
                end
            end
        end
    end
end



%% Show differences
sd = [ 1 2 3 4];
d1 = dat(sd(1));
d2 = dat(sd(2));
d3 = dat(sd(3));
d4 = dat(sd(4));
em = permute(errorMat,sd);


clf

prop = meanError;
prop = prop-min(prop);
prop = round(prop * 99/max(prop)+1);

propMat = em;
propMat = propMat-min(propMat(:));
propMat = round(propMat .* 99/max(propMat(:))+1);

pCol = jet(100);
[mX mY] = meshgrid(d2.ux,d1.ux);


minError = min(meanError);
bInd = find(meanError==minError,1);
bX1 = d1.x(bInd);
bX2 = d2.x(bInd);
bX3 = d3.x(bInd);
bX4 = d4.x(bInd);


[b1 b2 b3 b4] = ind2sub(size(em),bInd)
disp(sprintf('min error %0.4f',minError))
bestStr = sprintf('best %s = %0.1f, %s = %0.1f, %s = %0.1f, %s = %0.1f,',...
    d1.lab,bX1,d2.lab,bX2,d3.lab,bX3,d4.lab,bX4);
disp(bestStr)

bestPlane = mX * 0 + minError;

if 1 % show each plane

    for dC = 1:d4.num
        for dR = 1:d3.num
            sp = subplot(d3.num,d4.num,(dR-1)*d4.num+dC); cla(sp); hold on
            %sp = subplot(1,1,1); cla(sp); hold on

            isCond = find((d4.x == d4.ux(dC) ) & ( d3.x == d3.ux(dR)));
            isCondX1 = d1.x(isCond);
            isCondX2 = d2.x(isCond);
            isCondY = meanError(isCond);


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
            %sc = scatter3(d2.x,d1.x,meanError,'k','filled');
            %sc.CData =pCol(prop,:);
            zlim([0 max(meanError)])
            ylabel(d1.lab)
            xlabel(d2.lab)
            zlabel('error')
            view(-8,20)
            title(sprintf('%s = %0.1f, %s = %0.1f',d3.lab,d3.ux(dR),d4.lab,d4.ux(dC)));

            drawnow
            if (d4.ux(dC) == bX4) & (d3.ux(dR) == bX3)
                scatter3(d2.x(bInd),d1.x(bInd),minError,100,'o','r','linewidth',3)
            end
        end
    end

    return
end

%% Match errors for each cell
errorMeasureV = meanErrorV;

sp = subplot(1,2,1); cla(sp),hold on

plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
plot([-1 1],[-1 1],'k')
axis 'equal'

vs = unique(uV);
vCol = jet(length(vs));
for v = 1:length(useV)
    isV = (uV == useV(v));
    [sortError idxE] = sort(errorMeasureV{v},'ascend');
    bestScat = uPN(isV,idxE(1));
    scB = scatter(uCN(isV),bestScat,50,'markerfacecol',vCol(v,:),'markeredgealph',0,...
        'markerfacealpha',1)
end



sp = subplot(1,2,2); cla(sp),hold on

for c = 1:10
    cla, hold on
    for v = 1:length(useV)

        isV = (uV == useV(v));
        bestScat = uPN(isV,idxE(c));
        scB = scatter(uCN(isV),bestScat,50,'markerfacecol',vCol(v,:),'markeredgealph',0,...
            'markerfacealpha',1)
        otherStr = sprintf('other %s = %0.1f, %s = %0.1f, %s = %0.1f, %s = %0.1f,',...
            d1.lab,d1.x(idxE(c)),d2.lab,d2.x(idxE(c)),d3.lab,d3.x(idxE(c)),d4.lab,d4.x(idxE(c)));
        disp(otherStr)

    end
    pause(.1)
end

return

%% error measures = meanError errorOfMeans errorOfVars rms predCor errorOfBin binError

meanPol = mean(nPol,1);
clf


errorMeasure = meanError;
[sortError idxE] = sort(errorMeasure,'ascend');

bestScat = uPN(:,idxE(1));

sp = subplot(1,2,1); cla(sp),hold on

plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
plot([-1 1],[-1 1],'k')
axis 'equal'

vs = unique(uV);
vCol = jet(length(vs));
for v = 1:length(vs)
    isV = (uV == vs(v));
    scB = scatter(uCN(isV),bestScat(isV),50,'markerfacecol',vCol(v,:),'markeredgealph',0,...
        'markerfacealpha',1)
end


if 1 %% Show runners up
    sp = subplot(1,2,2); cla(sp),hold on

    plot([-1 1],[0 0],'k')
    plot([0 0],[-1 1],'k')
    plot([-1 1],[-1 1],'k')
    axis 'equal'

    scB = scatter(uCN,bestScat,50,'markeredgecol','k','markeredgealph',1,...
        'markerfacealpha',1)
    scO = scatter(uCN,bestScat,20,'filled','markeredgecol','r','markeredgealph',1,...
        'markerfacealpha',1)

    bestStr = sprintf('best %s = %0.1f, %s = %0.1f, %s = %0.1f, %s = %0.1f,',...
        d1.lab,bX1,d2.lab,bX2,d3.lab,bX3,d4.lab,bX4);
    disp(bestStr)
    disp(sprintf('min error %0.4f',minError))

    for c = 1:1;%length(idxE)
        c
        uPNs = uPN(:,idxE(c));
        delete(scO);
        scO = scatter(uCN,uPNs,20,'filled','markerfacecol','r','markeredgealph',0,...
            'markerfacealpha',1);

        otherStr = sprintf('other %s = %0.1f, %s = %0.1f, %s = %0.1f, %s = %0.1f,',...
            d1.lab,d1.x(idxE(c)),d2.lab,d2.x(idxE(c)),d3.lab,d3.x(idxE(c)),d4.lab,d4.x(idxE(c)));
        disp(otherStr)
        pause(.1)
    end
end


return
%% get length constants

%subplot(2,2,3),cla,hold on

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
logLc = log10(lc)

%
% plot(plotLogGm,logLc,'k')
% scatter(plotLogGm,logLc,'r','filled')
% %yticklabels(10.^[0:5])
% ylim([0 5])
%



return

%{


%% Find error for all
predList = nPol';
usePred = predList(:,useRoi);

%%Set error weights
if weightErrors == 1
    ew = SNR;% * 0 + 1;
else %dont weight errors
    ew = SNR * 0 + 1;
end
ew = ew(goodRoi>0);
ew = ew(:)';
[sew ewRank] = sort(ew,'ascend');
ewRank = ewRank/max(ewRank);
ewM = repmat(ew,[size(usePred,1) 1]);





%% Standardize polarities
numCond = size(predList,1);
caPolN = caPol(useRoi);% 
if standardize
    caPolN = caPolN - mean(caPolN);
    caPolN = caPolN/std(caPolN);
end
usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
if standardize
    usePredN = usePredN - repmat(mean(usePredN,2),[1 numUse]);
    usePredN = usePredN ./ repmat(std(usePredN,2),[1 numUse]);
end
caPolMatN = repmat(caPolN',[numCond 1]);

%% Polarity distributions
%{
meanPred = mean(usePred(~isnan(usePred)));
meanCa = mean(caPolN);

numOn = numOn(goodCid>0);
numOff = numOff(goodCid>0);
cidPol = (numOff-numOn)./(numOff+numOn);
%}



%% Get best result for each roi
condList = conductances';
testLengths = conductances';

numCondType = 1;
difMatN = 0-abs(caPolMatN - usePredN);
bestCondR = zeros(numCondType,numUse);
maxMatchR = zeros(1,numUse);
for r = 1:numUse
    
    maxMatchR(r) = max(difMatN(:,r));
    matchIdx = find(difMatN(:,r) == maxMatchR(r));
    bestCondRaw = condList(matchIdx,:);
    bestCondR(:,r) = mean(bestCondRaw,1);

end

bestLengthsR = bestCondR(1,:);
inBounds = ((bestLengthsR > min(testLengths)) & (bestLengthsR<max(testLengths)));


filterErrors = [-1:.01:0];
meanFiltL = filterErrors * 0;
for i = 1:length(filterErrors)
    goodR = (maxMatchR>filterErrors(i)) & (inBounds);
    meanFiltL(i) = mean(bestLengthsR(goodR));
end

matchThresh = -.2;
goodMatches = (maxMatchR>=matchThresh)&inBounds;
meanBestCondR = mean(bestCondR(:,goodMatches),2)';
medianBestCondR = median(bestCondR(:,goodMatches),2)';


lengthsR = bestCondR(1,goodMatches);
stdE_R = std(lengthsR)/sqrt(length(lengthsR))
mean_R = mean(bestCondR(1,goodMatches))

sortR = sort(bestCondR(1,goodMatches));
r95 = [sortR(round(.05*length(sortR))) sortR(round(.95*length(sortR))) ]
r50 = [sortR(round(.25*length(sortR))) sortR(round(.75*length(sortR))) ]

for c = 1:numCondType
    subplot(numCondType,2,c*2)
    hist(bestCondR(c,inBounds),unique(bestCondR(c,inBounds)))
end

% disp(sprintf('Mean for Rs: length = %0.2f, noise = %0.4f, on = %0.2f, off = %0.2f',...
%     meanBestCondR(1),meanBestCondR(2),meanBestCondR(3),meanBestCondR(4)));
% disp(sprintf('Median for Rs: length = %0.2f, noise = %0.4f, on = %0.2f, off = %0.2f',...
%     medianBestCondR(1),medianBestCondR(2),medianBestCondR(3),medianBestCondR(4)));
meanMatchR = mean(maxMatchR)
medianMatchR = median(maxMatchR)

subplot(4,2,6)
scatter(maxMatchR,bestCondR(1,:))
hold on
plot(filterErrors,meanFiltL)
hold off

subplot(4,2,8)
hist(maxMatchR,[-1:.01:0])
xlim([-.2 0])


%% Measure error for population
caPolMatM2 = caPolMatN - repmat(mean(caPolMatN,2),[1 numUse]);
usePredM2 = usePredN - repmat(mean(usePredN,2),[1 numUse]);
covN = mean(caPolMatM2 .* usePredM2 .* ewM,1);
std1 = std(caPolMatM2,1);
std2 = std(usePredM2,1);
cc = covN./(std1.*std2);

rmseN = zeros(numCond,1);
if 0
    %%Calculate error with rmse
    difMatN = caPolMatN - usePredN;
    meanErrN = mean(abs(difMatN));
    rmseN = sqrt(mean(difMatN.^2,2));
    rmseN = rmseN/ mean(abs(caPolN));
    rmseN = (std1*3)-rmseN./std1; %Adjust for positive error measure
elseif 0
    difMatN = caPolMatN - usePredN;
    meanDif = mean(abs(difMatN),2);
    rmseN = 0 - meanDif;
elseif 1
    difMatN = caPolMatN - usePredN;
    meanDif = mean(abs(difMatN),2);
    rmseN = 0 - meanDif; 
elseif 0
    sameSign = (caPolMatN./abs(caPolMatN)) == (usePredN./abs(usePredN));
    mean(sameSign,2);
    rmseN = mean(sameSign,2);
else
    %%Correlation coef
    for i = 1:numCond
        CC = corrcoef(caPolMatN(i,:),usePredN(i,:));
        rmseN(i) = CC(1,2);
    end
    rmseN(isnan(rmseN)) = 0;
end


matchVal = rmseN;



%%Find peak match (average) in case of multiple hits)
minE = max(matchVal(:));
indE = find(matchVal==minE);
numBest = length(indE);
bestCond = mean(condList(indE,:),1);

bestLength = bestCond(1);
% bestNoise = bestCond(2);
% bestOnScale = bestCond(3);
% bestOffScale = bestCond(4);
% bestStr = sprintf('For all: length = %0.2f, noise = %0.2f, onScale = %0.2f, offScale = %0.2f',...
%     bestCond(1),bestCond(2),bestCond(3),bestCond(4));
% disp(bestStr)




%% find change of each parameter at best of others
clf
hold off
isBestMat = condList==repmat(bestCond,[numCond 1]);
condNames = {'length','noise','scale ON', 'scale OFF'};
for i = 1:numCondType
    filtWithCond = setdiff([1:numCondType],i);
    isBest = find(sum(isBestMat(:,filtWithCond),2) == length(filtWithCond));
    cForBest = condList(isBest,i);
    mForBest = matchVal(isBest);
    subplot(numCondType,2,(i-1)*2+1)
    plot(cForBest,mForBest)
    %title(sprintf('best %s = %0.4f',condNames{i},bestCond(i)))
end


pause(1)


return




%{

clf

subplot(3,2,1)
showErr = rmseN(condList());
plot(testLengths,showErr)
title({sprintf('best length = %0.1f, OnScale = %0.1f,OffScale = %0.1f,', bestLength,bestOnScale,bestOffScale),...
    sprintf('noise = %0.4f, SNR filter = %0.2f',bestNoise,filterBySNR)})


subplot(3,2,2)
showErr = cc(:,:,round(bs1m),round(bs2m),round(bnm));
plot(testLengths,showErr)
title(sprintf('cc'))



showPred = usePredN(:,:,round(bs1m),round(bs2m),round(bnm));
%showPred = squeeze(usePredN(:,round(bxm),:,round(bnm))); %% show all on scaling
%showPred = squeeze(usePredN(:,round(bxm),round(bzm),:)); %% show all noise



%%Show shift in matching with length constants
if 1


    subplot(3,2,3)

    stdDif = abs(testLengths-15);
    std15 = find(stdDif==min(stdDif),1);
    scatS = scatter(usePredN(:,std15,round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    scatS.CData = SNRcol(useRoi,:);
    xlim([-1 1])
    ylim([-1 1])
    title(sprintf('distribution at 15 um standard'))
    set(gca,'color','k')

    drawnow


    subplot(3,2,4)

    if 1
        for r = 1
            for s = 1:size(showPred,2)
                scat1 = scatter(showPred(:,s),caPolN,15,'filled');
                hold on
                plot([-1 1],[-1 1],'w')
                scat1.CData = SNRcol(useRoi,:);
                title(sprintf('length constant %0.2f',testLengths(s)))
                xlim([-1 1])
                ylim([-1 1])
                set(gca,'color','k')
                drawnow
                pause(.1)
                hold off
            end
            for s = size(showPred,2):-1:1
                scat1 = scatter(showPred(:,s),caPolN,15,'filled');
                hold on
                scat1.CData = SNRcol(useRoi,:);
                title(sprintf('length constant %0.2f',testLengths(s)))
                xlim([-1 1])
                ylim([-1 1])
                set(gca,'color','k')
                drawnow
                hold off
                %pause(.1)
            end

        end
    end

    scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    scat1.CData = SNRcol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(bxm)))
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'color','k')

    drawnow


end

if 0
    clf
    jCol = jet(100);
    snrInd = SNR;
    snrInd = snrInd-min(snrInd);
    snrInd = round(snrInd * 299/max(snrInd(:))) + 1;
    snrInd(isnan(snrInd)) = 1;
    snrInd(snrInd<1) = 1;
    snrInd(snrInd>100) = 100;
    snrCol = jCol(snrInd,:)


    scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,20,'o','MarkerFaceColor',...
        'none','LineWidth',1);
    scat1.CData = snrCol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(bxm)))
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'color','w')

    hold on
    plot([-1 1],[-1 1])

    return

end

%% Compare errror to SNR

if 0
    bestErrors = difMatN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
    scatter(SNR(useRoi),bestErrors)
    caN = caPolMatN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
    uPN = usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
    scatter3(caN,uPN,SNR(useRoi))
end


%% Check each roi
useRoi = find(goodRoi>0);
subplot(3,2,5)
cla, hold on
for r = 1:length(useRoi);

    usePredR = usePredN(r,:,:,:,:);
    caPolNR = caPolN(r);% - mean(caPol(useRoi));
    dif = abs(usePredR-caPolNR);

    minE = min(dif(:));
    indE = find(dif==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(dif),indE);
    bymR = mean(by);
    bxmR = mean(bx);
    bs1mR = mean(bs1);
    bs2mR = mean(bs2);
    bnmR = mean(bn);

    bestOnScaleR(r) = onScale(round(bs1m));
    bestOffScaleR(r) = offScale(round(bs2m));
    bestLengthR(r) = testLengths(round(bxmR));
    bestNoiseR(r) = noise(round(bnm));

    %     subplot(2,2,3)
    %     hold off
    %     scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    %     scat1.CData = SNRcol(useRoi,:);
    %     title(sprintf('length constant %0.2f',testLengths(s)))
    %     %                 xlim([-1 1])
    %     %                 ylim([-1 1])
    %     set(gca,'color','k')
    %     hold on
    %     scat2 = scatter(usePredN(r,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN(r),150,'filled');
    %
    if 0
        subplot(3,2,5)
        showDif = dif(:,:,round(bs1m),round(bs2m),round(bnm))
        hold on
        plot(showDif)
        drawnow
    end
end

useBest = (bestLengths>testLengths(1)) & (bestLengths<testLengths(end));
useBestLengths = bestLengths(useBest);

percentLengthsUsed = mean(useBest)*100

subplot(3,2,5)
hold off
hist(useBestLengths,testLengths)

medBestLengthsR = median(useBestLengths)
sortBL = sort(useBestLengths,'ascend');
bl95 = [sortBL(round(length(sortBL) * .025))      sortBL(round(length(sortBL) * .975))]
SE = std(useBestLengths)/sqrt(length(useBestLengths))
title(sprintf('L = %0.2f, 95%% = %0.2f - %0.2f, SE = %0.2f',medBestLengthsR, bl95(1),bl95(2)),SE)



%% Show each cell

useCid = vCids(goodCid>0);
pCol  = [1 0 0; 0 1 0; 0 0 1; .75 .75 0; .75 0 .75; 0 .75 .75; .3 0 0; 0 .3 0; 0 0 .3; .3 .3 0; .3 0 .3; 0 .3 .3];

bestOnScales = zeros(1,length(useCid));
bestOffScales = zeros(1,length(useCid));
bestLengths = zeros(1,length(useCid));
bestNoises = zeros(1,length(useCid));


for v = 1:length(useCid)
    vCid = useCid(v);
    useRoi = find((goodRoi>0) & (GOI.roiCids==vCid));

    usePred = allPred(useRoi,:,:,:,:);
    caPolMat = repmat(caPol(useRoi),1, L, size(usePred,3), size(usePred,4), size(usePred,5));
    difMat = caPolMat - usePred;
    rms = sqrt(mean(difMat.^2,1));
    meanErr = mean(abs(difMat),1);


    %Set error weights
    if weightErrors == 1
        ew = SNR;% * 0 + 1;
    else %dont weight errors
        ew = SNR * 0 + 1;
    end
    ew = ew(useRoi);
    ew = ew/mean(ew);
    ewM = repmat(ew,[1 size(usePred,2) size(usePred,3) size(usePred,4) size(usePred,5)]);


    caPolN = caPol(useRoi);% - mean(caPol(useRoi));
    if standardize
        caPolN = caPolN/std(caPolN);
    end
    usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
    if standardize
        stdUsePred = std(usePredN,1);
        usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1 1]);
    end
    caPolMatN = repmat(caPolN,[1 L size(usePred,3) size(usePred,4) size(usePred,5)]);



    %%Get correlation coefficient
    caPolMatM2 = caPolMatN - repmat(mean(caPolMatN,1),[size(caPolMatN,1) 1 1 1 1]);
    usePredM2 = usePredN - repmat(mean(usePredN,1),[size(usePredN,1) 1 1 1 1]);
    covN = mean(caPolMatM2 .* usePredM2 .* ewM,1);
    std1 = std(caPolMatM2,1);
    std2 = std(usePredM2,1);
    cc = covN./(std1.*std2);


    if 0
        %%Calculate error with rmse
        difMatN = caPolMatN - usePredN;
        meanErrN = mean(abs(difMatN));
        rmseN = sqrt(mean(difMatN.^2,1));
        rmseN = rmseN/ mean(abs(caPolN));
        rmseN = (std1*3)-rmseN./std1; %Adjust for positive error measure
    elseif 1
        difMatN = caPolMatN - usePredN;
        meanErrN = mean(abs(difMatN));
        meanDif = mean(abs(difMatN),1);
        rmseN = 0 - meanDif;
    else
        sameSign = (caPolMatN./abs(caPolMatN)) == (usePredN./abs(usePredN));
        mean(sameSign,1);
        rmseN = mean(sameSign,1);
    end

    if matchType == 1
        matchVal = rmseN;
    else
        matchVal = cc;
    end

    minE = max(matchVal(:))
    indE = find(matchVal==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(matchVal),indE);
    bym = mean(by);
    bxm = mean(bx);
    bs1m = mean(bs1);
    bs2m = mean(bs2);
    bnm = mean(bn);

    bestOnScales(v) = onScale(round(bs1m))
    bestOffScales(v) = offScale(round(bs2m))
    bestLengths(v) = testLengths(round(bxm))
    bestNoises(v) = noise(round(bnm))

    showErr = matchVal(:,:,round(bs1m),round(bs2m),round(bnm));
    showMeanErr = meanErrN(:,:,round(bs1m),round(bs2m),round(bnm));

    showPred = usePredN(:,:,round(bs1m),round(bs2m),round(bnm));


    subplot(3,2,6)
    hold on
    plot(testLengths,showErr,'color',pCol(v,:))
    title(' ')

    %     ylim([0 1])
    %     xlim([0 40])
    title(sprintf('%d rois for cell %d',length(useRoi),vCid))

    drawnow
end
useCid

bestOnScales
bestOffScales
bestLengths
bestNoises

%}

%}
