
%%Check for unassigned bipolar cells
%%Determine if polarity of Ca responses are consistent with other cells. Is
%%it SNR dependent.
%%See if changing masks changes bimodality or polarity spread.
%%Consolidate redundant ROIs
clear all


%% Set options for what to compare
vCids = [2 3 4 5 13 14];
%vCids = [20];

readAll = 0;
doSmooth = 1;

%% Set options for filtering functional ROIs
filterBySNR = 0.96; %remove rois with low signal
filterByEdge = 1;
weightErrors = 0;
standardize = 0;

%% Experiments
expNames = {'test_exp_01_c'};
expNamesWithInhib = {'exp_01_c' 'exp_01_d' 'exp_01_e'  'exp_01_g' 'exp_01_h' 'exp_01_j' 'exp_01_k' 'exp_01_l'};
expNamesNoInhib = {'exp_02_a' 'exp_02_b' 'exp_02_c' 'exp_02_d' 'exp_02_e' 'exp_02_g' 'exp_02_i' 'exp_02_j' 'exp_02_k' 'exp_02_k'};
expNamesMore = {'exp2_03_a'};
expNames = cat(2,expNamesNoInhib,expNamesWithInhib,expNamesMore);
nVcName = 'nVc.mat';

if doSmooth
    expNames = {'expSm_2_a' 'expSm_2_b'};
    nVcName = 'nVcSm.mat';

end

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
if readAll
    readNames = {};
else
    load([SPN nVcName]);
    readNames = {nVc.exRes(:).expName};
end


for f = 1:length(expNames)
    expName = expNames{f};
    if sum(strcmp(readNames,expName)) & ~readAll
        disp(sprintf('alread read %s',expName));
    else
        disp(sprintf('reading %s',expName));

    aC = 0; % Count ALL Conditions
    clear nn gnn r nn2Roi trackC
    r.expName = expName;
    % r.conductances = [];
    % r.inhibWeight  = [];
    % r.exciteWeight  = [];
    % r.offInhibScale = [];
    % r.nn2Roi = {};

    for v = 1 : length(vCids)
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
            r.startV = expInfo.nneuron.params.v_init;

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
                    r.nnON{aC} = b.nnRes.maxV(gnn.nnNode(isCid))-r.startV;
                    r.nnOFF{aC} = a.nnRes.maxV(gnn.nnNode(isCid))-r.startV;

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
    r.offInhibScale(isnan(r.offInhibScale)) = 1;
    r.aC = aC;
    nVc.exRes(f) = r;
    end
end
save([SPN nVcName],'nVc');



%% pack variables
conductances = cat(2,nVc.exRes.conductances);
inhibWeight = cat(2,nVc.exRes.inhibWeight);
exciteWeight = cat(2,nVc.exRes.exciteWeight);
offInhibScale = cat(2,nVc.exRes.offInhibScale);
nn2Roi = cat(2,nVc.exRes.nn2Roi);
nnON = cat(2,nVc.exRes(:).nnON);
nnOFF = cat(2,nVc.exRes(:).nnOFF);


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
nVc.nn2Roi = nn2Roi;
nVc.caPol = caPol;
nVc.dat = dat;
nVc.cPol = cPol;
nVc.nPol = nPol;
nVc.roiCids = roiCids;
nVc.goodRoiCid = goodRoiCid;
nVc.gChecked = gChecked;
save([SPN nVcName],'nVc');
disp('finished saving nVc.mat')
