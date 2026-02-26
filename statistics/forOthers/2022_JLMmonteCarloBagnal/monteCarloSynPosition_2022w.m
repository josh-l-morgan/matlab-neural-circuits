%%Monte carlo synapses across nodes

%%Note current nodes are not evenly distributed so that conditioning of
%%skeleton will be required
%%should then reassign synapses to new nodes



%% define data locations
clear all
connectorFile = ['E:\Bagnal\2022\connectors_table.csv'];
skelCSVFolder = 'E:\Bagnal\2022\vestibular_neurons_CSV_files\'
affCSVFolder = 'E:\Bagnal\2022\afferent_files_SWC_CSV\'

%% Parse connections
cons = readtable(connectorFile,'Delimiter',',');
conLoc = cons{:,3};
conNode = cons{:,9};
conTarget = cons{:,10};
conSource = cons{:,5};
conSourceNode = cons{:,4};

uTargets = unique(conTarget);

%% Parse skeletons
dFold = dir([skelCSVFolder '*.csv']);
sNams = {dFold.name};
clear skel skelIDs skelType
for i = 1:length(sNams)
    
    nam = sNams{i};
    
    if regexp(nam,'VN')
        skelType(i) = 1;
    elseif regexp(nam,'VS')
        skelType(i) = 2;
    elseif regexp(nam,'tangential')
        skelType(i) = 3;
    else
        skelType(i) = 0;
    end
    
    
    skelFile = [skelCSVFolder nam];
    skelT = readtable(skelFile,'Delimiter',',');
    parent = skelT{:,3};
    
    skelID =  skelT{1,11}; %get ID of cell
    skelIDs(i) = skelID;
    skels(i).pos = skelT{:,[ 5 6 7]}; %get possition of nodes
    
    
    nid = skelT{:,2}; % get IDs of nodes
    lookupNid = zeros(1,max(nid));
    lookupNid(nid) = 1:length(nid); % make up lookup table for conversion of node ids to node indexes
    skels(i).pred = parent * 0 -1; %set all to -1
    skels(i).pred(parent>0) = lookupNid(parent(parent>0)); % look up correct values for non-root nodes
    
    skels(i).nid = nid;
    
    
    inp = find(conTarget == skelID);
    skels(i).inputID = inp;
    skels(i).inputNodeRaw = conNode(inp);
    skels(i).sNid = lookupNid(conNode(inp)); %node index for synapse
    
    
end

useID = skelIDs(skelType>0);


%% Parse afferents


dFold = dir([affCSVFolder '*.csv']);
aNams = {dFold.name};
clear aff affIDs affType
for i = 1:length(aNams)
    
    nam = aNams{i};
    
    if regexp(nam,'AnteriorMacula')
        affType(i) = 1;
    else
        affType(i) = 0;
    end
    
    
    skelFile = [affCSVFolder nam];
    skelT = readtable(skelFile,'Delimiter',',');
    parent = skelT{:,3};
    
    skelID =  skelT{1,11}; %get ID of cell
    affIDs(i) = skelID;
    affs(i).pos = skelT{:,[ 5 6 7]}; %get possition of nodes
    
    
    nid = skelT{:,2}; % get IDs of nodes
    lookupNid = zeros(1,max(nid));
    lookupNid(nid) = 1:length(nid); % make up lookup table for conversion of node ids to node indexes
    affs(i).pred = parent * 0 -1; %set all to -1
    affs(i).pred(parent>0) = lookupNid(parent(parent>0)); % look up correct values for non-root nodes
    
    affs(i).nid = nid;
    
        
    
    inp = find(conSource == skelID);
    affs(i).outputID = inp;
    affs(i).outputNodeRaw = conSourceNode(inp);
    affs(i).sNid = lookupNid(conSourceNode(inp)); %node index for synapse
    
    
end

useAffID = affIDs(affType>0);

%% Find nodes in proximity to afferents. 



allAffPos = [];
allAffSynPos = [];
for i = 1:length(useAffID)
    targ = find(affIDs == useAffID(i),1);
    s = affs(targ);
    pos = s.pos;
    sPos = pos(s.sNid,:);
    allAffPos = cat(1,allAffPos,pos);
    allAffSynPos = cat(1,allAffSynPos,sPos);
end
    


allPostPos = [];
allPostSynPos = [];
for i = 1:length(useID)
    targ = find(skelIDs == useID(i),1);
    s = skels(i);
    pos = s.pos;
    sPos = pos(s.sNid,:);
    allPostPos = cat(1,allPostPos,pos);
    allPostSynPos = cat(1,allPostSynPos,sPos);
end


clf
scatter3(allAffPos(:,1),allAffPos(:,2),allAffPos(:,3),'k','.')
hold on
scatter3(allAffSynPos(:,1),allAffSynPos(:,2),allAffSynPos(:,3),'g','o')
scatter3(allPostPos(:,1),allPostPos(:,2),allPostPos(:,3),'r','.')
scatter3(allPostSynPos(:,1),allPostSynPos(:,2),allPostSynPos(:,3),'m','o')
hold off
pause(.1)


%%Find real proximity weights of synapses
dif1 = allAffSynPos(:,1) - allPostSynPos(:,1)';
dif2 = allAffSynPos(:,2) - allPostSynPos(:,2)';
dif3 = allAffSynPos(:,3) - allPostSynPos(:,3)';
distsSyn2Syn = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
% 
% 
% %%Find real proximity weights of post synaptict cell nodes
% dif1 = allAffPos(:,1) - allPostSynPos(:,1)';
% dif2 = allAffPos(:,2) - allPostSynPos(:,2)';
% dif3 = allAffPos(:,3) - allPostSynPos(:,3)';
% distsAff2Syn = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;

syn2Aff = distsSyn2Syn;

distThresh = 10;
threshSyn2Aff = distThresh-syn2Aff;
threshSyn2Aff(threshSyn2Aff<0) = 0;
sumWeightSyn = sum(threshSyn2Aff,1);
hist(sumWeightSyn)

%%Find real proximity presynaptic syn nodes to all post synaptic cell nodes
dif1 = allAffSynPos(:,1) - allPostPos(:,1)';
dif2 = allAffSynPos(:,2) - allPostPos(:,2)';
dif3 = allAffSynPos(:,3) - allPostPos(:,3)';
distsAff2Post = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;

threshDistsAff2Syn = distThresh-distsAff2Post;
threshDistsAff2Syn(threshDistsAff2Syn<0) = 0;
sumWeightAll = sum(threshDistsAff2Syn,1);
hist(sumWeightAll)


maxWeight = max(sumWeightAll(:));
sumWeightSynN = sumWeightSyn/ maxWeight * 100;
sumWeightAllN = sumWeightAll/ maxWeight * 100;

%%Find fractions of synapses
checkRange = [1:max(sumWeightAllN(:))];
bin = 10;

numAll = zeros(length(checkRange),1);
numSyn = numAll;
for i = 1:length(checkRange)
    numAll(i) = sum((sumWeightAllN> (checkRange(i)-bin)) & (sumWeightAllN<=checkRange(i)+bin));
    numSyn(i) = sum((sumWeightSynN> (checkRange(i)-bin)) & (sumWeightSynN<=checkRange(i)+bin));
end
plot(numAll/max(numAll),'k');
hold on
plot(numSyn/max(numSyn),'g');
plot(numSyn./numAll,'r')
hold off

synProb = [0; numSyn./numAll]; %shift one to the right to allow for easy 0 values


%%Assign synapse probabilities to all nodes

for i = 1:length(skels)
    
    s = skels(i);
    pos = s.pos;
        
    %%Find real proximity presynaptic syn nodes to all post synaptic cell nodes
    dif1 = allAffSynPos(:,1) - pos(:,1)';
    dif2 = allAffSynPos(:,2) - pos(:,2)';
    dif3 = allAffSynPos(:,3) - pos(:,3)';
    distsAff2Post = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
    
    threshDistsAff2Syn = distThresh-distsAff2Post;
    threshDistsAff2Syn(threshDistsAff2Syn<0) = 0;
    sumWeightAll = sum(threshDistsAff2Syn,1);
    sumWeightAllN = sumWeightAll/ maxWeight * 100;
    weightInd = ceil(sumWeightAllN)+1;
    %weightInd(weightInd>100) = 100; %!!!!!!!!!!!
    
    
    skels(i).nodeSynProb = synProb(weightInd);
       
    
end

%% Get real values


meanDistReal = zeros(length(useID),1);
synNums = meanDistReal;
for i = 1:length(useID)
    targ = find(skelIDs==useID(i),1);
    s = skels(targ);
    useNodes = 1:length(s.nid);
    nodeNum = length(useNodes);
    synNum = length(s.sNid);
    synNums(i) = synNum;
    
    sPos = s.pos(s.sNid,:);
    
    %%Get distances
    dif1 = sPos(:,1) - sPos(:,1)';
    dif2 = sPos(:,2) - sPos(:,2)';
    dif3 = sPos(:,3) - sPos(:,3)';
    dists = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
    meanDist = sum(dists(:))/(synNum^2-synNum);
    meanDistReal(i) = meanDist;
    
    if 1
        clf
        scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
        hold on
        scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
        hold off
        title(sprintf('cell %d, %d synapses, mean distance = %0.1f',useID(i),synNum,meanDist))
        drawnow
    end
    
    
end


%% Run randomization, do each cell separately for speed
tic
reps = 1000
meanDistRand = zeros(length(useID),reps);
for i = 1:length(useID)
    disp(sprintf('running randomization for %d, %d of %d',useID(i),i,length(useID)))
    targ = find(skelIDs==useID(i),1);
    s = skels(targ);
    useNodes = 1:length(s.nid);
    useWeights = s.nodeSynProb(useNodes);
    nodeNum = length(useNodes);
    synNum = length(s.sNid);
    
    if 1
        clf
        sPos = s.pos(s.sNid,:);
        
        scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
        hold on
        scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
        scatter3(allPostSynPos(:,1),allPostSynPos(:,2),allPostSynPos(:,3),'b','o')

        drawnow
    end
        
    for r = 1:reps
        
        %sNidR = useNodes(randperm(nodeNum,synNum));
        
        sNidR = randsample(useNodes,synNum,true,useWeights);
        sPos = s.pos(sNidR,:);
        
        
        if r<10
            scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','g')
            drawnow
        end
        
        if 0
            scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
            hold on
            scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
            hold off
            pause(.1)
            
        end
        
        %%Get distances
        dif1 = sPos(:,1) - sPos(:,1)';
        dif2 = sPos(:,2) - sPos(:,2)';
        dif3 = sPos(:,3) - sPos(:,3)';
        dists = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
        meanDist = sum(dists(:))/(synNum^2-synNum);
        meanDistRand(i,r) = meanDist;
        
    end
    pause(.1)
end

toc
%% analyze results
mMeanDistReal = mean(meanDistReal)
mMeanDistRandom = mean(meanDistRand,1);

sVals = sort(mMeanDistRandom,'ascend');

P = mean(sVals<=mMeanDistReal)
ci95 = [sVals(round(0.025*reps)) sVals(round(0.975*reps))]















