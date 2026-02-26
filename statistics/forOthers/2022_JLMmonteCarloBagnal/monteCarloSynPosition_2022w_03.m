%%Monte carlo synapses across nodes

%%Note current nodes are not evenly distributed so that conditioning of
%%skeleton will be required
%%should then reassign synapses to new nodes



%% define data locations
clear all
datFolder = 'E:\Bagnal\2022\';

connectorFile = [datFolder 'connectors_table.csv'];
skelCSVFolder = [datFolder 'vestibular_neurons_CSV_files\']
affCSVFolder = [datFolder 'afferent_files_SWC_CSV\']

%% Pick variables
useAffSyn = 0;
distThresh = 10;
skelGap = 1000; %nanometers between skeleton nodes
affGap = 1000; %same for afferents


mcd.var.useAffSyn = useAffSyn;
mcd.var.distThresh = distThresh;
mcd.var.skelGap = skelGap;
mcd.var.affGap = affGap;

%% Parse connections
cons = readtable(connectorFile,'Delimiter',',');
conLoc = cons{:,3};
conNode = cons{:,9};
conTarget = cons{:,10};
conSource = cons{:,5};
conSourceNode = cons{:,4};
conTargetNode = cons{:,9};

uTargets = unique(conTarget);

%%Parse synapse locations
cPos = zeros(length(conLoc),3);
for i = 1:length(conLoc)
    nam = conLoc{i};
    b1 = regexp(nam,'[');
    b2 = regexp(nam,']');
    coms = regexp(nam,',');
    X = str2num(nam(b1+1:coms(1)-1));
    Y = str2num(nam(coms(1)+2:coms(2)-1));
    Z = str2num(nam(coms(2)+2:b2-1));
    cPos(i,:) = [X Y Z];
    
end

mcd.con.cons = cons;
mcd.con.cPos = cPos;
mcd.con.conTarget = conTarget;
mcd.con.conSource = conSource;
mcd.con.conNode = conNode;
save([datFolder 'mcd.mat'],'mcd');

%% Parse skeletons
datDir = dir([skelCSVFolder '*.csv']);
sNams = {datDir.name};
clear skel skelIDs skelType
for i = 1:length(sNams)
    
    nam = sNams{i};
    
    disp(sprintf('reading %s,%d of %d posts',nam,i,length(sNams)));
    
    if regexp(nam,'VN')
        skels(i).skelType = 1;
    elseif regexp(nam,'VS')
        skels(i).skelType  = 2;
    elseif regexp(nam,'tangential')
        skels(i).skelType = 3;
    else
        skels(i).skelType = 0;
    end
    
    
    skelFile = [skelCSVFolder nam];
    skelT = readtable(skelFile,'Delimiter',',');
    
    skelID =  skelT{1,11}; %get ID of cell
    skelIDs(i) = skelID;   
    
    parent = skelT{:,3};
    pos = skelT{:,[ 5 6 7]}; %get possition of nodes

    %%Asign IDs and pred
    nid = skelT{:,2}; % get IDs of nodes
    lookupNid = zeros(1,max(nid));
    lookupNid(nid) = 1:length(nid); % make up lookup table for conversion of node ids to node indexes
    root = find(parent<0);
    parent(root) = 1;
    pred = lookupNid(parent); % look up correct values for non-root nodes
    pred(root) = -1;
    
    %%resample arbor 
    [nPos nPred] = regularEdges(pos,pred,skelGap,1);
    newNid = 1:length(nPred);
    
    skels(i).pos = nPos;
    skels(i).pred = nPred;
    skels(i).oldPos = pos;
    skels(i).nid = nid;
   
    
    tags = skelT{:,13};
    
    %%Get synaptic inputs
    inp = find(conTarget == skelID);
    skels(i).inputID = inp;
    skels(i).inputNodeRaw = conNode(inp);
        
    %%Find closest node
    skels(i).sNidOld = lookupNid(conNode(inp)); %node index for synapse    
    oldSynPos = pos(skels(i).sNidOld,:);
    conPos = cPos(inp,:);
    skels(i).conPos = conPos;
    sNid = zeros(size(conPos,1),1);
    for c = 1:size(conPos,1);
        dists = sqrt((nPos(:,1)-conPos(c,1)).^2 +  (nPos(:,2)-conPos(c,2)).^2 +  ...
            (nPos(:,3)-conPos(c,3)).^2);
        sNid(c) = find(dists==min(dists),1);
    end
    skels(i).sNid = sNid;
     sPos = nPos(sNid,:);
    skels(i).sPos = sPos;
    
    
    %%Make edges (nx2)
    clear edges
    e1 = find(nPred>0);
    e2 = nPred(e1);
    skels(i).edges = [e1 e2];
    
end

useID = skelIDs([skels.skelType]>0);


%%Test skelT

mcd.skels = skels;
mcd.skelIDs = skelIDs;
save([datFolder 'mcd.mat'],'mcd');


%% Get distances between nodes
for i = 1:length(skels)
    edges = skels(i).edges;
    pos = skels(i).pos;
    nodes = 1:size(pos,1);
    e1 = edges(:,1);
    e2 = edges(:,2);
    lengths = sqrt((pos(e1,1)-pos(e2,1)).^2 + (pos(e1,2)-pos(e2,2)).^2 + (pos(e1,3)-pos(e2,3)).^2);
    % scatter(pos(:,1),pos(:,2),'.')
    
    
    
    num = length(nodes);
    linDist = zeros(num);
    for y = 1:num;
        pp = node2nodeDist(edges,lengths,nodes(y));
        linDist(y,:) = pp.dists;
        
        if 0 %show distances
            colormap jet(100)
            plot([pos(edges(:,1),1) pos(edges(:,2),1)]', [pos(edges(:,1),2) pos(edges(:,2),2)]','k')
            hold on
            s = scatter(pos(:,1),pos(:,2),'o','filled');
            scatter(pos(nodes(y),1),pos(nodes(y),2),200,'o','k')
            prop = linDist(y,:);
            s.CData = prop/max(prop) * 100;
            drawnow
            hold off
        end
        
    end
    skels(i).skel2skelLinDist = linDist;
    
    
    
end
mcd.skels = skels;
save([datFolder 'mcd.mat'],'mcd');


%% Parse afferents


dFold = dir([affCSVFolder '*.csv']);
aNams = {dFold.name};
clear aff affIDs affType
for i = 1:length(aNams)
    
    nam = aNams{i};
        disp(sprintf('reading %s,%d of %d afferents',nam,i,length(sNams)));

    if regexp(nam,'AnteriorMacula')
        affs(i).affType = 1;
    else
        affs(i).affType = 0;
    end
    
    
    skelFile = [affCSVFolder nam];
    skelT = readtable(skelFile,'Delimiter',',');
    parent = skelT{:,3};
    
    skelID =  skelT{1,11}; %get ID of cell
    affIDs(i) = skelID;
    pos = skelT{:,[ 5 6 7]}; %get possition of nodes
    
    
    nid = skelT{:,2}; % get IDs of nodes
    lookupNid = zeros(1,max(nid));
    lookupNid(nid) = 1:length(nid); % make up lookup table for conversion of node ids to node indexes
    pred = parent * 0 -1; %set all to -1
    pred(parent>0) = lookupNid(parent(parent>0)); % look up correct values for non-root nodes
    
   
    
    %%resample arbor 
    [nPos nPred] = regularEdges(pos,pred,affGap,1);
    newNid = 1:length(nPred);
    
    affs(i).pos = nPos;
    affs(i).pred = nPred;
    affs(i).oldPos = pos;
    affs(i).nid = nid;
   
        
    
    inp = find(conSource == skelID);
    affs(i).outputID = inp;
    affs(i).outputNodeRaw = conSourceNode(inp);
     
    %%Find closest node
%     affs(i).sNidOld = lookupNid(conSourceNode(inp)); %node index for synapse    
%     oldSynPos = pos(affs(i).sNidOld,:);
    conPos = cPos(inp,:);
    affs(i).conPos = conPos;
    sNid = zeros(size(conPos,1),1);
    for c = 1:size(conPos,1);
        dists = sqrt((nPos(:,1)-conPos(c,1)).^2 +  (nPos(:,2)-conPos(c,2)).^2 +  ...
            (nPos(:,3)-conPos(c,3)).^2);
        sNid(c) = find(dists==min(dists),1);
    end
    affs(i).sNid = sNid;
    sPos = nPos(sNid,:);
    affs(i).sPos = sPos;
    
    if 1 
        clf,hold on
        e1 = find(nPred>0);
        e2 = nPred(e1);
        plot3([nPos(e1,1) nPos(e2,1)]',[nPos(e1,2),nPos(e2,2)]',[nPos(e1,3),nPos(e2,3)]','r');
        scatter3(conPos(:,1),conPos(:,2),conPos(:,3),'o');
        scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'.','k');
        drawnow
    end
    
    
    
end

useAffID = affIDs([affs.affType]>0);

mcd.affs = affs;
mcd.affIDs = affIDs;
save([datFolder 'mcd.mat'],'mcd');

%% Find synaptic probability for nodes based on proximity to afferents


%%Collect
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

if useAffSyn
    prePos = allAffSynPos;
else
    prePos = allAffPos;
end



%%Assign synapse probabilities to all nodes

for i = 1:length(skels)
    disp(sprintf('measuring distances for %d of %d',i,length(skels)))
    s = skels(i);
    pos = s.pos;
        
    %%Find real proximity presynaptic syn nodes to all post synaptic cell nodes
    dif1 = prePos(:,1) - pos(:,1)';
    dif2 = prePos(:,2) - pos(:,2)';
    dif3 = prePos(:,3) - pos(:,3)';
    distsAff2Post = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
    
    threshDistsAff2Syn = distThresh-distsAff2Post;
    threshDistsAff2Syn(threshDistsAff2Syn<0) = 0;
    sumWeightAll = sum(threshDistsAff2Syn,1);
    skels(i).nodeWeights = sumWeightAll;
    
    
    %%Find real proximity presynaptic syn nodes to all post synaptic cell nodes
    sPos = s.pos(s.sNid,:);
    dif1 = prePos(:,1) - sPos(:,1)';
    dif2 = prePos(:,2) - sPos(:,2)';
    dif3 = prePos(:,3) - sPos(:,3)';
    distsAff2Post = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
    
    threshDistsAff2Syn = distThresh-distsAff2Post;
    threshDistsAff2Syn(threshDistsAff2Syn<0) = 0;
    sumWeightAll = sum(threshDistsAff2Syn,1);
    skels(i).synWeights = sumWeightAll;
        
end

%%Compare all weights
allSynWeights = [];
allNodeWeights = [];
for i = 1:length(skels)
    allSynWeights = cat(2,allSynWeights,skels(i).synWeights);
    allNodeWeights = cat(2,allNodeWeights,skels(i).nodeWeights);
end


%%normalize weights
maxWeight = max(allNodeWeights(:));
sumWeightSynN = allSynWeights/ maxWeight * 100;
sumWeightAllN = allNodeWeights/ maxWeight * 100;

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
pause(2)
synProb = [0; numSyn./numAll]; %shift one to the right to allow for easy 0 values


%%Assign synapse probabilities to all nodes

for i = 1:length(skels)
    
    s = skels(i);
    sumWeightAllN = skels(i).nodeWeights/ maxWeight * 100;
    weightInd = ceil(sumWeightAllN)+1;
    %weightInd(weightInd>100) = 100; %!!!!!!!!!!!
    skels(i).nodeSynProb = synProb(weightInd);      
    
end

mcd.skels = skels;
save([datFolder 'mcd.mat'],'mcd');

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

mcd.res.meanDistReal = meanDistReal;
save([datFolder 'mcd.mat'],'mcd');

%% Run randomization, do each cell separately for speed
tic
reps = 10000
meanDistRand = zeros(length(useID),reps);
for i = 1:length(useID)
    disp(sprintf('running randomization for %d, %d of %d',useID(i),i,length(useID)))
    targ = find(skelIDs==useID(i),1);
    s = skels(targ);
    useNodes = 1:size(s.pos,1);
    useWeights = s.nodeSynProb(useNodes);
    nodeNum = length(useNodes);
    synNum = length(s.sNid);
    
    if 1
        clf
        sPos = s.pos(s.sNid,:);
        
        scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
        hold on
        scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
        scatter3(allAffPos(:,1),allAffPos(:,2),allAffPos(:,3),'b','.')
        scatter3(allAffSynPos(:,1),allAffSynPos(:,2),allAffSynPos(:,3),'m','o')
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

mcd.res.meanDistRand = meanDistRand;
save([datFolder 'mcd.mat'],'mcd');

%% analyze results
mMeanDistReal = mean(meanDistReal)
mMeanDistRandom = mean(meanDistRand,1);

sVals = sort(mMeanDistRandom,'ascend');

P = mean(sVals<=mMeanDistReal)
ci95 = [sVals(round(0.025*reps)) sVals(round(0.975*reps))]

mcd.res.P = P;
mcd.res.ci95 = ci95;
mcd.res.sVals = sVals;
mcd.res.mMeanDistReal = mMeanDistReal;
mcd.res.testStat = mMeanDistReal;

save([datFolder 'mcd.mat'],'mcd');














