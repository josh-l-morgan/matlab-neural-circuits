

%% Plot axon to tcr
clear all
%MPN = GetMyDir;
load('MPN.mat')
load([MPN 'obI.mat'])
shouldPause = 10;
anchorScale = [-.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];


%% variables

shiftList = [-20 -10 -5 -4 -3 -2 -1.5 -1 -.5 -.25 0 .25 .5 1 1.5 2 3 4 5 10 20]

skelOverlapPred.predictionType = 'skelOverlap';
skelOverlapPred.shiftList = shiftList;

firstPredThresh = 0;
binWidth = .1;
skelOverlapPred.binWidth = binWidth;
seedList = [201]
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

% muchTraced = [106 107 108 109 111 112 117 120 123 129 133 134 148 156 159 162 163 169 ...
%     170 201 203 205 206 207 210 212 213 215 216 218];
% skelOverlapPred.muchTraced = muchTraced;

%% Get cells
useList =  obI2cellList_seedInput_RGC_TCR(obI,seedList);
axList = useList.preList;
cellList = useList.postList;
synMat = useList.con;
axList = setdiff(axList,noSkel);


useList.preList = axList;
skelOverlapPred.useList = useList;

synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


disp(sprintf('Results calculated without %d',noSkel));

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end

skelOverlapPred.con = con;

%%  Get Skeletons

for i = 1 : length(axList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,axList(i));
    load(fileName)
    axSkel(i).skel = cellStruct.skel;
    axSkel(i).name = axList(i);
end


for i = 1 : length(cellList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,cellList(i));
    load(fileName)
    cellSkel(i).skel = cellStruct.skel;
    cellSkel(i).name = cellList(i);
end

%% Create filtered node list
minLength = 3; %minimum branch length to be counted
clear axNodes cellNodes

axLength = zeros(length(axList),1);
cellLength = zeros(length(cellList),1);
for i = 1:length(axSkel)
    
    skel = axSkel(i).skel;
    nodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keepNodes = [];
    for b = 1:length(skel.bones);
        bone = skel.bones(b);
        L = 0;
        for e = 1:size(bone.edges,1);
            edgeNodes = nodeSubs(bone.edges(e,:),:);
            edgeDif = diff(edgeNodes,1);
            edgeLength = sqrt(sum(edgeDif.^2));
            L = L + edgeLength;
        end
        disp(L)
        boneLengths(b) = L;
        if L>minLength
            keepNodes = cat(2,keepNodes,bone.nodes);
            axLength(i) = axLength(i) + L;
            
        end
    end
    keepNodes = unique(keepNodes);
    
    
    allNodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keptNodeSubs = allNodeSubs(keepNodes,:);
    scatter(allNodeSubs(:,3),allNodeSubs(:,2),'.','k')
    hold on
    scatter(keptNodeSubs(:,3),keptNodeSubs(:,2),'.','r')
    hold off
    pause(.01)
    
    axNodes{i} = keptNodeSubs;
end

for i = 1:length(cellSkel)
    
    skel = cellSkel(i).skel;
    nodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keepNodes = [];
    for b = 1:length(skel.bones);
        bone = skel.bones(b);
        L = 0;
        for e = 1:size(bone.edges,1);
            edgeNodes = nodeSubs(bone.edges(e,:),:);
            edgeDif = diff(edgeNodes,1);
            edgeLength = sqrt(sum(edgeDif.^2));
            L = L + edgeLength;
        end
        disp(L)
        boneLengths(b) = L;
        if L>minLength
            keepNodes = cat(2,keepNodes,bone.nodes);
            cellLength(i) = cellLength(i) + L;
        end
    end
    keepNodes = unique(keepNodes);
    
    
    allNodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keptNodeSubs = allNodeSubs(keepNodes,:);
    scatter(allNodeSubs(:,3),allNodeSubs(:,2),'.','k')
    hold on
    scatter(keptNodeSubs(:,3),keptNodeSubs(:,2),'.','r')
    hold off
    pause(.01)
    
    cellNodes{i} = keptNodeSubs;
    
end




%% Map existing synapses
rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors = scaleSubs(double(rawSynAnchors),anchorScale);

goodSyn = zeros(size(synAnchors,1),1)>0;

for i = 1:size(con,1);
    for p = 1:size(con,2)
        conDists{i,p} = [];
    end
end
goodCon = con * 0;

for i = 1:size(synAnchors,1)
    if sum(cellList == edges(i,1)) & sum(axList == edges(i,2))
        %         preSkel = axSkel(find(axList == edges(i,2))).skel;
        %         postSkel = cellSkel(find(cellList == edges(i,1))).skel;
        %
        %         preNodes = scaleSubs(double(preSkel.node2subs),voxelScale);
        %         postNodes = scaleSubs(double(postSkel.node2subs),voxelScale);
        
        %     preNodes = preSkel.node2subs;
        %     preNodes  = [preNodes(:,1) * voxelScale(1) preNodes(:,2) * voxelScale(2) preNodes(:,3) * voxelScale(3)];
        %     postNodes = postSkel.node2subs;
        %     postNodes  = [postNodes(:,1) * voxelScale(1) postNodes(:,2) * voxelScale(2) postNodes(:,3) * voxelScale(3)];
        %
        preTarg = find(axList == edges(i,2));
        preNodes = axNodes{preTarg};
        postTarg = find(cellList == edges(i,1));
        postNodes = cellNodes{postTarg};
        postDists = sqrt((postNodes(:,1)-synAnchors(i,1)).^2 + (postNodes(:,2)-synAnchors(i,2)).^2 + ...
            (postNodes(:,3)-synAnchors(i,3)).^2);
        
        preDists = sqrt((preNodes(:,1)-synAnchors(i,1)).^2 + (preNodes(:,2)-synAnchors(i,2)).^2 + ...
            (preNodes(:,3)-synAnchors(i,3)).^2);
        
        notEmpty = ~isempty(preDists) & ~isempty(postDists);
        
        if (sum([min(preDists) min(postDists)])<10) & notEmpty
            goodSyn(i) = 1;
            
            synNodeId = [find(postDists == min(postDists),1) find(preDists == min(preDists),1)];
            synNodeSub = [postNodes(synNodeId(1),:); preNodes(synNodeId(2),:)];
            synDif = diff(synNodeSub,1);
            synDist = sqrt(sum(synDif.^2));
            
            conDists{preTarg,postTarg} = cat(1,conDists{preTarg,postTarg},synDist);
            goodCon(preTarg,postTarg) = goodCon(preTarg,postTarg) + 1;
        else
            synNodeSub = [0 0 0;0 0 0]
            
        end
    else
        synNodeSub = [0 0 0; 0 0 0];
    end
    
    
    
end

skelOverlapPred.conDists = conDists;
skelOverlapPred.goodCon = goodCon;

%% axon expanse

targCell = seedList;
%%Get distribution of distances for relevant cells
allDists  = zeros(10000,1);
lastDist = 0;
binWidth = 10
histBin = [0:binWidth:100];
tic
for i = 1:length(axList)
    disp(sprintf('checking axon %d of %d',i,length(axList)))
    preNodes = axNodes{i};
    if isempty(preNodes)
        maxDist(i) = 0;
        minDist(i) = 0;
    else
    
    useSyn = find((edges(:,1) == targCell) & (edges(:,2) == axList(i)));
    synNodes = synAnchors(useSyn,:);
    maxSyn = useSyn*0;
    for p = 1: length(useSyn)
        pairDists = sqrt((preNodes(:,1)-synNodes(p,1)).^2 + (preNodes(:,2)-synNodes(p,2)).^2 + ...
            (preNodes(:,3)-synNodes(p,3)).^2);
        %         histDists = histc(pairDists(:),histBin);
        %         axCellHists{i,p} = histDists;
        maxSyn(p) = max(pairDists);
    end
    if isempty(maxSyn)
        maxDist(i) = 0;
        minDist(i) = 0;
    else
        maxDist(i) = max(maxSyn);
        minDist(i) = min(maxSyn);
    end
    
    end
end
toc
axList(maxDist>100)

sum(maxDist>0)
sum(maxDist>100)
sum(minDist>100)


%% Show axons that are minimum distance to syn with seed
for i = 1:length(axList)
    preNodes = axNodes{i};
    maxDist(i)
    axLength(i)
    if maxDist(i)>100
        scatter(preNodes(:,3),preNodes(:,2),'.','k')
    else
        scatter(preNodes(:,3),preNodes(:,2),'.','r')
        
    end
    pause(.1)
    
end

plot(maxDist,'k')
hold on
plot(minDist,'r')
hold off


%% Length attatched to 108;

%pre108 = unique(edges(edges(:,1) == 108,2));
targSeed = find(cellList == seedList);
preSeed = find(con(:,targ108)>0);
lengthDist = axLength(preSeed);
sum(lengthDist)


%%  % Show cells
%{

for i = 1:length(axList)
    disp(sprintf('checking axon %d of %d',i,length(axList)))
    preSkel = axSkel(i).skel;
    preNodes = scaleSubs(preSkel.node2subs,voxelScale);
    %preNodes  = [preNodes(:,1) * voxelScale(1) preNodes(:,2) * voxelScale(2) preNodes(:,3) * voxelScale(3)];
    tic
    for p = 1: length(cellList)
        
        postSkel = cellSkel(p).skel;
        postNodes = scaleSubs(postSkel.node2subs,voxelScale);
        %postNodes  = [postNodes(:,1) * voxelScale(1) postNodes(:,2) * voxelScale(2) postNodes(:,3) * voxelScale(3)];
        
        scatter(preNodes(:,3),preNodes(:,2),'.','k')
        hold on
        scatter(postNodes(:,3),postNodes(:,2),'.','r')
        hold off
        disp(sprintf('ax %d and cell %d',axList(i),cellList(p)))
        ylim([0 400])
        xlim([0 300])
        
        pause(.01)
        
        
    end
end

%}
%% Find length without primaries and secondaries
minLength = 3; %minimum branch length to be counted
clear axNodes cellNodes

axLength = zeros(length(axList),1);
cellLength = zeros(length(cellList),1);
for i = 1:length(axSkel)
    
    skel = axSkel(i).skel;
    nodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keepNodes = [];
    for b = 1:length(skel.bones);
        bone = skel.bones(b);
        L = 0;
        for e = 1:size(bone.edges,1);
            edgeNodes = nodeSubs(bone.edges(e,:),:);
            edgeDif = diff(edgeNodes,1);
            edgeLength = sqrt(sum(edgeDif.^2));
            L = L + edgeLength;
        end
        disp(L)
        boneLengths(b) = L;
        if L>minLength
            keepNodes = cat(2,keepNodes,bone.nodes);
            axLength(i) = axLength(i) + L;
            
        end
    end
    keepNodes = unique(keepNodes);
    
    
    allNodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keptNodeSubs = allNodeSubs(keepNodes,:);
    scatter(allNodeSubs(:,3),allNodeSubs(:,2),'.','k')
    hold on
    scatter(keptNodeSubs(:,3),keptNodeSubs(:,2),'.','r')
    hold off
    pause(.01)
    
    axNodes{i} = keptNodeSubs;
end

%% compare length to synapse number

useLength = tertAxLength; % from editNodes
lengthMat = repmat(useLength, [ 1 size(synMat,2)]);

[gotSeed notSeed] = setdiff(cellList,seedList);

lengthList = lengthMat(:,notSeed);
synList = synMat(:,notSeed);


filtList = (synList>0) & ( lengthList >100);
lengthList = lengthList(filtList>0);
synList = synList(filtList>0);


scatter(lengthList,synList)


synPerMM = synList./lengthList * 1000

hist(synPerMM,[0:2:60])
mean(synPerMM)


synCount = sum(synMat,2)
scatter(tertAxLength,synCount)

sum(tertAxLength)










