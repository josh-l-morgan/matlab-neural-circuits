%%Weight synapses according to their distance to the cell body

%% Plot axon to tcr
clear all
MPN = GetMyDir;
load([MPN 'obI.mat'])
shouldPause = 1;
anchorScale = [-.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];

%% variables
seedList = [108]
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

%% Get cells
tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
linList = obI.nameProps.cellNum(obI.nameProps.lin);
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

preList = [];
for i = 1:length(seedList)
    isPost = edges(:,1) == seedList(i);
    preList{i} = unique(edges(isPost,2));
end
%axList =    setxor(preList{1},preList{2});  %dont use a
axList =    cat(1,preList{:});  %dont use a

axList = intersect(unique(axList),rgcList);
axList = axList((axList>=1000) & (axList<5000));
axList = setdiff(axList,noSkel);
disp(sprintf('Results calculated without %d',noSkel));

% axList = setdiff(axList, crossoverAxons)
% axList = axList((axList>=1000) & (axList<5000));
% 
% axList = axList(axList<2000);
% axList = [axList crossoverAxons]

postList = [];
for i = 1:length(axList)
    isPre = edges(:,2) == axList(i);
    postList = [postList; edges(isPre,:)];
end
cellList = intersect(unique(postList),tcrList);
cellList = cellList((cellList>0) & (cellList<1000));

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end



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


%% Map existing synapses
rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors(:,1) = rawSynAnchors(:,1).*anchorScale(1);
synAnchors(:,2) = rawSynAnchors(:,2).*anchorScale(2);
synAnchors(:,3) = rawSynAnchors(:,3).*anchorScale(3);

goodSyn = zeros(size(synAnchors,1),1)>0;
for i = 1:size(synAnchors,1)
    if sum(cellList == edges(i,1)) & sum(axList == edges(i,2))
    preSkel = axSkel(find(axList == edges(i,2))).skel;
    postSkel = cellSkel(find(cellList == edges(i,1))).skel;
    
    preNodes = preSkel.node2subs;
    preNodes  = [preNodes(:,1) * voxelScale(1) preNodes(:,2) * voxelScale(2) preNodes(:,3) * voxelScale(3)];
    postNodes = postSkel.node2subs;
    postNodes  = [postNodes(:,1) * voxelScale(1) postNodes(:,2) * voxelScale(2) postNodes(:,3) * voxelScale(3)];
   
    postDists = sqrt((postNodes(:,1)-synAnchors(i,1)).^2 + (postNodes(:,2)-synAnchors(i,2)).^2 + ...
        (postNodes(:,3)-synAnchors(i,3)).^2);
    
    preDists = sqrt((preNodes(:,1)-synAnchors(i,1)).^2 + (preNodes(:,2)-synAnchors(i,2)).^2 + ...
        (preNodes(:,3)-synAnchors(i,3)).^2);
    
    if sum([min(preDists) min(postDists)])<10
        goodSyn(i) = 1;        
    end
    
    synNodeId(i,:) = [find(postDists == min(postDists),1) find(preDists == min(preDists),1)];
    synNodeSub(i,:,:) = [postNodes(synNodeId(i,1),:); preNodes(synNodeId(i,2),:)];
    end
end

synDists = sqrt((synNodeSub(goodSyn,1,1)-synNodeSub(goodSyn,2,1)).^2 + ...
    (synNodeSub(goodSyn,1,2)-synNodeSub(goodSyn,2,2)).^2 + ...
    (synNodeSub(goodSyn,1,3)-synNodeSub(goodSyn,2,3)).^2);



%% count synapses with distance to cb
seedAnchor = obI.cell.anchors(obI.cell.name == seedList(1),:);
seedAnchor = seedAnchor.*anchorScale;

%%%%%%%   exp(1)^ (-1* dists/distConstant)
distConstant = 100;
for i = 1:length(axList)
    isBoth = find((edges(:,1) == seedList(1)) & (edges(:,2) == axList(i)));
    synNum(i) = length(isBoth);
    synPositions = synAnchors(isBoth,:);
    synDists = sqrt((synPositions(:,1) - seedAnchor(1)).^2 + ...
        (synPositions(:,2) - seedAnchor(2)).^2 + ...
        (synPositions(:,3) - seedAnchor(3)).^2);
    synWeights = exp(-1 * synDists/distConstant);
    meanDist(i) = mean(synDists);
    synW(i) = sum(synWeights);
    axWeights{i} = synWeights;
    if mean(synDists)>80
        break
    end
   
end


scatter(synNum,synW)
hold on
plot([0 max(synNum)],[0 max(synNum)])
hold off

scatter(synNum,meanDist)

%% Hist connectivity

preHits = synW;

binHit = 0:max(preHits)
clear survivalC
for i = 1:length(binHit)
    
survivalC(i) = sum(preHits >=binHit(i));
end
plot(binHit,survivalC,'r','LineWidth',5)
hold on
histSyn = histc(preHits,binHit)
bar(binHit,histSyn,'b')
set(gca,'XTick',binHit)
set(gca,'YTick',[0:max(survivalC)])
xlim([0 binHit(end)+1])
hold off
