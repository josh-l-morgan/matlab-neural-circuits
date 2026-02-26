%%Make predictions of how many synapses each pair of cells should form
%%based of the distribution of their skeletons. -
%-Skeletons are filtered



%% Plot axon to tcr
clear all
%MPN = GetMyDir;
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

seedList = [108 201 109 903 907]
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
preList = useList.preList;

%%

anchorScale = [-.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];


%% Get cells

axList = preList;
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


%%  Get Skeletons

for i = 1 : length(axList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,axList(i));
    load(fileName)
    axSkel(i).skel = cellStruct.skel;
    axSkel(i).name = axList(i);
    allStruct(i) = cellStruct;
end



for a = 1:length(axList)

    cellStruct = allStruct(i);
    arbor = cellStruct.arbor;
    nearestNode = arbor.vox.nearestNode
    subs = scaleSubs(arbor.vox.subs,voxelScale);
    scatter(subs(:,1),subs(:,2),'.')
    I = showSubs(ceil(abs(subs)));
    
 preDists = sqrt((preNodes(:,1)-synAnchors(i,1)).^2 + (preNodes(:,2)-synAnchors(i,2)).^2 + ...
            (preNodes(:,3)-synAnchors(i,3)).^2);
        
        if sum([min(preDists) min(postDists)])<10
            goodSyn(i) = 1;
        
        synNodeId = [find(postDists == min(postDists),1) find(preDists == min(preDists),1)];
        synNodeSub = [postNodes(synNodeId(1),:); preNodes(synNodeId(2),:)];

        end
end

