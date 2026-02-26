%%Find cell positions


%% Plot axon to tcr
clear all
load('MPN.mat')
load([MPN 'obI.mat'])
% anchorScale =  ;
% voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];
voxelScale = obI.em.dsRes * 2; %adjusted for skeleton down samp

%% variables

seedList = [108 201 903 907]
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);

axList = useList.preList;
cellList = useList.postList;
synMat = useList.con;

synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end

skelOverlapPred.con = con;


%% Map existing synapses
rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
%synAnchors = scaleSubs(double(rawSynAnchors),anchorScale);

dSamp =  (obI.em.res .* [4 4 1])./1000;%./obI.em.dsRes;

synAnchors(:,1) = rawSynAnchors(:,1)*dSamp(1);
synAnchors(:,2) = rawSynAnchors(:,2)*dSamp(2);
synAnchors(:,3) = rawSynAnchors(:,3)*dSamp(3);

%% Find ax positions
axPos = [];
for a = 1:length(axList)
   isPre = edges(:,2) == axList(a);
   axSyn = synAnchors(isPre,:);
   badSyn = sum(axSyn<=0,2);
   axSyn = axSyn(~badSyn,:);
   
   axPos(a,:) = mean(axSyn,1);
end

%% find cell positions
cellPos = [];
for c = 1:length(cellList)
    cellPos(c,:) = obI.cell.anchors(obI.cell.name == cellList(c),:);
end


cellPos(:,1) = cellPos(:,1)*dSamp(1);
cellPos(:,2) = cellPos(:,2)*dSamp(2);
cellPos(:,3) = cellPos(:,3)*dSamp(3);


cellIDs = [axList'; cellList'];
cellProps = [axPos; cellPos];








