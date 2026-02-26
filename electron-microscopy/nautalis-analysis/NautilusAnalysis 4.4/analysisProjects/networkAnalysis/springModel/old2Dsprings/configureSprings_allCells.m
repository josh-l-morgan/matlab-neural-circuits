
% %%
% F = -kX
% k = stiffness
% X = proportional to distance;

%% Get data
loadData = 1;
if loadData
    clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    %seedList = [108 201 601];
    %useList = obI2cellList_seedInput(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
end

springDir = 'D:\LGNs1\Analysis\springDat\allCellsOnlyTCRRGC_2\'
if ~exist(springDir,'dir'), mkdir(springDir), end



%% shape attributes

cellIDs = obI.cell.name;

isTCR = obI.nameProps.tcr( obI.cell.mainObID)
isRGC = obI.nameProps.rgc( obI.cell.mainObID)
isLIN = obI.nameProps.lin( obI.cell.mainObID)
isLDM = obI.nameProps.ldm( obI.cell.mainObID)

useCell = (cellIDs>0) & (isTCR | isRGC);
nodeIDs = cellIDs(useCell);
nodeNum = length(nodeIDs);


nodes.type = isRGC(useCell) * 1 + isTCR(useCell) * 2 + isLIN(useCell) * 3;

con = zeros(length(nodeIDs),length(nodeIDs));
useCon = con * 0;
for i = 1:length(nodeIDs)
    for p = 1:length(nodeIDs)
        con(i,p) = con(i,p) + sum( (allEdges(:,1) == nodeIDs(i)) & (allEdges(:,2) == nodeIDs(p)));
        %con(i,p) = con(i,p) + sum( (allEdges(:,2) == nodeIDs(p)) & (allEdges(:,1) == nodeIDs(i)));
           useCon(i,p) = (p>i);

    end
end

seedList = [108 201 170];
clear direct indirect
for s = 1:length(seedList)
    targ = find(nodeIDs == seedList(s));
    direct(:,s) = con(:,targ)>0;
    repSeed = repmat(con(:,targ),[1,nodeNum]);
    shared = (repSeed>0).*con;
    sumShared = sum(shared,1);
    indirect(:,s) = sumShared/sum(sumShared);
end

nodeCol = direct + indirect;%(direct + indirect);
% nodeCol = nodeCol/max(nodeCol(:));
% nodeCol = (nodeCol>0) * .8 + .2;
nodeCol = double(nodeCol>0)*.7;
nodeCol = nodeCol(:,[1 3 2]);
nodeCol(nodeCol==0) = .2;
nodeCol((nodeIDs>=900 ) & (nodeIDs<1000),2) = 1


%%Get new edges
[e1 e2] = find((con.*useCon)>0);
ei = sub2ind(size(con),e1,e2);
ew = con(ei);
edgeCol = nodeCol(e1,:) + nodeCol(e2,:);
edgeCol(ew==0,:) = edgeCol(ew==0,:)*.2;
edgeCol(edgeCol<0) = 0;
edgeCol(edgeCol>1) = 1;

edges.all = allEdges;
edges.e1 = e1;
edges.e2 = e2;
edges.ew = ew;
edges.col = edgeCol;
edges.symmetry = ew*0+1;
%edges.width = edgeWidth;














%
% tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
% rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
% linList = obI.nameProps.cellNum(obI.nameProps.lin);
% synapses = obI.nameProps.edges;
% edges = synapses(:,1:2);


%
%
%     useList.preList = nodeIDs;
%     useList.postList = nodeIDs;
%

% useList.con = con;
% seedPref = seedPreferences(seedList,useList);
%
%     postSynPref = seedPref.sharedSynNorm(1,:)./sum(seedPref.sharedSynNorm,1);
%     preSynPref = seedPref.ax2seed(1,:)./sum(seedPref.ax2seed,1);
%     synPref = nodeIDs*0;
%     isPref = logical(nodeIDs*0);
%     for i = 1:length(useList.postList)
%        targ = find(nodeIDs==useList.postList(i));
%         if nodes.type(targ) == 1
%             synPref(targ) = preSynPref(i);
%         else
%             synPref(targ) = postSynPref(i);
%         end
%
%         if ~isnan(synPref(targ))
%             isPref(targ) = 1>0;
%         end
%     end
% %
%     nodeCol = zeros(nodeNum,3);
%     nodeCol(isPref,1) = synPref(isPref);
%     nodeCol(isPref,3) = 1-synPref(isPref);
%     nodeCol(~isPref,2) = 1;
%



%% Create springDat
% 
seed.list = seedList;
% seed.seedPref = synPref;
% seed.isPref = isPref;

edges.all = allEdges;

nodes.labelIDs = nodeIDs;
nodes.color = nodeCol;


%%set param
param.zeroSeeds = 0;
param.k = .5;
param.noise = [10 1 .1 0];
param.damp = .5;
param.disperse = 20;
param.reps = 2000;
param.fsize = 200;
param.speedMax = 1;
param.centerSpring = 80;
param.imageFreq = 10;

param.startRepulse = 500;%
param.repulse = 5; %1
param.minDist = 5;

%%Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');

group(1).ind = isRGC(useCell)>0;
group(2).ind = isTCR(useCell)>0;
group(3).ind = isLIN(useCell)>0;
group(1).marker = '^';
group(2).marker = 'o';
group(3).marker = 's';
group(1).size = 10;
grooup(2).size = 40;
group(3).size = 20;
groupNum = length(group);


springDat.group = group;
springDat.param = param;
springDat.nodes = nodes;
springDat.edges = edges;
springDat.seed = seed;



for rerun = 1: 10
    %allResults{rerun} = runSprings(springDat);
        allResults{rerun} = runSpringMass(springDat);

    
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    saveas(gcf,imageName,'png')
    %print(gcf, imageName, '-depsc2','-painters')
    set(gcf, 'InvertHardCopy', 'off');
    %print(gcf, imageName, '-dpng')
    
    %imageName = sprintf('%stestRun_%s.0f.png',springDir,'openGL');
    print(gcf,imageName,'-dpng','-r1600','-opengl') 
    
    %datName
    
end

springDat.results = allResults;
save([springDir 'springDat.mat'],'springDat')




