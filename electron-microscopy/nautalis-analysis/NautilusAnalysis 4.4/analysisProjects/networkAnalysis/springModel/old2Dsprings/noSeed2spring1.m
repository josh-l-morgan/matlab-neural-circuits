
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
    seedList = [108 201];
    useList = obI2cellList_seedInput(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end

springDir = 'D:\LGNs1\Analysis\springDat\noSeedCrossGlom\'
if ~exist(springDir,'dir'), mkdir(springDir), end



%% get attributes

%cellIDs = [useList.preList useList.postList];

%
% tcrList = cellIDs(obI.nameProps.tcr(obI.cell.mainObID));
% rgcList = cellIDs(obI.nameProps.rgc(obI.cell.mainObID));
% linList = cellIDs(obI.nameProps.lin(obI.cell.mainObID));
% ldmList = cellIDs(obI.nameProps.ldm(obI.cell.mainObID));
% sdmList = cellIDs(obI.nameProps.sdm(obI.cell.mainObID));
% %unkList = obI.nameProps.cellNum(obI.nameProps.unk);
% allLab = unique([tcrList rgcList linList ldmList sdmList]);
% allIDs = unique(obI.nameProps.cellNum);
% allIDs = allIDs(allIDs>9);
% noLab = setdiff(allIDs,allLab);
%
%
% tcrList = intersect(tcrList,allIDs);
% rgcList = intersect(rgcList,allIDs);
% linList = intersect(linList,allIDs);
% ldmList = intersect(ldmList,allIDs);
% sdmList = intersect(sdmList,allIDs);
% allLab = intersect(allLab,allIDs);
% noLab = intersect(noLab,allIDs);

%% Filter for convergence

con = useList.con;
numIn = sum(con>0,1);
useIn = numIn>1;
numOn = sum((con>0).*repmat(useIn,[size(con,1),1]),2);
useOn = numOn>1;

nodeIDs = [useList.preList(useOn) setdiff(useList.postList(useIn),seedList)]
nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList(useIn),seedList)*0+2]
nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);

nodeCol = zeros(nodeNum,3);
for s = 1:length(seedList)
    for n = 1:nodeNum
        nodeCol(n,s) = sum((allEdges(:,1)==nodeIDs(n)) & (allEdges(:,2)) == seedList(s))>0;
    end
end

crossoverAxons = [2032	2033	2034	2035]

labelCells = crossoverAxons;
for i = 1:length(labelCells)
   nodeCol(find(nodeIDs==labelCells(i)),3) = 1;    
end


nodeCol(:,2) = nodeCol(:,2)*.75;
nodeCol(nodeCol==0) = 0;
nodeCol(nodeCol>1) = 1;
nodeCol(nodeCol<0) = 0;




%% Create springDat (need nodeCol, nodeIDs,

springIn.nodeIDs = nodeIDs;
springIn.allEdges = allEdges;
springIn.allWeights = allEdges * 0+1;
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = seedList;
springDat = springParameters_zeroSeed(springIn);

%% Run springs

for rerun = 1: 10
    
    allResults{rerun} = runSprings(springDat);
    set(gcf, 'InvertHardCopy', 'off');
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    
    epsName = sprintf('%sspringRun_%03.0f.eps',springDir,rerun);
    print(gcf, epsName, '-depsc2','-painters')
    
end
