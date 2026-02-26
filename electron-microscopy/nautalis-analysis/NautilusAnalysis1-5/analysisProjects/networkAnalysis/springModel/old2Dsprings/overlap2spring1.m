
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
    oldEdges = obI.nameProps.edges(:,[2 1]);
    
    load([MPN 'skelOverlapPred.mat'])
    
    useList = skelOverlapPred.useList;
    con = useList.con;
    nodeIDs = [useList.preList setdiff(useList.postList,useList.seedList)]
    nodeNum = length(nodeIDs);
    nodeType(1:length(nodeIDs)) = 2;
    nodeType(1:length(useList.preList))=1;
    lookUpID(nodeIDs+1) = 1:length(nodeIDs);
    
    predCon = skelOverlapPred.axCellLength(:,:,5);
    predCon = predCon/sum(predCon(:));
    predCon = predCon * sum(con(:));
    predCon = round(predCon);
    image(predCon * 50)
    
    %%Choose edges
    %eList = con2syn(predCon);
    eList = con2syn(con);
    allEdges = [useList.preList(eList(:,1))' useList.postList(eList(:,2))'];
    
    
end

springDir = 'D:\LGNs1\Analysis\springDat\overlap_real2\'
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
%
% con = useList.con;
% nodeIDs = [1:size(con,1) (size(con,1)+1):(size(con,1)+size(con,2))]
%
% nodeType = nodeIDs*0+2;





%% Node colors
seedList = [108 201];
nodeCol = zeros(nodeNum,3);
for s = 1:length(seedList)
    for n = 1:nodeNum
        nodeCol(n,s) = sum((oldEdges(:,1)==nodeIDs(n)) & (oldEdges(:,2)) == seedList(s))>0;
    end
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
springDat = springParameters_overlap(springIn);

%% Run springs

for rerun = 1: 10
    
    allResults{rerun} = runSprings(springDat);
    set(gcf, 'InvertHardCopy', 'off');
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    
    epsName = sprintf('%sspringRun_%03.0f.eps',springDir,rerun);
    print(gcf, epsName, '-depsc2','-painters')
    
end
