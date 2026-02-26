
% %%
% F = -kX
% k = stiffness
% X = proportional to distance;


%% Get data
loadData = 1;
if loadData
    %clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [ 108  201 109 907 903];
    %seedList = [ 108  201 109 ];

    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end

springDir = 'D:\LGNs1\Analysis\springDat\highlightGiantsText\'
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

con2 = con;
con2  = con2;
for i = 1:length(seedList)
    con2(:,find(useList.postList == seedList(i))) = 0;
end


minEdge = 2;
minSyn = 2;
minCon = 2;
binaryMat = 0;
con2(con2<minCon) = 0;


numIn = sum(con2>0,1);
synIn = sum(con2,1);
useIn = (numIn>=minEdge) & (synIn>=minSyn);

numOn = sum((con2>0).*repmat(useIn,[size(con2,1),1]),2);
synOn = sum((con2).*repmat(useIn,[size(con2,1),1]),2);
useOn = (numOn>=minEdge) & (synOn>=minSyn);

nodeIDs = [useList.preList(useOn) setdiff(useList.postList(useIn),seedList)]
nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList(useIn),seedList)*0+2]


nodeIDs = [useList.preList(useOn) setdiff(useList.postList(useIn),10000000)]
nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList(useIn),10000000)*0+2]



nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);


seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .3 1; 0 .3 1] * .6; 

nodeCol = zeros(nodeNum,3);
for s = 1:length(seedList)
    useCol = seedCol(s,:);
    for n = 1:nodeNum
        nodeCol(n,:) = nodeCol(n,:) + useCol *  ...
            double(sum((allEdges(:,1)==nodeIDs(n)) & (allEdges(:,2)) == seedList(s))>0);
    end
end



crossoverAxons = [2032	2033	2034	2035]
gotList = getList_giantBoutons(MPN);
labelCells = gotList;
for i = 1:length(labelCells)
   nodeCol(find(nodeIDs==labelCells(i)),:) = nodeCol(find(nodeIDs==labelCells(i)),:) + .5;    
end


nodeCol(nodeCol>1) = 1;
nodeCol(nodeCol<0) = 0;




%% Create springDat (need nodeCol, nodeIDs,

springIn.nodeIDs = nodeIDs;
springIn.allEdges = allEdges;
springIn.allWeights = allEdges * 0+1;
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = seedList;
springDat = springParameters_fiveSeedsText(springIn);

if binaryMat
    springDat.edges.ew = springDat.edges.ew>0;
end

%% Run springs

for rerun = 1: 1
    
    allResults{rerun} = runSprings(springDat);
    set(gcf, 'InvertHardCopy', 'off');
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    
    epsName = sprintf('%sspringRun_%03.0f.eps',springDir,rerun);
    print(gcf, epsName, '-depsc2','-painters')
    
end
