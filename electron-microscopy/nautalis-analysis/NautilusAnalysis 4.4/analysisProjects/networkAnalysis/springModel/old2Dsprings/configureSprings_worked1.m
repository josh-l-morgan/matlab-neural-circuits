
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

springDir = 'D:\LGNs1\Analysis\springDat\testSpringsSeedAx3\'
if ~exist(springDir,'dir'), mkdir(springDir), end



%% get attributes

cellIDs = [useList.preList useList.postList];


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

nodeCol(:,2) = nodeCol(:,2)*.75;
nodeCol(nodeCol==0) = 0;
nodeCol(nodeCol>1) = 1;
nodeCol(nodeCol<0) = 0;




%% Create springDat (need nodeCol, nodeIDs, 
% 

nodes.type =  [useList.preList*0+1 useList.postList*0+2];

seed.list = seedList;
% seed.seedPref = synPref;
% seed.isPref = isPref;

% 
% rawCon = zeros(nodeNum,nodeNum);
% useCon = rawCon * 0;
% for i = 1:length(nodeIDs)
%     for p = 1:length(nodeIDs)
%         rawCon(i,p) = rawCon(i,p) + sum((allEdges(:,1) == nodeIDs(i)) & (allEdges(:,2) == nodeIDs(p)));
%         %rawCon(i,p) = rawCon(i,p) + sum((allEdges(:,1) == nodeIDs(p)) & (allEdges(:,2) == nodeIDs(i)));
%         %useCon(i,p) = (p>i);
%     end
% end
% con  = rawCon;


[someEdges ew] = uniqueEdges(allEdges,nodeIDs);

%%Get new edges
[e1 e2] = find((con)>0);
ei = sub2ind(size(con),e1,e2);
ew = con(ei);
edgeCol = nodeCol(e1,:);
edgeCol(ew==0,:) = edgeCol(ew==0,:)*.2;
edgeCol(edgeCol<0) = 0;
edgeCol(edgeCol>1) = 1;

edges.all = allEdges;
edges.e1 = someEdges(:,1);
edges.e2 = someEdges(:,2);
edges.ew = ew;
edges.col = edgeCol;
edges.symmetry = ew*0;
edges.width = ew*0+.7;


nodes.labelIDs = nodeIDs;
nodes.color = nodeCol;


%%set param
param.k = 1;
param.disperse = 15;
param.centerSpring = 1;
param.centerSpring = .01;


param.reps = 2000;
param.zeroSeeds = 0;
param.noise = (param.reps:-1:1)/param.reps*10;
param.damp = .5;
param.fsize = 200;
param.speedMax = 1;
param.imageFreq = 10;

param.startRepulse = 500;%
param.repulse = 100; %1
param.minDist = 7;

%%Make groups
%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');


springDat.groupNum = 2;
group(1).ind = find(nodeType==1);
group(2).ind = find(nodeType==2);
%group(3).ind = lookUpID([linList sdmList ldmList noLab ]+1);
group(1).marker = '^';
group(2).marker = 'o';
%group(3).marker = 's';
group(1).size = 100;
grooup(2).size = 5;
%group(3).size = 60;
groupNum = length(group);
group(1).lineWidth = .5;
group(2).lineWidth = .7;

springDat.group = group;
springDat.param = param;
springDat.nodes = nodes;
springDat.edges = edges;
springDat.seed = seed;



for rerun = 1: 10
    allResults{rerun} = runSprings(springDat);
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    %saveas(gcf,imageName,'png')
    %print(gcf, imageName, '-depsc2','-painters')
    set(gcf, 'InvertHardCopy', 'off');
    %print(gcf, imageName, '-dpng')
    
    %imageName = sprintf('%stestRun_%s.0f.png',springDir,'openGL');
   
    print(gcf,imageName,'-dpng','-r1024','-opengl','-noui') 
        
    epsName = sprintf('%sspringRun_%03.0f.eps',springDir,rerun);
    print(gcf, epsName, '-depsc2','-painters')  
    
    %datName
    
end

return

imageName = sprintf('%sspringRun_%03.0f.pdf',springDir,rerun);
print(gcf,imageName,'-dpdf','-painter') 
print(gcf,imageName,'-dpdf') 


imageName = sprintf('%sspringRun_%03.0f.svg',springDir,rerun);
print(gcf,imageName,'-dsvg') 

-dformat 

imageName = sprintf('%sspringRun_%03.0f.emf',springDir,rerun);
print(gcf,imageName,'-dmeta') 

springDat.results = allResults;
save([springDir 'springDat.mat'],'springDat')



imageName = sprintf('%stestRun_%s.0f.png',springDir,'openGL');
print(gcf,imageName,'-dpng','-r1600','-opengl') 
    
    
opengl 'info'
opengl('save','hardware')
opengl('save','software')
imageName = sprintf('%stestRun_%s.0f.eps',springDir,'openGL');
print(gcf, imageName, '-depsc2','-painters','-cmyk')  
    
    
    %}
    
    
    
    








