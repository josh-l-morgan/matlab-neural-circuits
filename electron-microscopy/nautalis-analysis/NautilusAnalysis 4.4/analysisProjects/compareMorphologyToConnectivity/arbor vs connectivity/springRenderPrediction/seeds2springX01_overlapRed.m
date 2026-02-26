



%% Get data
loadData = 1;
if loadData
    clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [ 108 201 109 907 903];
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    %allEdges = obI.nameProps.edges(:,[2 1]);
    
    fileName = sprintf('%sskelOverlapPred.mat',MPN);
    load(fileName,'skelOverlapPred')
    useList = skelOverlapPred.useList;
    
end


springDir = 'D:\LGNs1\Analysis\springDat\test\';
springRes = 'D:\LGNs1\Analysis\springDat\results\';
if ~exist(springDir,'dir'), mkdir(springDir), end
if ~exist(springRes,'dir'), mkdir(springRes), end


%% Filter for convergence
conRaw = useList.con;

setEdge = 1; % 1

if setEdge == 1; % arrange based on overlap
    mask = conRaw*0+1;
    if 1 %zero seeds
        for i = 1:length(seedList)
            mask(:,find(useList.postList == seedList(i))) = 0;
        end
    end
    
    binRes = .1;
    [predSyn nodeFrac] = overlapPredictSynOfSubset(skelOverlapPred,binRes,mask);
    %con = round(predSyn);
    %con = conRaw;
    [con] = predNewCon(conRaw,predSyn,mask)
    con(~mask) = conRaw(~mask);

elseif setEdge ==2
    
     mask = conRaw*0+1;
    if 1 %zero seeds
        for i = 1:length(seedList)
            mask(:,find(useList.postList == seedList(i))) = 0;
        end
    end
    
    binRes = .1;
    [predSyn nodeFrac] = overlapPredictSynOfSubset(skelOverlapPred,binRes,mask);
    %con = round(predSyn);
    %con = conRaw;
    [con] = randNewCon(conRaw,rand(size(conRaw)),mask)
    con(~mask) = conRaw(~mask);

    
else
    con = conRaw;
end


[y x] = find(con);
allEdges = [];
for i = 1:length(y)
    newEdge = [useList.preList(y(i)) useList.postList(x(i))];
    allEdges = cat(1,allEdges,repmat(newEdge,[con(y(i),x(i)),1]));
end

%allEdges = con2syn(con>0)

con2 = con;
if 0 %zero seeds
    for i = 1:length(seedList)
        con2(:,find(useList.postList == seedList(i))) = 0;
    end
end

minEdge = 1;
minSyn = 1;
minCon = 1;
binaryMat = 0;
con2(con2<minCon) = 0;


numIn = sum(con2>0,1);
synIn = sum(con2,1);
useIn = (numIn>=minEdge) & (synIn>=minSyn);

numOn = sum((con2>0).*repmat(useIn,[size(con2,1),1]),2);
synOn = sum((con2).*repmat(useIn,[size(con2,1),1]),2);
useOn = (numOn>=minEdge) & (synOn>=minSyn);

nodeIDs = [useList.preList(useOn) setdiff(useList.postList,seedList)];
nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList,seedList)*0+2];


nodeIDs = [useList.preList(useOn) (useList.postList(useIn))];
nodeType = [useList.preList(useOn)*0+1 (useList.postList(useIn))*0+2];



nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);

%% Set color
nodeCol = zeros(nodeNum,3);

%%Set color according to attribute
if 0
    %[colorList cellCol] = getAttributes(obI);
    
    load([MPN 'cb2d.mat'])
    colorList = cb2d.IDs;
    nodeCol = nodeCol + .2;
    
    colPropRaw = cb2d.areaUM;
    
    nodeProp = [];nodePropRef = [];
    for i = 1:length(nodeIDs)
        targ = find(colorList == nodeIDs(i));
        if ~isempty(targ)
            if length(targ)>1
                'too many targets'
                colorList(targ)
            end
            nodeProp= [nodeProp colPropRaw(targ(1))];
            nodePropRef = [nodePropRef i];
            
        end
    end
    
    colProp = nodeProp-min(nodeProp)+1;
    colProp = ceil(colProp/max(colProp) * 100);
    
    colTable = jet(100);
    cellCol = colTable(colProp,:);
    
    nodeCol(nodePropRef,:) = cellCol;
    nodeType(nodePropRef) = nodeType(nodePropRef) + 2;
    
    
    cat(2,nodeIDs', nodeCol);
end

%%Set axons according to seed
if 1
    seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .3 1; 0 .3 1] * .6;
    for s = 1:length(seedList)
        useCol = seedCol(s,:);
        for n = 1:nodeNum
            nodeCol(n,:) = nodeCol(n,:) + useCol *  ...
                double(sum((allEdges(:,1)==nodeIDs(n)) & (allEdges(:,2)) == seedList(s))>0);
        end
    end
end

%%Label specific cells
if 0
    crossoverAxons = [2032	2033	2034	2035]
    gotList = getList_giantBoutons(MPN);
    labelCells = gotList;
    for i = 1:length(labelCells)
        nodeCol(find(nodeIDs==labelCells(i)),:) = nodeCol(find(nodeIDs==labelCells(i)),:) + .5;
    end
end


nodeCol(nodeCol>1) = 1;
nodeCol(nodeCol<0) = 0;

%% Create springDat (need nodeCol, nodeIDs,

springIn.nodeIDs = nodeIDs;
springIn.allEdges = allEdges;
springIn.allWeights = allEdges ;
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = seedList;
springDat = springParameters_X01(springIn);

if binaryMat
    springDat.edges.ew = springDat.edges.ew>0;
end

%% Run springs

for rerun = 1: 1
    
    allResults{rerun} = runSprings(springDat);
    
    %runSprings(springDat,allResults{1})
    set(gcf, 'InvertHardCopy', 'off');
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    
    epsName = sprintf('%sspringRun_%03.0f.eps',springDir,rerun);
    print(gcf, epsName, '-depsc2','-painters')
    
end

%% Save
%{


result = allResults{rerun};
save([springRes 'res_all_07.mat'],'result')


%}
