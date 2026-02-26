



%% Get data
loadData = 1;
if loadData
    %clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
    tempList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    useList = tempList;
    %useList = obI2nodes_rtl_plus1(obI,plusOne);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end
binaryMat = 0;

springDir = 'D:\LGNs1\Analysis\springDat\figDraft1b\';
springRes = 'D:\LGNs1\Analysis\springDat\results\';
if ~exist(springDir,'dir'), mkdir(springDir), end
if ~exist(springRes,'dir'), mkdir(springRes), end
load([springRes 'noSeedHouse4_edit2.mat'])

load([springRes 'noSeedHouse8_edit3.mat'])
load([springRes 'turkey2_edit2.mat'])
%load([springRes 'withFourSeeds_edit1.mat'])
load([springRes 'turkey_plusOne125.mat'])


%% Filter nodes
nodes = useList.nodes;
seedFilt = 108;
seedPos = find(nodes==seedFilt);
con = useList.con;
symCon = con + con';
[cellGroup] = segmentCon(symCon);
seedGroup = cellGroup(seedPos);
nodes = nodes(cellGroup == seedGroup);

nodeType = nodes*0;
for i = 1:length(nodes)
    nodeType(i) = useList.nodeType(useList.nodes==nodes(i));
end
nodeIDs = nodes;




nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);


allEdges(:,3) = 1;
allEdges(allEdges(:,1) == plusOne,3) = 1;
allEdges(allEdges(:,2) == plusOne,3) = 1;

%% Set color
nodeCol = zeros(nodeNum,3);

if 1
   
    [nodeCol useCol] = getList_seedColor(seedList,nodeIDs) ;
    nodeCol = nodeCol * .8;
    nodeCol(sum(nodeCol,2) == 0,:) = .3;
    
    seedCol = [1 0 0; 0 1 0; .5 0 1 ; 0 .5 1];
     for s = 1:length(seedList)
        targ = find(nodeIDs == seedList(s));
        nodeType(targ) = nodeType(targ)+2;
        nodeCol(targ,:) = seedCol(s,:);
     end
     
     showCol = permute(useCol.col,[1 3 2])
    
    image(uint8(showCol*256))
    useCol.comb{10}
    
    targ = find(nodeIDs == plusOne);
    nodeCol(targ,:) = [0 0 0]
    
end


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
if 0
    seedCol = [1 0 0; 0 1 0;  0 0 1; 0 0 1; 0 0 1] * 10;
    for s = 1:length(seedList)
        useCol = seedCol(s,:);
        for n = 1:nodeNum
            nodeCol(n,:) = nodeCol(n,:) + useCol *  ...
                double(sum((allEdges(:,1)==nodeIDs(n)) & (allEdges(:,2)) == seedList(s))>0);
        end
    end
end

if 0  %Color TCR according to seed preference
    prefID = seedPref.cellList;
    usePref = seedPref.sharedSyn;
    usePref(3,:) = sum(usePref(3:end,:),1);
    maxPref = max(usePref,[],1);
    usePref = usePref./repmat(maxPref,[size(usePref,1),1]);
    
    seedCol = [1 .2 .2; .2 1 .2;  .2 .2 1; .2 .2 1; .2 .2 1] * 1;
    for i = 1:length(prefID);
        
        targ = find(nodeIDs == prefID(i));
        if ~isempty(targ)
            nodeCol(targ,:) = sum([seedCol(:,1) * usePref(1,i)   seedCol(:,2) * usePref(2,i)...
                seedCol(:,3) * usePref(3,i)],1) ;
        end
    end
    
    nodeCol(nodeCol>.7) = 0.7;
end




%%Set axons according to seed indexed combinations
if 0
    seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .3 1; 0 .3 1] * .9;
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
%springDat = springParameters_X01_figDraftNoProp_WithSeeds(springIn);
%springDat = springParameters_X01_testSprings(springIn);
springDat = springParameters_LIN01(springIn);


if binaryMat
    springDat.edges.ew = springDat.edges.ew>0;
end

%% Run springs
movDir  = [springDir 'seedMovie_WithSeeds3\'];

for rerun = 1: 1
    
    allResults{rerun} = runSprings(springDat);
    
   
    
end


if 0
  %%  
     set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])

     %runSprings(springDat,allResults{1})
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'net125_white';
%     imageName = sprintf('%sspringRun_%s.png',springDir,tag);
%     print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
%     
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
    
    
    
    
end



%% Save
%{






result = allResults{rerun};
save([springRes 'turkey_plusOne125.mat'],'result')


%}






