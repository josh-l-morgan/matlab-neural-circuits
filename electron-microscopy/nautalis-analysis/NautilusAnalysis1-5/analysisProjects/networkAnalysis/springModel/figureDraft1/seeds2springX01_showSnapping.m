



%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903];
    %seedList = [ 108  201 109 ];
    
    tempList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
%     tracedList = obI2cellList_tracedCells(obI);
%     axList = getList_glomID([1,2]);
%     cellList = tempList.postList; 
%     useList = obI2cellList_prePostList2Use(obI,axList,cellList);
    useList = tempList;
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end


springDir = 'D:\LGNs1\Analysis\springDat\figDraft1b\';
springRes = 'D:\LGNs1\Analysis\springDat\results\';
if ~exist(springDir,'dir'), mkdir(springDir), end
if ~exist(springRes,'dir'), mkdir(springRes), end
load([springRes 'noSeedHouse4_edit2.mat'])

load([springRes 'noSeedHouse8_edit3.mat'])
load([springRes 'turkey2_edit2.mat'])
load([springRes 'fourSeeds_edit1.mat'])

%% Filter for convergence

con = useList.con;
%allEdges = con2syn(con>0)

con2 = con;
if 0 %zero seeds
    for i = 1:length(seedList)
        con2(:,find(useList.postList == seedList(i))) = .000000;
    end
end

minEdge = 2;
minSyn = 1;
minCon = 1;
binaryMat = 0;
con2(con2<minCon) = 0;


numIn = sum(con2>0,1);
synIn = sum(con2,1);
useIn = (numIn>=minEdge) & (synIn>=minSyn);

numOn = sum((con2>0).*repmat(useIn,[size(con2,1),1]),2);
synOn = sum((con2).*repmat(useIn,[size(con2,1),1]),2);
useOn = ((numOn>=minEdge) & (synOn>=minSyn));

% 
% numIn = sum(con2>0.* repmat(useOn,[1, size(con2,2)]),1);
% synIn = sum(con2 .* repmat(useOn ,[1, size(con2,2)]),1);
% useIn = (numIn>=minEdge) & (synIn>=minSyn);

% [a idx] = intersect(useList.preList,tracedList.preList);
% useOn = intersect(useOn, idx);
%useOn(useList.preList == 1125) = 0;
% 
% nodeIDs = [useList.preList(useOn) setdiff(useList.postList(useIn),[seedList])];
% nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList(useIn),seedList)*0+2];


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
    seedCol = [1 0 0; 0 1 0;  0 .3 1; .3 0 1; 0.5 0 1];
    for s = 1:length(seedList)
        useCol = seedCol(s,:);
        for n = 1:nodeNum
            nodeCol(n,:) = nodeCol(n,:) + useCol *  ...
                double(sum((allEdges(:,1)==nodeIDs(n)) & (allEdges(:,2)) == seedList(s))>0);
        end
    end
end

if 1  %Color TCR according to seed preference
    prefID = seedPref.cellList;
    usePref = seedPref.sharedSyn;
    usePref(3,:) = sum(usePref(3:end,:),1);
    maxPref = max(usePref,[],1);
    usePref = usePref./repmat(maxPref,[size(usePref,1),1]);
    
    seedCol = [1 0 0; 0 1 0; .3 0 1; 0 .3 1; 0 0 1];
    for i = 1:length(prefID);
        
        targ = find(nodeIDs == prefID(i));
        if ~isempty(targ)
            nodeCol(targ,:) = sum([seedCol(:,1) * usePref(1,i)>0   seedCol(:,2) * usePref(2,i)>0 ...
                seedCol(:,3) * usePref(3,i)>0],1) ;
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
springDat = springParameters_X01_figDraftNoProp_snap(springIn);
%springDat = springParameters_X01_testSprings(springIn);


if binaryMat
    springDat.edges.ew = springDat.edges.ew>0;
end

%% Run springs
movDir  = [springDir 'seedSnap13_edge2\'];
for rerun = 1: 1
    
    allResults{rerun} = runSprings(springDat,[],movDir);
    
   
    
end


if 0
  %%  
     set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])

     %runSprings(springDat,allResults{1})
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'snapState1';
%     imageName = sprintf('%sspringRun_%s.png',springDir,tag);
%     print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
%     
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
    
    
    
    
end



%% Save
%{






result = allResults{rerun};
save([springRes 'erode_edge2.mat'],'result')


%}
