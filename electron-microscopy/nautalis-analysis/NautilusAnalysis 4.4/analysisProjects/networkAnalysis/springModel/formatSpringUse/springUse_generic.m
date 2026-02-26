function[springIn] = springUse_generic()


    load('MPN.mat')
    %MPN = GetMyDir;
    load([MPN 'obI.mat']);
    
    seedList = [ 108  201 907 903];
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);

conTo = makeConTo(obI,seedList);

        itag  = 'giant';

springDir = 'D:\LGNs1\Analysis\springDat\figDraft1b\';
springRes = 'D:\LGNs1\Analysis\springDat\results\';
if ~exist(springDir,'dir'), mkdir(springDir), end
if ~exist(springRes,'dir'), mkdir(springRes), end

load([springRes 'res_all_Phage_edit3.mat'])
load([springRes 'fourSeeds_edit1.mat'])
load([springRes 'fourSeeds_edit1.mat'])


%% Filter for convergence

con = useList.con;
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

if 1 %define clade colors
    load('.\data\clade\cladePick_six2.mat')
    nodeCol = zeros(nodeNum,3)+.2;

    
    useCells = nodeIDs;%[conTo([1:4]).tcrList];
        useCells = [conTo(1).rgcList];

    colorTable = [ .3 .3 1; 1 1 0; 1 0 0; 1 1 1; 1 0 1 ; 0 1 0]
      colorTable = [ .3 .3 1; 1 1 0; 1 0 0; .3 .3 .3; .3 .3 .3 ;  .3 .3 .3]
 
    col = [];
    allMembers = [];
    for i = 1:length(cladePick.members)
        members = cladePick.members{i};
        members = intersect(useCells,members);
        allMembers = [allMembers members];
        col = cat(1,col,repmat(colorTable(i,:),[length(members),1]));
        
    end
    
    for i = 1:length(allMembers);
        targ = find(nodeIDs == allMembers(i));
        nodeCol(targ,:) = col(i,:);
        nodeType(targ) =  nodeType(targ)+2;
    end
        
end



%%Set color according to attribute
if 0
    %[colorList cellCol] = getAttributes(obI);
    nodeCol = nodeCol + .2;

    load([MPN 'cb2d.mat'])
    colorList = cb2d.IDs;
    colPropRaw = cb2d.areaUM;
        
    [colorList colPropRaw] = getList_mtaCount();
    
    [colorList colPropRaw] = getList_giantBoutonsCounts(MPN);
    [colorList colPropRaw] = getList_inSpineCount;

    
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
    
    propRange =[ 0 30]
    showRange = round([propRange(1) : diff(propRange)/10 : propRange(2)]*10)/10
    ticPos = [0:10:100];
    
    colProp = nodeProp-propRange(1);
    colProp = ceil(colProp/((propRange(2)-propRange(1))) * 100);
    colProp(colProp<1) = 1;
    colProp(colProp>100) = 100;
    
    colTable = bluered(100);
    cellCol = colTable(colProp,:);
    
    %cellCol(nodeProp==0,:) = 0;
    
    nodeCol(nodePropRef,:) = cellCol;
    nodeType(nodePropRef) = nodeType(nodePropRef) + 2;
    
   
    
    
    
    cat(2,nodeIDs', nodeCol);
end

%%Set axons according to seed
if 0
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
springIn.allWeights = cat(2,allEdges , ones(size(allEdges,1),1));
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = seedList;