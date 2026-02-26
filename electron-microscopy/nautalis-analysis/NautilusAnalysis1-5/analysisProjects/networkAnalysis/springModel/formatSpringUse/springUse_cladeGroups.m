function[springIn] = springUse_clades()


load('MPN.mat')
load([MPN 'obI.mat']);

seedList = [ 108  201 907 903];

useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
%seedPref = seedPreferences(seedList,useList);
allEdges = obI.nameProps.edges(:,[2 1]);
conTo = makeConTo(obI,seedList);

%%
   colorTable = [ .3 .3 1; 1 1 0; 1 0 0; 1 1 1; 1 0 1 ; 0 1 0]
      colorTable = [ .3 .3 1; 1 1 0; 1 0 0; .3 .3 .3; .3 .3 .3 ;  .3 .3 .3]
 
con = useList.con;

 load('.\data\clade\cladePick_six2.mat')
    
    useRGC =  [conTo(:).rgcList];
    useTCR = [conTo(:).tcrList];

    cladeCon = zeros(length(cladePick.members))
    preCon = zeros(length(cladePick.members),size(con,2));
    for i = 1:length(cladePick.members)
        members = cladePick.members{i};
        members = intersect(useRGC,members);
        
        for m = 1:length(members)
            targ = find(useList.preList == members(m));
            if ~isempty(targ)
               preCon(i,:) = preCon(i,:) + con(targ,:); 
            end
            
        end
        
    end
    
      for i = 1:length(cladePick.members)
        members = cladePick.members{i};
        members = intersect(useTCR,members);
        
        for m = 1:length(members)
            targ = find(useList.postList == members(m));
            if ~isempty(targ)
               cladeCon(:,i) = cladeCon(:,i) + preCon(:,targ); 
            end
            
        end
        
    end
    


useList = con2use(cladeCon);



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

    colorTable = [ .3 .3 1; 1 1 0; 1 0 0; 1 1 1; 1 0 1 ; 0 1 0]

nodeCol = cat(1,colorTable,colorTable);
nodeType = nodeType + 2;





%% Create springDat (need nodeCol, nodeIDs,

springIn.nodeIDs = nodeIDs;
springIn.allEdges = useList.allEdges;
springIn.allWeights = useList.allWeights;
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = [];







