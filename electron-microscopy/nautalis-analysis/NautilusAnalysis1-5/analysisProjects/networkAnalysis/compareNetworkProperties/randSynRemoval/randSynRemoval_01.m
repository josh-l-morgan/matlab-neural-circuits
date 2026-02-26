


%% Get data
loadData = 1;
if loadData
    %clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [ 108 201 109 907 903];
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end


randSynResDir = 'D:\LGNs1\Analysis\springDat\results\';
if ~exist(randSynResDir,'dir'), mkdir(randSynResDir), end

colormap bluered(256)
%% Filter for convergence

con2 = useList.con;
if 1 %zero seeds
    for i = 1:length(seedList)
        con2(:,find(useList.postList == seedList(i))) = 0;
    end
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

nodeIDs = [useList.preList(useOn) setdiff(useList.postList,seedList)];
nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList,seedList)*0+2];


nodeIDs = [useList.preList(useOn) (useList.postList(useIn))];
nodeType = [useList.preList(useOn)*0+1 (useList.postList(useIn))*0+2];



nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);



%%
realCon = zeros(length(nodeIDs));
for i = 1:length(allEdges)
   preTarg = find(nodeIDs == allEdges(i,1));
   postTarg = find(nodeIDs == allEdges(i,2));
   if ~isempty(preTarg) & ~isempty(postTarg)
      realCon(preTarg,postTarg) =  realCon(preTarg,postTarg) + 1;
   end
end
realCon = realCon;
%%
nn = length(nodeIDs);

transX = ones(nn);
transY = ones(nn);
for y = 1:nn;
    for x = 1:nn
        transY(y,x) = y;
        transX(y,x) = mod(y -x,nn)+1;
    end
end
        

forTrans = sub2ind(size(transX),transX,transY);
%% Remove synapses and find groups
 reps = 48;
allReps = zeros(nn,nn,reps);
parfor rep = 1:reps;
    rep
    con = realCon;
    sumGroup = groupRandRemovalFun(con,forTrans);
    allReps(:,:,rep) = sumGroup;
end


meanReps = mean(allReps,3);
maxMean = ceil(max(meanReps(:)));

randSynRes.allReps = allReps;
randSynRes.nodeIDs = nodeIDs;
randSynRes.con = realCon;
randSynRes.meanReps = meanReps;
randSynRes.useList = useList;

%save('D:\LGNs1\Analysis\groupNodes\randSynRemoval\randSynRes_24min2.mat','randSynRes')
    

%%

for r = 1:100
   
    use1 = rand(size(allReps,3),1)>.5;
    use2 = ~use1;
    subplot(2,1,1)
    image(mean(allReps(:,:,use1),3)*200/maxMean)
    subplot(2,1,2)
    image(mean(allReps(:,:,use2),3)*200/maxMean)
    pause(.1)
end



%%
[shuffleGraph,newOrder] = reshuffleSimilarity(meanReps);
newCon = realCon(newOrder,:);
newCon = newCon(:,newOrder);
image(newCon*100)


%%
randGroups = zeros(maxMean+1,nn);
snapTimes = 0:maxMean;
for s = 1:length(snapTimes)
        [randGroup] = segmentCon(meanReps>snapTimes(s));
        randGroups(s,:)  = randGroup;
end

image(randGroups)

results.cellGroups = randGroups;
results.snapTimes = snapTimes;
results.nodeIDs = nodeIDs;

%save('D:\LGNs1\Analysis\groupNodes\randSynRemoval\result_randSynRes_24min2.mat','results')








