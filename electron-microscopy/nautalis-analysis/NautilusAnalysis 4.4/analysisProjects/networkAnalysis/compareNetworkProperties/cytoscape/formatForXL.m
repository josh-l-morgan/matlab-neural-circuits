



%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 109 907 903];
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    seedPref = seedPreferences(seedList,useList);
    %allEdges = obI.nameProps.edges(:,[2 1]);
    
    
    [y x] = find(useList.con>0);
    v = pi * useList.con(useList.con>0)/8;
    allEdges = [useList.preList(y)' useList.postList(x)' v];
end

%% Filter for convergence

con = useList.con;
%allEdges = con2syn(con>0)

con2 = con;
if 1 %zero seeds
    for i = 1:length(seedList)
        con2(:,find(useList.postList == seedList(i))) = .000000;
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
useOn = find((numOn>=minEdge) & (synOn>=minSyn));
% [a idx] = intersect(useList.preList,tracedList.preList);
% useOn = intersect(useOn, idx);


nodeIDs = [useList.preList(useOn) setdiff(useList.postList,seedList)];
nodeType = [useList.preList(useOn)*0+1 setdiff(useList.postList,seedList)*0+2];


nodeIDs = [useList.preList(useOn) (useList.postList(useIn))];
nodeType = [useList.preList(useOn)*0+1 (useList.postList(useIn))*0+2];

nodeNum = length(nodeIDs);


%%  Create weighted edge list

[y x] = find(con2);
v = con2(con2>0);
preEdge = useList.preList(y);
postEdge = useList.postList(x);













