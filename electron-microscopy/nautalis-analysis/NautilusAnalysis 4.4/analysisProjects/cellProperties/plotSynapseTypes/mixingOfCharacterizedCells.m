

preNames = [1028 1032 2033 2034 2035 2003 2007]


seedList = [108 201 109 903 907];
seedNum = length(seedList);
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);
allEdges = obI.nameProps.edges(:,[2,1]);


postList = useList.postList;
edges = butSize.edges;
con = edge2con(edges);
axList = butSize.axList;

conPref =seedPreferences(seedList,useList);
isMixed = sum(conPref.sharedSyn([1 2],:)>0,1)==2;
mixList = conPref.cellList(isMixed);
%mixList = postList(randperm(length(postList),26));
unMixList = setdiff(postList,mixList);


preMix = zeros(1,length(preNames));
preUnMix = preMix;
for p = 1:length(preNames)
    
    isPre = allEdges(:,1) == preNames(p);
    isPost = allEdges(isPre>0,2);
    
    for i = 1:length(isPost);
        
        if sum(unMixList == isPost(i))>0
            preUnMix(p) = preUnMix(p) + 1;
        elseif sum(mixList == isPost(i))>0;
            preMix(p) = preMix(p) + 1;
        end
    end
end


mixNums = [preNames' preMix' preMix'+preUnMix']







