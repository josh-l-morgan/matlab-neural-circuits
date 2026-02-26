
%%Takes ax variable from displayInSpineTypes2.m

%% Choose counts

minSyn = 1;
seedList = [ 201];
useList108 = obI2cellList_seedInput_RGC_TCR(obI,seedList);

preNames = intersect(preID,useList108.preList);
postNames = intersect(postID,useList108.postList);




listIDs = []; listProp = [];
c = 0;
for a = 1:length(ax.axID)
    inList = ax.inList{a};
    if sum(preNames == ax.axID(a)) & (length(inList)>minSyn)
       
        c = c+1;
        listIDs = cat(2,listIDs,repmat(c,[1 length(inList)]));
        listProp = cat(2,listProp,inList);
    end
end
axNum = c;
synNum = length(listProp);

%%
realMeans = zeros(1,axNum);
for i = 1:axNum
    newSyn = listProp(listIDs == i);
    realMeans(i) = mean(newSyn);
    sort(newSyn)
end
bar(sort(realMeans))

realTest = cvrmse(realMeans);


reps = 1000;
sortTests = zeros(reps,axNum);
for r = 1:reps
    
    newIDs = listIDs(randperm(synNum));
    meanSyn = zeros(1,axNum);
    for i = 1:axNum
        newSyn = listProp(newIDs == i);
        meanSyn(i) = mean(newSyn);
    end
    sortTests(r,:) = sort(meanSyn);
    randTest(r) = cvrmse(meanSyn);
end

subplot(2,1,1)
bar([sort(realMeans)' mean(sortTests,1)'],1)


subplot(2,1,2)
hist(randTest,[0:.1:2])
hold on
scatter(realTest,1,'r')
hold off













