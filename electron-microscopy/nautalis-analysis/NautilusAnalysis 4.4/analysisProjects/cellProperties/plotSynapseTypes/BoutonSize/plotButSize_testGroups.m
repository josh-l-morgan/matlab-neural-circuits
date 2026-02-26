clear all
load('MPN.mat')
load([MPN 'butSize2.mat'])
load([MPN 'obI.mat'])

%%
seedList = [108 201 109 903 907];
seedNum = length(seedList);
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);

postList = useList.postList;
edges = butSize.edges;
con = edge2con(edges);
axList = butSize.axList;

conPref =seedPreferences(seedList,useList);
isMixed = sum(conPref.sharedSyn([1 2],:)>0,1)==2;
mixList = conPref.cellList(isMixed);
%mixList = postList(randperm(length(postList),26));
unMixList = setdiff(postList,mixList);


%%
allSynEdge = [];
allbutDiam = [];
for a = 1:length(butSize.synDat)
    synEdge = butSize.synDat(a).edges;
    butVols = butSize.butVols{a};
    butDiam = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2
    allSynEdge = [allSynEdge; synEdge];
    allbutDiam = [allbutDiam butDiam];
end
    
filtList = obI2cellList_seedInput_RGC_TCR(obI,[108 201]);

testPre = intersect(unique(allSynEdge(:,2)),filtList.preList);
testPost = intersect(unique(allSynEdge(:,1)),filtList.postList);

preMean = testPre * 0;
preCount = preMean;
for i = 1:length(testPre)
    isTarg = find(allSynEdge(:,2) == testPre(i));
    butDiam = allbutDiam(isTarg);
    
preMean(i) = mean(butDiam);
preCount(i) = length(butDiam);
end


postMean = testPost * 0;
postCount = postMean;
for i = 1:length(testPost)
    isTarg = find(allSynEdge(:,1) == testPost(i));
    butDiam = allbutDiam(isTarg);
    
postMean(i) = mean(butDiam);
postCount(i) = length(butDiam);
end

minSyn = 10;
meanBinWidth = [0:.2:3];
subplot(2,1,1) 
histPreMean = hist(preMean(preCount>=minSyn),meanBinWidth);
bar(meanBinWidth,histPreMean)
subplot(2,1,2)
histPostMean = hist(postMean(postCount>=minSyn),meanBinWidth);
bar(meanBinWidth,histPostMean)


[sortMean idx] = sort(postMean);
testPost = testPost(idx);
postCount = postCount(idx);
testPost = testPost(postCount>=minSyn);

[sortPre idx] = sort(preMean);
testPre = testPre(idx);
preCount = preCount(idx);
testPre = testPre(preCount>=minSyn);

synMat = zeros(length(testPre),length(testPost));
sumMat = synMat;
for i = 1:length(testPre)
    for p = 1:length(testPost)
        isTarg = find((allSynEdge(:,2) == testPre(i)) & ...
            (allSynEdge(:,1) == testPost(p)));
        butDiam = allbutDiam(isTarg);
        sumMat(i,p) = sum(butDiam);
        synMat(i,p) = length(butDiam);
    end
end
sizeMat(synMat>0) = sumMat(synMat>0)./synMat(synMat>0);


meanX = sum(sumMat,1)./sum(synMat,1);
meanY = sum(sumMat,2)./sum(synMat,2);

xgroup = kmeans(meanX,2);
ygroup = kmeans(meanY,2);
for y = 1:2
    for x = 1:2
        mask = repmat(ygroup == y,[1 length(xgroup)]) &...
            repmat(xgroup' == x,[length(ygroup) 1]);
        quads(y,x) = sum(sum(sumMat.*mask)) / sum(sum(synMat.*mask));
        quadCount(y,x) = sum(sum(synMat.*mask));
    end
end

quads
quadCount


cmap = bluered(99);
cmap = cat(1,[0 0 0],cmap);
colormap(cmap)
subplot(1,1,1)
image(sizeMat*100/max(sizeMat(:)))






