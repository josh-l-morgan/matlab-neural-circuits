function[newEdges] = ranomizeNonSeed(nodeIDs,seedList,allEdges);




for e = 1:size(allEdges,1)
    
    
    isNode(e) = (sum(nodeIDs == allEdges(e,1))>0) + (sum(nodeIDs == allEdges(e,2))>0);
    isSeed(e) = (sum(seedList == allEdges(e,1))>0) + (sum(seedList == allEdges(e,2))>0);
    
end

seedEdges = allEdges(isSeed>0,:);

nodeEdges = allEdges( (isNode == 2) & (isSeed == 0),:);


maxE = max(nodeEdges);

edgeInd = sub2ind(maxE,nodeEdges(:,1),nodeEdges(:,2));
uEdge = unique(edgeInd);
histEdge = hist(edgeInd,uEdge);

[uPre uPost] = ind2sub(maxE,uEdge);

newOrder = randperm(length(uPre));

newPre = uPre(newOrder);
histPre = round(histEdge(newOrder) * .501);
histPost = round(histEdge * .499);
newHist = histPre+histPost;

randEdges = [];
for h = 1:length(newHist)

    randEdges = cat(1,randEdges, repmat([newPre(h) uPost(h)],[newHist(h),1]));
    
end
% randEdges = [nodeEdges(randperm(size(nodeEdges,1)),1) nodeEdges(:,2)];

newEdges = cat(1,seedEdges,randEdges);









