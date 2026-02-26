function[bridgeEdge] = bridgeGaps(skel,seedPath)
%%Adds bridge edges to skeleton node list. 

%%
nodes = skel.node2surf;
nodeSeg = seedPath.segID(skel.node2surf);

seeds = unique(nodeSeg);
numSeeds = hist(nodeSeg,seeds);
firstSeed = seeds(numSeeds == max(numSeeds));

fixedNodes = nodeSeg == firstSeed;
mainNodes = find(fixedNodes);
lostNodes = find(~fixedNodes);

bridgeEdge = [];
while 1
    mainSubs = skel.node2subs(mainNodes,:);
    lostSubs = skel.node2subs(lostNodes,:);
    distMat = zeros(size(mainSubs,1),size(lostSubs,1));
    for m = 1:size(mainSubs,1)
        distMat(m,:) = sqrt((lostSubs(:,1)-mainSubs(m,1)).^2 + ...
            (lostSubs(:,2)-mainSubs(m,2)).^2 + (lostSubs(:,3)-mainSubs(m,3)).^2) ;
    end
    image(distMat/10),pause(.1)
    [m l] = find(distMat == min(distMat(:)),1);
    bridgeEdge = cat(1,bridgeEdge,[mainNodes(m) lostNodes(l)]);
    newSeg = nodeSeg(lostNodes(l));
    fixedNodes(nodeSeg == newSeg) = 1;
    
    if ~sum(fixedNodes == 0), break, end
    
    mainNodes = find(fixedNodes);
    lostNodes = find(~fixedNodes);
    
end

bridgeEdge
