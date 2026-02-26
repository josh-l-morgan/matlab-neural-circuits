function outputDat=getLocalHood(sm, synInd, hoodSize)
curSM=sm;
%synInd=find(curSM.syn.edges(:,3)==synID);
curSynDists=curSM.syn2Skel.syn2SynDist(:,synInd);
curSynDists(synInd)=inf;
[distSrtd pathInds]=sort(curSynDists,'ascend');
closeSynIDs=find(curSynDists<hoodSize);
outputDat=zeros(find(distSrtd>hoodSize,1),4);
for j=1:find(distSrtd>hoodSize,1)
    curNearID=pathInds(j);
    curNearDist=curSynDists(curNearID);
    curEdges=curSM.syn.edges(curNearID,:);
%     curPreType=allPreTypes{1}(curNearID);
%     curPreSubtype=allPreTypes{3}(curNearID);
%     curPreType=allPreTypes{1}(curNearID);
%     curPreSubtype=allPreTypes{3}(curNearID);
    isInput=curEdges(1)==curSM.cid;
    isOutput=curEdges(2)==curSM.cid;
    if isInput
    end
    outputDat(j,:)=[curNearID curNearDist curEdges(1) curEdges(2)];
    
    
end






end