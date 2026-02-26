%% check the syn to appo dists
allRGCcids=type2cid({'rgc'},{'all'},curTis);
vgcCidList=[2 3 4 5 10 11 13 14 20];
allRGCcids=allRGCcids{1}';
rgcIDs=find(ismember(curTis.syn.edges(:,1),allRGCcids));
vgcIDs=find(ismember(curTis.syn.edges(:,2),vgcCidList));

synIDs=intersect(rgcIDs,vgcIDs);
synLocs=curTis.syn.pos(synIDs,:);
%synLocs=synLocs.*[250 250 25];
[gcl,inl,synDepths]=getIPLdepth(synLocs(:,3),synLocs(:,2),synLocs(:,1),[],[]);
%distMat=pdist2(appoLocs,synLocs);
%figure(); histogram(distMat(:));
%appoSynMinDist=min(distMat,[],2);
%figure(); histogram(appoSynMinDist,115);



%% get all outputs for a vgc
vgcCidList=2;
for i=1:length(vgcCidList)
    curVGCcid=vgcCidList(i);
    curV2RIDS=intersect(synIDs,find(curTis.syn.edges(:,2)==curVGCcid));
    curVtargCids=curTis.syn.edges(curV2RIDS,1);
    targetCounts=arrayfun(@(x)length(find(curVtargCids==x)),unique(curVtargCids),'Uniform',false);
    targetCountMat=cell2mat(targetCounts);
    results=horzcat(unique(curVtargCids),targetCountMat);
    resultsSorted=sortrows(results,2,'descend');
    
end

%% as a test, bpc appositions close to the cid2-cid2002 synapses.
synList=find(curTis.syn.edges(:,1)==2002 & curTis.syn.edges(:,2)==2);
synLocs=curTis.syn.pos(synList,:);
outSynList=find(curTis.syn.edges(:,2)==2);
randSynList=randsample(outSynList,length(synList));
randSynLocs=curTis.syn.pos(randSynList,:);
% manually created the apposition lists from cellNav. - working on the
% getting the function to work. Uses the globs from the GUI, so I'll have
% to hack those to use the function independently.

distSynApp=pdist2(bpc_cid2_appo,synLocs);
distRandApp=pdist2(bpc_cid2_appo,randSynLocs);

histcounts(distRandApp(:));


