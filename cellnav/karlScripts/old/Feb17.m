%% check the syn to appo dists
targetRGCcidCell=type2cid({'rgc'},{'4ow'},curTis);
vgcCidList=[2 3 4 5 10 11 13 14 20];
targetRGCcids=targetRGCcidCell{1};
rgcIDs=find(ismember(curTis.syn.edges(:,1),targetRGCcids));
vgcIDs=find(ismember(curTis.syn.edges(:,2),vgcCidList));
synIDs=intersect(rgcIDs,vgcIDs);
synLocs=curTis.syn.pos(synIDs,[2 1 3]);
synLocs=synLocs.*[250 250 25];
distMat=pdist2(appoLocs,synLocs);
%figure(); histogram(distMat(:));
appoSynMinDist=min(distMat,[],2);
figure(); histogram(appoSynMinDist,115);