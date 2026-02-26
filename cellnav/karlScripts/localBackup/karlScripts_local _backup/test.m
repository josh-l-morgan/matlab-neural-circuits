function sourceSeg=getSourceSeg(postCid,preCid,VASTxyz,curTis,obI)
targPost=3;
targPre=0;
synPtRaw=[39216, 28272, 719];
synPt=synPtRaw([2 1 3])./[250 250 25];
goodSynIDs=find(curTis.syn.edges(:,1)==targPost&curTis.syn.edges(:,2)==targPre);
goodSynPos=curTis.syn.pos(goodSynIDs,:);
posDiff=goodSynPos-synPt;
absPosDiff=abs(posDiff);
closeSynIDid=find(sum(absPosDiff,2)==min(sum(absPosDiff,2)));
matchingSynID=goodSynIDs(closeSynIDid);
goodObId=curTis.syn.obID(matchingSynID);
joeID=obI.fuse.obSource(goodObId);
joe=obI.fuse.exportDir(joeID);
sourceSeg=joe;