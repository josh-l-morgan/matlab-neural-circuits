%get the locations of the syns from 1026 to vg3s
t=highTis;
allVG=[2 3 4 5 10 11 13 14 20];
vgish=[6735 3723 6724];
new2vg=find(t.syn.edges(:,2)==1026 & ...
    ismember(t.syn.edges(:,1),allVG));
new2vgish=find(t.syn.edges(:,2)==1026 & ...
    ismember(t.syn.edges(:,1),vgish));
allPreInfo=cid2type(t.syn.edges(:,2),t);
amcIn=find((allPreInfo{1}'==0 | allPreInfo{1}'==8) & ...
    ismember(t.syn.edges(:,1),allVG));
bpcIn=find(allPreInfo{1}'==7 & ...
    ismember(t.syn.edges(:,1),allVG));
loxB=t.syn.pos(amcIn,:);
edgesB=t.syn.edges(amcIn,:);
[a,b,depsB]=getIPLdepth(loxB(:,3),loxB(:,1),loxB(:,2),[],[]);

loxA=t.syn.pos(new2vg,:);
[a,b,depsA]=getIPLdepth(loxA(:,3),loxA(:,1),loxA(:,2),[],[]);

loxC=t.syn.pos(bpcIn,:);
[a,b,depsC]=getIPLdepth(loxC(:,3),loxC(:,1),loxC(:,2),[],[]);


binedges=[0:0.05:1];
histA=histcounts(depsA,binedges);
histB=histcounts(depsB,binedges);
histC=histcounts(depsC,binedges);


figure();
plot(histA*10);
hold on
plot(histB);
%plot(histC);

legend({"1026 to VG3","AMC to VG3"})



%finding out if 3008 is another one of the important ones
out1026=find(highTis.syn.edges(:,2)==1026);
lowDown=find(highTis.syn.pos(:,2)>150);
targetSynIds=intersect(out1026,lowDown);
targetEdges=highTis.syn.edges(targetSynIds,:);



%% SEARCHING
%search for more of these
%make a spreadsheet
SSheet={};
checkInds=find(depsB>0.40&depsB<0.53);
for chind=1:length(checkInds)
    chind
    curInd=checkInds(chind);
    curLoc=loxB(curInd,:);
    copDat=uint16(curLoc([2 1 3]).*[250 250 25])
    SSheet{chind,1}=copDat;
    SSheet{chind,2}=edgesB(curInd,:);
    %clipboard('copy',copDat);
    %pause();
end


