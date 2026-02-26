fvLib='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
vgcSkelList=[2 3 4 5 13 14];
skelFVs=cell(length(vgcSkelList),1);
for i=1:length(vgcSkelList)
    curSkel=load([fvLib num2str(vgcSkelList(i)) '.mat']);
    skelFVs{i}=curSkel;
end

synPos=curTis.syn.pos;
cid2synIDs=find(curTis.syn.edges(:,1)==2|curTis.syn.edges(:,2)==2);

f1=figure();
hold on
synScat=scatter3(synPos(cid2synIDs,3),synPos(cid2synIDs,1),synPos(cid2synIDs,2),1,'.');

vgcPatches=cell(length(skelFVs),1);
for i=1:length(skelFVs)
    curFV=skelFVs{i};
    vgcPatches{i}=patch(curFV.fv);
    
end

peaks;explore3d(double(gca))