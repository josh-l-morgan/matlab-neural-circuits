analDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\';
mergeDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\';

load([analDir 'tis.mat']);
load([mergeDir 'dsObj.mat']);

curTis=tis;
vgcCidList=[2 3 4 5 10 11 13 14 20];
figs=1;
sampSize=20;

%get the type information for all the cids for all the synapses
allPreTypes=cid2type(tis.syn.edges(:,2),curTis);
allPostTypes=cid2type(tis.syn.edges(:,1),curTis);

%get the IDs of all the bpc2vgc synapses
bpcPre=find(allPreTypes{1}==7)';
vgcPost=find(ismember(tis.syn.edges(:,1),vgcCidList));
bpc2vgcSynIDs=intersect(bpcPre,vgcPost);

if figs
    figure(); 
    histogram(tis.syn.edges(bpc2vgcSynIDs,1));
    title('target vgc of bpc inputs');
end

vgcCidListShort=[2 3 5 13 14];

posList=[];
sampIDs=[];
for i=1:length(vgcCidListShort)
    curV=vgcCidListShort(i);
    allIDs=bpc2vgcSynIDs(find(tis.syn.edges(bpc2vgcSynIDs,1)==curV));
    curIDs=randsample(allIDs,sampSize);
    sampIDs=[sampIDs;curIDs];
    curPos=tis.syn.pos(curIDs,:);
    posList=[posList;curPos];
end

posListVast=posList.*[250 250 25];
posListVast=horzcat(posListVast(:,2),posListVast(:,1),posListVast(:,3));


figure();
tesDS=tes./[250 250 25];
scatter3(tesDS(:,1),tesDS(:,2),tesDS(:,3));
hold on
