

if 0

tis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat');
curTis=tis.tis;
typeDat=cid2type(curTis.cells.cids,curTis);
rgcCidList=curTis.cells.cids(find(typeDat{1}==1));
vgcCidList=[2 3 4 5 10 11 13 14 20];

output=length(find(ismember(curTis.syn.edges(:,2),vgcCidList)));
out2RGC=length(find(ismember(curTis.syn.edges(:,1),rgcCidList) & ismember(curTis.syn.edges(:,2),vgcCidList)));

end

rgcSubTypes=typeDat{3}(typeDat{1}==1);
rgcSubTypes(rgcSubTypes==0)=56;
figure(); histogram(rgcSubTypes,[0.5:1:56.5]);
rgcSubHistDat=histcounts(rgcSubTypes,[0.5:1:56.5]);
xticks([0.5:1:56.5]);
xticklabels(curTis.cells.type.subTypeNames{1});
ylim([0 8]);

significantList=find(rgcSubHistDat>0);
significantNames=curTis.cells.type.subTypeNames{1}(significantList);
significantDat=rgcSubHistDat(rgcSubHistDat>0);

knownRGCcids=
