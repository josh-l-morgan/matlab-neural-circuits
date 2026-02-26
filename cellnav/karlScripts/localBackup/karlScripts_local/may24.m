%partners=[];
partTypes=cid2type(partners,curTis);

rgcPartTypes=partTypes{3}(partTypes{1}==1);

figure()
histogram(rgcPartTypes,[-0.5:1:length(curTis.cells.type.subTypeNames{1})-0.5]);
xticks([0:length(curTis.cells.type.subTypeNames{1})]);
labs=curTis.cells.type.subTypeNames{1};
labs={'none' labs{:}};
xticklabels(labs);

allRGCtargs=cid2type(curTis.syn.edges(find(ismember(curTis.syn.edges(:,2),vgcCidList)),1),curTis);
allRGCtargs=allRGCtargs{3}(allRGCtargs{1}==1);

figure();
histogram(allRGCtargs,[-0.5:1:length(curTis.cells.type.subTypeNames{1})-0.5]);