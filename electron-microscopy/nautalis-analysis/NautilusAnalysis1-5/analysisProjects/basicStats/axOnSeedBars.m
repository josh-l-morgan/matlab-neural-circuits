
load(['MPN.mat'])
load([MPN 'obI.mat'])

targCell = 108;
useList = obI2cellList_seedInput_RGC_TCR(obI,targCell);
useList = useListBoutSize(useList);

size(useList.con);

allEdges = obI.nameProps.edges(:,1:2);






%%
clf


con = useList.con;
targ = find(useList.postList == targCell);

isPre = con(:,targ)>0;
rawPre = con(:,targ);
linkedSyn = con(isPre,:);
linkedIDs = useList.preList(isPre);

onSeed = linkedSyn(:,targ);

[onSeed idx] = sort(onSeed,'descend');
sortPre = linkedIDs(idx);
linkedSyn = linkedSyn(idx,:);
onAll = sum(linkedSyn,2);

axNum = length(onSeed);
maxOn = max(onAll);

bar(onAll,'FaceColor',[.8 .8 .8])
hold on
bar(onSeed,'FaceColor',[0 0 0])
hold off

ylim([0 70]);
return

%%
goingUp = sort(onSeed,'ascend');
sumOn = cumsum(goingUp);
plot(sumOn)

hist(onSeed,[1:1:max(onSeed)]);
;
tabPreSyn = [sortPre' onSeed];

%%
useAx = onAll> 10;

hist(onSeed(useAx),[1:1:max(onSeed)]);

bar(sort(onSeed(useAx),'descend'));

median(onSeed(useAx));


%%
sum(onSeed);

sum(onSeed(onSeed>=10));
sum(onSeed(onSeed<=10));

sum(onSeed(onSeed<=5));

%%
allEdges
preList = useList.preList;
synCount = preList * 0;
for i = 1:length(preList)
    synCount(i) = sum((allEdges(:,1) == targCell) & (allEdges(:,2) == preList(i)));

end

[sortSyn idx] = sort(synCount,'descend')
sortPre = preList(idx);
preSyn = [sortPre' sortSyn'];












