



clf
useList = obI2cellList_seedInput(obI,[108 201])
conPref = seedPreferences([108 201],useList);
cellNum = length(conPref.cellList);

con = useList.con

targ = find(useList.postList == 201);

isPre = con(:,targ)>0;

linkedSyn = con(isPre,:)

sum(linkedSyn(:))-sum(con(:,targ))

sum(sum(linkedSyn,1)>0)