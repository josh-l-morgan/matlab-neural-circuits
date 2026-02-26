


seedCells = [201];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedCells);

conRaw = useList.con;
isSeed = find(useList.postList == seedCells);
notSeed = setdiff(1:size(con,2),isSeed);
con = conRaw(:,notSeed);

axNum = sum(con>0,1)
synNum = sum(con,1);
meanSyn = synNum./axNum;

scatter(axNum,meanSyn)
%kruskalwallis(axNum,meanSyn)