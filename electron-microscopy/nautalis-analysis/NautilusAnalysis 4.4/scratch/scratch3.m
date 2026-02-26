
clf
seedList = [108 201 109 903 907]
%seedList = 903
useList = obI2cellList_seedInput(obI,seedList);

con = useList.con;
image(con*10)

%%
gotList = getList_giantBoutons(MPN);

conTo = makeConTo(obI,seedList);
mix1 = intersect(conTo(1).tcrList, conTo(2).tcrList);
mix2 = intersect(mix1, conTo(3).tcrList);
mix3 = intersect(mix2,gotList)