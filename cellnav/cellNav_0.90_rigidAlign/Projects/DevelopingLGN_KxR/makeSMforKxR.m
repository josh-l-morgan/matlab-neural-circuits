

global tis


%%makeCOI ?

tcrs = tis.cids(find(tis.cells.type.typeID==2));
lins = tis.cids(find(tis.cells.type.typeID==3));

tcsWithDendritesCharacterized = [3004 3107 3210 3211 3213 3018 3097 3032 3098 3200 3202 3203 3204 3206 3209 3219];

cidsForSMs = setdiff(tcsWithDendritesCharacterized,[3210 3204])

fig_SM = figure;
for i = 1:length(cidsForSMs)
    sm = makeSM(cidsForSMs(i));
end













