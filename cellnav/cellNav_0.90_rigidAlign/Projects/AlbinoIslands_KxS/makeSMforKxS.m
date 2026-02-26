

global tis


%%makeCOI ?

tcrs = tis.cids(find(tis.cells.type.typeID==2));
lins = tis.cids(find(tis.cells.type.typeID==3));

cidsForSMs = [tcrs lins];
cidsForSMs = [1005]

fig_SM = figure;
for i = 1:length(cidsForSMs)
    sm = makeSM(cidsForSMs(i));
end













