global glob tis



tis.syn.post == 



isVGC = find((tis.cells.type.typeID == 8) & (tis.cells.type.subTypeID == 1));
vgcCids = tis.cids(isVGC);

allCids = tis.cids;

notPre = allCids;

for i = 1 : length(isVGC)
    
    vPost = tis.syn.post == isVGC(i);
    isPre = tis.syn.pre(vPost);
    notPre = setdiff(notPre,isPre);

end

notPre = notPre'






