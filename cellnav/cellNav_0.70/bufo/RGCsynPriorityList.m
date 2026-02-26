 global tis
 
vgcs = tis.cells.cids(((tis.cells.type.typeID == 8 )...
& (tis.cells.type.subTypeID == 1)));

v2rSyn = [];
for i = 1:length(vgcs)
   v2rSyn = cat(1,v2rSyn, find(( tis.syn.pre == vgcs(i) ) & ( tis.syn.postClass == 1)));
end

rgcSynIDs = tis.syn.post(v2rSyn);
rgcs = unique(rgcSynIDs);

subNames = tis.cells.type.subTypeNames{1};
clear rgcNames
for i = 1:length(rgcs)
    targ = find(tis.cells.cids==rgcs(i));
    subID = tis.cells.type.subTypeID(targ);
    if subID
    rgcNames{i,1} = subNames{subID};
    else
        rgcNames{i} = '      ';
    end
end

rgcSynCount = histc(rgcSynIDs,rgcs);

[rgcSynCount idx] = sort(rgcSynCount,'descend');
rgcs = rgcs(idx);
rgcNames = rgcNames(idx);





