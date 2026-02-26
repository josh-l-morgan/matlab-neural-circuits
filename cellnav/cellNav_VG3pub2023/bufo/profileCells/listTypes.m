

syn = tis.syn;

post = COI.postToVGC

allSubTypes = [];

for i = 1:length(post)

    targ = find(tis.cells.cids == post(i));
    cType = tis.cells.type.typeID(targ);
    sType = tis.cells.type.subTypeID(targ);
    if cType == 1
        allSubTypes = [allSubTypes sType];
    end
end

knownTypes = allSubTypes(allSubTypes>0);

uKnownTypes = unique(knownTypes);

tis.cells.type.subTypeNames{1}(uKnownTypes);




useRNames = { '2an' '4i'  '4on'  '4ow' '3i' '37'  '5si' '5to' ...
    '25' '28'  '5so' '5ti' '63' '85'  '7i' '6sw' '6sn' '6t' '8w' };





