global tis





typeIDs = tis.cells.type.typeID;
tis.cells.type.typeNames
subTypeIDs = tis.cells.type.subTypeID;
tis.cells.type.subTypeNames

vgcList = tis.cids((typeIDs==8) & (subTypeIDs == 1));

vgcSyn = [];
for i = 1:length(vgcList)
    vgcSyn = cat(1,vgcSyn, find(tis.syn.post == vgcList(i)));
end
vgcSyn = unique(vgcSyn);

pre2vgc = tis.syn.pre(vgcSyn)
rib2vgc = vgcSyn(find(tis.syn.synType(vgcSyn) == 2));
ribCids = tis.syn.pre(rib2vgc);

noIDrib = sum(ribCids==0)


num7 = 0; numSub = 0;
for i = 1:length(ribCids)
   targ = find(tis.cids == ribCids(i));
    typeID = typeIDs(targ);
    subTypeID = subTypeIDs(targ);
    hasSub(i) = 0;
    if typeID == 7
        num7 = num7+1;
        if subTypeID >0
            numSub = numSub+1;
            hasSub(i) = 1;
        end
    end
    
end

numberOfRibs = length(ribCids)

ribNoSub = rib2vgc(hasSub==0)
length(ribNoSub)



tis.syn.post(ribNoSub);
badSynPos = tis.syn.synPosRaw(ribNoSub,[2 1 3])

notBip = find(tis.syn.preClass(ribNoSub)==0)
notBipPos = tis.syn.synPosRaw(ribNoSub(notBip),[2 1 3])








