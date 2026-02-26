
syn = tis.syn
targ = 4;
targPre = find(syn.pre==targ);
targPost = find(syn.post==targ);

tis.obI.nameProps.names{syn.obID(syn.synType==1)}

isRibIn = find((syn.post==targ) & (syn.synType == 2))
isNonRibIn = find((syn.post==targ) & (syn.synType == 1))
isOutRGC = find((syn.pre==targ) & (syn.postClass == 1))
isOutAMC = find((syn.pre==targ) & (syn.postClass == 8))
isOutUNK = find((syn.pre==targ) & (syn.postClass == 4))

length(isRibIn)
length(isNonRibIn)
length(isOutRGC)
length(isOutAMC)
length(isOutUNK)
