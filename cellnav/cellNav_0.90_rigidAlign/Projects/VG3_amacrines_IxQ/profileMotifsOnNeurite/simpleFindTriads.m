


vgcCids = [2 3 4 5 10 11 13 14 20 21 25 6720 6735 6738];
pre = tis.syn.pre;
post = tis.syn.post;
isBip = tis.syn.preClass == 7;
targCid = 3119;
    


pre2targ = pre((post== targCid) & isBip);

pre2vgc = [];
for v = 1:length(vgcCids)
    pre2vgc = [pre2vgc; pre((post== vgcCids(v)) & isBip)];
end
pre2vgc = unique(pre2vgc);

c = 0
for i = 1:length(pre2targ)
    if sum(pre2vgc == pre2targ(i)); 
    c = c+1;
    end
end
c

sharedCids = intersect(pre2vgc,pre2targ)
nonSharedCids = setdiff(pre2targ,pre2vgc)






