

if 0
    dat = [];
end

depth = dat(:,1);
cal = dat(:,2);
scatter(depth,cal)

uDepth = unique(depth);

clear gMeans
for u = 1:length(uDepth)
    gMeans(u) = mean(cal(depth==uDepth(u)));
end

g1 = cal(depth==45);
g2 = cal(depth==50);

[h,p,ci,stats] = ttest2(g1,g2)












