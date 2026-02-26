function[r] = standardDifCI(dat1,dat2,ciFrac)

if ~exist('ciFrac','var')
    ciFrac = 0.95;
end

r.numSEs = norminv(1-(1-ciFrac)/2);

r.n1 = length(dat1);
r.m1 = mean(dat1);
r.std1 = std(dat1);
r.se1 = r.std1/sqrt(r.n1);
r.ci1 = [r.m1-r.se1*r.numSEs r.m1+r.se1*r.numSEs];

r.n2 = length(dat2);
r.m2 = mean(dat2);
r.std2 = std(dat2);
r.se2 = r.std2/sqrt(r.n2);
r.ci2 = [r.m2-r.se2*r.numSEs r.m2+r.se2*r.numSEs];

%%Calculate difference range
r.md = r.m2-r.m1;
r.sed = sqrt(r.se1^2+r.se2^2);
r.ciD = [r.md-r.sed*r.numSEs r.md+r.sed*r.numSEs];
