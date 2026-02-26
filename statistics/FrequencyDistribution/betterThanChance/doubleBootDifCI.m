function[res] = bootDifCI(dat1,dat2,ciFrac,reps)

%%Generate two bootstrap distributions with replacement and find
%%Confidence interval for their median difference using a confidence
%%interval that covers ciFrac portion of results. Itteration number = reps
dat1 = dat1(:);
dat2 = dat2(:);

if ~exist('ciFrac','var')
    ciFrac = .95;
end


if ~exist('reps','var')
    reps = 10000;
end

obMedDif = median(dat2) - median(dat1);
obMeanDif = mean(dat2) - mean(dat1);

res.ciFrac = ciFrac;
res.reps = reps;
res.observedMedianDifference = obMedDif;
res.observedMeanDifference = obMeanDif;



n1 = length(dat1);
n2 = length(dat2);

p1 = ceil(rand(n1,reps)*n1);
p2 = ceil(rand(n2,reps)*n2);

res.r1 = dat1(p1);
res.r2 = dat2(p2);

m1 = median(res.r1,1);
m2 = median(res.r2,1);

mean1 = mean(res.r1,1);
mean2 = mean(res.r2,1);

d = sort(m2-m1);
dMean = sort(mean2-mean1);

a = 1-ciFrac;

lowB = round(reps*a);
highB = length(d) - round(reps * a);
lowB = max(lowB,1);
highB = min(highB,length(d));

ciRaw = [d(lowB) d(highB)];
ciMeanRaw = [dMean(lowB) dMean(highB)];

%%pivot
ci = 2*obMedDif - fliplr(ciRaw);
ciMean = 2* obMeanDif - fliplr(ciMeanRaw);

res.ciMedian = ci;
res.ciMean = ciMean;
res.ciMedianNoPivot = ciRaw;
res.ciMeanNoPivot = ciMeanRaw;





