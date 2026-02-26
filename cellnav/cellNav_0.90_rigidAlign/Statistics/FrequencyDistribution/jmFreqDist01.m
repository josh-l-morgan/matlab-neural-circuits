
%{
datC1 = []
datC2 = []
datI1 = []
datI2 = []
%}

%% what is the frequency of a given difference producing some result?
%% what is the range that contains 95% of results at some observation?


dat1 = datI1;
dat2 = datI2;
testDif = -.2:.001:.2;
reps = 100000
hitNum = 100000;



n1 = length(dat1)
n2 = length(dat2)

m1 = median(dat1)
m2 = median(dat2)


md = m2-m1;


std1 = std(dat1)

testNum = length(testDif)

rmd = zeros(testNum,reps);

for td = 1:testNum
    disp(sprintf('testing %d of %d',td,testNum))

    r1 = randn(n1,reps)*std1+m1;
    r2 = randn(n2,reps)*std1+m1+m1*testDif(td);

    rm1 = median(r1,1);
    rm2 = median(r2,1);

    rmd(td,:) = rm2-rm1;

end


%% Find best fit
rmTd = mean(rmd,2);
difRvO = abs(rmTd-md);
minDifRvO = min(difRvO)
bestTarg = find(difRvO==minDifRvO);
bestFit = testDif(bestTarg)


%% Find range by fraction

clf
subplot(2,1,1)
allDifs = abs(rmd-md);

sortDifs = sort(abs(allDifs(:)),'ascend');
hitThresh = sortDifs(hitNum);
inRange = abs(allDifs)<=hitThresh;



%inRange =  abs(allDifs)<=(std1/sqrt(n2));
testHits = sum(inRange,2)/hitNum;
plot(testDif,testHits)
cumHit = cumsum(testHits)/sum(testHits(:));
subplot(2,1,2)
plot(testDif,cumHit)
lowTarg = max(find(cumHit<.025));
highTarg = min(find(cumHit>.975));
bestRange = [testDif(lowTarg) testDif(highTarg)]
drawnow









