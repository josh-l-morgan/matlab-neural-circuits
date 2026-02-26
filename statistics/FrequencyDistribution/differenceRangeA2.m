function[res] = differenceRangeJM(datIn1,datIn2,param)
%%Code by Josh L. Morgan 2023
%%Code for generating a range of effect magnitudes that are likely to
%%generate a given difference between two groups
%%Note: useMean = 1 depends on standard deviations and randn to model data.
%%Input parameters as structure with fields confidenceFraction,
%%convertToPercent, useMean, testResolution, reps, hitNum, showProcess
%%what is the frequency of a given difference producing some result?
%%what is the range that contains 95% of results at some observation?

dat1 = datIn1;
dat2 = datIn2;


%% Define variables

if ~exist('param','var')
    param.default = 1;
end

if ~isfield(param,'confidenceFraction')
    param.confidenceFraction = 0.95;
end

if ~isfield(param,'useMean')
    param.useMean = 0;
end

if ~isfield(param,'testResolution')
    param.testResolution = 40; % Number of tests values per standard error of the differences
end

if ~isfield(param,'reps')
    param.reps = 100000; % number of repetitions per search value
end

if ~isfield(param,'hitNum')
    param.hitNum = 100000; % number of hits included in final range
end

if ~isfield(param,'showProcess')
    param.showProcess = 1;
end

if ~isfield(param,'testDifferences')
    param.testDifferences = 1;
end

reps = param.reps;
confidenceFraction = param.confidenceFraction;

numSEs = norminv(1-(1-confidenceFraction)/2);
testBuffer = numSEs*1.5; %Number of standard errors of difference to cover with test values


%% Get basic data statistics
n1 = length(dat1);
n2 = length(dat2);

if 1
    m1 = mean(dat1);
    m2 = mean(dat2);
else
    m1 = median(dat1);
    m2 = median(dat2);
end

std1 = std(dat1);
std2 = std(dat2);
se1 = std1/sqrt(n1);
se2 = std2/sqrt(n2);


%% Define test statistic
md = m2-m1;

%% Define search range

%%Standard error of difference
seDif = sqrt(se1^2+se2^2);


res.observedDifference = md;
res.standardErrorOfDifference = seDif;
res.standardDifferenceRange = [md-seDif*numSEs md+seDif*numSEs];
res.parameters = param;


if param.testDifferences

    %%Generate range of difference values to test. If test range is not user
    %%defined (param.testDif), then cover the observed mean difference plus or
    %%minus the standard error of the difference times testBuffer (3)
    if ~isfield(param,'testDif')
        if seDif == 0 %no error
            testRange = [md-m1 md+m1]; % no sense
            testStep = m1/param.testResolution/3;
        else
            testRange = [md-seDif*testBuffer md+seDif*testBuffer];
            testStep = seDif/param.testResolution;
        end
        testDif = [testRange(1) : testStep : testRange(2)];
    else
        testDif = param.testDif;
    end

    %% Run tests
    testNum = length(testDif);

    rmd = zeros(testNum,reps); % matrix for holding all test results
    tic
    disp('running randomizations')
    r1 = zeros(n1,reps);
    r2 = zeros(n2,reps);
    for td = 1:testNum

        %print progress every 10 seconds
        t = toc;
        if t>10
            tic
            disp(sprintf('testing %d of %d',td,testNum))
        end

        if param.useMean % Generate distribution using randn normal distribution
            r1 = randn(n1,reps)*std1  + m1; % Make randomized group 1
            r2 = randn(n2,reps)*std2 + m1 + testDif(td); % make randomized group 2
        else % Generate distribution using random picks from group1
           
            rPool = [dat1(:)-m1; dat2(:)-m2]; %pool groups on assumption that testDif is only difference
            np = length(rPool);
            pickMat = ceil(rand(np,reps)*np);
            r1 = rPool(pickMat) + m1 ;

            pickMat = ceil(rand(np,reps)*np);
            r2 = rPool(pickMat) + m1 + testDif(td);
        end


        if 1
            rm1 = mean(r1,1);
            rm2 = mean(r2,1);
        else
            rm1 = median(r1,1);
            rm2 = median(r2,1);
        end

        rmd(td,:) = rm2-rm1;

    end
    disp('finished randomizations')


    %% Find range by fraction

    allDifs = abs(rmd-md); % Find difference between random and real results
    sortDifs = sort(abs(allDifs(:)),'ascend');
    hitThresh = sortDifs(param.hitNum); % Find absolute difference threshold that includes n (hitNum) random results
    inRange = abs(allDifs)<=hitThresh; % Select random results within range

    testHits = sum(inRange,2)/param.hitNum; % Get distribution of hits for each test value

    %%Find 95% range for results
    cumHit = cumsum(testHits)/sum(testHits(:));
    lowTarg = max(find(cumHit<= ((1-confidenceFraction)/2)));
    highTarg = min(find(cumHit>= ((1+confidenceFraction)/2)));
    bestRange = [testDif(lowTarg) testDif(highTarg)];


else

    bestRange = [0 0];
    testDif = 0;
    testHits = 0;

end

res.bestRange = bestRange;
res.testedDifferences = testDif;
res.hitsForDifferences = testHits;

%% Show results
disp(sprintf(['\nMean difference = %f \n' ...
    'Standard error of differences = \n%f to %f \n'...
    'Magnitudes composing 95%% of observed results = \n'...
    '%f to %f'],md,md-seDif,md+seDif,bestRange(1), bestRange(2)));

if param.showProcess
    try
        close(figDifRange);
    end
    figDifRange = figure;
    clf

    subplot(1,5,1)
    scatter(dat1*0+1,dat1,'k')
    hold on
    scatter(dat2*0+2,dat2,'k')
    hold off
    xlim([0 3])
    %ylim([0 max([dat1(:);dat2(:)])])


    subplot(1,5,2:5)
    plot(testDif,testHits)
    hold on
    maxY = max(testHits);
    plot([bestRange(1) bestRange(1)],[0 maxY],'r')
    plot([bestRange(2) bestRange(2)],[0 maxY],'r')
    hold off
    %     subplot(2,1,2)
    %     plot(testDif,cumHit)
    drawnow
end




%% CUT
%{

%% Find best fit
if 0
    rmTd = mean(rmd,2);
    difRvO = abs(rmTd-md);
    minDifRvO = min(difRvO)
    bestTarg = find(difRvO==minDifRvO);
    bestFit = testDif(bestTarg)
end

%inRange =  abs(allDifs)<=(std1/sqrt(n2));






%}




