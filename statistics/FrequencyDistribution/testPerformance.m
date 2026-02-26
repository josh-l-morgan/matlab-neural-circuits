
%% Set test parameters
param.confidenceFraction = 0.95;
param.convertToPercent = 1;
param.useMean = 1;
param.testResolution = 40; % Number of tests values per standard error of the differences
param.reps = 100000; % number of repetitions per search value
param.hitNum = 100000; % number of hits included in final range
param.showProcess = 0;

param1 = param;
param2 = param;
param2.useMean = 0;



%% Set experiment conditions
m2s = [120];%[ 100 105 110 120 140 180];
m1s = 100;
n1s = [3 5 10 20 40];
n2s = [3 5 10 20 40];
std1s = [10];
std2s = [1 5 10 20 40];

sameNs = 1;
sameSTDs = 1;

if sameNs
    n2sL = 1;
else
    n2sL = length(n2s);
end

if sameSTDs
    std2L = 1;
else
    std2L = length(std2s)
end

allRes = zeros(5,length(m2s),length(n1s),n2sL,length(std1s),2);

expNum = length(allRes(:))/10;

c = 0;
m1 = m1s(1);
for im2 = 1:length(m2s)
    m2 = m2s(im2);
    for in1 = 1:length(n1s)
            n1 = n1s(in1);
        

        for in2 = 1:n2sL;
            if sameNs
                n2 = n1;
            else
                n2 = n2s(in2);
            end

            for istd1 = 1:length(std1s)
                std1 = std1s(istd1);


                for istd2 = 1:length(std2s)
                    if sameSTDs
                        std2 = std1;
                    else
                        std2 = std2s(istd2);
                    end
                
                    c = c+1;
                    disp(sprintf('running experiment %d of %d',c,expNum))
                    %% Run test
                     dat1 = randn(n1,1)*std1 + m1;
                     dat2 = randn(n2,1)*std2 + m2;
                
                     res1 = differenceRangeA1(dat1,dat2,param1);
                     res2 = differenceRangeA1(dat1,dat2,param2);

                     allRes(1,im2,in1,in2,istd1,:) = res1.bestRange;
                     allRes(2,im2,in1,in2,istd1,:) = res2.bestRange;
                     allRes(3,im2,in1,in2,istd1,:) = res2.standardDifferenceRange;
                     allRes(4,im2,in1,in2,istd1,:) = ((m2-m1)/m1)*100;
                     allRes(5,im2,in1,in2,istd1,:) = res1.observedDifference;
                
                    
                
                end
            end
        end
    end
end


disp('Finished running tests')

%% Show results
%allRes = zeros(5,length(m2s),length(n1s),n2sL,length(std1s),2);

plotN1 = 1;
plotS1 = 1;

numVar = plotN1 + plotS1;

classCol = [1 0 0; 1 0 1; 0 1 0; 0 0 0; 0 0 1];
classNames = {'NormModel','NonParaModel','StandardDifference','TargetDifference','RealDifference'};

if plotN1

    for s = 1:length(std1s)

        subplot(length(std1s),numVar,(s-1)*numVar+1)
        vals = squeeze(allRes(:,1,:,1,s,:));

        cla
        hold on
        for p = 1:size(vals,1)
            plot(n1s,vals(p,:,1),'color',classCol(p,:));
        end
        for p = 1:size(vals,1)
            plot(n1s,vals(p,:,2),'color',classCol(p,:));
        end
        title(sprintf('ns for std = %d',std1s(s)))

    end
end




if plotS1

    for s = 1:length(n1s)

        subplot(length(n1s),numVar,(s)*numVar)
        vals = squeeze(allRes(:,1,s,1,:,:));

        cla
        hold on
        for p = 1:size(vals,1)
            plot(std1s,vals(p,:,1),'color',classCol(p,:));
        end
        for p = 1:size(vals,1)
            plot(std1s,vals(p,:,2),'color',classCol(p,:));
        end
        title(sprintf('stds for n = %d',n1s(s)))
    end
end











