function[e] = getEdge(colMean)

global watProps


crossWin1 = round(length(colMean) * watProps.win1);
crossWin2 = round(length(colMean) * watProps.win2);

%% dummy

fSize = 11;
filt1 = ones(fSize,1);
filt1 = filt1/sum(filt1);
%fMean = filter(filt1,1,colMean)
fMean = imfilter(colMean,filt1,'same');


if 0
    
    plot(colMean)
    hold on
    plot(ones(length(colMean),1) * watProps.thresh1)
    plot(fMean)
    plot([crossWin1 crossWin1],[0 70],'b')
    plot([crossWin2 crossWin2],[0 70],'c')

    plot(ones(length(colMean),1) * watProps.thresh1)

    hold off
end

bMean = fMean > watProps.thresh1;

cross = bMean(2:end) -  bMean(1:end-1);
crosses = find(cross);

e.mask = bMean;
e.croses = crosses;

goodCross = find(cross(crossWin1 :crossWin2));

if isempty(crosses)
    firstCross = length(colMean);
else
    firstCross = crosses(1);
end

darkWindow = sum(~bMean(crossWin1:crossWin2));
tooHigh = darkWindow == 0;

tooLow = (darkWindow>0) & isempty(goodCross);

e.tooLow = tooLow;
e.tooHigh = tooHigh;
e.darkWindow = darkWindow;
e.firstCross = firstCross;




