function[sI] = erodePeaks(sI,tooHigh,dilateNum)

if ~exist('dilateNum')
    dilateNum = 5;
end

image(255-sI)
pause(1)
%% expand maxes
for i = 1:5
maxed = find(sI>=tooHigh);
shiftInd = getShiftSubFun(size(sI),maxed);
sI(shiftInd) = tooHigh;
image(255-sI)
pause(1)
end
%% erode
for i = 1:100
maxed = find(sI>=tooHigh);
size(maxed)
if isempty(maxed)
    break
end
shiftInd = getShiftSubFun(size(sI),maxed);
shiftVals = sI(shiftInd);

shiftVals(shiftVals >= tooHigh) = 0;
hasVal = shiftVals>0;
numVals = sum(hasVal,2);
passVals = numVals>2;
sumVals = sum(shiftVals(find(passVals),:),2);
meanVals = sumVals./numVals(passVals);

sI(maxed(passVals)) = meanVals;

image(255 - sI)
pause(.1)
end

%% 
function[shiftInd] = getShiftSubFun(sizeI,indList)


[y x] = ind2sub(sizeI, indList);

yshift = [-1 0 1 -1  1 -1 0 1];
xshift = [-1 -1 -1 0  0 1 1 1];

yshifted = repmat(y,[1 8]) + repmat(yshift,[length(y) 1]);
xshifted = repmat(x,[1 8]) + repmat(xshift,[length(x) 1]);

ywall = (yshifted<1) | (yshifted>sizeI(1));
xwall = (xshifted<1) | (xshifted>sizeI(2));

yshifted(ywall) = y(1);
xshifted(xwall) = x(1);

shiftInd = sub2ind(sizeI, yshifted,xshifted);











