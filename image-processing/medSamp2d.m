

SPN = 'D:\KarlsRetina\HxQ\tempStack\'
TPN = 'D:\KarlsRetina\HxQ\tempStack_med\'
binNum = 16;


inam = dir([SPN '*.png'])


I = imread([SPN inam(1).name]);
[ys xs a] = size(I);

secNum = length(inam);

numPlan = floor(secNum/binNum);
medStack = zeros(ys, xs, numPlan);
for i = 1:numPlan
    sprintf('%d of %d',i,numPlan)
    subStack = zeros(ys, xs, binNum);
    for p = 1:binNum
        targ = (i-1)*binNum+p;
        subStack(:,:,p) = 
    end
    I = imread([SPN inam(targ).name]);
    medPlan = median(subStack,3);
    image(medPlan)
    filename = sprintf('%s%05.0f.png',TPN,i);
    imwrite(uint8(medPlan),filename)
    pause(.01)
end