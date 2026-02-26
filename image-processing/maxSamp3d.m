

SPN = 'E:\Pratyush\EMprocStacks\Overview1\FullStack\'
TPN = 'E:\Pratyush\EMprocStacks\Overview1\FullStack_3Dmax16\'
mkdir(TPN)
binNum = 8;


inam = dir([SPN '*.tif'])


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
        subStack(:,:,p) = imread([SPN inam(targ).name]);
    end
    
    medPlan = max(subStack,[],3);
    image(medPlan)
    filename = sprintf('%s%05.0f.png',TPN,i);
    imwrite(uint8(medPlan),filename)
    medStack(:,:,i) = medPlan;
    pause(.01)
end