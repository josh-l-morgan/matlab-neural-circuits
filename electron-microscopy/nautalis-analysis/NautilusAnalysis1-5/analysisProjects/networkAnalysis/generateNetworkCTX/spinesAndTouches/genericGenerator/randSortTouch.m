function[synMat] = randSortTouch(synMat,touchMat);


axNum = size(synMat,1);
spNum = size(synMat,2);
aSynNum = sum(synMat,2);
sSynNum = sum(synMat,1);
synNum = sum(synMat(:));
touchNum = sum(touchMat(:));


%% Make matrix with correct number of synapses
makeSyn = synMat * 0;
randMat = rand(size(touchMat));

touchInd = find(touchMat>0);
[sortedRand randOrd] = sort(rand(size(touchInd)));

makeSyn(touchInd(randOrd(1:synNum))) = 1;


%%  find differences to constraints 

aRandSyn = sum(makeSyn,2);
sRandSyn = sum(makeSyn,1);

difAx = aSynNum - aRandSyn;
difSp = sSynNum - sRandSyn;


%% run rearrangement

c = c + 1;
while 1
    tI = mod(c,touchNum);
    touchNow = touchInd(randOrd(tI)); %Get index of current touch.
    
    
    
    if ~sum([difAx' difSp]) %if there are no differences
        break
    end
    
end





















