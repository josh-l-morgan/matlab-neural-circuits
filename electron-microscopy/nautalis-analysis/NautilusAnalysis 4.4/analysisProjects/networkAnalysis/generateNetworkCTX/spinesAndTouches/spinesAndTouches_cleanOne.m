%why dont we just scramble the identity of the spines

clear all

%% Dummy Data
axNum = 100;
spineNum = 700;
touchNum = 10;

%%make touchMat
while 1
touchMat = rand(axNum,spineNum);
touchMat = touchMat<(touchNum/axNum);
sumTouch = sum(sum(touchMat,1)==0);
    if ~sumTouch
        break
    end
end

%%Make Synapses
synMat = touchMat*0;
for i = 1:spineNum
    touches = find(touchMat(:,i));
    pickAx = touches(fix(rand*length(touches))+1);
    synMat(pickAx,i) = 1;
end
aSynNum = sum(synMat,2)';


%% Run model

syn2go = aSynNum;
touch2go = touchMat;
touchSyn = touchMat * 0;
c = 0;
while 1
    c = c+1;
    pickA = find((syn2go>0));% & (syn2go > max(syn2go)-30));
    pickA = pickA(randi(length(pickA),1));
    touches = find(touch2go(pickA,:));
    
    if ~isempty(touches)
        pickTouch = touches(randi(length(touches),1));
        touch2go(:,pickTouch) = 0;
        touchSyn(pickA,pickTouch) = 1;
        syn2go(pickA) = syn2go(pickA)-1;
    else
        oldTouches = find(touchMat(pickA,:));
        pickTouch =  oldTouches(fix(randi(length(oldTouches),1)));
        %% whipe old
        oldAx = find(touchSyn(:,pickTouch));
        
        touchSyn(oldAx,pickTouch) = 0;
        syn2go(oldAx) = syn2go(oldAx)+1;
        touch2go(oldAx,pickTouch) = 1; %might be faster without
        
        %% add new
        touchSyn(pickA,pickTouch) = 1;
        syn2go(pickA) = syn2go(pickA)-1;
        touch2go(:,pickTouch) = 0;
        
    end
        
    if ~sum(syn2go) % if no synapses left to place
        'finished'
        break
    end
    
    if ~mod(c,100001)
        remainingSyn= sum(syn2go);
        disp(sprintf('%d syn remaining on cycle %d.',remainingSyn,c));
        image(touchSyn*10000);
        pause(.01)
    end
end

remainingSyn= sum(syn2go);
randSynNum = sum(touchSyn,2)';
mistakeNum = sum(abs(aSynNum-randSynNum));
spineSynNum = sum(touchSyn,1);
badSpine = sum(spineSynNum ~=1);

disp(sprintf('%d remaining syn and %d mistakes and %d bad spines',remainingSyn,mistakeNum,badSpine))

col = synMat;
col(:,:,2) = touchSyn;
col(:,:,3) = touchMat;

image(uint8(col*1000))
pause(.010)