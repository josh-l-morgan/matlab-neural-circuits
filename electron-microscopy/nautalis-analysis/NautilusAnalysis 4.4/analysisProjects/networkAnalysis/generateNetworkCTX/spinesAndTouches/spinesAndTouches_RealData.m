%why dont we just scramble the identity of the spines

clear all
colormap gray(256)
reps = 100;
displayInterval = 10000;
reportOn = 1;

load('singlespinematrices2.mat')
synMat2 = singlespinesynapsematrix>0;
touchMat2 = singlespinetouchmatrix>0;

c = 0;
clear synMat touchMat
for i = 1:size(synMat2,2)
    if sum(synMat2(:,i)) == 1
        c = c+1;
        synMat(:,c) = synMat2(:,i);
        touchMat(:,c) = touchMat2(:,i);
    end
end



%% Dummy Data
%{
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
%}

aSynNum = sum(synMat,2)';

%% Run model
tic

for r = 1:reps
    
    syn2go = aSynNum;
    touch2go = touchMat;
    makeSyn = touchMat * 0;
    trackSyn = zeros(1,4000);
    
    
    c = 0;
    while 1
        c = c+1;
        pickA = find((syn2go>0));% & (syn2go > max(syn2go)-30));
        pickA = pickA(randi(length(pickA),1));
        touches = find(touch2go(pickA,:));  %% look for remaining touches
        return
        dends = dendMat(pickA,(touch2go(pickA,:)>0)) 
        
        
        if ~isempty(touches)
            pickSpine = touches(randi(length(touches),1));
            dends = dendMat
            touch2go(:,pickSpine) = 0;
            makeSyn(pickA,pickSpine) = 1;
            syn2go(pickA) = syn2go(pickA)-1;
        else
            oldTouches = find(touchMat(pickA,:));
            pickSpine =  oldTouches(fix(randi(length(oldTouches),1)));
            
            %% whipe old
            oldAx = find(makeSyn(:,pickSpine));
            
            makeSyn(oldAx,pickSpine) = 0;
            syn2go(oldAx) = syn2go(oldAx)+1;
            
            %% add new
            makeSyn(pickA,pickSpine) = 1;
            syn2go(pickA) = syn2go(pickA)-1;
            
        end
        
        remainingSyn= sum(syn2go);
        %trackSyn(c) = remainingSyn;
        
        if ~sum(syn2go) % if no synapses left to place
            break
        end
        
        %     if ~mod(c,displayInterval)
        %         subplot(2,1,1)
        %         disp(sprintf('%d syn remaining on cycle %d.',remainingSyn,c));
        %         image(makeSyn*10000);
        %         pause(.01)
        %     end
    end
    
    
    %% Report
    if reportOn
        
        checkTouch = makeSyn - touchMat;
        badSyn = sum(checkTouch(:)>0);
        remainingSyn= sum(syn2go);
        randSynNum = sum(makeSyn,2)';
        mistakeNum = sum(abs(aSynNum-randSynNum));
        randSpineNum = sum(makeSyn,1);
        badSpine = sum(randSpineNum~=1);
        disp(sprintf('%d bad syn, %d axon mistakes, %d spine mistakes after %d cycles', badSyn,mistakeNum,badSpine,c))
    end
    
    subplot(2,1,1)
    
    col = synMat;
    col(:,:,2) = makeSyn;
    col(:,:,3) = touchMat;
    
    image(uint8(col*1000))
    
    
    subplot(2,1,2)
    plot(trackSyn)
    pause(.001)
    
    
end
toc