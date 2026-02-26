function[makeSyn] = randBiasSteal(synMat,touchMat,dendMat,prefRat);


%% Run model
tic
displayInterval = 10000;

axNum = size(synMat,1);
dendNum = max(dendMat(:));
aSynNum = sum(synMat,2);
sSynNum = sum(synMat,1);
synNum = sum(synMat(:));


aSyn2go = aSynNum; %how many synapses left for each axon to deliver
touch2go = touchMat.* repmat(sSynNum,[size(touchMat,1) 1]); % which touches are available for synapsing (win spine limited number)
makeSyn = touchMat * 0; %new synaptic matrix
trackSyn = zeros(1,4000); % keep track of how many synapses are left for the sake of understanding the process
ax2Dend = zeros(axNum,dendNum); %track how many of each dendrite each axon synapses with.
axPref = rand(axNum,dendNum);

for i = 1:axNum  %build spine sized preference matrix out of ax Pref
    prefMat(i,:) = axPref(i,dendMat(i,:));
end

c = 0; %start counter
prefCount = [0 0; 0 0];
while 1
    c = c+1;
    
    
    %% eliminate extra synapse
    newSynNum = sum(makeSyn(:));
    if newSynNum >= synNum
        [y x] = find(makeSyn>0);
        pickSyn = y(ceil(rand*length(y)));
        makeSyn(y(pickSyn),x(pickSyn))  = makeSyn(y(pickSyn),x(pickSyn)) -1;
        touch2go(:,x(pickSyn)) = touch2go(:,x(pickSyn)) + touchMat(:,x(pickSyn))>0;
    end
    
    currentAxSyn = sum(makeSyn,2);
    noSyn = find(currentAxSyn == 0);
    if isempty(noSyn)
        pickA = ceil(rand*axNum); % pick an axon
    else
         pickA = noSyn(ceil(rand* length(noSyn)));
   
    end
    touches = find(touch2go(pickA,:)>0);  %% look for remaining touches
    
    
    if sum(touches)
        
        if prefRat> rand
            hitPref = prefMat(pickA,touches);
            prefHits = find((hitPref== max(hitPref)));
            pickSpine = touches(prefHits(ceil(rand*length(prefHits)))); %random pick from prefered
            prefCount(1,1) = prefCount(1,1) + 1;
        else
            pickSpine =  touches(fix(ceil(rand * length(touches))));  %pick occupied spine to steal
            prefCount(1,2) = prefCount(1,2) + 1;
        end
        
        touch2go(:,pickSpine) = touch2go(:,pickSpine) - 1; % reduce available synapses for that spine
        makeSyn(pickA,pickSpine) = makeSyn(pickA,pickSpine) +1; % mark the new synapse
        aSyn2go(pickA) = aSyn2go(pickA)-1; % reduce synapses left to distribute
        
        
    else
        oldTouches = find(touchMat(pickA,:) & ~(makeSyn(pickA,:))); % find all touches (occupied) (but not by self)
        if ~isempty(oldTouches)
            if prefRat > rand
                hitPref = prefMat(pickA,oldTouches);
                prefHits = find((hitPref== max(hitPref)));
                pickSpine = oldTouches(prefHits(ceil(rand*length(prefHits)))); %random pick from prefered
                prefCount(2,1) = prefCount(2,1) + 1;
            else
                pickSpine =  oldTouches(fix(ceil(rand * length(oldTouches))));  %pick occupied spine to steal
                prefCount(2,2) = prefCount(2,2) + 1;
            end


            %% wipe old
            oldAx = find(makeSyn(:,pickSpine)); % find the axon being stolen from
            pickOldAx = oldAx(ceil(rand * length(oldAx)));

            makeSyn(pickOldAx,pickSpine) = makeSyn(pickOldAx,pickSpine) -1; % clear the old synapse
            aSyn2go(pickOldAx) = aSyn2go(pickOldAx)+1; % ad synapse to distribution list for axon that was stolen from

            %% add new
            makeSyn(pickA,pickSpine) = makeSyn(pickA,pickSpine) + 1; % mark new synapse
            aSyn2go(pickA) = aSyn2go(pickA)-1; % subtract synapse from distribution list of active axon
        end
    end
    
    remainingSyn= sum(aSyn2go); %check distribution list to see if process is finished
    trackSyn(c) = remainingSyn;
    
    if ~sum(sum(makeSyn,2)==0) & (sum(makeSyn(:)) == synNum)% if every axon has a synapse
        break
    end
    
    if ~mod(c,displayInterval) %disply process state
        disp(sprintf('%d syn remaining on cycle %d.',remainingSyn,c));
        pause(.0001)
    end
end



