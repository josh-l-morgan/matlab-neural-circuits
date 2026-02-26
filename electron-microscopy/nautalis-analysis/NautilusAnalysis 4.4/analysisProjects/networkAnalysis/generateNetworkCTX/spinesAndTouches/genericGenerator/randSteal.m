function[makeSyn] = randSteal(synMat,touchMat);

%%Potential concern s
%%if axons are selected randomly

reportOn = 0;
imageOn = 0;
displayInterval = 1000000;  % number of iterations per randomaization before displaying status


%% get per cell synapse numbers (constraints)
aSynNum = sum(synMat,2)';
sSynNum = sum(synMat,1);



% %% reformat dends
% sortDend = sort(unique(spineparentids));
% lookupDend(sortDend) = 1:length(sortDend);
% parentDend = lookupDend(spineparentids);
% dendMat = repmat(parentDend,[size(touchMat,1) 1]);
% histDend = [.5:max(dendMat(:))]
% dendNum = max(dendMat(:));







% hitDend = double(dendMat).* synMat;
% for i = 1:size(dendMat,1)
%     grabHits = hitDend(i,:);
%     grabHits = grabHits(grabHits>0);
%     hGrab = hist(grabHits,histDend);
%     maxHit(i) = max(hGrab);
% end
% realHits = sum(maxHit);

%% Filter
% for i = 1:max(dendMat(:));
%     dSynNum(i) = length(find((dendMat==i).* synMat));
% end
% dendSynNumMat = dSynNum(dendMat);
% 
% numDend = max(dendMat(:));
% dendTouch = dendMat.*touchMat;
% for i = 1:size(dendMat,1)
%     histCon(i,:) =  hist(dendMat(i,touchMat(i,:)>0),.5:numDend);
% end
% multCon = histCon-1;
% multOpportunities = sum(sum(multCon(multCon>0)));
% 
% 
% useMat = repmat(aSynNum'>=minAxSyn,[1 size(touchMat,2)]);
% useMat = useMat .* dendSynNumMat>minDendSyn;
% 
% 
% %% Analyze
% 
% hitInd = find(synMat);
% [y x] = ind2sub(size(synMat),hitInd);
% dendVal = dendMat(hitInd);
% combIDs = y* max(dendVal*2) + dendVal;
% realUniques = length(unique(combIDs));
% realMults = length(hitInd)-realUniques;
% realSynLeft = length(hitInd);
% %
% % for i = 1:max(dendVal);
%     hitAx = y(dendVal==i);
%     uDend(i) = length(unique(hitAx));
% end
% realUniques = sum(uDend);
%
% 
% 
% 
% %% Get Cohorts
% commonSyn = matchDend(synMat,dendMat);
% randCommonSyn = zeros(dendNum,dendNum,reps);
% 
% %% Run model
% tic
% % sumSyn = double(synMat*0);
% for r = 1:reps
%     
    
    
    aSyn2go = aSynNum; %how many synapses left for each axon to deliver
    touch2go = touchMat.* repmat(sSynNum,[size(touchMat,1) 1]); % which touches are available for synapsing (win spine limited number)
    makeSyn = touchMat * 0; %new synaptic matrix
    trackSyn = zeros(1,4000); % keep track of how many synapses are left for the sake of understanding the process
    
    c = 0; %start counter
    while 1
        c = c+1;
        pickA = find((aSyn2go>0));% & (aSyn2go > max(aSyn2go)-30)); % find axons with synapses to deliver
        pickA = pickA(randi(length(pickA),1)); % pick an axon
        touches = find(touch2go(pickA,:)>0);  %% look for remaining touches
        
        if ~isempty(touches)
            pickSpine = touches(ceil(rand * length(touches))); %pick spine to synapse with
            touch2go(:,pickSpine) = touch2go(:,pickSpine) - 1; % reduce available synapses for that spine
            makeSyn(pickA,pickSpine) = makeSyn(pickA,pickSpine) +1; % mark the new synapse
            aSyn2go(pickA) = aSyn2go(pickA)-1; % reduce synapses left to distribute
        else
            oldTouches = find(touchMat(pickA,:) & ~(makeSyn(pickA,:))); % find all touches (occupied) (but not by self)
            pickSpine =  oldTouches(fix(ceil(rand * length(oldTouches))));  %pick occupied spine to steal
            
            %% wipe old
            oldAx = find(makeSyn(:,pickSpine)); % find the axon being stolen from
            pickOldAx = oldAx(ceil(rand * length(oldAx)));
            
            makeSyn(pickOldAx,pickSpine) = makeSyn(pickOldAx,pickSpine) -1; % clear the old synapse
            aSyn2go(pickOldAx) = aSyn2go(pickOldAx)+1; % ad synapse to distribution list for axon that was stolen from
            
            %% add new
            makeSyn(pickA,pickSpine) = makeSyn(pickA,pickSpine) + 1; % mark new synapse
            aSyn2go(pickA) = aSyn2go(pickA)-1; % subtract synapse from distribution list of active axon
            
        end
        
        remainingSyn= sum(aSyn2go); %check distribution list to see if process is finished
        trackSyn(c) = remainingSyn;
        
        if ~sum(aSyn2go) % if no synapses left to place exit loop
            break
        end
        
        if ~mod(c,displayInterval) %disply process state
            subplot(2,1,1)
            disp(sprintf('%d syn remaining on cycle %d.',remainingSyn,c));
            image(makeSyn);
            pause(.01)
        end
    end
    
    
    %% Report
    if reportOn
        
        checkTouch = (makeSyn>0) - (touchMat>0);
        badSyn = sum(checkTouch(:)>0);
        remainingSyn= sum(aSyn2go);
        randSynNum = sum(makeSyn,2)';
        mistakeNum = sum(abs(aSynNum-randSynNum));
        randSpineNum = sum(makeSyn,1);
        badSpine = sum(randSpineNum ~= sSynNum);
        disp(sprintf('%02.1f%% done,%d bad syn, %d axon mistakes, %d spine mistakes after %d cycles', r/reps*100,badSyn,mistakeNum,badSpine,c))
        pause(.0001)
        
        hitInd = find(makeSyn);
        [y x] = ind2sub(size(synMat),hitInd);
        dendVal = dendMat(hitInd);
        combIDs = y* max(dendVal*2) + dendVal;
        numUniques = length(unique(combIDs));
        randUniques(r) = numUniques;
        randMults(r) = length(hitInd)-numUniques;
        randSynLeft(r) = length(hitInd);
        %ndUDend(r,:));
        
        sumSyn = sumSyn + makeSyn;
         randCommonSyn(:,:,r) = matchDend(makeSyn,dendMat);
        
    end
    
    if imageOn
        
        col = uint8(synMat*1000);
        col(:,:,2) = synMat*0;
        col(:,:,3) = touchMat*1000;
        
        %%Show before and after matrix
        subplot(2,1,1)
        col(:,:,2) = makeSyn;
        image(col)
        
        %%Show synapse distribution history
        subplot(2,1,2)
        plot(trackSyn)
        pause(.0001)
        
    end
    
end

