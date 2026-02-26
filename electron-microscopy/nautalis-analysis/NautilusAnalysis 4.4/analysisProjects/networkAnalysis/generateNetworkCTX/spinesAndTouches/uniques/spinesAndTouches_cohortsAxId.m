%%Potential concern s
%%if axons are selected randomly

%% Required input
% Every axon and spine must be able to make at least one synapse

clear all
colormap gray(256)

reps = 1000;  % number of randomizations
displayInterval = 1000000;  % number of iterations per randomaization before displaying status
reportOn = 1; % display confirmation of proper matrix after each randomization
imageOn = 0;

minAxSyn = 2;
minDendSyn = 2;

%% get touch and synapse matrix
load('singlespinematrices2.mat')
synMat = singlespinesynapsematrix>0;
touchMat = singlespinetouchmatrix>0;
spineparentids;

totSyn= sum(synMat(:));

%% get per cell synapse numbers (constraints)
aSynNum = sum(synMat,2)';
sSynNum = sum(synMat,1);



%% reformat dends
sortDend = sort(unique(spineparentids));
lookupDend(sortDend) = 1:length(sortDend);
parentDend = lookupDend(spineparentids);
dendMat = repmat(parentDend,[size(touchMat,1) 1]);
histDend = [.5:max(dendMat(:))];




% hitDend = double(dendMat).* synMat;
% for i = 1:size(dendMat,1)
%     grabHits = hitDend(i,:);
%     grabHits = grabHits(grabHits>0);
%     hGrab = hist(grabHits,histDend);
%     maxHit(i) = max(hGrab);
% end
% realHits = sum(maxHit);

%% Filter
for i = 1:max(dendMat(:));
    dSynNum(i) = length(find((dendMat==i).* synMat));
end
dendSynNumMat = dSynNum(dendMat);

numDend = max(dendMat(:));
dendTouch = dendMat.*touchMat;
for i = 1:size(dendMat,1)
    histCon(i,:) =  hist(dendMat(i,touchMat(i,:)>0),.5:numDend);
end
multCon = histCon-1;
multOpportunities = sum(sum(multCon(multCon>0)));


useMat = repmat(aSynNum'>=minAxSyn,[1 size(touchMat,2)]);
useMat = useMat .* dendSynNumMat>minDendSyn;


%% Analyze

hitInd = find(synMat);
[y x] = ind2sub(size(synMat),hitInd);
dendVal = dendMat(hitInd);
combIDs = y* max(dendVal*2) + dendVal;
realUniques = length(unique(combIDs));
realMults = length(hitInd)-realUniques;
realSynLeft = length(hitInd);

uComb = unique(combIDs);
clear axId
for i = 1:length(uComb)
   sameComb =  find(combIDs==uComb(i));
   countComb(i) = length(sameComb);
   axId(i) = y(sameComb(1));
end

multCombs = countComb(countComb>1);
sum(multCombs)-length(multCombs)
synNumWithMults = sum(multCombs);
axMults = axId(countComb>1);
axNumWithMults = length(unique(axMults));
pause

%% Count axon



%
% for i = 1:max(dendVal);
%     hitAx = y(dendVal==i);
%     uDend(i) = length(unique(hitAx));
% end
% realUniques = sum(uDend);
%



%% Get Cohorts
[cohortSum] = findCohorts(synMat,touchMat,dendMat)

%% Run model
tic
sumSyn = double(synMat*0);
for r = 1:reps
    
    
    
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
            
            %% whipe old
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
        badSpine = sum(randSpineNum~= sSynNum);
        disp(sprintf('%02.1f%% done,%d bad syn, %d axon mistakes, %d spine mistakes after %d cycles', r/reps*100,badSyn,mistakeNum,badSpine,c))
        pause(.0001)
        
        
        %         hitDend = double(dendMat).* makeSyn;
        %         for i = 1:size(dendMat,1)
        %             grabHits = hitDend(i,:);
        %             grabHits = grabHits(grabHits>0);
        %             hGrab = hist(grabHits,histDend);
        %             maxHit(i) = max(hGrab);
        %         end
        %         randHits(r) = sum(maxHit);
        
        
        
        
        hitInd = find(makeSyn);
        [y x] = ind2sub(size(synMat),hitInd);
        dendVal = dendMat(hitInd);
        combIDs = y* max(dendVal*2) + dendVal;
        numUniques = length(unique(combIDs));
        randUniques(r) = numUniques;
        randMults(r) = length(hitInd)-numUniques;
        randSynLeft(r) = length(hitInd);
        %
        % hitInd = find(makeSyn);
        % [y x] = ind2sub(size(synMat),hitInd);
        % dendVal = dendMat(hitInd);
        % for i = 1:max(dendVal);
        %     hitAx = y(dendVal==i);
        %     randUDend(r,i) = length(unique(hitAx));
        % end
        % randUniques(r) = sum(randUDend(r,:));
        %
        
        sumSyn = sumSyn + makeSyn;
        
        [randCohortSum(r) memDifs(:,:,r)] = findCohorts(makeSyn,touchMat,dendMat);
        
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
toc

% hist(randHits)
% Pmaxes = sum(randHits>=realHits)/reps;
% disp(sprintf('P = %1.5f',Pmaxes))

Puniques = sum((randUniques./randSynLeft)<=realUniques/realSynLeft)/reps;
disp(sprintf('Total uniques P = %0.05f',Puniques))

Pmults = sum((randMults)>=realMults)/reps;
disp(sprintf('Total mults P = %0.05f',Pmults))


Pcohort = sum(randCohortSum>=cohortSum)/reps;
disp(sprintf('Cohort Sum P = %0.05f',Pcohort))
hist(randCohortSum)


sortedUs = sort(randMults);
threshMults = sortedUs(round(length(sortedUs)*.95));
difMults = threshMults - realMults;
percentChangeMults = difMults/realMults*100;
meanRandMults = mean(randMults)

hist(randMults,[0:1:100])
hold on
scatter(46,1,'r')
hold off



col = uint8(synMat*1000);
col(:,:,2) = synMat*0;
col(:,:,3) = touchMat*1000;

%%Show before and after matrix
subplot(1,1,1)
col(:,:,2) = sumSyn*10000;
image(col)

onOld = sum(sum(sumSyn(synMat>0)));
onNew = sum(sum(sumSyn(synMat <1)));

newRat = onNew/(onOld+onNew)*100;
binSynReps = [0:.01:1];
synReps = (hist(sumSyn(touchMat>0)/reps,binSynReps))/sum(synMat(:));
bar(synReps)
%
% uMean = mean(randUDend,1);
% uSTD = std(randUDend,1);
% bar(uMean,'b','barWidth',.8)
% hold on
% bar(uDend,'r','barWidth',.3)
% hold off
%


image(mean(memDifs,3)*100+100)
