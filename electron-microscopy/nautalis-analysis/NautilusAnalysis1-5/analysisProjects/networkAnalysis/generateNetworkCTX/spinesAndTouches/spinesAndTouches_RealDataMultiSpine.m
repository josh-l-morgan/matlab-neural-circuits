%why dont we just scramble the identity of the spines

clear all
colormap gray(256)

reps = 100;  % number of randomizations
displayInterval = 10000;  % number of iterations per randomaization before displaying status
reportOn = 1; % display confirmation of proper matrix after each randomization
imageOn = 1;

%% get touch and synapse matrix
load('singlespinematrices2.mat')
synMat = singlespinesynapsematrix>0;
touchMat = singlespinetouchmatrix>0;


%% get per cell synapse numbers (constraints)
aSynNum = sum(synMat,2)';
sSynNum = sum(synMat,1);

col = uint8(synMat*1000);
col(:,:,2) = synMat*0;
col(:,:,3) = touchMat*1000;

%% Run model
tic

for r = 1:reps
    
    aSyn2go = aSynNum; %how many synapses left for each axon to deliver
    sSyn2go = sSynNum; % how many synapse left for each spine
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
            pickSpine = touches(randi(length(touches),1)); %pick spine to synapse with
            
            
            touch2go(:,pickSpine) = touch2go(:,pickSpine) - 1; % reduce available synapses for that spine
            
            
            makeSyn(pickA,pickSpine) = makeSyn(pickA,pickSpine) +1; % mark the new synapse
            
            aSyn2go(pickA) = aSyn2go(pickA)-1; % reduce synapses left to distribute
        else
            oldTouches = find(touchMat(pickA,:)); % find all touches (occupied)
            pickSpine =  oldTouches(fix(randi(length(oldTouches),1)));  %pick occupied spine to steal
            
            %% whipe old
            oldAx = find(makeSyn(:,pickSpine)); % find the axon being stolen from
            
            makeSyn(oldAx,pickSpine) = makeSyn(oldAx,pickSpine) -1; % clear the old synapse
            aSyn2go(oldAx) = aSyn2go(oldAx)+1; % ad synapse to distribution list for axon that was stolen from
            
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
            col(:,:,2) = makeSyn;
            image(col);
            pause(.01)
        end
    end
    
    
    %% Report
    if reportOn
        
        checkTouch = makeSyn - touchMat;
        badSyn = sum(checkTouch(:)>0);
        remainingSyn= sum(aSyn2go);
        randSynNum = sum(makeSyn,2)';
        mistakeNum = sum(abs(aSynNum-randSynNum));
        randSpineNum = sum(makeSyn,1);
        badSpine = sum(randSpineNum~=1);
        disp(sprintf('%d bad syn, %d axon mistakes, %d spine mistakes after %d cycles', badSyn,mistakeNum,badSpine,c))
        pause(.0001)
    end
    
    if imageOn
        
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