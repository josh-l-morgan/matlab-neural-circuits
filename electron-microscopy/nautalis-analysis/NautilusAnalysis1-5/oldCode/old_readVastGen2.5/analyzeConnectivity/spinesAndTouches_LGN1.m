%%Potential concern s
%%if axons are selected randomly

%% Required input
% Every axon and spine must be able to make at least one synapse

clear all
colormap gray(256)

reps = 100;  % number of randomizations
displayInterval = 100000;  % number of iterations per randomaization before displaying status
reportOn = 1; % display confirmation of proper matrix after each randomization
imageOn = 1;

%% get touch and synapse matrix
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
TPN = [MPN 'skel\'];
load([TPN 'gridArbor.mat'],'gridArbor')


synMat = synMat;
touchMat = touchMat;
spineparentids


%% get per cell synapse numbers (constraints)
aSynNum = sum(synMat,2)';
sSynNum = sum(synMat,1);



%% Run model
tic

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
        touches = find((touch2go(pickA,:)>0) & ~(makeSyn(pickA,:)));  %% look for remaining touches
        
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
        disp(sprintf('%d bad syn, %d axon mistakes, %d spine mistakes after %d cycles', badSyn,mistakeNum,badSpine,c))
        pause(.0001)
        
        tooManySyn = sum(makeSyn(:)>1);
        disp(sprintf('Too many syn = %d',tooManySyn))
    end
    
    if imageOn
        
        col = uint8(synMat*1000);
        col(:,:,2) = synMat*0;
        col(:,:,3) = touchMat*1000;
        
        %%Show before and after matrix
        subplot(1,1,1)
        col(:,:,2) = makeSyn*1000;
        image(col)
        
        %%Show synapse distribution history
%         subplot(2,1,2)
%         plot(trackSyn)
%         pause(.0001)
        
    end
    
end
toc