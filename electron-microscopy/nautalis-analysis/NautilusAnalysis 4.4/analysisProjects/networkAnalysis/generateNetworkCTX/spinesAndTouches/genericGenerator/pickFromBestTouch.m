function[makeSyn] = pickFromTouch(synMat,touchMat)

displayInterval = 1000;

axNum = size(synMat,1);
spNum = size(synMat,2);
aSynNum = sum(synMat,2);
sSynNum = sum(synMat,1);
synNum = sum(synMat(:));

aSyn2go = aSynNum; %how many synapses left for each axon to deliver
sSyn2go = sSynNum;
makeSyn = touchMat * 0; %new synaptic matrix


c = 0; %start counter
while 1
    c = c+1;
    
    
    spMatLeft = repmat(sSyn2go>0,[axNum,1])& touchMat;
    axMatLeft = repmat(aSyn2go>0,[1,spNum])  & touchMat;
    
    %%Pick needed for both
    openT = find(spMatLeft & axMatLeft);
    
    %%Pick needed for one
    if isempty(openT)
        openT = find(spMatLeft | axMatLeft);
    end
    
    %%Pick needed for none
    if isempty(openT)
        openT = find(touchMat & ~makeSyn);
    end
    
    pickT = openT(ceil(rand*length(openT)));
    [y x] = ind2sub(size(touchMat),pickT);
    
    
    if aSyn2go(y) & sSyn2go(x)
        makeSyn(y,x) = makeSyn(y,x) + 1;
        aSyn2go(y) = aSyn2go(y) - 1;
        sSyn2go(x) = sSyn2go(x) - 1;
        
    elseif ~aSyn2go(y) & sSyn2go(x) %steal
        
        oldSp = find(makeSyn(y,:)); % find the axon being stolen from
        pickoldSp= oldSp(ceil(rand * length(oldSp)));
        
        makeSyn(y,pickoldSp) = makeSyn(y,pickoldSp) -1; % clear the old synapse
        aSyn2go(y) = aSyn2go(y) + 1;
        sSyn2go(pickoldSp) = sSyn2go(pickoldSp) + 1;
        
        %% add new
        makeSyn(y,x) = makeSyn(y,x) + 1; % mark new synapse
        aSyn2go(y) = aSyn2go(y) - 1;
        sSyn2go(x) = sSyn2go(x) - 1;
        
    elseif ~sSyn2go(x) & aSyn2go(y)% spine conflict
        oldAx = find(makeSyn(:,x)); % find the axon being stolen from
        pickOldAx = oldAx(ceil(rand * length(oldAx)));
        
        makeSyn(pickOldAx,x) = makeSyn(pickOldAx,x) -1; % clear the old synapse
        aSyn2go(pickOldAx) = aSyn2go(pickOldAx) + 1;
        sSyn2go(x) = sSyn2go(x) + 1;
        
        
        %% add new
        makeSyn(y,x) = makeSyn(y,x) + 1; % mark new synapse
        aSyn2go(y) = aSyn2go(y) - 1;
        sSyn2go(x) = sSyn2go(x) - 1;
        openTouch(pickT) = 0;
        
    elseif ~aSyn2go(y) & ~sSyn2go(x)
        
        oldAx = find(makeSyn(:,x)); % find the axon being stolen from
        pickOldAx = oldAx(ceil(rand * length(oldAx)));
        
        makeSyn(pickOldAx,x) = makeSyn(pickOldAx,x) -1; % clear the old synapse
        aSyn2go(pickOldAx) = aSyn2go(pickOldAx) + 1;
        sSyn2go(x) = sSyn2go(x) + 1;
        openTouch(find(touchInd == sub2ind(size(touchMat),pickOldAx,x))) = 1;
        
        oldSp = find(makeSyn(y,:)); % find the axon being stolen from
        pickoldSp= oldSp(ceil(rand * length(oldSp)));
        
        makeSyn(y,pickoldSp) = makeSyn(y,pickoldSp) -1; % clear the old synapse
        aSyn2go(y) = aSyn2go(y) + 1;
        sSyn2go(pickoldSp) = sSyn2go(pickoldSp) + 1;
        openTouch(find(touchInd == sub2ind(size(touchMat),y,pickoldSp))) = 1;
        
        
        %% add new
        makeSyn(y,x) = makeSyn(y,x) + 1; % mark new synapse
        aSyn2go(y) = aSyn2go(y) - 1;
        sSyn2go(x) = sSyn2go(x) - 1;
        openTouch(pickT) = 0;
        
        
    else
        'something is wrong'
    end
    
    if ~sum(abs(aSyn2go)) & ~sum(abs(sSyn2go)) % if no synapses left to place exit loop
        break
    end
    
    if ~mod(c,displayInterval) %disply process state
        disp(sprintf('%d syn remaining on cycle %d.',synNum-sum(makeSyn(:)),c));
        pause(.0001)
    end
end



