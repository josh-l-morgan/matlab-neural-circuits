function[makeSyn] = driftSolver(makeSyn,synMat,touchMat,disruptReps)


if ~exist('disruptReps','var')
    disruptReps = 0;
end

axNum = size(synMat,1);
spNum = size(synMat,2);
aSynNum = sum(synMat,2);
sSynNum = sum(synMat,1);
synNum = sum(synMat(:));

aSyn2go = aSynNum - sum(makeSyn,2); %how many synapses left for each axon to deliver
sSyn2go = sSynNum - sum(makeSyn,1);


%% disrupt the matrix
for d = 1:disruptReps
    
    if rand>.5  %% add synapses
        
        openT = find(touchMat & ~makeSyn);
        pickT = openT(ceil(rand*length(openT)));
        [y x] = ind2sub(size(touchMat),pickT);
        makeSyn(y,x) = makeSyn(y,x) + 1;
        aSyn2go(y) = aSyn2go(y) - 1;
        sSyn2go(x) = sSyn2go(x) - 1;
        
    else
        openT = find(makeSyn);
        if ~isempty(openT)
            pickT = openT(ceil(rand*length(openT)));
            
            [y x] = ind2sub(size(touchMat),pickT);
            makeSyn(y,x) = makeSyn(y,x) - 1;
            aSyn2go(y) = aSyn2go(y) + 1;
            sSyn2go(x) = sSyn2go(x) + 1;
        end
    end
end

choiceMat = [sum(aSyn2go<0)  sum(sSyn2go<0) ; sum(aSyn2go>0)  sum(sSyn2go>0)];
fixRep = 0;
while sum(choiceMat(:))
    fixRep = fixRep + 1;
    
    choiceMat = [sum(aSyn2go<0)  sum(sSyn2go<0) ; sum(aSyn2go>0)  sum(sSyn2go>0)];
    
    if choiceMat(1,1) & ((~choiceMat(1,2) | (rand>.5)))
        badAx = find(aSyn2go<0);
        pickA = badAx(ceil(rand*length(badAx)));
        badS = find(makeSyn(pickA,:));
        pickS = badS(ceil(rand*length(badS)));
        aSyn2go(pickA) = aSyn2go(pickA) + 1;
        sSyn2go(pickS) = sSyn2go(pickS) + 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) -1;
        
    elseif  choiceMat(1,2)
        
        badS = find(sSyn2go<0);
        pickS = badS(ceil(rand*length(badS)));
        badA = find(makeSyn(:,pickS));
        pickA = badA(ceil(rand*length(badA)));
        aSyn2go(pickA) = aSyn2go(pickA) + 1;
        sSyn2go(pickS) = sSyn2go(pickS) + 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) -1;
        
    elseif choiceMat(2,1) & ((~choiceMat(2,2) | (rand>.5)))
        
        badAx = find(aSyn2go>0);
        pickA = badAx(ceil(rand*length(badAx)));
        badS = find((makeSyn(pickA,:)==0) & touchMat(pickA,:));
        pickS = badS(ceil(rand*length(badS)));
        aSyn2go(pickA) = aSyn2go(pickA) - 1;
        sSyn2go(pickS) = sSyn2go(pickS) - 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) +1;
        
    elseif choiceMat(2,2) %
        
        badS = find(sSyn2go>0);
        pickS = badS(ceil(rand*length(badS)));
        badA = find((makeSyn(:,pickS)==0) & touchMat(:,pickS));
        pickA = badA(ceil(rand*length(badA)));
        aSyn2go(pickA) = aSyn2go(pickA) - 1;
        sSyn2go(pickS) = sSyn2go(pickS) - 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) + 1;
        
    end
    
    
end %loop until matrix is fixed




