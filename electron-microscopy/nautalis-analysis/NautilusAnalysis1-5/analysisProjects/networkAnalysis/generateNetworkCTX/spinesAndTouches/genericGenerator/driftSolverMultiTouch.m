function[makeSyn] = driftSolverMultiTouch(makeSyn,synMat,touchMat,disruptReps)


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
    
    
    if rand>.5 | ~sum(makeSyn(:)) %% add synapses
        lessMat = touchMat - makeSyn;
        openT = find(lessMat>0);
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
    
    if choiceMat(1,1) & ((~choiceMat(1,2) | (rand>.5))) %subtract synapse by picking axon
        
        pickA = pickWeighted(-aSyn2go);
        pickS = pickWeighted(makeSyn(pickA,:));
       
        aSyn2go(pickA) = aSyn2go(pickA) + 1;
        sSyn2go(pickS) = sSyn2go(pickS) + 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) -1;
        
    elseif  choiceMat(1,2) %% subtract synapse by picking spine
        
        pickS = pickWeighted(-sSyn2go);
        pickA = pickWeighted(makeSyn(:,pickS));
        
        aSyn2go(pickA) = aSyn2go(pickA) + 1;
        sSyn2go(pickS) = sSyn2go(pickS) + 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) -1;
        
    elseif choiceMat(2,1) & ((~choiceMat(2,2) | (rand>.5)))
        
        pickA = pickWeighted(aSyn2go);
        pickS = pickWeighted(touchMat(pickA,:) - makeSyn(pickA,:));
        
        aSyn2go(pickA) = aSyn2go(pickA) - 1;
        sSyn2go(pickS) = sSyn2go(pickS) - 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) +1;
        
    elseif choiceMat(2,2) %
        
        pickS = pickWeighted(sSyn2go);
        pickA = pickWeighted(touchMat(:,pickS) - makeSyn(:,pickS));
        
        aSyn2go(pickA) = aSyn2go(pickA) - 1;
        sSyn2go(pickS) = sSyn2go(pickS) - 1;
        makeSyn(pickA,pickS) = makeSyn(pickA,pickS) + 1;
        
    end
    
    
end %loop until matrix is fixed




