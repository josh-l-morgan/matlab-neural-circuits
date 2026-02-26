function[cohortSum,sumSameSyn] = findCohorts(synMat,touchMat,dendMat)

binAx = 1; %Binarize touch and synapse data
%% Reformat into axon by dend touch and synapse matrix
dendNum = max(dendMat(:));
touchDend = zeros(size(dendMat,1),dendNum);
synDend = touchDend*0;
for i = 1:dendNum;
    dendFilt = dendMat == i;
    touchDend(:,i) = sum(touchMat.*dendFilt,2);
    synDend(:,i) = sum(synMat.*dendFilt,2);
end


if binAx
    touchDend = touchDend>0;
    synDend = synDend>0;
end

sumSyn = sum(synDend,1);
sumTouch = sum(touchDend,1);
fracDend = sumSyn./sumTouch;
crossDend = repmat(fracDend,[length(fracDend) 1]) .* repmat(fracDend',[1 length(fracDend)]);

%% Count similarities
sumSameTouch = zeros(size(touchDend,2));
sumSameSyn = sumSameTouch;
for i = 1:dendNum
    
    sameTouch = min(touchDend, repmat(touchDend(:,i),[1 dendNum]));
    sameSyn = min(synDend, repmat(synDend(:,i),[1 dendNum]));
    
    sumSameTouch(i,:) = sum(sameTouch,1);
    sumSameSyn(i,:) = sum(sameSyn,1);
end

fracSame = sumSameSyn./sumSameTouch;
predSame = sumSameTouch.*crossDend; %predicted number of shared synapses
synDifs = sumSameSyn-predSame;

selfs = sub2ind(size(synDifs),[1:dendNum],[1:dendNum]);
synDifs(selfs) = 0;
countDifs = floor(abs(synDifs));
cohortSum = sum(countDifs(:));
disp(sumSameSyn(selfs(1:10)))
%
% ratDifs = synDifs./sumSameTouch;
%
% hist(ratDifs(sumSameTouch>0))

%%
%
% image(synDifs*100+100)
%




