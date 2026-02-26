function[numUniques numMults] = findUniques(synMat,dendMat)


hitInd = find(synMat);
[y x] = ind2sub(size(synMat),hitInd);
dendVal = dendMat(hitInd);
combIDs = y* max(dendVal*2) + dendVal;
numUniques = length(unique(combIDs));
numMults = length(hitInd)-numUniques;




%         hitDend = double(dendMat).* makeSyn;
%         for i = 1:size(dendMat,1)
%             grabHits = hitDend(i,:);
%             grabHits = grabHits(grabHits>0);
%             hGrab = hist(grabHits,histDend);
%             maxHit(i) = max(hGrab);
%         end
%         randHits(r) = sum(maxHit);


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
