%function[sharedAx sumSameSyn] = matchDend(synMat,dendMat)


parentDend = dendMat(1,:);
dendNum = max(parentDend);

catSyn = [];
for i = 1:max(synMat(:))
    catSyn = cat(1,catSyn,find(synMat==1));
end

[y x] = ind2sub(size(synMat),catSyn);
d = parentDend(x);
newInd = sub2ind([size(synMat,1) dendNum],y,d');








tic
binAx = 1; %Binarize touch and synapse data
%% Reformat into axon by dend touch and synapse matrix
dendNum = max(dendMat(:));
touchDend = zeros(size(dendMat,1),dendNum);
synDend = touchDend*0;
toc
tic



for i = 1:dendNum;
    dendFilt = dendMat == i;
    synDend(:,i) = sum(synMat.*dendFilt,2);
end

if binAx
    synDend = synDend>0;
end
toc
tic
%% Count similarities
sumSameSyn = zeros(dendNum,dendNum);
for i = 1:dendNum
    sameSyn = min(synDend, repmat(synDend(:,i),[1 dendNum]));
    sumSameSyn(i,:) = sum(sameSyn,1);
end

selfs = sub2ind(size(sumSameSyn),[1:dendNum],[1:dendNum]);
sumSameSyn(selfs) = 0; 
sharedAx = sum(sumSameSyn(:))/2;
toc
