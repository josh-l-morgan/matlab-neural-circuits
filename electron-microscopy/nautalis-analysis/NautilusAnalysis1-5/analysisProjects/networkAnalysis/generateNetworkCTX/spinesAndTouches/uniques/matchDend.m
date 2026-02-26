function[sumSameSyn] = matchDend(synMat,dendMat)

binAx = 1; %Binarize touch and synapse data
%% Reformat into axon by dend touch and synapse matrix
dendNum = max(dendMat(:));
touchDend = zeros(size(dendMat,1),dendNum);
synDend = touchDend*0;

for i = 1:dendNum;
    dendFilt = dendMat == i;
    synDend(:,i) = sum(synMat.*dendFilt,2);
end

if binAx
    synDend = synDend>0;
end

%% Count similarities
sumSameSyn = zeros(dendNum,dendNum);
for i = 1:dendNum
    sameSyn = min(synDend, repmat(synDend(:,i),[1 dendNum]));
    sumSameSyn(i,:) = sum(sameSyn,1);
end

