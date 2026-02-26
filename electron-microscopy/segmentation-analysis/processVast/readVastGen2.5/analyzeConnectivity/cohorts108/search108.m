

refPost = 3;
refSynMat = synMat(:,3);
[refSynMat idx] = sort(refSynMat,'descend');
allBut = find(1:size(synMat,2)~=refPost);
subSynMat = synMat(idx,allBut);



%%

sharedMatRef = getMinSharedPre(refSynMat);

sharedMatSub = getMinSharedPre(subSynMat);

subplot(2,1,1)
image(sharedMatRef * 10)
subplot(2,1,2)
image(sum(sharedMatSub,3) * 10)





