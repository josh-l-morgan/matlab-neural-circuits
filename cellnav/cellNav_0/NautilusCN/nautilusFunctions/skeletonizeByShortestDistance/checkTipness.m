function[tipShape] = checkTipness(surfVox,longTip);

%%
checkVox = find(longTip.countVox(end,:)>(size(longTip.countVox,1)));
tipLength = zeros(size(longTip.vox));
tipRat = tipLength;

for i = 1:length(checkVox)
   
    childVox = find(longTip.vox == checkVox(i));
    childDists = surfVox.seedPath.dist(childVox);
    selfDist = surfVox.seedPath.dist(checkVox(i));
    
    tipLength(checkVox(i)) = selfDist-mean(childDists);
    tipRat(checkVox(i)) = tipLength(i)/length(childDists);
   
end

tipShape.tipLength = tipLength;
tipShape.tipRat = tipRat;

