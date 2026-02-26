
scatter(conPref.sharedAxNorm(1,:),conPref.sharedAxNorm(2,:)+ ...
    rand(size(conPref.sharedAxNorm(2,:)))*.01,'r')

%%
conPref

testCon = conPref.sharedAxNorm;
testCon = conPref.geoMean;
testCon = conPref.sharedSyn;

bothCon = find((testCon(1,:)>0) & (testCon(2,:)>0))
is1 = find(testCon(1,:)>0);
is2 = find(testCon(2,:)>0);
val1 = testCon(1,:)
val2 = testCon(2,:)

maxTest = max(testCon(:))
histBin = [0:maxTest/20:maxTest];
cross1 = hist(testCon(1,intersect(is1,bothCon)),histBin)
nocross1 = hist(testCon(1,setdiff(is1,bothCon)),histBin)
cross2 = hist(testCon(2,intersect(is2, bothCon)),histBin)
nocross2 = hist(testCon(2,setdiff(is2,bothCon)),histBin)


subplot(4,1,1)
bar(histBin,nocross1)
subplot(4,1,2)
bar(histBin,cross1)
subplot(4,1,3)
bar(histBin,nocross2)
subplot(4,1,4)
bar(histBin,cross2)


%%
scatter(val1,val2,'b')
hold on
scatter(val1(bothCon),val2(bothCon),'r')
hold off

%%

notSeed = ones(length(conPref.cellList),1);
for s = 1:length(conPref.seedCells)
    notSeed(find(conPref.cellList == conPref.seedCells(s))) = 0;
end

val1p = val1(notSeed>0);
val2p = val2(notSeed>0);

maxV = [max(val1p) max(val2p)]+1;
valInd = sub2ind(maxV,val1p+1,val2p+1);
uVal = unique(valInd);
hVal = hist(valInd,uVal);
[valY valX] = ind2sub(maxV,uVal);
valY = valY - 1;
valX = valX - 1;

max(hVal)

radiusHVal = sqrt(hVal/pi);
radiusHVal = radiusHVal/min(radiusHVal);

scatter(valY,valX,radiusHVal*15,'k','filled')


%% Find most linked


sumLinks = sum(testCon,1);
[sortLinks idx] = sort(sumLinks,'descend');
sortCells = conPref.cellList(idx);

mostLinked = [sortCells(1:10)' sortLinks(1:10)']




