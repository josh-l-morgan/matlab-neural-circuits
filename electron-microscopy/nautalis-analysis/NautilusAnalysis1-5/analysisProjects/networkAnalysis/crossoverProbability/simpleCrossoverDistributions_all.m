


%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])
seedCells = [108 201 109 903 907];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedCells);
conPref = seedPreferences(seedCells,useList);


sharedSyn = conPref.sharedSyn;

isA = sharedSyn(1,:)>0;
isB = sharedSyn(2,:)>0;
isBig = sum(sharedSyn([3 4 5],:),1)>0;
isBigMult = sum(sharedSyn([3 4 5],:),1)>1;


sum(isBig)
sum(isA & isBig)
sum(isB & isBig)





scatter(conPref.sharedAxNorm(1,:),conPref.sharedAxNorm(2,:)+ ...
    rand(size(conPref.sharedAxNorm(2,:)))*.01,'r')
hold on

scatter(conPref.sharedSynNorm(1,:),conPref.sharedSynNorm(2,:)+ ...
    rand(size(conPref.sharedSynNorm(2,:)))*.01,'g')

scatter(conPref.geoMeanNorm(1,:),conPref.geoMeanNorm(2,:)+ ...
    rand(size(conPref.geoMeanNorm(2,:)))*.01,'b')
hold off


ax1 = conPref.sharedAx(1,:)+1;
ax2 = conPref.sharedAx(2,:)+1;
fSize = [max(ax1) max(ax2)]
axInd = sub2ind(fSize,ax1,ax2);
uax = unique(axInd);
histAx = hist(axInd,uax);
axField = zeros(fSize);
axField(uax) = histAx;
subplot(1,2,1)
image(axField*3);


%%

sum1 = sum(axField,1);
sum1 = sum1/max(sum1);
sum2 = sum(axField,2);
sum2 = sum2/max(sum2);
for y = 1:length(sum1)
    for x = 1:length(sum2)
        crossSum(y,x) = sum1(y)* sum2(x);
    end
end
subplot(1,2,2)
image(crossSum*256)

%% 

sharedAx = conPref.sharedAx;
sum(sharedAx(1,:)>0)
sum(sharedAx(2,:)>0)
sum(sum(sharedAx,1)>0)
sum((sharedAx(1,:)>0) & (sharedAx(2,:)>0))


