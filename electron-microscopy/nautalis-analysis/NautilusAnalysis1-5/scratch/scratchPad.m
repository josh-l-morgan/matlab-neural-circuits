MPN = GetMyDir
load([MPN 'obI.mat'])
synapses = obI.nameProps.edges(:,1:2)


%%
giants = [1012 1013]

is108 = synapses(:,1)==108;

pre108 = synapses(is108,2);

upre108 = unique(pre108)
upre108 = upre108((upre108>1000) & (upre108<2000))


for i = 1:length(upre108)
    with108(i) = sum((synapses(:,2) == upre108(i)) & is108);
    withAny(i) = sum(synapses(:,2) == upre108(i));
end

scatter(with108,withAny);



%%

showCon = con;

targ = 273
showCon(:,find(cellList == targ)) = showCon(:,find(cellList == targ)) +1;

image(showCon * 10)

%% 

prefList = seedPref.cellPrefNoExclusion(2,:)
con = seedPref.con;

[sortPref idx] = sort(prefList,'descend');

sortCells = seedPref.cellList(idx);

uint16([sortCells(1:40)' sortPref(1:40)'])


%% get post
targ  = 2032;
isTarg = synapses(:,2) == targ;
isPost = synapses(isTarg,1)

%% Show new preferences
% 
% MPN = GetMyDir
% load([MPN 'obI.mat'])
useList = obI2cellList_seedInput(obI,[108 201])
conPref = seedPreferences([108 201],useList);

subplot(2,2,1)
scatter(conPref.sharedAxNorm(1,:) + rand(size(conPref.sharedAxNorm(2,:)))*.01 , ...
    conPref.sharedAxNorm(2,:) + rand(size(conPref.sharedAxNorm(2,:)))*.01,'k','.')
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 201')
title('Shared Axons')
hold on

xlim([-0.01  .41]);
ylim([-0.01  .41]);

subplot(2,2,2)
scatter(conPref.sharedSynNorm(1,:)  + rand(size(conPref.sharedAxNorm(2,:)))*.003 , ...
    conPref.sharedSynNorm(2,:) + rand(size(conPref.sharedSynNorm(2,:)))*.003,'k','.')
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 201')
title('Synapses from Axons linked to seed')

xlim([-0.01  .41]);
ylim([-0.01  .41]);


subplot(2,2,3)
scatter(conPref.geoMeanNorm(1,:) + rand(size(conPref.sharedAxNorm(2,:)))*.003 , ...
    conPref.geoMeanNorm(2,:) + rand(size(conPref.geoMeanNorm(2,:)))*.003,'k','.')
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 001')
title('geometric mean - sqrt(S1 * S2)')

hold off

xlim([-0.01  .41]);
ylim([-0.01  .41]);


subplot(2,2,4)
scatter(conPref.zcovNorm(1,:) + rand(size(conPref.zcovNorm(2,:)))*.0003 , ...
    conPref.zcovNorm(2,:) + rand(size(conPref.zcovNorm(2,:)))*.0003,'k','.')
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 001')
title('covariance without mean subtraction')

hold off

xlim([-0.01  .2]);
ylim([-0.01  .2]);




% 
% subplot(2,2,4)
% scatter(conPref.seedCov(1,:) + rand(size(conPref.seedCov(2,:)))*.003 , ...
%     conPref.seedCov(2,:) + rand(size(conPref.seedCov(2,:)))*.003,'k','.')
% xlabel('similartiy (0-1) to seed cell 108')
% ylabel('similarity (0-1) to seed cell 001')
% title('covariance')
% hold on
% minCov = min(conPref.seedCov,[],2);
% maxCov = max(conPref.seedCov,[],2);
% 
% plot([0 0 ],[minCov(1) maxCov(1)])
% plot([minCov(2) maxCov(2)],[0 0 ])
% 
% 
% hold off
% % 
% % xlim([-0.01  .41]);
% % ylim([-0.01  .41]);

%%
%% Show new preferences
% 
% MPN = GetMyDir
% load([MPN 'obI.mat'])
clf
useList = obI2cellList_seedInput(obI,[108 201])
conPref = seedPreferences([108 201],useList);
cellNum = length(conPref.cellList);

minAx = min(conPref.sharedAx,[],1);
PercentCrossover = sum(minAx>0)/length(minAx) * 100;
disp(sprintf('Percent Crossover %.1f%%',PercentCrossover))

sumAx = sum(conPref.sharedAx,1);




subplot(1,4,1)
minVal = min(conPref.sharedAxNorm);
histBin = [0:max(minVal)/30:max(minVal)];
hist(minVal,histBin)
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 201')
title('Shared Axons')
hold on

ylim([-0.01  cellNum]);
%xlim([-0.01  1]);

subplot(1,4,2)
minVal = min(conPref.sharedSynNorm);
histBin = [0:max(minVal)/30:max(minVal)];
hist(minVal,histBin)
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 201')
title('Synapses from Axons linked to seed')

ylim([-0.01  cellNum]);
%xlim([-0.01  1]);


subplot(1,4,3)
minVal = min(conPref.geoMeanNorm);
histBin = [0:max(minVal)/30:max(minVal)];
hist(minVal,histBin)
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 001')
title('geometric mean - sqrt(S1 * S2)')

hold off

ylim([-0.01  cellNum]);
%xlim([-0.01  1]);


subplot(1,4,4)
minVal = min(conPref.zcovNorm);
histBin = [0:max(minVal)/30:max(minVal)];
hist(minVal,histBin)
xlabel('similartiy (0-1) to seed cell 108')
ylabel('similarity (0-1) to seed cell 001')
title('covariance without mean subtraction')

hold off

ylim([-0.01  cellNum]);
%xlim([-0.01  1]);

%%
clf

minAx = min(conPref.sharedAx,[],1);
sumAx = sum(conPref.sharedAx,1);
multSums = sumAx(minAx>0);
histBin = [1:max(sumAx)];
histMin = hist(multSums,histBin);
histSum = hist(sumAx,histBin);
bar(histBin,histSum,'b')
hold on
bar(histBin,histMin,'r')
hold off


multNonSeed = sum(sumAx>1)-2;
crossoverNum = sum(minAx>0);
PercentCrossoverOfMults = crossoverNum/multNonSeed * 100

PercentCrossover = sum(minAx>0)/length(minAx) * 100;
disp(sprintf('Percent Crossover %.1f%%',PercentCrossover))




%%
ax1 = conPref.sharedAx(1,:)+1;
ax2 = conPref.sharedAx(2,:)+1;
fSize = [max(ax1) max(ax2)]
axInd = sub2ind(fSize,ax1,ax2);
uax = unique(axInd);
histAx = hist(axInd,uax);
axField = zeros(fSize);
axField(uax) = histAx;

image(axField*10);






