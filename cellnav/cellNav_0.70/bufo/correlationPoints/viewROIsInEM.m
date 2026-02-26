
clear all

SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasksSpread\'

%load([SPN 'maskDat.mat'])
load([SPN 'ptDat.mat']);
load([SPN 'ROI.mat']);

pos = ptDat(:,6:8);
pol1 = ROI.Polarity;




pol2 = mean(ROI.Polarity2,2);
onMean = mean(ROI.ON,2);
offMean = mean(ROI.OFF,2);
onMax = max(ROI.ON,[],2);
offMax = max(ROI.OFF,[],2);
pol3 = onMean ./ (onMean + offMean);


polRange = [-1.5:.1:1.5];
subplot(4,1,1)
hist(pol1,polRange)
xlim([-1.5 1.5])
subplot(4,1,2)
hist(pol2,polRange)
xlim([-1.5 1.5])

subplot(4,1,3)
hist(onMean)
subplot(4,1,4)
hist(onMax)



%% Show positions
showPol1 = (pol1 +1) * 50;
showPol2 = (pol2 + 1) * 50;
showOnMean = onMean * 50;
showOnMax = onMax * 50;

showOffMean = offMean * 50;
showOffMax = offMax * 50;
prop = showOnMean;

clf
subplot(1,1,1)
allScat = scatter3(pos(:,1),pos(:,2),pos(:,3),100,'.','k')
ax = gca
ax.Color = [0 0 0];
% prop = prop - min(prop);
% prop = prop * 100/max(prop);
colormap jet(100)

allScat.CData = prop;






