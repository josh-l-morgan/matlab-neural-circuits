global glob

%% List scripts for getting correlation

if 0
   glob.datDir = 'Y:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\';
end

if 0
    runCellNav %% if no COI exists
end

if 1
    makeMaskDatJm
    tweakMasks
    %SignalExtraction_FlashBar_CorrER_121721
    SignalExtraction_FlashBar_CorrER_091222
    SignalExtraction_jmAutoCorr
    %jmSignalExtraction_FlashBar_CorrER_100921
    makeCOI %% if no COI exists
    % makeMOI %% to get data from neuron model
    makeSOI %% makes SOI and NOI
    makeGOI %% Groups ROIs that are close together
    makePOI
    jmFindLengthConstant14
end
if 0
    showPolaritiesAndSyn
end

return

%SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasksSpread\'

SPN =  [glob.datDir 'Analysis\Data\preproc\'];

load([SPN 'ptDat.mat']);
load([SPN 'ROI.mat']);
load([SPN 'ROI2.mat']);
load([SPN 'SOI.mat']);
load([SPN 'NOI.mat']);
load([SPN 'MOI.mat']);



%%Get roiCids
roiCids = ptDat(:,3);
uCids = unique(roiCids);
lookup = zeros(max(uCids),1);
lookup(uCids) = 1:length(uCids);
roiCidIDs = lookup(roiCids);

emPos = ptDat(:,[6 7 8]);

useCids = uCids;
goodRoiCid = [];
for i = 1:length(useCids)
    goodRoiCid = [goodRoiCid; find(roiCids==useCids(i))];
end

useROI = intersect(find(SOI.closeNode>0), goodRoiCid) ;


%%Get roiExp
roiExp = ptDat(:,1);
uExp = unique(roiExp);
lookup = zeros(max(uExp),1);
lookup(uExp) = 1:length(uExp);
expIDs = lookup(roiExp);


%%Get zColors
z = SOI.pos(:,3);
showZ = z - min(z);
showZ = round(showZ * 99/max(showZ) + 1);

sPol = (SOI.off-SOI.on)./(SOI.off+SOI.on);

colormap jet(100)
% 
% roiON = ROI2.maxON(useROI);
% roiOFF = ROI2.maxOFF(useROI);
%roiPol = ROI2.Polarity5(useROI);
roiPolH = ROI.Polarity(useROI);
roiPolH2 = ROI.Polarity2(useROI);

%%SNR
SNR = ROI.SNR(useROI);
SNRc = SNR - min(SNR);
SNRc = floor(round(SNRc * 100/max(SNRc)))+1 ;
SNRc(SNRc>10) = 10;

scatSNR = scatter(SNR,roiPolH,'filled')
scatSNR.CData = SNRc;



%% Compare polarity measures
if 0
    clf
    scatPolPol = scatter(ROI.Polarity(useROI),ROI2.Polarity(useROI));
    scatPolPol.CData = showZ(useROI);
    ylim([0 1]);
    xlim([0 1])
    hold on
    plot([0 1],[0 1])
    hold off
    
    %corFig = figure;
    subplot(3,1,1)
    scatter(SOI.on(useROI),roiON,'.')
    subplot(3,1,2)
    scatter(SOI.off(useROI),roiOFF,'.')
    subplot(3,1,3)
end

%% Compare calcium polarity to predicted polarity
subplot(3,1,1)
scatPol = scatter(SOI.offBias(useROI),roiPolH,'o','filled');
scatPol.CData = roiCidIDs(useROI) * ceil(100/length(uExp));
subplot(3,1,2)
scatPol2 = scatter(SOI.offBias(useROI),roiPolH,'o','filled');
scatPol2.CData = roiCidIDs(useROI) * ceil(100/length(uExp));
subplot(3,1,3)
scatPol3 = scatter(SOI.offBias(useROI),roiPolH2,'o','filled');
scatPol3.CData = roiCidIDs(useROI) * ceil(100/length(uExp));
%scatPol3.CData = SNRc;


%%
clf
plot([-1 1],[-1 1],'r')
hold on
plot([-1 1],[ 0 0],'k')
plot([0 0],[ -1 1],'k')

xDat = sPol(useROI);
xDat = xDat(:);
yDat = roiPolH;
yDat = yDat(:);

scatPol4 = scatter(xDat,yDat,'o','filled');

scatPol4.CData = roiCidIDs(useROI) * ceil(100/length(uExp));


xlim([-1 1])
ylim([-1 1])

if 1 % show point information, run, click, press enter
    
    [x y] = ginput;
    dists = sqrt((yDat-y(end)).^2 + (xDat-x(end)).^2);
    targ = find(dists==min(dists),1);
    rTarg = useROI(targ);
    SOI.pos(rTarg,:)
    sPol(rTarg)
    roiPolH(targ)
    
    ptEmPos = emPos(rTarg,:);
    ptCid = roiCids(rTarg);
    
    cidTarg = find(SOI.vgcCids == ptCid);
    skelPos = SOI.cell(cidTarg).skelPos;
    scatter3(skelPos(:,1),skelPos(:,2),skelPos(:,3),'.','k')
    hold on
    cN = SOI.closeNode(rTarg);
    scatter3(skelPos(cN,1),skelPos(cN,2),skelPos(cN,3),200,'r','filled')
    hold off
    
    
    
end
hold off





%% Plot by Z
clf
subplot(3,1,1)
scatZ0 = scatter(SOI.pos(useROI,3),roiPolH,'o','filled','r');
scatZ0.CData = expIDs(useROI) * ceil(100/length(uExp));

ylim([-1 1])

subplot(3,1,2)
scatZ = scatter(SOI.pos(useROI,3),SOI.offBias(useROI),'o','filled','r');
scatZ.CData = expIDs(useROI) * ceil(100/length(uExp));
ylim([-1 1])

subplot(3,1,3)
scatZ2 = scatter(SOI.pos(useROI,3),roiPolH2,'o','filled','g');
scatZ2.CData = expIDs(useROI) * ceil(100/length(uExp));
ylim([-1 1])

%% where did the original distribution go?
clf
scatPol2 = scatter(SOI.offBias(useROI),roiPolH,10,'o','filled');
scatPol2.CData = expIDs(useROI) * ceil(100/length(uExp));
set(gca,'color',[0 0 0]);





