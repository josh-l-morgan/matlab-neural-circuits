
global glob tis


if 0
    
    SPN = [glob.datDir 'Analysis\Data\preproc\'];
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'GOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load([SPN 'COI.mat']);



    %%Load all sms
    smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
    clear sms
    roiCid = ptDat(:,3);
    runCids = unique(roiCid);%MOI.cids;
    for i = 1:length(runCids);
        cid = runCids(i);
        %%Get distances between nodes
        disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
        fileName = sprintf('sm_cid%d.mat',cid);
        useSM(i) = 1;
        %sm = load([smDir fileName],'skel2skel','nep');
        load([smDir fileName]);
        sms(i).sm = sm;
    end


end


roiCid = ptDat(:,3);
runCids = unique(roiCid);%MOI.cids;
numCid = length(runCids);

posCal = ptDat(:,[7 6 8]);
posCal(:,1) = posCal(:,1) * 0.004;
posCal(:,2) = posCal(:,2) * 0.004;
posCal(:,3) = posCal(:,3) * .04;

pol = ROI.Polarity;
snr = ROI.SNR;

smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];

allPolDif = [];
allRoiDists = [];
allSNRMins = [];
allPol1 = [];
allPol2 = [];

for i = 1:length(sms);
    
    
    %%Get distances between nodes
%     disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
    useSM(i) = 1;
    sm = sms(i).sm;
    cid = sm.cid;
    dAll = sm.skel2skel.linDist;
    eAll = sm.skel2skel.eucDist;

    pos = sm.nep.pos;
       
    %%Get Rois
    useRois = find(roiCid == cid);
    roiPos = SOI.pos(useRois,:);
    roiPol = ROI.Polarity(useRois);
    roiSN = ROI.SNR(useRois);
    closeNodes = SOI.closeNode(useRois);
    
    
    dRoi = dAll(closeNodes,closeNodes);
    eRoi = eAll(closeNodes,closeNodes);
    
    %%Find difference in polarities
    pRoi = abs(pol(useRois) - pol(useRois)');
    uRoi = 1:length(useRois) > [1:length(useRois)]';
    snrRoi = min(snr(useRois),snr(useRois)');

    allPolDif = cat(1,allPolDif,pRoi(uRoi));
    allRoiDists = cat(1,allRoiDists,dRoi(uRoi));
    allSNRMins = cat(1,allSNRMins,snrRoi(uRoi));

    polTest = repmat(pol(useRois)',[length(useRois) 1]);
    allPol1 = cat(1,allPol1,polTest(uRoi));

    polTest = repmat(pol(useRois),[1 length(useRois) ]);
    allPol2 = cat(1,allPol2,polTest(uRoi));

    scatter(dRoi(uRoi),pRoi(uRoi),'.')
    drawnow

%     for r = 1:length(closeNodes)
%         sRoi = sort(dRoi(r,:),'ascend');
%         
% 
%     end
       
    
end


pause(1)


%% show data
binD = {[0 3] [10 13] [30 33] [60 63]}
subP = length(binD) + 1;
usePair = allSNRMins> .0;


subplot(subP,1,1)
for i = 1:length(binD)
end
scatAll = scatter(allRoiDists(usePair),allPolDif(usePair),5,'o','MarkerFaceColor',[0 0 0],...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.002)
xlim([0 100])


pRange = 0:.2:3;
for i = 1:length(binD)

    pickPair = (usePair) & (allRoiDists>=binD{i}(1)) & (allRoiDists<binD{i}(2));

    subplot(subP,2,i*2+1)
    scatter(allPol1(pickPair),allPol2(pickPair),5,'o','MarkerFaceColor',[0 0 0],...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.002)
    xlim([-1 1])
    ylim([-1 1])

    subplot(subP,2,i*2+2)
    bPol = allPolDif(pickPair);
    hPol = hist(bPol,pRange);
    bar(pRange,hPol);
    title(sprintf('bin %d to %d',binD{i}(1),binD{i}(2)))

end










