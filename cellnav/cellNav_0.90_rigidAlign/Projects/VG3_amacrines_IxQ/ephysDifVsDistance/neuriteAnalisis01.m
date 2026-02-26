
%%Compair distance between points with difference in physiology

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
        fileName = sprintf('smx_cid%d.mat',cid);
        useSM(i) = 1;
        sm = load([smDir fileName],'skel2skel','nep');
        sms(i) = sm;
    end


end
jCol = jet(100);

roiCid = ptDat(:,3);
runCids = unique(roiCid);%MOI.cids;
numCid = length(runCids);

posCal = ptDat(:,[7 6 8]);
posCal(:,1) = posCal(:,1) * 0.004;
posCal(:,2) = posCal(:,2) * 0.004;
posCal(:,3) = posCal(:,3) * .04;

pol = ROI.Polarity;
snr = ROI.SNR;

%Get experiment list
exp = ptDat(:,1);
uExp = unique(exp);
lookup(uExp) = [1 4 7 10 2 5 8 11 3 6 9 12];
lookup(uExp) = [1 5 9 2 6 10 3 7 11 4 8 12];
exIds = lookup(exp);


%% show all polarities
clf

polInd = pol;
polInd = round((polInd + 1) * 100);
polInd(polInd<1) = 1;
polInd(polInd>100) = 100;
showPol = ~isnan(polInd);
polCol = jCol(polInd(showPol),:);
blankPol = isnan(polInd);


clf
jCol = jet(100);
snrInd = snr;
snrInd = round(snrInd * 100/0.05);
snrInd(snrInd<1) = 1;
snrInd(snrInd>100) = 100;
showPol = ~isnan(snrInd);
snrCol = jCol(snrInd(showPol),:);
blankSNR = isnan(snrInd);


expCol = jCol(ceil(exIds(showPol)*99/max(exIds)),:);

scatP = scatter3(posCal(showPol,1),posCal(showPol,2),posCal(showPol,3),'o','filled');
scatP.CData = expCol;
hold on
scatter3(posCal(blankPol,1),posCal(blankPol,2),posCal(blankPol,3),'o','k','filled')
hold off
pause(1)

%% Run each cell
allPolDif = [];
allRoiDists = [];
allSNRMins = [];
allExp1 = [];
allExp2 = [];
allPol1 = [];
allPol2 = [];
allCids = [];

for i = 1:length(runCids);
    
    cid = runCids(i);
    
    %%Get distances between nodes
    disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
    useSM(i) = 1;
    sm = sms(i);
    dAll = sm.skel2skel.linDist;
    eAll = sm.skel2skel.eucDist;

    pos = sm.nep.pos;
       
    %%Get Rois
    useRois = find(roiCid == cid);
    roiPos = SOI.pos(useRois,:);
    roiPol = ROI.Polarity(useRois);
    roiSN = ROI.SNR(useRois);
    closeNodes = SOI.closeNode(useRois);
    
    
    dRoi = dAll(closeNodes,:);
    dRoi = dRoi(:,closeNodes);
    eRoi = eAll(closeNodes,:);
    eRoi = eRoi(:,closeNodes);

    pos = sm.nep.pos;

    %%Test Mapping
    if 0
        clf
        for ry = 1:size(dRoi,1)
            for rx = 1:size(dRoi,2)


                if (dRoi(ry,rx)>0) & (dRoi(ry,rx)< 3)

                    scatter(pos(:,1),pos(:,2),'.')
                    hold on
                    scatter(posCal(useRois(ry),1),posCal(useRois(ry),2),2,'o','g','linewidth',3)
                    scatter(posCal(useRois(rx),1),posCal(useRois(rx),2),2,'o','r','linewidth',3)
                    scatter(pos(closeNodes(ry),1),pos(closeNodes(ry),2),100,'d','g','linewidth',2)
                    scatter(pos(closeNodes(rx),1),pos(closeNodes(rx),2),100,'d','r','linewidth',2)
                    title(sprintf('distance = %0.1f',dRoi(ry,rx)))
                    hold off
                    drawnow

                    pause
                end
            end
        end
    end


    if 1 % show roi polarities
            jCol = jet(100);
            polInd = pol(useRois);
            polInd = round((polInd + 1) * 100);
            polInd(polInd<1) = 1;
            polInd(polInd>100) = 100;
            showPol = ~isnan(polInd);
            polCol = jCol(polInd(showPol),:);
            
            scatter(pos(:,1),pos(:,3),'.','k')
            hold on
            scatR = scatter(posCal(useRois(showPol),1),posCal(useRois(showPol),3),200,'o','g','markerfacecolor','flat');
            scatR.CData = polCol;
            hold off
            pause(1)

    end

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

    polTest = repmat(exIds(useRois)',[1 length(useRois) ]);
    allExp1 = cat(1,allExp1,polTest(uRoi));
    polTest = repmat(exIds(useRois),[length(useRois) 1]);
    allExp2 = cat(1,allExp2,polTest(uRoi));

    allCids = cat(1,allCids,ones(sum(uRoi(:)),1) * i);

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


%% Figure out close points
clf
pickPair = (usePair) & (allRoiDists>=binD{1}(1)) & (allRoiDists<binD{1}(2));
cidCol = [1 0 0; 0 1 0; 0 0 1; .7 .7 0; .7 0 .7; 0 .7 .7 ;1 .5 0; .5 1 0; 0 1 .5; 0 .5 1; 1 0 .5; .5 0 1];
useCol = cidCol(allCids(pickPair>0),:);

jCol = jet(100);
snrInd = allSNRMins(pickPair);
snrInd = snrInd-min(snrInd);
snrInd = round(snrInd * 299/max(snrInd(:))) + 1;
snrInd(isnan(snrInd)) = 1;
snrInd(snrInd<1) = 1;
snrInd(snrInd>100) = 100;
snrCol = jCol(snrInd,:)

expCol2 = jCol(ceil(allExp2*99/max(allExp2)),:);
expCol1 = jCol(ceil(allExp1*99/max(allExp1)),:);

clf
% scatP = scatter(allPol1(pickPair),allPol2(pickPair),5,'o',...
%     'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.002)
scatP = scatter(allPol1(pickPair),allPol2(pickPair),40,'o','MarkerFaceColor',...
    'none','LineWidth',3);
scatP.CData = expCol1(pickPair,:);
hold on
scatP = scatter(allPol1(pickPair),allPol2(pickPair),20,'o','MarkerFaceColor',...
    'flat','LineWidth',0.1,'markeredgealpha',0);
scatP.CData = expCol2(pickPair,:);

xlim([-1 1])
ylim([-1 1])

hold on
plot([-1 1],[-1 1])
hold off




