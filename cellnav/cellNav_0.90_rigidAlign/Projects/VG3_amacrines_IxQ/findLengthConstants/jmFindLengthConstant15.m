
%%Check for unassigned bipolar cells
%%Determine if polarity of Ca responses are consistent with other cells. Is
%%it SNR dependent.
%%See if changing masks changes bimodality or polarity spread.
%%Consolidate redundant ROIs

%{

1) Karl finishes cell 14, bip identification, ROI correlation
2) Try filtering ROIs by edge proximity

4) In parrallel run model against population of vGlut

Responses to flashing Bar


%}



if 1
    global glob
    datFold = [glob.datDir 'Analysis\Data\preproc\'];
    SPN =datFold;
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'GOI.mat']);
end

realCal = GOI.Polarity1;
SNR = GOI.qual;

targetType = 'real'; %Keep original calcium data or change to test model
%% real = keep observed values
%% allMean = replace values with the mean
%% meanPlusNorm = use mean + normal noise * caNoise
%% evenRand = use random distribution evenly spread between -1 and 1
%% meanPlusNorm = use mean + normal noise * caNoise

caNoise = 3;

if strcmp(targetType,'real')
    caPol = realCal;
elseif strcmp(targetType,'allMean')
    caPol = randn(size(caPol))*0 + meanCaPol; % all values equal mean
elseif strcmp(targetType,'real')
    caPol = randn(size(caPol)) * caNoise + meanCaPol;
elseif strcmp(targetType,'evenRand')
    caPol = rand(size(caPol))*2-1; % even random distribution between -1 and 1
elseif strcmp(targetType,'real')


end
meanCaPol = mean(caPol(:));

%%Correct for weird errors that create outliers
caPol = caPol * -1;
caPol(caPol>1) = 1;
caPol(caPol<-1) = -1;




standardize = 0; %scale distributions by standard deviation before calculating error
matchType = 1; % 1 = rmse, 2 = cc
weightErrors = 0; % multiply covariance by weight (SNR)
filterBySNR = 0.9;

filterByEdge = 0;
onScale = 1;%[.6:.1:3];
offScale = [1];
noise = [0:.001:.005];

%Remove bad frames
frames = ptDat(:,1);
allFrames = unique(frames); %all
useFrames = allFrames;
%useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 2006]; %all
useFrames = [1005 1006 1007 1008 1009 1010 2001 2002 2003 2004 2005 ]; %all but worst
%useFrames = [1005 1006 1007 1008 1009 1010 ]; % only first stack
removeFrames = setdiff(allFrames,useFrames);

%% Remove ROIs near edge
roiNearEdge = [];[66 67 68 69 71 120 144 70 198 ]; % manually identified by jm 1/6/2022 (GOI?)




SNR(isnan(SNR)) = 0;
SNRcol = colorProp(SNR,'STD');
minBip = 5;
testLengths = 0:2:100;
L = length(testLengths);


roiCids = GOI.roiCids;%ptDat(:,3);
numRoi = length(roiCids);

allPred = zeros(numRoi,L,length(onScale),length(offScale),length(noise));
goodRoiCid = zeros(numRoi,1);
vCids = SOI.cids;
vCids = [2 3 4  13];
goodCid = zeros(length(vCids),1);

%%Get weight


for v = 1:length(vCids)
    vCid = vCids(v); %Get cid of VG3 to run
    vTarg = find(SOI.cids == vCid);
    isCid = find(GOI.roiCids==vCid);  % Find all grouped rois for the cid being run
    useNodes = GOI.closeNode(isCid); % regtrieve list of skeleton nodes for rois on cid

    d = SOI.cell(vTarg).d(:,useNodes); % Get matrix of distiances from all synapse nodes (Y) to all functional roi nodes (X)
    preSign = SOI.cell(vTarg).preSign; % Get list of signs (0=ama, 1 = on Bip, 2 = off Bip) for synapses onto VG3
    numBip = sum(preSign>0); % count bipolar cell inputs
    preSignMat = repmat(preSign,[1 size(d,2)]); % make matrix of synapse signs with the same shape as d
    onMat = preSignMat == 2; % make logical matrix of on synapses with shape of d
    offMat = preSignMat == 1; % as onMat
    if (sum(preSign==1)>=minBip) && (sum(preSign==2)>=minBip) % if there are at least minBip of on and off bip inputs (each)
        fprintf('%d is good\n',vCid)
        goodRoiCid(isCid) = 1;
        goodCid(v) = 1;
    else
        fprintf('%d is bad\n',vCid)
    end

    predPol = zeros(length(isCid),L,length(onScale),length(offScale)); %Make matrix for results of each test variation
    for c = 1:L
        lc = testLengths(c); % get length constant

        W = exp(-d/lc); % Apply length constant to d (synapse to roi matrix)
        for n = 1:length(noise)
            for s = 1:length(onScale)
                for s2 = 1:length(offScale)
                    Won = W .* onMat * onScale(s) + noise(n); %calculate weight of all ON inputs for each roi
                    Woff = W .* offMat * offScale(s2) + noise(n); % same as Won

                    sumOn = sum(Won,1) ; % sum all on influences across each ROI
                    sumOff = sum(Woff,1) ; % as sumOn

                    predPol(:,c,s,s2,n) = (sumOff-sumOn)./(sumOff+sumOn); % enter polarity estemat for vector of rois
                end
            end
        end

    end
    allPred(isCid,:,:,:,:) = predPol; %Enter predictions for cell into full matrix of experiments. First dimention ges back into full roi space
end




%% Filter Rois
'finish reading out noise!!!!'
%useRoi = find(~isnan(sum(allPred,2)));

goodRoi = goodRoiCid;
mean(goodRoi)

if filterBySNR
    goodRoi = goodRoi & (SNR>filterBySNR);
end
mean(goodRoi)

if filterByEdge
    for r = 1:length(GOI.roiID)
        hit = intersect(GOI.roiID{r},roiNearEdge);
        if ~isempty(hit)
            goodRoi(r) = 0;
        end
    end
end
mean(goodRoi)


if exist('useFrames')
    for i = 1:length(removeFrames)
        for g = 1:length(GOI.roiID)
            grIDs = GOI.roiID{g};
            gFrames = frames(grIDs);
            if sum(gFrames==removeFrames(i))
                goodRoi(g) = 0;
            end
        end
    end
end
mean(goodRoi)

useRoi = find(goodRoi>0);


%% Find error for all
usePred = allPred(useRoi,:,:,:,:);

%%Set error weights
if weightErrors == 1
    ew = SNR;% * 0 + 1;
else %dont weight errors
    ew = SNR * 0 + 1;
end
ew = ew(goodRoi>0);
ew = ew/mean(ew);
ewM = repmat(ew,[1 size(usePred,2) size(usePred,3) size(usePred,4) size(usePred,5)]);


caPolN = caPol(useRoi);% - mean(caPol(useRoi));
if standardize
    caPolN = caPolN/std(caPolN);
end
usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
if standardize
    stdUsePred = std(usePredN,1);
    usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1 1]);
end
caPolMatN = repmat(caPolN,[1 L size(usePred,3) size(usePred,4) size(usePred,5)]);



%%Get correlation coefficient
caPolMatM2 = caPolMatN - repmat(mean(caPolMatN,1),[size(caPolMatN,1) 1 1 1 1]);
usePredM2 = usePredN - repmat(mean(usePredN,1),[size(usePredN,1) 1 1 1 1]);
covN = mean(caPolMatM2 .* usePredM2 .* ewM,1);
std1 = std(caPolMatM2,1);
std2 = std(usePredM2,1);
cc = covN./(std1.*std2);


if 0
    %%Calculate error with rmse
    difMatN = caPolMatN - usePredN;
    meanErrN = mean(abs(difMatN));
    rmseN = sqrt(mean(difMatN.^2,1));
    rmseN = rmseN/ mean(abs(caPolN));
    rmseN = (std1*3)-rmseN./std1; %Adjust for positive error measure
elseif 1
    difMatN = caPolMatN - usePredN;
    meanDif = mean(abs(difMatN),1);
    rmseN = 0 - meanDif;
else
    sameSign = (caPolMatN./abs(caPolMatN)) == (usePredN./abs(usePredN));
    mean(sameSign,1);
    rmseN = mean(sameSign,1);
end


if matchType == 1
    matchVal = rmseN;
else
    matchVal = cc;
end

%%Apply weights




%%Find peak match (average) in case of multiple hits)
minE = max(matchVal(:))
indE = find(matchVal==minE);
[by bx bs1 bs2 bn] = ind2sub(size(matchVal),indE);
disp(sprintf('Found %d hits',length(by)))
bym = mean(by);
bxm = mean(bx);
bs1m = mean(bs1);
bs2m = mean(bs2);
bnm = mean(bn);


bestOnScale = onScale(round(bs1m))
bestOffScale = offScale(round(bs2m))
bestLength = testLengths(round(bxm))
bestNoise = noise(round(bnm))


clf

subplot(3,2,1)
showErr = rmseN(:,:,round(bs1m),round(bs2m),round(bnm));
plot(testLengths,showErr)
title({sprintf('best length = %0.1f, OnScale = %0.1f,OffScale = %0.1f,', bestLength,bestOnScale,bestOffScale),...
    sprintf('noise = %0.4f, SNR filter = %0.2f',bestNoise,filterBySNR)})


subplot(3,2,2)
showErr = cc(:,:,round(bs1m),round(bs2m),round(bnm));
plot(testLengths,showErr)
title(sprintf('cc'))



showPred = usePredN(:,:,round(bs1m),round(bs2m),round(bnm));
%showPred = squeeze(usePredN(:,round(bxm),:,round(bnm))); %% show all on scaling
%showPred = squeeze(usePredN(:,round(bxm),round(bzm),:)); %% show all noise



%%Show shift in matching with length constants
if 1


    subplot(3,2,3)

    stdDif = abs(testLengths-15);
    std15 = find(stdDif==min(stdDif),1);
    scatS = scatter(usePredN(:,std15,round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    scatS.CData = SNRcol(useRoi,:);
    xlim([-1 1])
    ylim([-1 1])
    title(sprintf('distribution at 15 um standard'))
    set(gca,'color','k')

    drawnow


    subplot(3,2,4)

    if 1
        for r = 1
            for s = 1:size(showPred,2)
                scat1 = scatter(showPred(:,s),caPolN,15,'filled');
                hold on
                plot([-1 1],[-1 1],'w')
                scat1.CData = SNRcol(useRoi,:);
                title(sprintf('length constant %0.2f',testLengths(s)))
                xlim([-1 1])
                ylim([-1 1])
                set(gca,'color','k')
                drawnow
                pause(.1)
                hold off
            end
            for s = size(showPred,2):-1:1
                scat1 = scatter(showPred(:,s),caPolN,15,'filled');
                hold on
                scat1.CData = SNRcol(useRoi,:);
                title(sprintf('length constant %0.2f',testLengths(s)))
                xlim([-1 1])
                ylim([-1 1])
                set(gca,'color','k')
                drawnow
                hold off
                %pause(.1)
            end

        end
    end

    scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    scat1.CData = SNRcol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(bxm)))
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'color','k')

    drawnow


end

if 0
    clf
    jCol = jet(100);
    snrInd = SNR;
    snrInd = snrInd-min(snrInd);
    snrInd = round(snrInd * 299/max(snrInd(:))) + 1;
    snrInd(isnan(snrInd)) = 1;
    snrInd(snrInd<1) = 1;
    snrInd(snrInd>100) = 100;
    snrCol = jCol(snrInd,:)


    scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,20,'o','MarkerFaceColor',...
        'none','LineWidth',1);
    scat1.CData = snrCol(useRoi,:);
    title(sprintf('length constant %0.2f',testLengths(bxm)))
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'color','w')

    hold on
    plot([-1 1],[-1 1])

    return

end

%% Compare errror to SNR

if 0
    bestErrors = difMatN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
    scatter(SNR(useRoi),bestErrors)
    caN = caPolMatN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
    uPN = usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm));
    scatter3(caN,uPN,SNR(useRoi))
end


%% Check each roi
useRoi = find(goodRoi>0);
subplot(3,2,5)
cla, hold on
for r = 1:length(useRoi);

    usePredR = usePredN(r,:,:,:,:);
    caPolNR = caPolN(r);% - mean(caPol(useRoi));
    dif = abs(usePredR-caPolNR);

    minE = min(dif(:));
    indE = find(dif==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(dif),indE);
    bymR = mean(by);
    bxmR = mean(bx);
    bs1mR = mean(bs1);
    bs2mR = mean(bs2);
    bnmR = mean(bn);

    bestOnScale = onScale(round(bs1m));
    bestOffScale = offScale(round(bs2m));
    bestLengths(r) = testLengths(round(bxmR));
    bestNoise = noise(round(bnm));

    %     subplot(2,2,3)
    %     hold off
    %     scat1 = scatter(usePredN(:,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN,15,'filled');
    %     scat1.CData = SNRcol(useRoi,:);
    %     title(sprintf('length constant %0.2f',testLengths(s)))
    %     %                 xlim([-1 1])
    %     %                 ylim([-1 1])
    %     set(gca,'color','k')
    %     hold on
    %     scat2 = scatter(usePredN(r,round(bxm),round(bs1m),round(bs2m),round(bnm)),caPolN(r),150,'filled');
    %
    if 0
        subplot(3,2,5)
        showDif = dif(:,:,round(bs1m),round(bs2m),round(bnm))
        hold on
        plot(showDif)
        drawnow
    end
end

useBest = (bestLengths>testLengths(1)) & (bestLengths<testLengths(end));
useBestLengths = bestLengths(useBest);

percentLengthsUsed = mean(useBest)*100

subplot(3,2,5)
hold off
hist(useBestLengths,testLengths)

medBestLengthsR = median(useBestLengths)
sortBL = sort(useBestLengths,'ascend');
bl95 = [sortBL(round(length(sortBL) * .025))      sortBL(round(length(sortBL) * .975))]
SE = std(useBestLengths)/sqrt(length(useBestLengths))
title(sprintf('L = %0.2f, 95%% = %0.2f - %0.2f, SE = %0.2f',medBestLengthsR, bl95(1),bl95(2)),SE)



%% Show each cell

useCid = vCids(goodCid>0);
pCol  = [1 0 0; 0 1 0; 0 0 1; .75 .75 0; .75 0 .75; 0 .75 .75; .3 0 0; 0 .3 0; 0 0 .3; .3 .3 0; .3 0 .3; 0 .3 .3];

bestOnScales = zeros(1,length(useCid));
bestOffScales = zeros(1,length(useCid));
bestLengths = zeros(1,length(useCid));
bestNoises = zeros(1,length(useCid));


for v = 1:length(useCid)
    vCid = useCid(v);
    useRoi = find((goodRoi>0) & (GOI.roiCids==vCid));

    usePred = allPred(useRoi,:,:,:,:);
    caPolMat = repmat(caPol(useRoi),1, L, size(usePred,3), size(usePred,4), size(usePred,5));
    difMat = caPolMat - usePred;
    rms = sqrt(mean(difMat.^2,1));
    meanErr = mean(abs(difMat),1);


    %Set error weights
    if weightErrors == 1
        ew = SNR;% * 0 + 1;
    else %dont weight errors
        ew = SNR * 0 + 1;
    end
    ew = ew(useRoi);
    ew = ew/mean(ew);
    ewM = repmat(ew,[1 size(usePred,2) size(usePred,3) size(usePred,4) size(usePred,5)]);


    caPolN = caPol(useRoi);% - mean(caPol(useRoi));
    if standardize
        caPolN = caPolN/std(caPolN);
    end
    usePredN = usePred;% - repmat(mean(usePred,1),[size(usePred,1) 1]);
    if standardize
        stdUsePred = std(usePredN,1);
        usePredN = usePredN./repmat(stdUsePred,[size(usePred,1) 1 1]);
    end
    caPolMatN = repmat(caPolN,[1 L size(usePred,3) size(usePred,4) size(usePred,5)]);



    %%Get correlation coefficient
    caPolMatM2 = caPolMatN - repmat(mean(caPolMatN,1),[size(caPolMatN,1) 1 1 1 1]);
    usePredM2 = usePredN - repmat(mean(usePredN,1),[size(usePredN,1) 1 1 1 1]);
    covN = mean(caPolMatM2 .* usePredM2 .* ewM,1);
    std1 = std(caPolMatM2,1);
    std2 = std(usePredM2,1);
    cc = covN./(std1.*std2);


    if 0
        %%Calculate error with rmse
        difMatN = caPolMatN - usePredN;
        meanErrN = mean(abs(difMatN));
        rmseN = sqrt(mean(difMatN.^2,1));
        rmseN = rmseN/ mean(abs(caPolN));
        rmseN = (std1*3)-rmseN./std1; %Adjust for positive error measure
    elseif 1
        difMatN = caPolMatN - usePredN;
        meanErrN = mean(abs(difMatN));
        meanDif = mean(abs(difMatN),1);
        rmseN = 0 - meanDif;
    else
        sameSign = (caPolMatN./abs(caPolMatN)) == (usePredN./abs(usePredN));
        mean(sameSign,1);
        rmseN = mean(sameSign,1);
    end

    if matchType == 1
        matchVal = rmseN;
    else
        matchVal = cc;
    end

    minE = max(matchVal(:))
    indE = find(matchVal==minE);
    [by bx bs1 bs2 bn] = ind2sub(size(matchVal),indE);
    bym = mean(by);
    bxm = mean(bx);
    bs1m = mean(bs1);
    bs2m = mean(bs2);
    bnm = mean(bn);

    bestOnScales(v) = onScale(round(bs1m))
    bestOffScales(v) = offScale(round(bs2m))
    bestLengths(v) = testLengths(round(bxm))
    bestNoises(v) = noise(round(bnm))

    showErr = matchVal(:,:,round(bs1m),round(bs2m),round(bnm));
    showMeanErr = meanErrN(:,:,round(bs1m),round(bs2m),round(bnm));

    showPred = usePredN(:,:,round(bs1m),round(bs2m),round(bnm));


    subplot(3,2,6)
    hold on
    plot(testLengths,showErr,'color',pCol(v,:))
    title(' ')

    %     ylim([0 1])
    %     xlim([0 40])
    title(sprintf('%d rois for cell %d',length(useRoi),vCid))

    drawnow
end
useCid

bestOnScales
bestOffScales
bestLengths
bestNoises




