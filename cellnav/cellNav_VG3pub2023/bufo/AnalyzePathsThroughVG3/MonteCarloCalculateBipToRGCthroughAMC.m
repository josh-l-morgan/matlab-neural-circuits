

%% load data
if 0

    global glob tis

    figure
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
        %sm = load([smDir fileName],'skel2skel','nep','syn2skel','syn');
        load([smDir fileName]);
        sm.nep.swcS = nep2swc(sm.nep);
        sms(i).sm = sm;
    end
    useVgc = runCids(useSM>0);

    s = makeSynapseClassifyer(COI); %make structure describing types of synapses

end

%% Choose cells

%useRGCs = {'1wt' '2aw' '2an' '4i' '4ow' '28' '37' '5si' '5ti' '6sw' '6sn' '63' '85' '8w'}
%useRGCs = { '2an'  '25' '3i' '37' '4on' '4ow' '4i' '5to' '5so' '5ti' '5si'  '63' '6t' '6sw' '6sn' '7i' '8w'}
useRGCs = COI.checkRGCLabels;
useBPCs = COI.checkBPCLabels;

%% Bip to RGC influence

%%load SM
smDir = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\SMs\'
influenceMat = zeros(length(COI.bpcGroupCids), length(COI.rgcGroupCids));
influenceMatEach = zeros(length(COI.bpcGroupCids), length(COI.rgcGroupCids),length(useVgc));
storeInfluenceOnRgcSyn = cell(length(COI.bpcGroupCids),length(COI.rgcGroupCids), length(useVgc));

lc = 16;
for v  = 1:length(useVgc)
    disp(sprintf('checking vgc cid%d. %d of %d.',useVgc(v),v,length(useVgc)))
    sm = sms(v).sm;
    d = sm.syn2Skel.syn2SynDist;
    %d = sm.skel2skel.linDist;

    pre = sm.syn.pre;
    post = sm.syn.post;

    %% get all synapse to node distances

    W = exp(-d/lc); % Apply length constant


    %% run groups
    verbose = 0;
    trackRgcCids = [];
    for rG = 1:length(COI.rgcGroupCids)
        rCids = COI.rgcGroupCids{rG}
        trackRgcCids = cat(1,trackRgcCids,rCids(:));
        %useBGroup


        for bG = 1:length(COI.bpcGroupCids)
            bCids = COI.bpcGroupCids{bG};
            trackPost = zeros(1,length(pre),1);

            storeEach = [];
            %%run cells
            for b = 1:length(bCids)
                bCid = bCids(b);
                for r = 1:length(rCids)
                    gCid = rCids(r);

                    isPre = find(pre == bCid);
                    isPost = find(post == gCid);

                    if ~isempty(isPre) & ~isempty(isPost)

                        cutEDist = W(isPre,isPost);
                        if verbose
                            rgcGroupLabel{rG}
                            COI.bpcGroupLabel{bG}
                            cutEDist
                            pause
                        end

                        influenceMatEach(bG,rG,v) = influenceMatEach(bG,rG,v) + sum(cutEDist(:));
                        trackPost(isPost) = trackPost(isPost) + sum(cutEDist,1);

                        pause(.01)
                    end

                end %run each rgc cid

            end % run each bip cid
            storeInfluenceOnRgcSyn{bG,rG,v} = trackPost;

        end % run each bpc type
    end %% run each rgc type


end
influenceMat = sum(influenceMatEach,3);




%% extract influence matrix
clf
% -1 = ON, Eyewire
% 1wt = -1, 2aw = -1, 4i = -1, 4on = -1, 4ow= -1,
%28 = 0, 37 = -.75, 5si = -.5, 5ti = 0, 5to = -1, 51, 6sw, 63, 7, 72 = -1,
% 8w
%-1 = ON, rgctypes
% 1wt = 1, 2aw = 1, 2an = 1, 4i = 1, 4on = 1, 4ow= 1,
%28 = 0.1, 37 = -0.5, 5si = -.5, 5ti = 0, 5so = -.1, 51 = -.2, 6sw = -1,
%63 = -.9, 7 = -1, 72 = -1,
% 8w


rgcTypesPol = {'1wt'  1; '2aw'  1; '2an'  1;  '4i'  1; '4on'  1; '4ow' 1
    ; '3i' .7;'37'  0.5; '5si'  .5;'5to' .3
    ;'25' 0.1 ;'28'  0.1; '5so' 0.1
    ;'5ti'  0;'63'   -.3
    ;'85' -0.9;   '7i'  -0.9; '72'   -0.9; 
    '6sw' -1;'6sn'   -1;'6t' -1;'8w' -1;  '51' -.1};

shiftPolAll = [-.15 -.1  -.05 0 .05 .1 0 0 .05 0 -.05 0 .05 0 0 -.05 0 .05 -.15 -.1 0.5 0 .05 ];

for i = 1:length(useRGCs)
    useRGroup(i) = find(strcmp(COI.rgcGroupLabel,useRGCs{i}));
end


%useBGroup = [3 4 5 6 7 8 9 20];


for i = 1:length(useBPCs)
    useBGroup(i) = find(strcmp(COI.bpcGroupLabel,useBPCs{i}));
end

bLab = COI.bpcGroupLabel(useBGroup)
rLab = COI.rgcGroupLabel(useRGroup)

clear physPol shiftPol
for i = 1:length(rLab)
    [y x] = find(strcmp(rgcTypesPol,rLab{i}));
    physPol(i) = rgcTypesPol{y,2};
    shiftPol(i) = shiftPolAll(y);
end


yPos = physPol + shiftPol;


useInfluenceMat = influenceMat(useBGroup,useRGroup');


%% get bip to RGC type matricies

clf
allW = [];
allPreTypeID = [];
allPostTypeID = [];
for v  = 1:length(useVgc)
    disp(sprintf('checking vgc cid%d. %d of %d.',useVgc(v),v,length(useVgc)))
    targ = find(runCids==useVgc(v),1);
    sm = sms(targ).sm;
    d = sm.syn2Skel.syn2SynDist;
    pre = sm.syn.pre;
    post = sm.syn.post;
    W = exp(-d/lc); % Apply length constant
    preTypeID = zeros(1,length(pre));
    postTypeID = zeros(1,length(post));
    for bG = 1:length(useBGroup)
        b = useBGroup(bG);
        bCids = COI.bpcGroupCids{b};
        for c = 1:length(bCids)
            preTypeID(pre==bCids(c)) = bG;
        end
    end
    for rG = 1:length(useRGroup)
        r = useRGroup(rG);
        rCids = COI.rgcGroupCids{r};
        for c = 1:length(rCids)
            postTypeID(post==rCids(c)) = rG;
        end
    end

    allW(size(allW,1)+1:size(allW,1)+size(W,1),...
        size(allW,1)+1:size(allW,1)+size(W,1)) = W;
    allPreTypeID = cat(2,allPreTypeID,preTypeID);
    allPostTypeID = cat(2,allPostTypeID,postTypeID);

end

image(allW*100)
sNum = size(allW,2);
pause(.01)

%% Squish matrix to bpc groups
bGW = zeros(length(useBGroup),sNum);
for bG = 1:length(useBGroup)
    cutBG = allW(allPreTypeID==bG,:);
    if ~isempty(cutBG)
        bGW(bG,:) = sum(cutBG,1);
    end
end
image(bGW*10)

%% Squish matrix into ON OFF bip polarity
offG = [1 2 3];
onG = [4:8];

onW = sum(bGW(onG,:),1);
offW = sum(bGW(offG,:),1);
pW = (offW-onW)./(offW+onW);



%% Get real RGC group matrix
clf

%yPos = physPol + [.2 .1 0 -.1 -.2 0 0 .1 0 0 0 -.1 -.1 -.2];
%yPos = physPol ;
rCol = hsv(length(yPos));

%subplot(1,2,1)
hold on
pR = zeros(1,length(useRGroup));
eR = zeros(size(pW));
for rG = 1:length(useRGroup);
    isRGCType = allPostTypeID==rG;
    eR(isRGCType) = physPol(rG); %make matching vector of post synaptic polarity based on ephys
    typeSynPol = pW(isRGCType);
    % scatter(typeSynPol,zeros(length(typeSynPol),1)+yPos(rG),'k')
    scatter(typeSynPol,zeros(length(typeSynPol),1)+yPos(rG),60,'MarkerEdgeColor','none',...
        'MarkerFaceColor','k',MarkerFaceAlpha=.2)

    pR(rG) = mean(typeSynPol);
    scatter(pR(rG),yPos(rG),400,'r','+','linewidth',2)

end
[sortYPos idx] = sort(yPos);
yticks(sortYPos)
yticklabels(rLab(idx))
xlim([-1.3 1.3])
ylim([-1.3 1.3])

plot([0 0],[-1 1],'k')
title('published polarity vs synapse predictions')

if 0
    epsName = 'Z:\Active\morganLab\PUBLICATIONS\VG3\Figure\Draft4\Pics3\shadedScat.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
end
%% Plot relative to polarity and get fit
subplot(1,2,2)
hold on
isRgc = allPostTypeID>0;
scatter(pW(isRgc),eR(isRgc),'k')

linFit = fit(eR(isRgc)',pW(isRgc)','poly1')
linX = [-1:.01:1];
linY = feval(linFit,linX);
plot(linY,linX,'r')

plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
xlim([-1.1 1.1])
ylim([-1.1 1.1])

% %%
% xticks(linX)
% xticklabels(rgcTypesPol(useRGroup))

%% Condition data for analysis

eRc = eR(isRgc);
pWc = pW(isRgc);
rTypeC = allPostTypeID(isRgc);

[CC P RL RU] = corrcoef(eRc,pWc);
scatter(eRc,pWc)
realCC = CC(2,1)

reps = 10000;
synNum = length(eRc);
for r = 1:reps
    randE = eRc(randsample(synNum,synNum,1));
    CC =  corrcoef(randE,pWc);
    randCC(r) = CC(1,2);
end

ci = bounds95(randCC)


%[p tbl sts] = kruskalwallis(pWc,rTypeC)


%%Compare 4s to 6





%%
subplot(4,1,1:2)
maxInfluence = max(useInfluenceMat(:));
colormap(jet(100))
image(useInfluenceMat * 100/maxInfluence);
title('Bipolar (y) to RGC (x) influence normalized for RGC')
normInfluence = useInfluenceMat./repmat(sum(useInfluenceMat,1),...
    [size(useInfluenceMat,1) 1]);
image(normInfluence*100./max(normInfluence(:)))
xticks([1:length(rLab)])
xticklabels(rLab)
yticks([1:length(bLab)])
yticklabels(bLab)

sumInfluenceRGC = sum(useInfluenceMat,1);
sumInfluenceBPC = sum(useInfluenceMat,2);

subplot(4,1,3)

bar(sumInfluenceRGC,'k','edgecolor','none')
title('total RGC')
subplot(4,1,4)

bar(sumInfluenceBPC,'k','edgecolor','none')
title('total Bipolar')

return

%% Recheck merges when changes

onInfluence = sum(normInfluence(4:end,:),1);
offInfluence = sum(normInfluence(1:3,:),1);
OFFbias = (offInfluence-onInfluence)./(offInfluence + onInfluence);
bar(OFFbias)

bar(onInfluence);
bar(offInfluence);

onInfluence = sum(useInfluenceMat(4:end,:),1);
offInfluence = sum(useInfluenceMat(1:3,:),1);
OFFbias = (offInfluence)./(offInfluence + onInfluence);
bar(OFFbias)
ylim([0 1])


%% Get influence stats by RGC type

totInfluence = zeros(length(useBGroup),length(useRGroup));
clear sMeanPol sN
clf
rangeP = [0:.1:1];
for r = 1:length(useRGroup)
    onStack = [];
    for b = 1:length(useBGroup)
        allV = [];
        for v = 1:length(useVgc)
            allV = cat(2,allV,storeInfluenceOnRgcSyn{useBGroup(b),useRGroup(r),v});
        end
        onStack(b,:) = allV;
    end

    totInfluence(:,r) = sum(onStack,2);

    sumB = sum(onStack,1);
    onStack = onStack(:,sumB>0);

    normInfluenceS = onStack./repmat(sum(onStack,1), [size(onStack,1) 1]);

    onInfluenceS = sum(normInfluenceS(4:end,:),1);
    offInfluenceS = sum(normInfluenceS(1:3,:),1);
    OFFbiasS = (offInfluenceS)./(offInfluenceS + onInfluenceS);
    sMeanPol(r) = mean(OFFbiasS);
    sN(r) = length(OFFbiasS);


    histP = hist(OFFbiasS,rangeP);
    subplot(length(useRGroup),1,r)
    bar(rangeP,histP)
    axis on
    ylim([0 40])

end

%% Stratification









