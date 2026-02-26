
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

%% set variables
useVgc = [2 3 4 5 13 14];
zJit = .5;
reps = 10;
lc = 16;


%% colect data on cells

cellSynCount = zeros(length(s),length(useVgc));
clear synPos cellSynPos cellSource sClose synIdsA
countG = zeros(length(s),1);
for i = 1:length(s)
    synPos{i} = [];
    cellSource{i} = [];
    sClose{i} = [];
    synIdsA{i} = [];
end
% allPos = [];
% for v = 1:length(useVgc)
%
%     cid = useVgc(v);
%     useSM(v) = 1;
%     disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
%     sm = sms(v).sm;%load([smDir fileName],'syn','syn2Skel','nep');
%
%     syn = sm.syn; % get syn information
%     synD = sm.syn2Skel.syn2SynDist;
%     skelD = sm.syn2Skel.syn2SkelDist;
%     nodeLengths = sm.nep.props.nodeLength;
%     closest = sm.syn2Skel.closest;
%
%
%     allPos = cat(1,allPos,sm.arbor.subs);
%
%     %%Create synapse groups filled with relevant IDs using COI and subtypes
%
%
%
%
% end
%
% allSynCount = sum(cellSynCount,2);


%%  pick synapses

sNames = {s.name};
clear sgR sgB sNum rWDAllV wDAllV


sgB(1).names = {'off'};
sgB(1).col = [0 1 0];
sgB(2).names = {'on'};
sgB(2).col = [1 0 0];

sgR(1).names = {'4on'};
sgR(1).col = [1 0 0];
sgR(2).names = {'4i'};
sgR(2).col = [.6 .6 0];
sgR(3).names = {'4ow'};
sgR(3).col = [.6 0 .6];
sgR(4).names = {'37'};
sgR(4).col = [.6 0 .6];
sgR(5).names = {'5ti'};
sgR(5).col = [.5 .5 .5];
sgR(6).names = {'63'};
sgR(6).col = [0 .6 .6];
sgR(7).names = {'6sw'};
sgR(7).col = [0 0 1];

for i = 1:length(sgR)
    rNames{i} = sgR(i).names{1};
end

%% pick s's for sgs
expNames = cat(1,{sgB(:).names});
for g = 1:length(sgB)
    nams = sgB(g).names;
    for n = 1:length(nams)
        m1 = find(strcmp(sNames,sgB(g).names{n}));
        if ~isempty(m1)
            sgB(g).sT(n) = m1;
        else
            sgB(g).sT(n) = [];
        end
    end
end

expNames = cat(1,{sgR(:).names});
for g = 1:length(sgR)
    nams = sgR(g).names;
    for n = 1:length(nams)
        m1 = find(strcmp(sNames,sgR(g).names{n}));
        if ~isempty(m1)
            sgR(g).sT(n) = m1;
        else
            sgR(g).sT(n) = [];
        end
    end
end



%% Find depths
%binDepths(sg,nDepth)

%% rand each cell

pol = zeros(length(sgR),length(useVgc));
rPol = zeros(length(sgR),length(useVgc),reps);
for v = 1:length(useVgc)
    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('running sm %d, %d of %d',cid,v,length(useVgc)))

    targ = find(runCids==cid);
    sm = sms(targ).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    D = sm.skel2skel.linDist; %use distance of every skeleton node to every other skeleton node
    nodeLengths = sm.nep.props.nodeLength;
    closest = sm.syn2Skel.closest;
    synIds = classifySynInCell(syn,s,cid);


    %allPos = cat(1,allPos,sm.arbor.subs);
    %%Get new node positions
    pos = sm.nep.pos;
    [zGCL,zINL,nDepth] = getIPLdepth(pos(:,3),pos(:,1),pos(:,2),[],[]);
    iplDepth = mean(zINL)-mean(zGCL);
    nDepth = nDepth * iplDepth;

    zDif = abs(nDepth-nDepth'); %matrix of every node to every node in depth
    zMask = zDif<zJit;

    %%Measure

    for g = 1:length(sgR)
        sgRIds = [synIds{sgR(g).sT}];
        sgRClose{g} = closest(sgRIds);
        sRNum(g,v) = length(sgRIds);
    end

    for g = 1:length(sgB)
        sgBIds = [synIds{sgB(g).sT}];
        sgBClose{g} = closest(sgBIds);
        sBNum(g,v) = length(sgBIds);
    end


    clear minD meanD wD wDAll%% wD is influence of every group on every group
    for g1 = 1:length(sgR)
        for g2 = 1:length(sgB)
            d = D(sgRClose{g1},sgBClose{g2});
            W = exp(-d/lc);
            if sum(d(:))
                minD(g1,g2) = mean(min(d,[],1));
                meanD(g1,g2) = mean(d(:));
                wD(g1,g2) = sum(W(:));
                wDAll{g1,g2} = sum(W,2);

            else
                minD(g1,g2) = 0;
                meanD(g1,g2) = 0;
                wD(g1,g2) = 0;
                wDAll{g1,g2} = 0;

            end
        end
    end




    %%Randomize
    rMinD = zeros(length(sgR),2,reps);
    rWD = rMinD;
    clear rWDAll;

    if 1 %show synapses
        clf
        scatter3(pos(:,1),pos(:,2),pos(:,3),'k','.');
        hold on
        axis equal
        set(gca,'clipping', 'off')
        view([0 0])
        %scatter3(pos(closest(isBip),1),pos(closest(isBip),2),pos(closest(isBip),3),'m','o','filled');
        synCol = hsv(length(sgR));
        for g = 1:length(sgR)
            scatter3(pos(sgRClose{g},1),pos(sgRClose{g},2),pos(sgRClose{g},3),'markerfacecolor',synCol(g,:));
        end
        pause (.1)
    end


    for r = 1:reps
        disp(sprintf('running sm %d, %d of %d, rep %d',cid,v,length(useVgc),r))

        %make random pointers from nodes to node with similar z

        closeMask = zMask(closest,:);
        if 0
            while 1
                zRand = rand(size(closeMask)) .* closeMask;
                zMax = repmat(max(zRand,[],2),[1 size(zRand,2)]);
                [y rClose] = find(zRand==zMax);
                if length(rClose) == size(zMax,1)
                    break
                end
            end
        else

            %pick random node that is within distance from close mask
            zRand = rand(size(closeMask)) .* closeMask;
            [sZRand idx] = sort(zRand,2,'descend');
            rClose = idx(:,1);

        end

        %%Measure Random
        for g = 1:length(sgR)
            sgRIds = synIds{sgR(g).sT}; % Get IDs of synapses for rgc subtype
            sgRClose{g} = rClose(sgRIds); % Get nodes closest to synapses
        end
        for g1 = 1:length(sgR) % For every cell type
            for g2 =  1:length(sgB)
                d = D(sgRClose{g1},sgBClose{g2});
                W = exp(-d/lc);
                if sum(d(:))
                    rMinD(g1,g2,r) = mean(min(d,[],1));
                    rWD(g1,g2,r) = sum(W(:));
                    rWDAll{g1,g2,r} = sum(W,2);
                else
                    rMinD(g1,g2) = 0;
                    rWD(g1,g2) = 0;
                    rWDAll{g1,g2,r} =0;
                end
            end
        end


        if 0 %show scramble
            clf
            isBip = find(syn.preClass == 7);
            scatter3(pos(:,1),pos(:,2),pos(:,3),'k','.');
            hold on
            %scatter3(pos(closest(isBip),1),pos(closest(isBip),2),pos(closest(isBip),3),'m','o','filled');
            scatter3(pos(oldClose,1),pos(oldClose,2),pos(oldClose,3),'r','filled');
            scatter3(pos(newClose,1),pos(newClose,2),pos(newClose,3),'g','filled');
            hold off
            axis equal
            set(gca,'clipping', 'off')
            view([0 0])
            pause (.1)
        end
    end

    %%Get polarities

    pol(:,v) = (wD(:,2)-wD(:,1))./(wD(:,2)+wD(:,1)); %Calculate polarity
    rPol(:,v,:) = (rWD(:,2,:)-rWD(:,1,:))./(rWD(:,2,:)+rWD(:,1,:)); % random polarity
    rWDAllV{v} = rWDAll;
    wDAllV{v} = wDAll;

    if 1 %%Show hist
        polRange = [-1:.01:1];
        clear hPol;
        clf
        for p = 1:size(pol,1)
            subplot(size(pol,1),1,p)
            hold on
            hPol(p,:) = hist(squeeze(rPol(p,v,:)),polRange);
            bar(polRange,hPol(p,:),'k');
            scatter(pol(p),0,'r','filled')
        end
        pause(.5)
    end

end
mrPol = squeeze(mean(rPol,3)); %mean random polarity

sPol = (sBNum(2,:)-sBNum(1,:))./(sBNum(2,:)+sBNum(1,:));


%% Collect results for cells
clear raPol mraPol rWDallV
for r = 1:reps
    for g1 = 1:length(sgR)
        for g2 = 1:length(sgB)
            vals = [];
            for v = 1:length(useVgc)
                vals = [vals; rWDAllV{v}{g1,g2,r}];
            end
            rWDallV{g1,g2,r} = vals;
        end
        raPol{g1,r} =  (rWDallV{g1,2,r} - rWDallV{g1,1,r})./(rWDallV{g1,2,r}+rWDallV{g1,1,r}); % random polarity
        mraPol(g1,r) = mean(raPol{g1,r});
    end
end

clear aPol wDallV maPol
for g1 = 1:length(sgR)
    for g2 = 1:length(sgB)
        vals = [];
        for v = 1:length(useVgc)
            vals = [vals; wDAllV{v}{g1,g2}];
        end
        wDallV{g1,g2} = vals;
    end
    aPol{g1} =  (wDallV{g1,2} - wDallV{g1,1})./(wDallV{g1,2}+wDallV{g1,1}); % random polarity
    maPol(g1) = mean(aPol{g1});
end

totON = sum(sBNum(2,:));
totOFF = sum(sBNum(1,:));
totPol = (totON-totOFF)/(totON+totOFF);


clf
statMraPol = quartileBars(mraPol');
scatter(1:length(maPol),maPol,300,'filled','r','markeredgecolor','k','markerfacealpha',.3)
plot([1:length(maPol)],zeros(length(maPol),1),'k')
plot([1:length(maPol)],zeros(length(maPol),1)+totPol,'r')

xticks(1:length(sgR))
xticklabels(rNames)
ylim([-1 1])

if 0
    epsName = 'Z:\Active\morganLab\PUBLICATIONS\VG3\Figure\Draft4\Pics2\MonteStratRGCPred.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
end


%% Compare 4s to 6s
clf
pickA = [1 2 3];
pickB = [7];

polA = cat(1,aPol{pickA});
polB = aPol{pickB};

wDallV{g1,2}


meanA = mean(polA);
seA = std(polA)/sqrt(length(polA));
[meanA-seA meanA+seA]


meanB = mean(polB);
seB = std(polB)/sqrt(length(polB));
[meanB - seB meanB+seB]

hRange = [-1:.2:1];
hA = hist(polA,hRange);
hB = hist(polB,hRange);
bar(hRange,[hA/sum(hA);hB/sum(hB)],'barwidth',1.8)


%% Linear regression
statMraPol
maPol

y = statMraPol(:,3);
p = polyfit(y,maPol,1);
yp =  polyval(p,y);
yresid = y- yp;
SSresid = sum(yresid.^2);
SStotal = sum((y-mean(y)).^2);
Rs = 1 - SSresid/SStotal;
Rs_adg = 1 - SSresid/SStotal * (length(y)-1)/(length(y) - length(p));
corrcoef(y,maPol)




return

%% Display each cell
clf
for g = 1:size(pol,1)

    subplot(1,size(pol,1),g)
    hold on
    plot([ones(1,length(useVgc));ones(1,length(useVgc))+1],[pol(g,:); mrPol(g,:)],'k','linewidth',2)
    plot([ones(1,length(useVgc))+1;ones(1,length(useVgc))+2],[mrPol(g,:); sPol],'color',[.3 .3 1],'linewidth',2)

    scatter(ones(1,length(useVgc)),pol(g,:),'markerfacecolor',[.8 .8 .8],'markeredgecolor',[0 0 0]);
    scatter(ones(1,length(useVgc))+1,mrPol(g,:),'markerfacecolor',[.8 .8 .8],'markeredgecolor',[0 0 0]);
    scatter(ones(1,sum(~isnan(pol(g,:))))+2,sPol(~isnan(pol(g,:))),'markerfacecolor',[.7 .7 1],'markeredgecolor',[.3 .3 1]);

    plot([ones(1,length(useVgc));ones(1,length(useVgc))+1],[pol(g,:); mrPol(g,:)],'k','linewidth',2)
    plot([ones(1,length(useVgc))+1;ones(1,length(useVgc))+2],[mrPol(g,:); sPol],'color',[.3 .3 1],'linewidth',2)

    scatter(ones(1,length(useVgc)),pol(g,:),'markerfacecolor',[.8 .8 .8],'markeredgecolor',[0 0 0]);
    scatter(ones(1,length(useVgc))+1,mrPol(g,:),'markerfacecolor',[.8 .8 .8],'markeredgecolor',[0 0 0]);
    scatter(ones(1,sum(~isnan(pol(g,:))))+2,sPol(~isnan(pol(g,:))),'markerfacecolor',[.7 .7 1],'markeredgecolor',[.3 .3 1]);

    ylim([-1 1])
    title(sgR(g).names)
end





meanSem(pol(1,:))
meanSem(pol(2,:))
meanSem(mrPol(1,:))
meanSem(mrPol(2,:))

%%

rgcTypesPol = {'1wt'  1; '2aw'  1; '2an'  1;  '4i'  1; '4on'  1; '4ow' 1
    ; '3i' .7;'37'  0.5; '5si'  .5;'5to' .3
    ;'25' 0.1 ;'28'  0.1; '5so' 0.1
    ;'5ti'  0;'63'   -.3
    ;'85' -0.9;   '7i'  -0.9; '72'   -0.9;
    '6sw' -1;'6sn'   -1;'6t' -1;'8w' -1;  '51' -.1};

shiftPolAll = [-.15 -.1  -.05 0 .05 .1 0 0 .05 0 -.05 0 .05 0 0 -.05 0 .05 -.1 -.05 0 .05 1];


for i = 1:length(useRGCs)
    useRGroup(i) = find(strcmp(COI.rgcGroupLabel,useRGCs{i}));
end

useBGroup = [ 8 9 11 15 12 13 14 17];
bLab = COI.bpcGroupLabel(useBGroup);
rLab = COI.rgcGroupLabel(useRGroup);

clear physPol shiftPol
for i = 1:length(rLab)
    [y x] = find(strcmp(rgcTypesPol,rLab{i}));
    physPol(i) = rgcTypesPol{y,2};
    shiftPol(i) = shiftPolAll(y);
end



hRange = [0:.2:2];
physPol2 = physPol * -1;
realSac = [];
randSac = [];
for g = 1:size(pol,1);
    realPol = pol(g,:);
    randPol = mrPol(g,:);

    isReal = ~isnan(realPol) & ~isnan(randPol);
    realSac = [realSac abs(physPol2(g) - realPol(isReal))];
    randSac = [randSac abs(physPol2(g) - randPol(isReal))];


end


hReal = hist(realSac,hRange);
hRand = hist(randSac,hRange);

if 0
    clf
    plot(hRange,hReal,'k')
    hold on
    plot(hRange,hRand,'r')
    hold off
end

%% Get correlation coefficent for randoms and real cells
isnum = find(~isnan(pol));

expected = repmat(physPol2', [1 size(rPol,2)]);
cc = corrcoef(pol(isnum),expected(isnum));
cco = cc(2,1); % observed cell to cell correlation coeficient
ccr = zeros(reps,1);
for r = 1:reps

    testR = rPol(:,:,r);
    cc = corrcoef(testR(isnum),expected(isnum));
    ccr(r) = cc(2,1);

end

cco
sCCR = sort(ccr,'ascend');
medianCCR = median(ccr)
ci95 = [sCCR(ceil(reps * 0.05)) sCCR(ceil(reps * 0.95))]
realFrac = mean(sCCR<cco)

%% Analyze trend



%{
function binDepths(sg,nDepth)

%% binDepths
clf
bin = .04;
dRange = -.1:.001:1.1;

maxCount = 0;

for i = 1:length(sg)

    dCount = dRange*0;
    for b = 1:length(dRange)
        d = dRange(b);
        dCount(b) = sum((sg(i).depth>(d-bin/2)) & (sg(i).depth<=(d+bin)));
    end
    sg(i).dCount = dCount;
    maxCount = max(maxCount,max(dCount));
end

clf
subplot(2,1,1)
hold on
for i = 1:length(sg)
    plot(dRange,sg(i).dCount/maxCount,'color',sg(i).col)
end

subplot(2,1,2)
hold on
nCount = dRange*0;
bin2 = .1;
for b = 1:length(dRange)
    d = dRange(b);
    nCount(b) = sum((nDepth>(d-bin2/2)) & (nDepth<=(d+bin2)));
end
nCount = nCount/max(nCount(:));
plot(dRange,nCount,'color',[.3 .3 .3])


ei = sg(1).dCount./sg(2).dCount;
plot(dRange,ei,'color',[0 0 1])



iDens = (sg(2).dCount/maxCount)./(nCount/max(nCount(:)));
plot(dRange,iDens,'color',[1 0 0])
% 
% 
% eDens = (sg(1).dCount/maxCount)./(nCount/max(nCount(:)));
% plot(dRange,eDens,'color',[0 1 0])

end

%}






