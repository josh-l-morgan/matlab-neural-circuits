%%Plot density of synapses relative to the tips of neurites

global glob tis

clf
if 0
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

end

captureLength = 20; % How long a length of neurite to collect (from tip to upstream nodes)
reuseNodes = 1; % allow the same nodes (and synapses) to be analyzed for multiple tips
stopAtBranch = 0; % stop collecting neurite when you hit a branch
normalizeToLastSyn = 0; % set synapse (bip) nearest to tip as zero position
mustBeLong = 0; % neurite must be as long as capture length to be analyzed
useGroup = 0; %only use neurites containing particular combinations of synapses
%tipL = 10; % length of tip used to define synapse groups
filterGroups = 0;

useVgc = [2 3 4 5 13 14];

clear n
c = 0; %count tips
for i  = 1:length(useVgc)
    targ = find(runCids==useVgc(i));
    sm = sms(targ).sm;
    nep = sm.nep;
    swc = sm.nep.swcS;
    pred = swc.pred + 1;
    wPos = swc.pos;
    pos = nep.pos;
    closest = sm.syn2Skel.closest;

    check = pred * 0 + 1;
    uNodes = 1:length(pred);
    hNodes = hist(pred,uNodes);
    tips = find(hNodes == 0);
    branches = hNodes >1;


    scatter(wPos(:,1),wPos(:,2),'.','k')
    hold on
    scatter(wPos(tips,1),wPos(tips,2),'o','g','filled')
    hold off
    pause(.01)

    tips = tips(randperm(length(tips)));
    for t = 1:length(tips)
        %disp(sprintf('tip %d of %d',t,length(tips)))
        nodeList = tips(t);
        nodeL = 0;
        L = 0;
        last = nodeList;
        check(last) = 0;
        for p = 1:length(pred)

            next = pred(last);
            if next<1
                break
            end

            dist = sqrt((wPos(last,1)-wPos(next,1)).^2 + ...
                (wPos(last,2)-wPos(next,2)).^2 + ...
                (wPos(last,3)-wPos(next,3)).^2);
            L = L + dist;

            shouldBreak = 0;
            if~reuseNodes & ~check(next)
                shouldBreak = 1;
            end
            if stopAtBranch & branches(next)
                shouldBreak = 1;
            end
            if L>captureLength
                shouldBreak = 1;
            end

            if shouldBreak
                break
            else
                check(next) = 0;
                nodeList = [nodeList next];
                last = next;
                nodeL = [nodeL L];
            end

        end
        nodes = swc.swc2arborID(nodeList);

        synList = [];
        synL = [];
        for s = 1:length(nodes)
            syns = find(closest == nodes(s));
            synList = [synList syns'];
            synL = [synL syns' * 0 + nodeL(s)];
        end

        cid = sm.cid;
        isOut = sm.syn.pre(synList) == cid;
        synClasses = [sm.syn.preClass(synList) sm.syn.postClass(synList)];
        synIDs =  [sm.syn.pre(synList) sm.syn.post(synList)];
        partID = synIDs(sub2ind(size(synIDs),1:size(synIDs,1),isOut'+1));
        partClass = synClasses(sub2ind(size(synIDs),1:size(synIDs,1),isOut'+1));
        synType = isOut * 100 + partClass';
        synType(synType == 0) = 8; %assume unknown ins are amacrines


        %%Get branches
        isBranch = find(branches(nodeList));

        %%Record neurite
        c = c+1;
        n(c).L = L;
        n(c).nodeList = nodeList;
        n(c).cid = sm.cid;
        n(c).synList = synList;
        n(c).synL = synL;
        n(c).isOut = isOut;
        n(c).partID = partID;
        n(c).partClass = partClass;
        n(c).nodes = nodes;
        n(c).nodeL = nodeL;
        n(c).nodeLength = nep.props.nodeLength(n(c).nodes);
        n(c).synType = synType;
        n(c).branches = nodeList(isBranch);
        n(c).branchL = nodeL(isBranch);

    end
    scatter(wPos(:,1),wPos(:,2),'.','k')
    hold on
    scatter(wPos(check==0,1),wPos(check==0,2),'.','r','filled')
    hold off
    pause(1)

end



if normalizeToLastSyn
    for i = 1:length(n)
        isBip = find(n(i).synType == 7);
        if ~isempty(isBip)
            minL = min(n(i).synL(isBip));
            n(i).synL = n(i).synL - minL;
            n(i).nodeL = n(i).nodeL - minL;
            n(i).branchL = n(i).branchL - minL;
        end
    end
end


allTypes = cat(1,n(:).synType);
%uTypes = unique(allTypes);

uTypes = [7 8 101  108];
hTypes = hist(allTypes,uTypes);
bar(hTypes)
clear typeLab
for i = 1:length(uTypes)
    typeLab{i} = num2str(uTypes(i));
end
xticklabels(typeLab)


%% Analyze synapse combinations
chunks = [0 .5; .5 1.5; 2 5; 10 20];
removeNoSyn = 0;
%chunks = [0 2; 8 10];

chunkNum = size(chunks,1);
clf
clear sg sgLab

sg{1} = [7 8 101];
sg{2} = [7 101];
sg{3} = [7 8];
sg{4} = [8 101];
sg{5} = 7;
sg{6} = 8;
sg{7} = 101;
sg{8} = [];
for sgC = 1:length(sg)
    sgLab{sgC} = num2str(sg{sgC});
    if isempty(sg{sgC})
        sgLab{sgC} = 'none';
    end
end

% sgCol = hsv(length(sg));
bg = 0.3;
sgCol = [.7 .7 .7; bg 1 1; 1 1 bg; 1 bg 1; bg 1 bg; 1 bg bg; bg bg 1; 0 0 0];



for ch = 1:length(chunks)
    checkL = chunks(ch,2)-chunks(ch,1);

    nGroup = zeros(length(n),1);
    isMultiple = nGroup;
    clear longEnough hasBip hasRgc hasAmcIn
    for i = 1:length(n)

        hitS = (n(i).synL>=chunks(ch,1)) & (n(i).synL<chunks(ch,2) );
        hitT = n(i).synType(hitS);
        uT = unique(hitT);

        if length(hitT)>length(uT)
            isMultiple(i) = 1;
        end

        for g = 1:length(sg)
            intT = unique(intersect(sg{g},uT));
            intT = setdiff(intT,[100 108]);
            if length(intT) == length(sg{g})
                nGroup(i) = g;
                break
            end
        end

        longEnough(i) = n(i).L>=chunks(ch,2);
        hasBip(i) = sum(hitT==7);
        hasRgc(i) = sum(hitT==101);
        hasAmcIn(i) = sum(hitT==8);


    end

    fractionMultiple = mean(isMultiple)
    hTypes = hist(nGroup(longEnough),1:length(sg))
    subplot(1,chunkNum,ch)
    hitNum = sum((nGroup>0) & (nGroup<=7));
    disp(sprintf('%d synaptic neurites from %d neurites',hitNum,sum(longEnough)))

    pieD = hTypes(1:end-removeNoSyn);
    pieD = pieD/sum(pieD);
    pChart = pie(pieD,sgLab(1:end-removeNoSyn))
    c = 0;
    for p = 1:length(pChart)
        if strcmp(class(pChart(p)),    'matlab.graphics.primitive.Patch')
            c = c+1;
            pChart(p).FaceColor = sgCol(c,:);
        end
    end


    clear br

    %     br{1} =  longEnough & ~hasBip;
    %     br{2} = longEnough & hasBip ;

    br{1} =  longEnough & ~hasBip & ~hasRgc;
    br{2} = longEnough & hasBip & ~hasRgc;
    br{3} = longEnough & ~hasBip & hasRgc;
    br{4} = longEnough & hasBip & hasRgc;

    bipFrac = mean(longEnough & hasBip);
    rgcFrac = mean(longEnough & hasRgc);
    amcFrac = mean(longEnough & hasAmcIn);
    bipRgcFrac = mean(longEnough & hasBip & hasRgc);
    bipAmcFrac = mean(longEnough & hasBip & hasAmcIn);
    amcRgcFrac = mean(longEnough & hasAmcIn & hasRgc);
    predBipRgc = bipFrac * rgcFrac;
    predAmcRgcFrac = amcFrac * rgcFrac;
    predBipAmcFrac = bipFrac * amcFrac;

    if 0% bootstrap combined probabilities
        num = sum(longEnough);
        fReps = 1000;
        rBipRgc = zeros(fReps,1);
        rBipAmc = zeros(fReps,1);
        rAmcRgc = zeros(fReps,1);
        for r = 1:fReps
            rBip = rand(num,1)<=bipFrac;
            rRgc = rand(num,1)<=rgcFrac;
            rAmc = rand(num,1)<=amcFrac;
            rBipRgc(r) = mean(rBip & rRgc);
            rBipAmc(r) = mean( rBip & rAmc);
            rAmcRgc(r) = mean( rAmc & rRgc);

        end

        subplot(2,1,2)
        cla
        hold on
        scatter(zeros(fReps,1)+1,rBipRgc,'k')
        scatter(1,bipRgcFrac,'r','filled')
        scatter(zeros(fReps,1)+2,rBipAmc,'k')
        scatter(2,bipAmcFrac,'r','filled')
        scatter(zeros(fReps,1)+3,rAmcRgc,'k')
        scatter(3,amcRgcFrac,'r','filled')


        sBipRgc = sort(rBipRgc);
        disp('bip frac')
        [chunks(ch,1) chunks(ch,2)]
        [sBipRgc(fReps*.025) sBipRgc(fReps * 0.975)]
        bipRgcFrac


        pause




    end



    averageAMCDensity = mean(hasAmcIn(longEnough)/checkL);
    averageRGCDensity = mean(hasRgc(longEnough)/checkL);
    if 0 %show dependent densities of amacrine cells
        subplot(2,1,2)
        cla
        hold on
        plot([0 4],[averageAMCDensity averageAMCDensity],'g')
        for b = 1:length(br)
            brDensity = hasAmcIn(br{b})/checkL;
            meanBr(b) = mean(brDensity);
            seBr(b) = std(brDensity)/sqrt(length(brDensity));
            swarmchart(zeros(length(brDensity),1)+b,brDensity,'k')
            scatter(b,meanBr(b),'r','filled')

        end

        subplot(2,1,2)
        cla
        hold on
        plot([0 4],[averageRGCDensity averageRGCDensity],'g')
        for b = 1:length(br)
            brDensity = hasRgc(br{b})/checkL;
            meanBr(b) = mean(brDensity);
            seBr(b) = std(brDensity)/sqrt(length(brDensity));
            swarmchart(zeros(length(brDensity),1)+b,brDensity,'k')
            scatter(b,meanBr(b),'r','filled')

        end


    end

end

legend(sgLab)
pause(1)


%% Pick tips
useT = 1:length(n);
if useGroup
    isGroup = find(nGroup == useGroup);
    useT = intersect(useT,isGroup);
end
if mustBeLong
    La = [n(:).L];
    longEnough = find(La>=captureLength);
    useT = intersect(useT,longEnough);
end

%% quantify by bins
subplot(2,1,1)
bin = .5;
binS = [0:.1:captureLength ];
clear nB

for uC = 1:length(useT)
    c = useT(uC);
    clear bL bS bT bB
    for b = 1:length(binS)
        nodeHit = find((n(c).nodeL >= (binS(b)-bin/2)) & (n(c).nodeL < (binS(b) + bin/2)));
        sHit = find((n(c).synL >= (binS(b)-bin/2)) & (n(c).synL < (binS(b) + bin/2)));
        hitType = n(c).synType(sHit);
        hType = hist(hitType,uTypes);
        bT(:,b) = hType;
        bL(b) = sum(n(c).nodeLength(nodeHit));
        bS(b) = length(sHit);
        bB(b) = length(find((n(c).branchL >= binS(b)) & (n(c).branchL < (binS(b) + bin))));

    end
    nB(uC).bL = bL; %total length at that position
    nB(uC).bS = bS; %Total synapses number
    nB(uC).bT = bT; %Synapse number for each type
    nB(uC).bB = bB; %Branch points

end

bLa  = cat(1,nB(:).bL);
bSa = cat(1,nB(:).bS);
bTa = cat(3,nB.bT);
bDa = sum(bSa,1)./sum(bLa,1);
bBa = cat(1,nB(:).bB);
bNs = sum(bLa>0,1);

cla
typeCol = [0 .8 0; .8 0 0 ;0 0 .8; .8 0 .8];
hold on
for i = 1:length(uTypes)
    typeS = squeeze(bTa(i,:,:))';
    clear meanDs seDs
    for b = 1:size(typeS,2) %calculate mean and standard error for each
        useN = bLa(:,b) >0;
        typeDs = typeS(useN,b)./bLa(useN,b);
        meanDs(b) = mean(typeDs);
        seDs(b) = std(typeDs)/sqrt(length(typeDs));
    end

    seLow = meanDs - seDs;
    seHigh = meanDs + seDs;
    Y = [seLow fliplr(seHigh)];
    X = [binS fliplr(binS)]
    fill(X,Y,typeCol(i,:),'FaceAlpha', .2,'linestyle','none')
    plot(binS,meanDs,'color', typeCol(i,:))
end
%legend({'seBip';'bip';'seAmcIn';'amcIn';'seRgcOut';'rgcOut';'seAmcOut';'amcOut'})

%plot(binS,sum(bBa,1)./sum(bLa,1),'color','k')
%%ylim([0 .3])
pause(1)

str = sprintf('N = %d, reuse = %d, stopAtBranch = %d, mustBeLong = %d, filterGroups = %d',...
    length(useT),reuseNodes, stopAtBranch, mustBeLong, useGroup);
title(str)


%% quantify by chunks
Y = get(gca,'YLim');
for i = 1:size(chunks,1)
    plot([chunks(i,1) chunks(i,1)],Y,'k')
    plot([chunks(i,2) chunks(i,2)],Y,'k')
end


subplot(2,1,2)

clear nC

for uC = 1:length(useT)
    c = useT(uC);
    clear bT bL bS bB
    for b = 1:size(chunks,1)
        nodeHit = find((n(c).nodeL >= (chunks(b,1))) & (n(c).nodeL < (chunks(b,2))));
        sHit = find((n(c).synL >= (chunks(b,1))) & (n(c).synL < (chunks(b,2))));
        hitType = n(c).synType(sHit);
        hType = hist(hitType,uTypes);
        bT(:,b) = hType;
        bL(b) = sum(n(c).nodeLength(nodeHit));
        bS(b) = length(sHit);
        bB(b) = length(find((n(c).branchL >= binS(b)) & (n(c).branchL < (binS(b) + bin))));


    end
    nC(uC).bL = bL; %total length at that position
    nC(uC).bS = bS; %Total synapses number
    nC(uC).bT = bT; %Synapse number for each type
    nC(uC).bB = bB; %Branch points

end

cLa  = cat(1,nC(:).bL);
cSa = cat(1,nC(:).bS);
cTa = cat(3,nC.bT);
cDa = sum(cSa,1)./sum(cLa,1);
cBa = cat(1,nC(:).bB);
cNs = sum(cLa>0,1);

cla
hold on
clear sub
for s = 1:size(chunks,1);
    sub(s) = subplot(2,size(chunks,1),size(chunks,1) + s);
    ax(s) = gca;
    title(sprintf('chunk %0.02f to %0.02f',chunks(s,1),chunks(s,2)))
    hold on
end

pieDs = []
for i = 1:length(uTypes)
    typeS = squeeze(cTa(i,:,:))';
    clear meanDs seDs allDs
    for b = 1:size(typeS,2) %calculate mean and standard error for each
        useN = cLa(:,b) >0;
        typeDs = typeS(useN,b)./cLa(useN,b);
        meanDs(b) = mean(typeDs);
        seDs(b) = std(typeDs)/sqrt(length(typeDs));
        allDs{b} = typeDs;

    end

    pieDs(i,:) = meanDs;

    seLow = meanDs - seDs;
    seHigh = meanDs + seDs;

    for s = 1:length(sub)
        Ds = allDs{s};
        %bar(sub(s),i,meanDs(s),'facecolor',typeCol(i,:))
        swarmchart(ax(s),zeros(length(Ds),1)+i,allDs{s},'markerfacealpha',.2,...
            'markerfacecolor',typeCol(i,:),'markeredgecolor','none')
        scatter(ax(s),i,mean(Ds),300,'+','k')
        %ylim(sub(s),[0 .5])
    end

    %     Y = [seLow fliplr(seHigh)];
    %     X = [1:size(chunks,1) size(chunks,1):-1:1]
    %     fill(X,Y,typeCol(i,:),'FaceAlpha', .2,'linestyle','none')
    %     plot(1:size(chunks,1) ,meanDs,'color', typeCol(i,:))
end
%legend({'seBip';'bip';'seAmcIn';'amcIn';'seRgcOut';'rgcOut';'seAmcOut';'amcOut'})

pause(1)

for i = 1:size(pieDs,2)
    hold(ax(i),'off')
    pieD = pieDs(:,i);
    pieD = pieD/sum(pieD);
    pChart = pie(ax(i),pieD,typeLab)
    c = 0;
    for p = 1:length(pChart)
        if strcmp(class(pChart(p)),    'matlab.graphics.primitive.Patch')
            c = c+1;
            pChart(p).FaceColor = typeCol(c,:);
        end
    end

end

return

%% display
subplot(2,1,2)

if 0
    Ls = [n(:).L];
    isLong = find(Ls>=captureLength);
    isLong = isLong(1:min(100,length(isLong)));
    useT = isLong;
elseif 0
    useT = [];
    for t = 1:length(n)
        isType = find(n(t).synType);
        typeL = n(t).synL(isType);
        if sum(typeL)
            useT = [useT t];
        end
    end
end

showT = useT(1:min(100,length(useT)));

clear scatSynType
scatSynType{length(uTypes)} = [];
scatSynTip = scatSynType;
scatBranchT = [];
scatBranchL = [];
allNodeL = [];
allNodeT = [];
for t = 1:length(showT)
    tip = showT(t);
    for i = 1:length(uTypes)
        isType = find(n(tip).synType == uTypes(i))';
        scatSynType{i} = [scatSynType{i}  n(tip).synL(isType)];
        scatSynTip{i} = [scatSynTip{i} isType*0+t];
    end
    scatBranchL = [scatBranchL n(tip).branchL];
    scatBranchT = [scatBranchT n(tip).branchL * 0 + t];
    allNodeL = [allNodeL n(tip).nodeL];
    allNodeT = [allNodeT n(tip).nodeL * 0 + t];
end

cla
scatter(allNodeT, allNodeL,'|','k')
hold on
%scatter(scatBranchT, scatBranchL,'x')
hold on
for i = 1:length(uTypes)
    scatter(scatSynTip{i},scatSynType{i},10,'o','markerfacecolor',typeCol(i,:),...
        'markeredgealpha',0)
end



a = gca;
% a.XGrid = 'on';
% a.XTick = [1:length(isLong)];
% a.YLim = [-0.1 captureLength];
a.XTickLabel = []




%% Print figure
if 0
    fDir = uigetdir;
    filename = [fDir '\tipProfile2mComboPiesNoSynOneUm2']
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])

end












