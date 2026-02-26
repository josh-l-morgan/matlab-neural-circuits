function[] = makeSOI()

global tis glob


polFig = figure;
showCells = 0;
lc = 16;

%% Get functiona data

datFold = [glob.datDir 'Analysis\Data\preproc\'];
SPN = datFold;%'\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasksSpread\'


%load([SPN 'maskDat.mat'])
load([SPN 'ptDat.mat']);
load([SPN 'ROI2.mat']);
load([SPN 'ROI.mat']);
load([SPN 'MOI.mat']);


load([datFold 'COI.mat']);
clear SOI;

jetMap = jet(100);

posCal = ptDat(:,[7 6 8]);
posCal(:,1) = posCal(:,1) * 0.004;
posCal(:,2) = posCal(:,2) * 0.004;
posCal(:,3) = posCal(:,3) * .04;
cidCal = ptDat(:,3);
vgcCids = unique(cidCal);
pol1 = ROI2.Polarity;
polX = ROI2.PolarityX;
pol2 = mean(ROI2.Polarity2,2);
onMean = mean(ROI2.ON,2);
offMean = mean(ROI2.OFF,2);
onMax = ROI2.maxON;
offMax = ROI2.maxOFF;
pol3 = ROI2.Polarity3;
pol4 = ROI2.Polarity4;

roiNum = size(ptDat,1);

%% Show different polarities
polRange = [-1.5:.1:1.5];
subplot(5,1,1)
hist(pol1,polRange)
xlim([-1.5 1.5])
subplot(5,1,2)
hist(polX,polRange)
xlim([-1.5 1.5])
subplot(5,1,3)
hist(pol2,polRange)
xlim([-1.5 1.5])
subplot(5,1,4)
hist(pol3,polRange)
xlim([-1.5 1.5])
subplot(5,1,5)
hist(pol4,polRange)
xlim([-1.5 1.5])

%%
subplot(4,1,3)
hist(onMean)
subplot(4,1,4)
hist(offMean)

drawnow

%%Show positions
showPol1 = round((pol1 +1) * 50);
showPol2 = round((pol2 + 1) * 50);
showPol3 = round(pol3 * 100);
showOnMean = round(onMean * 100);
showOnMax = round(onMax * 50);
showOffMean = round(round(offMean * 100));
showOffMax = round(offMax * 50);

showOffMean(showOffMean<1) = 1;
showOffMean(showOffMean>100) = 100;
showOnMean(showOnMean<1) = 1;
showOnMean(showOnMean>100) = 100;
showOffMax(showOffMax<1) = 1;
showOffMax(showOffMax>100) = 100;
showOnMax(showOnMax<1) = 1;
showOnMax(showOnMax>100) = 100;
showPol1(showPol1<1) = 1;
showPol1(showPol1>100) = 100;
showPol3(showPol1<1) = 1;
showPol3(showPol1>100) = 100;

prop = showOnMean;

clf
subplot(1,1,1)
allScat = scatter3(posCal(:,1),posCal(:,2),posCal(:,3),100,'.','k')
ax = gca
ax.Color = [0 0 0];
% prop = prop - min(prop);
% prop = prop * 100/max(prop);
colormap jet(100)

allScat.CData = prop;




%% Sort anatomical data
offTags = {'bc1' 'bc2' 'bc3a' 'bc3b' 'bc4' 'boff' 'bcoff'};
onTags =  {'bc5i' 'bc5o' 'bc5t' 'bc6' 'bc7' 'bc8' 'bc11' 'bon' '5' 'bcon' 'bc5' 'rbc' 'xbc'};
unkTags = {'bcunk' 'bunk' };

testTag = {unkTags offTags onTags}; %off == 1, on == 2

subTypeNames = tis.cells.type.subTypeNames{7};
unkOffOn = zeros(length(subTypeNames),1);
for i = 1:length(subTypeNames)

    nam = subTypeNames{i};
    for s = 1:length(testTag)
        for t = 1:length(testTag{s})
            if strcmp(testTag{s}{t},nam)
                unkOffOn(i) = s-1;
                break
            end
        end
    end
end



allR = [];
for i = 1:length(COI.rgcGroupCids)
    checkCids = COI.rgcGroupCids{i};
    allR = cat(1,allR,checkCids);
end

allB = [];
for i = 1:length(COI.bpcGroupCids)
    checkCids = COI.bpcGroupCids{i};
    allB = cat(1,allB,checkCids);
end


% 
% isAMC = find(tis.cells.type.typeID == 8);
% isVGC = intersect(isAMC,find(tis.cells.type.subTypeID == 1));
% vgcCids = tis.cids(isVGC)

%%load SM

smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];



for bG = 1:length(COI.bpcGroupCids)
    bCids = COI.bpcGroupCids{bG};
end

%%  Plot cells


SOI.cids = vgcCids;
SOI.closeNode= zeros(roiNum,1);
SOI.on = zeros(roiNum,1);
SOI.off = zeros(roiNum,1);
SOI.offBias = zeros(roiNum,1);
SOI.pos =zeros(roiNum,3);
numRoi = roiNum
allD = zeros(numRoi)+inf;

NOI.cids = vgcCids;




for v  = 1:length(vgcCids)

    disp(sprintf('running cell %d. %d of %d',vgcCids(v),v,length(vgcCids)))
    useCal = find(cidCal == vgcCids(v));


    if showCells
        figSkel{i} = figure;
    end
    fileName = sprintf('sm_cid%d.mat',vgcCids(v));
    if exist([smDir fileName],'file')
        useSM(v) = 1;
        if exist('sms','var')
            sm = sms(v).sm;
        else
            load([smDir fileName]);
        end

    

        d = sm.syn2Skel.syn2SkelDist;
        skel2skel = sm.skel2skel.linDist;

        pre = sm.syn.pre;
        post = sm.syn.post;
        preSign = zeros(length(pre),1);

        NOI.neps(v).nep = sm.nep;
        NOI.syns(v).syn = sm.syn;
        NOI.ds(v).d = d;
        NOI.preSigns(v).preSign = preSign;

        if 1;%~isempty(pre)


            %% get all synapse to node distances

            W = exp(-d/lc); % Apply length constant

            for p = 1:length(pre)
                targ = find(tis.cids == pre(p),1);
                if tis.cells.type.typeID(targ) == 7;
                    subT = tis.cells.type.subTypeID(targ);
                    if subT
                        preSign(p,1) = unkOffOn(subT);
                    end
                elseif tis.cells.type.typeID(targ) == 8
                    preSign(p,1) = -1;
                end
            end
            NOI.preSigns(v).preSign = preSign;

            pos = sm.nep.pos;
            Won = W .* repmat(preSign==1,[1 size(W,2)]);
            Woff = W .* repmat(preSign==2,[1 size(W,2)]);
            Winhib = W.* repmat(preSign==-1,[1 size(W,2)]);

            synPos = sm.syn.pos;
            onSynPos = synPos(preSign == 2,:);
            offSynPos = synPos(preSign == 1,:);
            inhibSynPos = synPos(preSign == -1,:);

            sumOn = sum(Won,1);
            sumOff = sum(Woff,1);
            sumInhib = sum(Winhib,1);

            colormap(jet(100))
            sumOnN = sumOn * 50/mean(sumOn(:));
            sumOffN = sumOff * 50/mean(sumOff(:));
            sumInhibN = sumInhib * 50/mean(sumInhib(:));

            sumOnN = sumOn * 50/mean(sumOn(:));

            maxI = max([sumOn(:); sumOff(:)]);
            col2 = [sumOff(:)/max(sumOff(:))  sumOn(:) * 0 sumOn(:)/max(sumOn(:))];

            offBias = (sumOff-sumOn)./(sumOn+sumOff);
            onBias = (sumOn - sumOff)./(sumOn+sumOff);

            pos = sm.nep.pos;
            markerSize = 150;


            %% Match lightEM
            clear nearCV
            for c = 1:length(useCal)
                cDists = sqrt((pos(:,1) - posCal(useCal(c),1)).^2 + (pos(:,2) - posCal(useCal(c),2)).^2  + ...
                    (pos(:,3) - posCal(useCal(c),3)).^2 );
                minDist = min(cDists);
                nearCV(c) = find(cDists == minDist,1);

            end

            SOI.cell(v).skelPos = pos;
            SOI.cell(v).syn = sm.syn;
            SOI.cell(v).d = d;
            SOI.cell(v).onBias = onBias;
            SOI.cell(v).sumOn = sumOn;
            SOI.cell(v).sumOff = sumOff;
            %SOI.cell(v).skel2skel = skel2skel;
            SOI.cell(v).preSign = preSign;

            SOI.closeNode(useCal) = nearCV;
            SOI.on(useCal) = sumOn(nearCV);
            SOI.off(useCal) = sumOff(nearCV);
            SOI.offBias(useCal) = offBias(nearCV);
            SOI.pos(useCal,:) = pos(nearCV,:);
            SOI.cid(useCal) = vgcCids(v);

            closeNodes = SOI.closeNode(useCal);
            dAll = sm.skel2skel.linDist;
            dRoi = dAll(closeNodes,closeNodes);
            allD(useCal,useCal) = dRoi;




            if showCells

                %% Plot
                subplot(2,2,1)
                scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),'.','clipping','off');
                axis equal
                view([ 0 90])
                scatSkel.CData = ceil(sumOnN);
                hold on
                scatCalOff = scatter3(posCal(useCal,1),posCal(useCal,2),posCal(useCal,3),100,'o','filled')
                scatCalOff.CData = showOffMean(useCal);
                %         hold on
                %         scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),300,'k','+','linewidth',2,'clipping','off');
                %         scatOff = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),300,'k','linewidth',2,'clipping','off');
                %         hold off
                ylim([min(pos(:,2)) max(pos(:,2))]);
                xlim([min(pos(:,1)) max(pos(:,1))]);
                zlim([min(pos(:,3)) max(pos(:,3))]);
                title('OFF influence')

                subplot(2,2,2)
                scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),'.','clipping','off');
                axis equal
                view([ 0 90])
                scatSkel.CData = ceil(sumOffN);
                hold on
                scatCalOn = scatter3(posCal(useCal,1),posCal(useCal,2),posCal(useCal,3),100,'o','filled')
                scatCalOn.CData = showOnMean(useCal);
                %         hold on
                %         scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),300,'k','+','linewidth',2,'clipping','off');
                %         scatOff = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),300,'k','linewidth',2,'clipping','off');
                %         hold off
                ylim([min(pos(:,2)) max(pos(:,2))]);
                xlim([min(pos(:,1)) max(pos(:,1))]);
                zlim([min(pos(:,3)) max(pos(:,3))]);
                title('ON influence')


                subplot(2,2,3)
                scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),'.','clipping','off');
                axis equal
                view([ 0 90])
                showOffBias = offBias;
                showOffBias = (showOffBias + 1) * 255;
                scatSkel.CData = ceil(showOffBias);
                hold on
                scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),markerSize,'k','+','linewidth',2,'clipping','off');
                scatOff = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),markerSize,'k','linewidth',2,'clipping','off');
                hold on
                scatCalOff = scatter3(posCal(useCal,1),posCal(useCal,2),posCal(useCal,3),200,'o','filled')
                scatCalOff.CData = showPol1(useCal);
                hold off
                ylim([min(pos(:,2)) max(pos(:,2))]);
                xlim([min(pos(:,1)) max(pos(:,1))]);
                zlim([min(pos(:,3)) max(pos(:,3))]);
                title('Polarity')

                subplot(2,2,4)
                scatSkel4 = scatter3(pos(:,1),pos(:,2),pos(:,3),'.','clipping','off');
                axis equal
                view([ 0 90])
                scatSkel4.CData = ceil(sumInhibN);
                %         hold on
                %         scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),300,'k','+','linewidth',2,'clipping','off');
                %         scatOff = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),300,'k','linewidth',2,'clipping','off');
                %         hold off
                ylim([min(pos(:,2)) max(pos(:,2))]);
                xlim([min(pos(:,1)) max(pos(:,1))]);
                zlim([min(pos(:,3)) max(pos(:,3))]);
                title('Amacrine input')

                pause(1)


            end



        end




        %
        %         subTypePre = zeros(length(pre),1);
        %         for p = 1:length(pre)
        %             targ = find(tis.cids == pre(p),1);
        %             if tis.cells.type.typeID(targ) == 7;
        %                 subTypePre(p) = tis.cells.type.subTypeID(targ);
        %             end
        %         end
        %
        %         preBipGroup = zeros(length(pre),1);
        %           for bG = 1:length(COI.bpcGroupCids)
        %                 bCids = COI.bpcGroupCids{bG};
        %
        %           end


        else
            disp('pre is empty')


            NOI.neps(v).nep = [];
            NOI.syns(v).syn = [];
            NOI.ds(v).d = [];
            NOI.preSigns(v).preSign = [];

    end
end

pos = SOI.pos;
allE = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2 + ...
    (pos(:,3)-pos(:,3)').^2);
SOI.r2rLin = allD;
SOI.r2rEuc = allE;

save([datFold 'NOI.mat'],'NOI');
save([datFold 'SOI.mat'],'SOI','-v7.3');
try close polFig, end









