function[] = SignalExtraction_jmAutoCorr()

global tis glob
SPN = [glob.datDir 'Analysis\Data\preproc\'];

clear R

KdatDir = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\';
clear RawOI

%%
% Pre-defined parameters
Topic = 'Ai148_129SVG3';
ListId = '122618';
Fz = 9.47;
UpSamplingFz = 100;
FigureFolder = 'CorrelatedEM-122618';
MaskResource = 'EM';% 'EM' 'AutoReg'
AmpEstMethod = 'mean';
% addpath('../../Multitude/Analysis');
% addpath('../../NaturalSceneVideo/Analysis');
% addpath('./Functions');
%%
IsBidirectional = 0;
upFz = UpSamplingFz; % in Hz
BaseT = 0.2; % in second
DriftT = 3; % in second

%%
StimulusType = 'FLBar';
MaskPath = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\Mask';
ROIPath = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\';

%% Load file list
FilNam = sprintf('%sROIs/ListRecordingFiles_%s.m',ROIPath,ListId);
run(FilNam);
NumFile = length(R);

StimType = 'FlashBar';
clear ROI maxResp PixResp
ROI.Polarity = [];
ROI.Polarity2 = [];
ROI.SNR      = [];
clear meanIs
for r = 1:NumFile
    Exps = R(r).Exps{1};
    if isempty(Exps), continue; end
    clc;
    fprintf('Flashbar... %d/%d \n', r, NumFile);

    for i = 1:length(Exps)
        close all;
        Day = R(r).Day;
        Cel = R(r).Cel;
        Exp = Exps(i);
        FileStamp = sprintf('%s_%d%03d', Day, Cel, Exp);
        fprintf('\t Exp %s... %d/%d \n', FileStamp, i, length(Exps));


        % Load qualified pixel mask
        %         MaskFileName =  sprintf('%s_AutoRregMask_%s_%d%03d.mat',Topic, Day, Cel, Exp);
        %         load([MaskPath '\' MaskFileName], 'Mask');
        %         Mask = Mask';

        % Load calcium signal after image registration
        FilNam = sprintf('%s_Translation_%s_%d%03d.mat',Topic, Day, Cel, Exp );
        load([KdatDir 'Translation/' FilNam], 'I');
        %         I = permute(I, [2 1 3]);
        I = medfilt3(I, [3, 3, 1]);
        meanIs{r,Exp} = mean(I,3);

        % Get correct timeline due to directional scanning
        TIds = reshape(1:numel(I), size(I));
        if IsBidirectional
            TIds = BidirecitonTimeCorrect3D(TIds);
        end

        % Load Time Table that aligned visual and scanning
        AlignName = sprintf('%s_Alignment_%s_%d%03d.mat',Topic, Day, Cel, Exp );
        load([KdatDir 'Alignment/' AlignName], 'TimeTable')

        % Loop for each qualified pixel
        % - x: recording Id (from the first pixel)
        % - v: pixel value
        % - xq: time Id and a fixed upsampling points
        % - vq (result): upsample epoching
        % (1) Spot (2) DGSpot


        %% Get the ROI signals
        switch lower(MaskResource)
            case 'em'
                ROITable = load([SPN 'ptDat.mat']);
                ROITable = ROITable.ptDat;
                ROIMask = load([SPN 'maskDat.mat']);
                ROIMask = ROIMask.maskDat;
                useRois = find(ROITable(:, 1) == (Cel*1000 + Exp));
                ROIMask = ROIMask(:, :, useRois);
            case 'autoreg'
                MaskName = sprintf('%s_AutoRregMask_%s_%d%03d.mat',Topic, Day, Cel, Exp );
                ROIMask = load([KdatDir 'Mask/' MaskName], 'Mask');
                ROIMask = ROIMask.Mask;
        end
        nROI = size(ROIMask, 3);


        Mask = ones(size(I, 1), size(I, 2));
        %% Create a 4D array for traces
        % [Repeats, Time, Condition, ROI]
        switch lower(StimulusType)
            case 'flbar'
                pixFz = mean(diff(TimeTable(:, 5))./diff(TimeTable(:, 3)));
                DispSize = unique(TimeTable(:, 1));
                SeqId = 2;
                TimeId = 5;
                nRepeat = max(TimeTable(:, 4));

                TLDrift = -BaseT:1/upFz:DriftT;
                PixResp = nan(nRepeat, length(TLDrift), length(DispSize), nROI);
                %% 1. Repeat, 2. Time, 3) Location, 4 ROI
                for e = 1:length(DispSize)
                    EpochId = TimeTable(:, SeqId) == 1 & TimeTable(:, 1) == DispSize(e);
                    for q = 1:nROI
                        Mask = squeeze(ROIMask(:, :, q));
                        [rId, cId, ~] = ind2sub(size(Mask), find(Mask));
                        cPixResp = nan(nRepeat, length(TLDrift), 1, length(rId));
                        for p = 1:length(rId)
                            x = squeeze(TIds(rId(p), cId(p), :));
                            v = squeeze(I(rId(p), cId(p), :));
                            xq = TimeTable(EpochId, TimeId)+round(pixFz*TLDrift);
                            vq = interp1(x, v, xq, 'linear');
                            cPixResp(:, :, 1, p) = vq;
                        end
                        PixResp(:, :, e, q) = nanmean(cPixResp, 4);
                    end
                end
            otherwise
                error('No such stimulus type');
        end

        %% Check Drift response overall
        close all
        %         figure('visible','off')
        switch lower(StimulusType)
            case 'flbar'
                for e = 1:length(DispSize)
                    subplot(2, ceil(length(DispSize)/2), e);
                    for q = 1:nRepeat
                        Trace = squeeze(nanmean(PixResp(q, :, e, :), 4));
                        Trace = Trace - mean(Trace(TLDrift<0));
                        plot(TLDrift, Trace); hold on
                    end
                    Trace = squeeze(nanmean(mean(PixResp(:, :, e, :), 4), 1));
                    Trace = Trace - mean(Trace(TLDrift<0));
                    plot(TLDrift, Trace, 'k', 'LineWidth', 2); hold on
                    xlabel('Time (s)');
                    ylabel('Responses');
                    title(sprintf('Display position: %d (um)', DispSize(e)));
                end
            otherwise
                error('No such stimulus type');
        end

        %% view individual exps
        if 0
            clf
            for r = 1:nROI
                for e = 1:length(DispSize) %% run each position
                    Traces = squeeze(PixResp(:,:,e,r));
                    maxResp(:,e,r) = max(Traces,[],1);
                    meanResp(:,e,r) = mean(Traces,1);
                    medResp(:,e,r) = median(Traces,1);
                    clf
                    plot(mean(Traces,1),'linewidth',3,'color','g')
                    hold on
                    title(sprintf('position %d, roi %d',e,r))
                    plot(median(Traces,1),'linewidth',3,'color','r')
                    plot(max(Traces,[],1),'linewidth',3,'color','c')
                    plot(Traces','k')
                    hold off

                    pause
                end
            end
        end

        %% Estimate the calcium response value
        % normalization by individual ROI
        PixRespNrm = PixResp./reshape(std(reshape(PixResp, [], nROI), [], 1), 1, 1, 1, nROI);
        % quick check
        figure; subplot(2, 3, 1);
        imagesc(squeeze(mean(mean(PixResp, 1), 2)));colorbar
        subplot(2, 3, 4);
        imagesc(squeeze(mean(mean(PixRespNrm, 1), 2)));colorbar
        % baseline substraction
        PixRespNrm = PixRespNrm-mean(PixRespNrm(:, TLDrift < 0, :, :), 2);
        mPixRespNrm = squeeze(median(PixRespNrm, 1));
        % Amplitude estimation
        Amplitude = squeeze(quantile(mPixRespNrm(TLDrift>0, :, :), 0.75, 1));
        wAmplitude = Amplitude./sum(Amplitude, 1);
        subplot(2, 3, 3);
        imagesc(wAmplitude);colorbar
        % polarity estimation
        TON  = TLDrift>0.05 & TLDrift<1.2;
        TOFF = TLDrift>1.55 & TLDrift<2.7;

        %%Polarity 1
        Polarity = mean(mPixRespNrm(TOFF, :, :) - mPixRespNrm(TON, :, :), 1);
        Polarity = squeeze(Polarity ./ mean(mPixRespNrm(TOFF, :, :) + mPixRespNrm(TON, :, :) + eps, 1));
        Polarity(Polarity>1) = 1;
        Polarity(Polarity<-1) = -1;
        subplot(2, 3, 2);
        imagesc(Polarity);colorbar

        %%Polarity 2
        Polarity2 = mean(mPixRespNrm(TOFF, :, :) , 1);
        Polarity2 = squeeze(Polarity2 ./ mean(mPixRespNrm(TOFF, :, :) + mPixRespNrm(TON, :, :) + eps, 1));
        Polarity2(Polarity2>1) = 1;
        Polarity2(Polarity2<0) = 0;
        %         subplot(2, 3, 2);
        %         imagesc(Polarity2);colorbar



        %         minVal = min(mPixRespNrm(:));
        %         maxVal = max(mPixRespNrm(:));
        %         onMean = squeeze(mean(mPixRespNrm(TON, :, :)))-minVal+maxVal/100;
        %         offMean = squeeze(mean(mPixRespNrm(TOFF,:,:)))-minVal+maxVal/100;
        %         offMean(offMean<0) = 0;
        %         pPolarity = offMean./(onMean + offMean);
        %         Polarity2 = mean(pPolarity,1);
        %         subplot(2, 3, 3);
        %         imagesc(Polarity2);colorbar


        % SNR estimation
        nRept = size(PixRespNrm, 1);
        Pairs = perms(1:nRept);
        Pairs = unique(Pairs(:, 1:2), 'rows');
        nPair = size(Pairs, 1);
        nCond = size(PixRespNrm, 3);
        SNRM = nan(nCond, nROI, nPair);
        for p = 1:nPair
            for c = 1:nCond
                for q = 1:nROI
                    CorrR = corrcoef(reshape(PixRespNrm(Pairs(p, 1), TLDrift>0, c, q), [], 1),...
                        reshape(PixRespNrm(Pairs(p, 2), TLDrift>0, c, q), [], 1));
                    SNRM(c, q, p) = CorrR(1, 2);
                end
            end
        end
        SNRM = median(SNRM, 3).^2;
        subplot(2, 3, 5);
        imagesc(SNRM);colorbar;

        ROI.Polarity = [ROI.Polarity; mean(Polarity'*wAmplitude)'];
        ROI.Polarity2 = [ROI.Polarity2; mean(Polarity2' * wAmplitude)'];
        ROI.SNR = [ROI.SNR; median(SNRM, 1)'];

        %% record raw values
        RawOI.PixResp{r,i} = PixResp;
        RawOI.useRois{r,i} = useRois;
        rawPixResp = zeros(nROI,size(I,3));
        for q = 1:nROI
            Mask = squeeze(ROIMask(:, :, q));
            [rId, cId, ~] = ind2sub(size(Mask), find(Mask));
            rawPixRec = zeros(length(rId),size(I,3));
            for p = 1:length(rId)
                v = squeeze(I(rId(p), cId(p), :));
                rawPixRec(p,:) = v;
            end
            rawPixResp(q,:) = mean(rawPixRec,1);
        end
        RawOI.rawPixResp{r,i} = rawPixResp;


        clf
        plot(rawPixResp')
        pause(1)

    end
end
ROI.ROITable = [ROITable, ROI.Polarity, ROI.SNR];
%save('../individualPointMasksSpread\ROI.mat','ROI')
RawOI.ROI = ROI;
RawOI.R = R;
RawOI.ROIMask = ROIMask;
RawOI.meanIs = meanIs;

global glob
try
    SPN =  [glob.datDir 'Analysis\Data\preproc\'];
    save([SPN 'RawOI.mat'],'RawOI');
end




