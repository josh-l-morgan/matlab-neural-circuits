global glob
SPN = [glob.datDir 'Analysis\Data\preproc\'];
KdatDir = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\';
f = figure
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

                ROITable = load([SPN 'ptDat.mat']);
ROITable = ROITable.ptDat;
                ROIMasks = load([SPN 'maskDat.mat']);
ROIMasks = ROIMasks.maskDat;

%% Load file list
FilNam = sprintf('%sROIs/ListRecordingFiles_%s.m',ROIPath,ListId);
run(FilNam);
NumFile = length(R);

StimType = 'FlashBar';
clear ROI
ROI.Polarity = [];
ROI.SNR      = [];
ROI.ON = [];
ROI.OFF = [];
ROI.roiIDs = [];
ROI.Polarity2 = [];
ROI.maxON = [];
ROI.maxOFF = [];
ROI.Polarity3 = [];
ROI.Polarity4 = [];
ROI.Polarity5 = [];


for r = 1:NumFile
    Exps = R(r).Exps{1};
    if isempty(Exps), continue; end
    fprintf('Flashbar... %d/%d \n', r, NumFile);
    
    for i = 1:length(Exps)
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
        rawMean = mean(I(:));
        rawMode = mode(I(:));
        targetMean = 5;
        
        if 1 %% running std scaling
            disp('normalizing image stack')
            clf, colormap jet(100)
            binI = 3;
            colormap gray(256)
            image(squeeze(sum(I,1)))
            for s = 1:size(I,3) %scale I
                %!!!!!!!!
                startI = max(1,s-binI);
                stopI = min(size(I,3),s+binI);
                Ib = I(:,:,setdiff(startI:stopI,s));
                Ip = I(:,:,s);
                Imode = mode(Ib(:));
                Ib = Ib-Imode;
                Ip = Ip-Imode;
                %I(:,:,s) = Ip/std(Ib(:))*10;
                I(:,:,s) = Ip * targetMean/mean(Ib(:));%std(Ib(:))*10;

                if ~mod(s,10), s,
                    subplot(2,1,1)
                    image(I(:,:,s)*10)
                    subplot(2,1,2)
                    image(squeeze(sum(I,1))), drawnow;
                end
            end
        elseif 1 %% scale whole image by std
            
            I = I - mode(I(:));
            I = I/std(I(:));
            
        end
        
        
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
        grabROIs = find(ROITable(:, 1) == (Cel*1000 + Exp));
        ROIMask = ROIMasks(:, :, grabROIs);
        nROI = size(ROIMask, 3);
        colPick = randperm(nROI);
        %%Display mask
        roiCol = hsv(nROI);
        IcRoi = zeros(size(ROIMask,1),size(ROIMask,2),3);
        for m = 1:size(ROIMask,3)
            for c = 1:3
                IcRoi(:,:,c) = IcRoi(:,:,c) + ROIMask(:,:,m) * roiCol(colPick(m),c);
            end
        end
        clf
        image(uint8(IcRoi*256))
        imshow(IcRoi)


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
                timeMask = I * 0;
                
                for q = 1:nROI
                    Mask = squeeze(ROIMask(:, :, q));
                    %%JM signal extraction
                    maskInd = find(Mask);
                    cutI = I.* repmat(Mask,[1 1 size(I,3)]);
                    meanI = squeeze(sum(sum(cutI,1),2)/length(maskInd));
                    
                    %                         fun = @(x) max(x(:));
                    %                         meanIfilt = nlfilter(meanI,[3 1],fun);
                    if 0 %dilate
                        SE = strel('rectangle',[1 1]);
                        meanIfilt = imdilate(meanI,SE);
                    else
                        meanIfilt = meanI;
                    end
                    
                    %% Mexican Hat
                    if 0
                        sigma1 = 1;
                        sz = 300;    % length of gaussFilter vector
                        x = linspace(-sz / 2, sz / 2, sz);
                        gaussFilter = exp(-x .^ 2 / (2 * sigma1 ^ 2));
                        gaussFilter1 = gaussFilter / sum (gaussFilter); % normalize
                        sigma2 = 10;
                        gaussFilter = exp(-x .^ 2 / (2 * sigma2 ^ 2));
                        gaussFilter2 = gaussFilter / sum (gaussFilter); % normalize
                        mexHat = gaussFilter1 - gaussFilter2;

                        buff = length(gaussFilter);
                        buffMeanI = zeros(buff * 2 + length(meanI),1);
                        buffMeanI(buff + 1: buff  + length(meanI)) = meanI;

                        meanIfilt = filter(gaussFilter1,1,buffMeanI);
                        meanIfilt = meanIfilt(buff * 2 + 1:buff  + length(meanI));
                        
                        plot(meanI)
                        hold on
                        plot(meanIfilt)
                        hold off
                    else
                        meanIfilt = meanI;
                    end
                    
                    for e = 1:length(DispSize)
                        EpochId = TimeTable(:, SeqId) == 1 & TimeTable(:, 1) == DispSize(e);
                        
                        [rId, cId, ~] = ind2sub(size(Mask), find(Mask));
                        p = 1;
                        x = squeeze(TIds(rId(p), cId(p), :));
                        v = meanIfilt;%squeeze(I(rId(p), cId(p), :));
                        xq = TimeTable(EpochId, TimeId)+round(pixFz*TLDrift);
                        vq = interp1(x, v, xq, 'linear');
                        PixResp(:, :, e, q) = vq;
                        
                    end
                end
            otherwise
                error('No such stimulus type');
        end
        
        %% Check Drift response overall
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
        pause(.1)
        
        if 0 % Emily normalization
            %% Estimate the calcium response value
            % normalization by individual ROI
            stdRoi = std(reshape(PixResp, [], nROI), [], 1);
            stdRoi = reshape(stdRoi, 1, 1, 1, nROI);
            %stdRoi = repmat(stdRoi,[size(PixResp,1) size(PixResp,2) size(PixResp,3) 1]);
            PixRespNrm = PixResp./stdRoi;
            
            
            % quick check
            subplot(2, 3, 1);
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
            
            
            
            
        else % JM normalization
            
            
            if 1  %Show individual traces
                clf
                for f = 1:nROI
                    for e = 1:length(DispSize) %% run each position
                        Traces = squeeze(PixResp(:,:,e,f));
                        maxResp(:,e,f) = max(Traces,[],1);
                        meanResp(:,e,f) = mean(Traces,1);
                        medResp(:,e,f) = median(Traces,1);
                        clf
                        plot(mean(Traces,1),'linewidth',3,'color','g')
                        hold on
                        title(sprintf('experiment %d, position %d, roi %d',Exps(i),e,f))
                        plot(median(Traces,1),'linewidth',3,'color','r')
                        plot(max(Traces,[],1),'linewidth',3,'color','c')
                        plot(Traces','k')
                        hold off
                        
                        pause(.01)
                    end
                end
            end
            
            %% Show each display position position
            clf
            for e = 1:length(DispSize) %% run each position
                subplot( ceil(length(DispSize)),1, e);
                for q = 1:nRepeat  %% run each repeat
                    Trace = squeeze(nanmean(PixResp(q, :, e, :), 4));
                    Trace = Trace - mean(Trace(TLDrift<0));
                    plot(TLDrift, Trace,'k'); hold on
                    ylim([0 70])
                end
                Trace = squeeze(nanmean(mean(PixResp(:, :, e, :), 4), 1));
                Trace = Trace - mean(Trace(TLDrift<0));
                plot(TLDrift, Trace, 'r', 'LineWidth', 1); hold on
                xlabel('Time (s)');
                ylabel('Responses');
                title(sprintf('Display position: %d (um)', DispSize(e)));
                pause
            end

        %% Show each ROI position
            clf
            roiNumForFrame = size(PixResp,4);
            for e = 1: roiNumForFrame%% run each position
                subplot( ceil(roiNumForFrame/4),4, e);
                for q = 1:nRepeat  %% run each repeat
                    Trace = squeeze(nanmean(PixResp(q, :, :, e), 3));
                   % Trace = squeeze(max(PixResp(q, :, :, e),[], 3));
                    Trace = Trace - mean(Trace(TLDrift<0));
                    plot(TLDrift, Trace,'r'); hold on
                end
                Trace = squeeze(nanmean(mean(PixResp(:, :, :, e), 3), 1));
                %Trace = squeeze(nanmean(max(PixResp(:, :, :, e),[], 3), 1));
                Trace = Trace - mean(Trace(TLDrift<0));
                plot(TLDrift, Trace, 'k', 'LineWidth', 2); hold on
                                    ylim([-10 50])
                                    xlim([-.3 3.1])
                xticks([])
                yticks([])
%                 xlabel('Time (s)');
%                 ylabel('Responses');
                %title(sprintf('Display position: %d (um)', DispSize(e)));
                
            end

            
            stdRoi = std(reshape(PixResp, [], nROI), [], 1);
            stdRoi = reshape(stdRoi, 1, 1, 1, nROI);
            %stdRoi = repmat(stdRoi,[size(PixResp,1) size(PixResp,2) size(PixResp,3) 1]);
            %PixRespNrm = PixResp./stdRoi;
            PixRespNrm = PixResp;  %%!!!!! denormalize
            
            % quick check
            subplot(2, 3, 1);
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
            
        end
        
        TON  = TLDrift>0.05 & TLDrift<1.2;
        TOFF = TLDrift>1.55 & TLDrift<2.7;
        
        %%Polarity 1
        Polarity = mean(mPixRespNrm(TOFF, :, :) - mPixRespNrm(TON, :, :), 1);
        Polarity = squeeze(Polarity ./ mean(mPixRespNrm(TOFF, :, :) + mPixRespNrm(TON, :, :) + eps, 1));
        Polarity(Polarity>1) = 1;
        Polarity(Polarity<-1) = -1;
        subplot(2, 3, 2);
        imagesc(Polarity);colorbar
        
        %%Polarity 3
        maxON = squeeze(max(mPixRespNrm(TON, :, :),[],1));
        maxOFF =  squeeze(max(mPixRespNrm(TOFF, :, :),[],1));
        maxON = squeeze(max(maxON,[],1));
        maxOFF = squeeze(max(maxOFF,[],1));
        Polarity3 = maxOFF./(maxOFF + maxON);
        Polarity4 = (maxOFF-maxON)./(maxOFF+maxON);
        
        %%Plarity 5
        onVals = mPixRespNrm(TON, :, :);
        sortOn = sort(onVals,1,'descend');
        top = ceil(size(onVals,1) * .01);
        bottom = floor(size(onVals,1) * .3);
        onTop = squeeze(sortOn(top,:,:));
        onBottom = squeeze(sortOn(bottom,:,:));
        onRange = onTop - onBottom;
        onAll = squeeze(sqrt(sum(onRange.^2,1)));
        
        offVals = mPixRespNrm(TOFF, :, :);
        sortOff = sort(offVals,1,'descend');
        top = ceil(size(offVals,1) * .01);
        bottom = floor(size(offVals,1) * .3);
        offTop = squeeze(sortOff(top,:,:));
        offBottom = squeeze(sortOff(bottom,:,:));
        offRange = offTop - offBottom;
        offAll = squeeze(sqrt(sum(offRange.^2,1)));
        
        Polarity5 = offAll./(offAll+onAll);
        
        
        %%ONOFF
        ON = squeeze(mean(mPixRespNrm(TON,:,:),1));
        OFF = squeeze(mean(mPixRespNrm(TOFF,:,:),1));
        
        if 0
            roiNum = size(mPixRespNrm,3);
            posNum = size(mPixRespNrm,2);
            for g = 1:roiNum
                %clf
                subplot(roiNum,1,g)
                for h = 1:posNum
                    plot(mPixRespNrm(:,h,g));
                    hold on
                end
                ylim([-2 160])
                pause
            end
        end
        
        
        %%Polarity 2
        onMinusOff = mean(mPixRespNrm(TOFF, :, :)) - mean(mPixRespNrm(TON, :, :), 1);
        onPlusOff =  mean(mPixRespNrm(TOFF, :, :),1) + mean(mPixRespNrm(TON, :, :),1);
        Polarity2 = squeeze( onMinusOff ./ onPlusOff);
        Polarity2(Polarity2>1) = 1;
        Polarity2(Polarity2<-1) = -1;
        subplot(2, 3, 6);
        imagesc(Polarity2);colorbar
        
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
        ROI.PolarityX = [ROI.Polarity; mean(Polarity')'];
        ROI.SNR = [ROI.SNR; median(SNRM, 1)'];
        ROI.Polarity2 = [ROI.Polarity2; Polarity2'];
        ROI.ON = [ROI.ON; ON'];
        ROI.OFF = [ROI.OFF; OFF'];
        ROI.roiIDs = [ROI.roiIDs; grabROIs];
        ROI.maxON = [ROI.maxON; maxON'];
        ROI.maxOFF = [ROI.maxOFF; maxOFF'];
        ROI.Polarity3 = [ROI.Polarity3; Polarity3'];
        ROI.Polarity4 = [ROI.Polarity4; Polarity4'];
        ROI.Polarity5 = [ROI.Polarity5; Polarity5'];

        
        
    end
end
ROI.ROITable = [ROITable, ROI.Polarity, ROI.SNR];
ROI2 = ROI;
%save('../individualPointMasksSpread\ROI2.mat','ROI2')

save([SPN 'ROI2.mat'],'ROI2');
