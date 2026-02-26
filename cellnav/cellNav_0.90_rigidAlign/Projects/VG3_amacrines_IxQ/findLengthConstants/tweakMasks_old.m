

global tis glob
SPN = [glob.datDir 'Analysis\Data\preproc\'];
clf
SE = strel('square',3);
distThresh = 5; % was 5
valThresh = .2; % was 0.2, Minimum allowable threshold of new value of pixel relative to startting value






%Define source directory for masks
% 
KdatDir = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\';
maskDir = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasks\';
TPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\individualPointMasksSpread\';

ROIMasks = load([SPN 'maskDatPts.mat']);
ROIMasks = ROIMasks.maskDat;

ROITable = load([SPN 'ptDat.mat']);
ROITable = ROITable.ptDat;


% 
% if 0 % Copy raw mask dir to Mask dir that is used
%     rawMaskDir = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\Mask_3x3_11+1+2021'
%     MaskPath = '\\storage1.ris.wustl.edu\jlmorgan\Active\kerschensteinerLab\Emily\Mask';
%     delete(MaskPath);
%     copyfile(rawMaskDir,MaskPath)
%     
% end

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


%% Load file list
FilNam = sprintf('%sROIs/ListRecordingFiles_%s.m',KdatDir,ListId);
run(FilNam);
NumFile = length(R);

StimType = 'FlashBar';
clear ROI
ROI.Polarity = [];
ROI.SNR      = [];
for r = 1:NumFile
    Exps = R(r).Exps{1};
    if isempty(Exps), continue; end
    clc;
    fprintf('Flashbar... %d/%d \n', r, NumFile);
    

    for i = 1:length(Exps)
%         close all;
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
        load([KdatDir '/Translation/' FilNam], 'I');
        %         I = permute(I, [2 1 3]);
        
        if 0 %should median filter I
            I = medfilt3(I, [3, 3, 1]);
        end
        
        % Get correct timeline due to directional scanning
        TIds = reshape(1:numel(I), size(I));
        if IsBidirectional
            TIds = BidirecitonTimeCorrect3D(TIds);
        end
        
        % Load Time Table that aligned visual and scanning
        AlignName = sprintf('%s_Alignment_%s_%d%03d.mat',Topic, Day, Cel, Exp );
        load([KdatDir '/Alignment/' AlignName], 'TimeTable')
        
        % Loop for each qualified pixel
        % - x: recording Id (from the first pixel)
        % - v: pixel value
        % - xq: time Id and a fixed upsampling points
        % - vq (result): upsample epoching
        % (1) Spot (2) DGSpot
        
        
        %% Get the ROI signals
            
            grabMasks = find(ROITable(:, 1) == (Cel*1000 + Exp));
            
            ROIMask = ROIMasks(:, :, grabMasks);
            
            
            
            nROI = size(ROIMask, 3);
            
            %Expand ROIs to edges of signal
            Imean = mean(I,3);
            Imean = Imean - mode(Imean(:));
            
                 [ay ax] = find(Imean*0+1);
            showEach = 1;
            roiCent = [];
            for m = 1:nROI
                oldMask = ROIMask(:,:,m);
                [y x] = find(oldMask);
                
                my = mean(y);
                mx = mean(x);
                roiCent(m,:) = [my mx];

                allDists = sqrt((my-ay).^2 + (mx - ax).^2);
                distMask = Imean*0;
                distMask(allDists<=distThresh) = 1;
                
                %startVal = Imean(round(my),round(mx));
                startVal = mean(Imean(oldMask>0));
                bug = Imean*0;
                bug(round(my),round(mx)) = 1;
                Icol = zeros(size(Imean,1),size(Imean,2),3);
                Icol(:,:,1) = oldMask*1000;
                Icol(:,:,3) = Imean*60/mean(Imean(:));
                Icol(:,:,2) = bug * 1000;
                Iclose = Imean .* distMask;
                closeVals = Imean(distMask>0);
                
                %determine intensity threshold
                relativeThresh = (startVal * valThresh);
                medThresh = median(closeVals);
                useThresh = max(medThresh,relativeThresh);
                
                Ibright = Iclose .* (Iclose>useThresh);
                Ilab = bwlabel(Ibright,8);
                pickLab = Ilab(round(my),round(mx));
                if pickLab == 0 % If picked is out of range
                    if 0
                        [py px] = find(Ilab>0);
                        pDists = sqrt((py-my).^2 + (px - mx).^2);
                        closeP = find(pDists == min(pDists),1);
                        pickLab = Ilab(py(closeP),px(closeP));
                    else
                        %%Make mask that is disk around point of half
                        %%radius of mask
                        bug = Imean*0;
                        bug(allDists<=(distThresh/2)) = 1;
                    end
                else
                    bug = Ilab == pickLab;
                end
                
                
                
                
                
                %                 if showEach,image(uint8(Icol)),drawnow, end
                %                 for w = 1:1000
%                     bug2 = imdilate(bug,SE); %dilate bug;
%                     [by bx] = find(bug2-bug);% get new pixels
%                     inds = sub2ind(size(bug),by,bx);
%                     vals = Imean(inds); % get mean image values
%                     dists = allDists(inds); % get distance
%                     closeEnough = dists<= distThresh; % check distance threshold
%                     if sum(closeEnough)
%                         maxVal =  max(vals(closeEnough));
%                         if maxVal >= (startVal * valThresh)
%                             bestPick = find((vals ==maxVal) & closeEnough); %find high values
%                             newInds = inds(bestPick);
%                         else
%                             break
%                         end
%                     else
%                         break
%                     end
%                     if isempty(newInds)
%                         break
%                     else
%                         bug(newInds) = 1;
%                         Icol(:,:,2) = bug * 1000;
%                         if showEach,image(uint8(Icol)),drawnow, end
%                     end
%                     pause(.1)
%                 end
                Icol(:,:,2) = bug*1000;
                if 1,image(uint8(Icol)),drawnow, pause(.1),end
                ROIMask(:,:,m) = bug;

            end

            %% Make sure no pixels are shared between ROIs.
            roiCent

            sumM = sum(ROIMask,3);
            ROIMask2 = ROIMask;
            [fixY fixX] = find(sumM>1);
            for p = 1:length(fixY)
                dists = sqrt((roiCent(:,1)-fixY(p)).^2 + (roiCent(:,2)-fixX(p)).^2); 
                closest = find(dists == min(dists),1);
                ROIMask2(fixY(p),fixX(p),:) = 0;
                ROIMask2(fixY(p),fixX(p),closest) = 1;
            end
            %%Force center pixel
            colormap gray(256)
            for p = 1:size(roiCent,1)
                ROIMask2(roiCent(p,1),roiCent(p,2),p) = 1;
                image(ROIMask2(:,:,p)*100);
                hold on
                scatter(roiCent(p,2),roiCent(p,1),'r','filled')
                hold off
                drawnow
            end
            sumM2 = sum(ROIMask2,3);
            image(sumM2*50)
%             for i = 1:size(ROIMask,3)
%                 image(ROIMask2(:,:,i)*100)
%                 pause
%             end

            %%Record masks
            ROIMasks(:, :, grabMasks) = ROIMask2;
            
            
        end
end


s = sort(squeeze(sum(sum(ROIMasks,1),2)),'descend')';
disp(s)


maskDat = ROIMasks;
ptDat = ROITable;
%if ~exist(TPN,'dir'), mkdir(TPN),end
save([SPN 'maskDat.mat'],'maskDat');
%save([TPN 'ptDat.mat'],'ptDat');

disp('finished')

    
    
    
    
    