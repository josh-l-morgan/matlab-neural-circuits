
%% Define folder
%SPN = 'Y:\Active\morganLab\DATA\LGN_Regeneration\Export\SM26-0_1_LR_TDE\';
%SPN = 'Y:\Active\morganLab\DATA\LGN_Regeneration\Export\SM26-0_1_Wide\';
SPN = 'Y:\Active\morganLab\DATA\LGN_Regeneration\Export\SM26-0_2_Wide\'


%% Parse image names
dSPN = dir([SPN '*.tif']);
nams = {dSPN.name};
for i = 1:length(nams)
    nam = nams{i};
    und = regexp(nam,'_');
    ch = regexp(nam,'ch');
    ch = ch(end);
    und = und(end-1:end);
    chan(i) = str2num(nam(ch+2:ch+3));
    iCount(i) = str2num(nam(und(1)-3:und(1)-1)); %count in file name
end
uICount = unique(iCount); % list of counts in file names
numImages = length(uICount)
uChan = unique(chan);


%% Assign images
global globSort
globSort.ic = uICount(1);
globSort.uICount = uICount;
globSort.numImages = numImages;
globSort.iCount = iCount;
globSort.chan = chan;
globSort.uChan = uChan;
globSort.nams = nams;
globSort.SPN = SPN;
globSort.checkChannel = 0;
globSort.groupAssignment = zeros(length(uICount),1);
save([SPN 'globSort.mat'],'globSort')

f = SortImages
return
% checkChannel = 0;
% colormap gray(256)
% for i = 1:numImages
%     ic = uICount(i)
%     imIdx = find((iCount==ic) & (chan == checkChannel));
%     I = double(imread([SPN nams{imIdx}]));
%     clf
%     imshow(uint8(I * 100/mean(I(:))))
%     pause
% end



%% set image groups
load([SPN 'globSort.mat'])
imageTags = {'L' 'R' 'W'};
clear imageOrder
for i = 1:3
imageOrder{i} =  find(globSort.groupAssignment==i); % order of images acquired 
end
numChan = length(uChan);
ys = 2048;
xs = 2048;



%% Convert images
for i = 1:length(imageTags)
    disp(sprintf('group %d',i))
    useCs = uICount(imageOrder{i});
    if ~isempty(useCs)

        groupDir = [SPN imageTags{i} '\'];
        if ~exist(groupDir,'dir'),mkdir(groupDir);end

        for z = 1:length(useCs)
            disp(sprintf('z %d of %d',z,length(useCs)))
            useCs(z);
            imIdxs = find(iCount==useCs(z));
            imChan = chan(imIdxs);
            I = zeros(ys,xs,numChan);
            for t = 1:length(imIdxs);
                I(:,:,imChan(t)+1) = imread([SPN nams{imIdxs(t)}]);
            end
            Is = I(:,:,[1 2 5]);
            Is2 = I(:,:,[3 4 6]);
            subplot(1,2,1)
            image(uint16(Is))
            subplot(1,2,2);
            image(uint16(Is2))
            drawnow
            fileName = sprintf('%ssec%03.0f.mat',groupDir,z);
            save(fileName,'I')
        end
    end
end

%% %%%%%%%  %run alinment script

alignMatSections

















