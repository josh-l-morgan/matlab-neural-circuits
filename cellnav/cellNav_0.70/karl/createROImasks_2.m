%% load points
dataDir='W:\individualPointMasks\';
destDir='W:\individualPointMasks\';
ptDat=load([dataDir 'ptDat.mat']);
ptDat=ptDat.pts.ptDat;
roiList=unique(ptDat(:,1));

vizOutput=0;

%% define the kernel for masking each point
ptMask=ones(3,3);

%% create masks
maskDat=zeros(32,256,length(ptDat(:,1)));
ptIterator=1;
for i=1:length(ptDat(:,1))
    curMask=zeros(32,256);
    curMask(ptDat(i,5),ptDat(i,4))=1;
    outputMask=conv2(curMask,ptMask,'same');
    maskDat(:,:,i)=outputMask;
end


if vizOutput==1
    figure();
    hold on
    for i=1:length(ptDat(:,1))
        curImg=maskDat(:,:,i);
        imshow(curImg);
        pause(0.2);
    end
end

% for i=1:length(roiList)
%     curMask=blankMask;
%     fieldPtDat=ptDat(ptDat(:,1)==roiList(i),:);
%     curROIlocs=fieldPtDat(:,[4 5]);
%     for j=1:length(curROIlocs(:,1))
%         curMask(curROIlocs(j,2),curROIlocs(j,1))=1;
%         outputMask=conv2(curMask,ptMask,'same');
%         outputFilename
%     end
% end

%% save masks