%% run planeFitPhys.m first

% This is going to heavily use the deplanes structure created in previous 

ROIimgMult=25;

zRad=25;

yMidTop=1150;
yMidBot=1350;
yRad=80;
xMid=1400;
xRad=500;

fvDir = 'Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\fvLibrary\';
mergeDir = 'Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\';
load([mergeDir 'obI.mat']);
load([mergeDir 'dsObj.mat']);
vgcOBIDs=[];
for i=1:10
    curIDs=obI.cell.obIDs{i};
    vgcOBIDs=[vgcOBIDs curIDs];    
end

allVGCvoxels=[];
for curObIt=1:length(vgcOBIDs)
    curObID=vgcOBIDs(curObIt);
    curObVox=dsObj(curObID).subs;
    allVGCvoxels=[allVGCvoxels; curObVox];
end


figure()
hold on
scatter3(allVGCvoxels(1:10:end,2),allVGCvoxels(1:10:end,1),allVGCvoxels(1:10:end,3),2,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.8 0.8 0.8]);

for i=1:length(roiList)
    bbox=deplanes(i).bbox;
    bbds=bbox./[25 25 2.5];
    curCent=[bbds(2,1)-bbds(1,1)/2+bbds(1,1) ...
        bbds(3,2)-bbds(1,2)/2+bbds(1,2) ...
        bbds(1,3)];
    
    if i<7
        yMid=yMidTop;
    else
        yMid=yMidBot;
    end

    vgcsInBB=abs(double(allVGCvoxels(:,2))-xMid)<xRad & ...
        abs(double(allVGCvoxels(:,1))-yMid)<yRad & ...
        abs(double(allVGCvoxels(:,3))-curCent(3))<zRad ;
    %vgcsInBB= ...
    %    abs(curCent(1)-double(allVGCvoxels(:,1)))<(bbds(2,1)-bbds(1,1)) & ...
    %    abs(curCent(2)-double(allVGCvoxels(:,2)))<(bbds(3,2)-bbds(1,2)) & ...
    %    abs(curCent(3)-double(allVGCvoxels(:,3)))<10 ;
    
    vgcInBBVox=allVGCvoxels(vgcsInBB,:);
    scatter3(vgcInBBVox(:,2),vgcInBBVox(:,1),vgcInBBVox(:,3),10,colmap(mod(i,6)+1,:));
    curROIimg=zeros(yRad*2,xRad*2);
    vgcVoxCrp=double(vgcInBBVox)-[yMid-yRad xMid-xRad curCent(3)-10];
    for j=1:length(vgcVoxCrp)
        curPx=vgcVoxCrp(j,:);
        curROIimg(curPx(1),curPx(2))=curROIimg(curPx(1),curPx(2))+1;
    end
    deplanes(i).curROIEMZ=uint8(curROIimg.*ROIimgMult);
    outputImg=uint8(curROIimg.*ROIimgMult);
    outputFilename=[num2str(roiList(i)) '_emz.png'];
    %imwrite(outputImg,outputFilename);
end



%% compare summed images

for i=1:length(roiList)
    figure();
    hold on
    imgA=deplanes(i).imDat;
    imgB=deplanes(i).curROIEMZ;
    subplot(2,1,1);
    imshow(imgA*5);
    subplot(2,1,2);
    imshow(imgB);
end