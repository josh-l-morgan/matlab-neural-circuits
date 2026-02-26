%
%Plan
% take a list of N cids. Return an N x N matrix showing the size of overlap
% in DS voxels.
% Input: cidList, dsObj, obi or tis
% steps:
%  for each cell, get a list of object IDs
%  get the voxels for each cid
%  do a convex hull on a zprojection
%  compare via same code as overlap
% test params
if 0
    load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\dsObj.mat');
    load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\obI.mat');
    cids=[1006 1075 1117 1121 6151 1160 1228 2008 2009 2033 1033 ...
        1061 1065 1074 1083 1084 1105 1113 1137 2015 3333 5314 6012 ...
        1095 1111 1138 1159 1179 1195 2030 3116 3310 3311 3329 6035 6113];
    DS=dsObj;
    curtis=tis;
    emptyImg=zeros(2048,2048);
    emptyImg=uint16(emptyImg);
    totImg=emptyImg;
end

grp3a=[1006 1075 1117 1121 1160 1228 2008 2009 2033 3334 6137 6149];
grp3b=[1033 1061 1065 1074 1083 1084 1105 1113 1137 2015 3333 5314 6012];
grp4=[1095 1111 1138 1159 1179 1195 2013 2030 3116 3310 3311 3329 6035 6113];

function outputMat=getMosaicOverlap(cids,DS,obI,debugBool)

%% TEST CHUNK 1


kl=[0 0 1 0 0; ...
    0 1 1 1 0; ...
    1 1 1 1 1; ...
    0 1 1 1 0; ...
    0 0 1 0 0];
km=kl;
kl=uint8(kl);
for i=1:8
    kn=conv2(km,kl);
    km=kn;
end
kl=uint16((km).^0.25);
%Get all of the voxels for the cells
cidList=cids;
dsObj=DS;
obI=obI;
allVox={};
allHulls={};
allImgs={};
if debugBool
    figure();
    hold on
end
for cidIt=1:length(cidList)
    curCidID=find(obI.cell.name==cidList(cidIt));
    obIdList=obI.cell.obIDs{curCidID};
    cidVoxList=[];
    cidImgMask=emptyImg;
    for obIt=1:length(obIdList)
        cidVoxList=[cidVoxList;dsObj(obIdList(obIt)).subs];
    end
    allVox{cidIt}=cidVoxList;
    cidZproj=unique(cidVoxList(:,1:2),'rows');
    for voxIt=1:length(cidZproj)
        cidImgMask(cidZproj(voxIt,1),cidZproj(voxIt,2))=1;
    end
    cidImgDial=conv2(cidImgMask,kl,'same');
    allImgs{cidIt}=cidImgDial;
    if debugBool
        totImg=totImg+uint16(cidImgDial);
    end
    if 0
        figure();
        scatter3(cidVoxList(:,2),cidVoxList(:,1),cidVoxList(:,3));
    end
    
    
    %let's make some convex hulls to visualize the territory
    STDxy=std(double(cidZproj));
    center=mean(cidZproj,1);
    conv=convhull(double(cidZproj));
    convPts=cidZproj(conv,:);
    %convPtsSparse=[];
    dists=sqrt(sum((double(convPts([1:end-1],:)) ...
        - double(convPts([2:end],:))).^2,2));
    convPtsSparse=convPts(dists>1,:);
    allHulls{cidIt}=convPtsSparse;
    
    if debugBool
        %figure();
        %hold on
        %scatter(cidZproj(1:10:end,1),cidZproj(1:10:end,2),1,'ko');
        shap=polyshape(convPtsSparse(:,1),convPtsSparse(:,2));
        text(center(1),center(2),num2str(cidList(cidIt)),'Color','red','FontSize',14,'HorizontalAlignment','center');
        %'Color','red','FontSize',14
%         t = linspace(0,2*pi);
%         x = 4 + 4*cos(t);
%         y = 5 + 5*sin(t);
%         plot(x,y)
        %shap=polyshape(convPts(:,1),convPts(:,2));
        plot(shap,'FaceColor',faceCol,'FaceAlpha',0.1);
    end
    
    %make the comparisons among the different cells

end

%% NEXT
outputMat=zeros(length(cids),length(cids));
comboList=nchoosek(1:length(cidList),2);
for compIt=1:length(comboList)
    compAID=comboList(compIt,1);
    compBID=comboList(compIt,2);
    compAcid=cidList(compAID);
    compBcid=cidList(compBID);
    compAVox=allImgs{compAID};
    compBVox=allImgs{compBID};
    corrMat=compAVox.*compBVox;
    corrVal=sum(corrMat>100,'all');
    outputMat(compAID,compBID)=corrVal;
    outputMat(compBID,compAID)=corrVal;
    
    if debugBool2
        figure();
        imMat=zeros(length(emptyImg),length(emptyImg),3);
        imMat=uint16(imMat);
        imMat(:,:,1)=compAVox/max(compAVox(:))*65000;
        imMat(:,:,2)=compBVox/max(compBVox(:))*65000;
        imMean=mean(imMat(imMat>0));
        %imStd=std(imMat(:));
        %imMat(imMat>imMean*2)=imMean*2;
        imshow(imMat);
    end
end

figure();
imshow(outputMat,[0 max(outputMat(:))]);
colormap('jet');

%% 





end