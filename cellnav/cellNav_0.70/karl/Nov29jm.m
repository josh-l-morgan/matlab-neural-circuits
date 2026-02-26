%Nov11
%% setup
%directory with Emily's corrected images
physDir='Y:\Active\kerschensteinerLab\CorrectRegistratedFiles\CorrectRegistratedFiles\';
%directory with the corr pts
ptsDir='Y:\Active\kerschensteinerLab\individualPointMasks\';
%dir with analysis
analDir='Y:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\';
%get a struct for the results
physCorrDat=struct;

%% parameters
slabThick=100;
loadAll=1;


%% get necessary data for IPL depth calcs
fvDir = [analDir 'Analysis\fvLibrary\'];
if exist([fvDir 'ref_gcl nucEdge.mat'],'file')
    ipl_bord_GCL = load([fvDir 'ref_gcl nucEdge.mat']);
    ipl_bord_INL = load([fvDir 'ref_inl nucEdge.mat']);
    GCLbord=ipl_bord_GCL.fv.vertices(:,:);
    INLbord=ipl_bord_INL.fv.vertices(:,:);
else
    GCLbord = [0 0 0; 100 0 0; 0 100 0; 100 100 0];
    INLbord = [0 0 100; 100 0 100; 0 100 100; 100 100 100];
end

Locs={GCLbord;INLbord}; %These are in z,x,y I think.
for i=1:2
    P=Locs{i};
    B(:,i) = [P(:,3), P(:,2), ones(size(P,1),1)] \ P(:,1);
end
GCLplane=struct();
INLplane=struct();
GCLplane.Parameters=[-1 B(2,1) B(1,1) B(3,1)];
INLplane.Parameters=[-1 B(2,2) B(1,2) B(3,2)];

%% loading

load([analDir 'Merge\dsObj.mat']);
load([analDir 'Analysis\tis.mat']);
load([ptsDir 'ptDat_newest.mat']);
regionNames=[1005:1010 2001:2006];
physImgs=cell(12,1);
corrPts=cell(12,1);
for i=1:length(regionNames)
    %get the data and make an average image of the scan region
    if loadAll
        curPhysImgFilename=strcat(['Ai148_129SVG3_Translation_122618_' num2str(regionNames(i)) '.mat']);
        curPhysImgDat=load([physDir curPhysImgFilename]);
    end
    curPhysMean=mean(curPhysImgDat.I,3);
    curPhysMeanImg=curPhysMean*(256/prctile(curPhysMean(:),95));
    physImgs{i}=curPhysMeanImg;
    curPhysImFilename=['Physimage_' num2str(regionNames(i)) '.png'];
    imwrite(curPhysMeanImg,['Y:\Active\morganLab\MATLAB\cellNav\cellNav_0.70\karl\images\emPhys\' curPhysImFilename]);
    %get the pt data into a parallel format
    curPts=pts.ptDat(pts.ptDat(:,1)==regionNames(i),:);
    corrPts{i}=curPts;
end

redo = 0;
if redo==1
    for i=1:length(regionNames)
        curPts=newDat(newDat(:,1)==regionNames(i),:);
        corrPts{i}=curPts;
    end
end

%% get the VG3 data from the plane fitted from the pts
%get the voxdata for VG3 cells
vgcCidList=[2 3 4 5 10 11 13 14 20 ];
cidVoxs=getCidVox(vgcCidList,1,dsObj,tis);

emptyMat=zeros(3200,3200,500,'uint8');
vgcMat=emptyMat;
%emptyMat=uint8(emptyMat);
for i=1:length(cidVoxs)
    curVoxList=cidVoxs{i};
    curVoxList=curVoxList(sum(curVoxList,2)>2,:);
    vgcMat(sub2ind(size(vgcMat),curVoxList(:,2),curVoxList(:,1),curVoxList(:,3)))=i;
end
%fit a plane to the phys corr pts
planeFits=cell(12,1);
emImages=cell(12,1);
emDat=cell(12,1);
emImageOffsets=cell(12,1);
tfs=cell(12,1);
for i=1:length(planeFits)
    curPts=corrPts{i};
    vastCoords=curPts(:,6:8);
    vastCoordsDS=vastCoords./[25 25 2.5];
    %B = [vastCoordsDS(:,2), vastCoordsDS(:,1), ones(size(vastCoordsDS,1),1)] \ ...
    %    vastCoordsDS(:,3);
    physPts=horzcat(curPts(:,4:5),repmat(1,[size(curPts,1),1]));
    %     emCorners=[-1,-1,-1;-1,-1,-1;-1,-1,-1;-1,-1,-1];
    %     while length(find(emCorners<0))>0
    %     try
    %         tf3d=estimateGeometricTransform3D(physPts,vastCoordsDS,'similarity');
    %         corners=[0,0,1;0,32,1;256,0,1;256,32,1];
    %         emCorners=transformPointsForward(tf3d,corners);
    %     catch
    fitIterations=49;
    tftr=zeros(3,3,fitIterations);
    goodTFsubsht=zeros(fitIterations,4);
    for q=1:fitIterations
        tf3d=estimateGeometricTransform(physPts(:,[1 2]),vastCoordsDS(:,[1 2]),'affine');
        tftr(:,:,q)=tf3d.T;
    end
    
    finalTF=affine2d();
    stdTF=std(tftr,0,3);
    meanTF=mean(tftr,3);
    highTF=meanTF+(stdTF);
    lowTF=meanTF-(stdTF);
    goodTFsubsht(:,[1:2])=squeeze(tftr(3,[2:3],:)>=lowTF(3,[2:3]))';
    goodTFsubsht(:,[3:4])=squeeze(tftr(3,[2:3],:)<=highTF(3,[2:3]))';
    goodTFs=tftr(:,:,sum(goodTFsubsht,2)==4);
    [sortedTX,TxIdx]=sort(goodTFs(3,1,:));
    
    finalTF.T=mean(goodTFs,3);
    
    debug=0;
    if debug
        figure();
        plot(0,0);
        hold on
        for w=1:fitIterations
            scatter(tftr(3,1,w),tftr(3,2,w));
        end
        scatter(mean(goodTFs(3,1,:),3),mean(goodTFs(3,2,:),3),40,'r+');
        title(string(regionNames(i)));
        %legend();
    end
    
    tfs{i}=finalTF;
    corners=[-8,-8;264,-8;-8,40;264,40];
    corners2=[0,0;256,0;0,32;256,32];
    emCorners=transformPointsForward(finalTF,corners);
    emCorners2=transformPointsForward(finalTF,corners2);
    %physEMsp=imwarp(physImgs{i},tf3d);
    %emCorners=horzcat(emCorners, [min(vastCoordsDS(:,3));repmat(mean(vastCoordsDS(:,3)), ...
    %    [size(emCorners,1)-2,1]);max(vastCoordsDS(:,3))]);
    emCorners=horzcat(emCorners, [(max(vastCoordsDS(:,3))-min(vastCoordsDS(:,3)))/2+min(vastCoordsDS(:,3)) - slabThick/2 ...
        ;repmat(mean(vastCoordsDS(:,3)),[size(emCorners,1)-2,1]); ...
        (max(vastCoordsDS(:,3))-min(vastCoordsDS(:,3)))/2+min(vastCoordsDS(:,3)) + slabThick/2]);
    
    %
    %     %corners=[0,0,1;0,32,1;256,0,1;256,32,1];
    %     %emCorners=transformPointsForward(tf3d,corners);
    %     end
    %curBoundaryMat=emptyMat;
    %curBoundaryMat(int32(min(emCorners(:,1))):int32(max(emCorners(:,1))), ...
    %    int32(min(emCorners(:,2))):int32(max(emCorners(:,2))), ...
    %    int32(min(emCorners(:,3))):int32(max(emCorners(:,3))))=1;
    %curPhysEMMat=emptyMat;
    %curPhysEMMat(curBoundaryMat>0)=vgcMat(curBoundaryMat>0);
    curEMROI=vgcMat(int32(min(emCorners(:,1))):int32(max(emCorners(:,1))), ...
        int32(min(emCorners(:,2))):int32(max(emCorners(:,2))), ...
        int32(min(emCorners(:,3))):int32(max(emCorners(:,3))));
    emImageOffsets{i}=[int32(min(emCorners(:,1))),int32(min(emCorners(:,2))),int32(min(emCorners(:,3))), ...
        int32(max(emCorners(:,1))),int32(max(emCorners(:,2))),int32(max(emCorners(:,3)))];
    test3dplot=0;
    if test3dplot
        figure();
        [x y z] = ind2sub(size(curEMROI), find(curEMROI>0));
        plot3(x, y, z, 'k.');
    end
    imageOut=max(curEMROI,[],3);
    emImages{i}=imageOut;
    curEMImFilename=['EMimage_' num2str(regionNames(i)) '.png'];
    imwrite(imageOut, ['Y:\Active\morganLab\MATLAB\cellNav\cellNav_0.70\karl\images\emPhys\' curEMImFilename]);
    
    
    
end



% PlaneFitCode
% Locs={GCLbord;INLbord}; %These are in z,x,y I think.
% for i=1:2
%     P=Locs{i};
%     B(:,i) = [P(:,3), P(:,2), ones(size(P,1),1)] \ P(:,1);
% end
% GCLplane=struct();
% INLplane=struct();
% GCLplane.Parameters=[-1 B(2,1) B(1,1) B(3,1)];
% INLplane.Parameters=[-1 B(2,2) B(1,2) B(3,2)];
%% set up the figure quick
firstFigure=0;
%% asdf
if firstFigure
    for i=1:length(physImgs)
        curPts=corrPts{i};
        curf=figure();
        %subplot(4,length(physImgs)/2,i+(floor((i-1)/6)*6));
        subplot(3,1,1);
        hold on
        imshow(physImgs{i},[0 255]);
        scatter(curPts(:,4),curPts(:,5),25,'ro','filled');
        %subplot(4,length(physImgs)/2,i+(floor((i-1)/6)*6)+1);
        subplot(3,1,2);
        newImage=imrotate(emImages{i},90);
        newImage2=flipdim(newImage,1);
        imshow(newImage2,[0,1]); %length(vgcCidList)]);
        hold on
        %colormap('turbo');
        scatter(curPts(:,6)/25-double(emImageOffsets{i}(1)), ...
            curPts(:,7)/25-double(emImageOffsets{i}(2)),50,'ro','filled');
        
        subplot(3,1,3);
        hold on
        invtf=invert(tfs{i});
        [emPhys, p]=imwarp(physImgs{i},tfs{i});
        %outputImage=max(emPhys(:,:,emImageOffsets{i}(3):emImageOffsets{i}(6)),[],3);
        %newImage=imrotate(physEMimg,90);
        %newImage2=flipdim(newImage,1);
        %imshow(outputImage,[0 1]);
        physEMoffsets=[p.XWorldLimits(1)-emImageOffsets{i}(1), ...
            p.YWorldLimits(1)-emImageOffsets{i}(2)];
        emptyImg=zeros([size(newImage2) 3]);
        emptyImg(:,:,1)=newImage2;
        %Need to accomodate situations where the phys goes off the edge
        if physEMoffsets(1)<0
            %chop the side off the phys to accomodate EM
            emptyImg(1:physEMoffsets(1)+p.ImageSize(1), ...
                physEMoffsets(2)+1:physEMoffsets(2)+p.ImageSize(2),2)=emPhys(-1*physEMoffsets(1)+1:end,:)/255;
        else
            emptyImg(physEMoffsets(1)+1:physEMoffsets(1)+p.ImageSize(1), ...
                physEMoffsets(2)+1:physEMoffsets(2)+p.ImageSize(2),2)=emPhys/255;
        end
        
        %% stuck on this paragraph
        if physEMoffsets(2) < 0
            %chop the side off the phys to accomodate EM
            emptyImg(1:physEMoffsets(1)+1, ...
                physEMoffsets(2)+1+p.ImageSize(2):physEMoffsets(2),2)=emPhys(-1*physEMoffsets(1)+1:end,:)/255;
        else
            emptyImg(physEMoffsets(1)+1:physEMoffsets(1)+p.ImageSize(1), ...
                physEMoffsets(2)+1:physEMoffsets(2)+p.ImageSize(2),2)=emPhys/255;
        end
        
        imshow(emptyImg);
        scatter(curPts(:,6)/25-double(emImageOffsets{i}(1)), ...
            curPts(:,7)/25-double(emImageOffsets{i}(2)),50,'co','filled');
    end
end
%% new plot
newCorrPts={};
secondFigure=0;
if secondFigure
for i=1:length(physImgs)
    figure();
    curPts=corrPts{i};
    %subplot(4,length(physImgs)/2,i+(floor((i-1)/6)*6));
    %subplot(3,1,1);
    %hold on
    %imshow(physImgs{i},[0 255]);
    %scatter(curPts(:,4),curPts(:,5),25,'ro','filled');
    %subplot(4,length(physImgs)/2,i+(floor((i-1)/6)*6)+1);
    %subplot(3,1,2);
    %hold on
    newImage=imrotate(emImages{i},90);
    newImage2=flipdim(newImage,1);
    %imshow(newImage2,[0,1]); %length(vgcCidList)]);
    %colormap('turbo');
    %scatter(curPts(:,6)/25-double(emImageOffsets{i}(1)), ...
    %    curPts(:,7)/25-double(emImageOffsets{i}(2)),50,'ro','filled');
    
    %subplot(4,3,i);
    
    invtf=invert(tfs{i});
    [emPhys, p]=imwarp(physImgs{i},tfs{i});
    %outputImage=max(emPhys(:,:,emImageOffsets{i}(3):emImageOffsets{i}(6)),[],3);
    %newImage=imrotate(physEMimg,90);
    %newImage2=flipdim(newImage,1);
    %imshow(outputImage,[0 1]);
    physEMoffsets=[p.XWorldLimits(1)-emImageOffsets{i}(1), ...
        p.YWorldLimits(1)-emImageOffsets{i}(2)];
    emptyImg=zeros([size(newImage2) 3]);
    emptyImg(:,:,1)=newImage2;
    %emImageOffsets{i}
    emptyImg(physEMoffsets(1)+1:physEMoffsets(1)+p.ImageSize(1), ...
        physEMoffsets(2)+1:physEMoffsets(2)+p.ImageSize(2),2)=emPhys/255;
    subplot(3,2,[1:4]);
    imshow(emptyImg);
    hold on
    scatter(curPts(:,6)/25-double(emImageOffsets{i}(1)), ...
        curPts(:,7)/25-double(emImageOffsets{i}(2)),50,'co','filled');
    title('Left Click Green; Right Click Red');
    
    [clickX,clickY,clickButt]=ginput();
    physCoords=horzcat(clickX(clickButt==1),clickY(clickButt==1));
    emCoords=horzcat(clickX(clickButt==3),clickY(clickButt==3));
    if ~isempty(emCoords)
    scatter(emCoords(:,1),emCoords(:,2),50,'cx');
    scatter(physCoords(:,1),physCoords(:,2),50,'mx');
    
    subplot(3,2,5);
    newImage=imrotate(emImages{i},90);
    newImage2=flipdim(newImage,1);
    imshow(newImage2,[0,1]); %length(vgcCidList)]);
    hold on
    %colormap('turbo');
    scatter(emCoords(:,1),emCoords(:,2),75,'ro','filled');
    
    subplot(3,2,6);
    imshow(physImgs{i},[0 255]);
    curtf=tfs{i};
    
    physCoordsNoOffset=(physCoords)-repmat([double(physEMoffsets(2)) double(physEMoffsets(1))],size(physCoords,1),1);
    physPtsTfd=curtf.transformPointsInverse(physCoordsNoOffset);
    curPhysOffsets=emImageOffsets{i};
    physCoordsAdj=physCoords+double(curPhysOffsets(1:2));
    physPts=curtf.transformPointsInverse(physCoordsAdj);
    hold on
    scatter(physPts(:,1),physPts(:,2),75,'ro','filled');
    
    newCorrPts{i}= horzcat(physPts,emCoords);
    end
end
newDat=pts.ptDat;
for j=1:length(newCorrPts)
    curPtList=newCorrPts{j};
    if ~isempty(curPtList)
        curEmImageOffsets=emImageOffsets{j};
        corrEMCoords=curPtList(:,[3:4]);
        corrEMCoords2=corrEMCoords+double(curEmImageOffsets(1:2));
        corrEMCoords3=horzcat(corrEMCoords2,repmat(mean(curEmImageOffsets([3,6])),size(curPtList,1),1));
        corrEMCoords4=corrEMCoords3.*[25 25 2.5];
        padDat=horzcat(repmat([regionNames(j) 0 0],size(curPtList,1),1), ...
            curPtList(:,[1:2]),corrEMCoords4);
        newDat=[newDat;padDat];
    end
end
end
%% get new points from the images.

for i=1:length(physImgs)
    figure();
    curPts=corrPts{i};
    %subplot(4,length(physImgs)/2,i+(floor((i-1)/6)*6));
    %subplot(3,1,1);
    %hold on
    %imshow(physImgs{i},[0 255]);
    %scatter(curPts(:,4),curPts(:,5),25,'ro','filled');
    %subplot(4,length(physImgs)/2,i+(floor((i-1)/6)*6)+1);
    %subplot(3,1,2);
    %hold on
    newImage=imrotate(emImages{i},90);
    newImage2=flipdim(newImage,1);
    %imshow(newImage2,[0,1]); %length(vgcCidList)]);
    %colormap('turbo');
    %scatter(curPts(:,6)/25-double(emImageOffsets{i}(1)), ...
    %    curPts(:,7)/25-double(emImageOffsets{i}(2)),50,'ro','filled');
    
    %subplot(3,1,1);
    
    invtf=invert(tfs{i});
    [emPhys, p]=imwarp(physImgs{i},tfs{i});
    %outputImage=max(emPhys(:,:,emImageOffsets{i}(3):emImageOffsets{i}(6)),[],3);
    %newImage=imrotate(physEMimg,90);
    %newImage2=flipdim(newImage,1);
    %imshow(outputImage,[0 1]);
    physEMoffsets=[p.XWorldLimits(1)-emImageOffsets{i}(1), ...
        p.YWorldLimits(1)-emImageOffsets{i}(2)];
    emptyImg=zeros([size(newImage2) 3]);
    emptyImg(:,:,1)=newImage2;
    emptyImg(physEMoffsets(1)+1:physEMoffsets(1)+p.ImageSize(1), ...
        physEMoffsets(2)+1:physEMoffsets(1)+p.ImageSize(2),2)=emPhys/255;
    imshow(emptyImg);
    hold on
    scatter(curPts(:,6)/25-double(emImageOffsets{i}(1)), ...
        curPts(:,7)/25-double(emImageOffsets{i}(2)),50,'co','filled');
    
    
    
    %% extra code
    
    allvgImage=max(vgcMat(:,:,180:260),[],3);
    figure(); imshow(allvgImage,[0 1])
    
    outputSynIDs=find(ismember(tis.syn.edges(:,2),[4]));
    targetSynDat=cid2type(tis.syn.edges(outputSynIDs,1),tis);
    rgcTargetIDs=find(targetSynDat{1}==1);
    rgcTargetSubtypes=targetSynDat{3}(rgcTargetIDs);
    figure();
    h1=histogram(rgcTargetSubtypes,'BinEdges',[0:55]);
    subTypeList=tis.cells.type.subTypeNames{1};
    subTypeList(2:56)=subTypeList(1:55);
    subTypeList{1}='na';
    xticks([0:55]);
    xticklabels(subTypeList);
    
end








