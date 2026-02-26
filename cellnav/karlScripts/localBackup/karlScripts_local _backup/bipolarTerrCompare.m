%bipolar overlap calculation

testFigs=1;
loadAll=0;

%size of dilation in pixels (at mip4)
dilatorMat=zeros(11,11);
dilatorMatRGB=insertShape(dilatorMat,'FilledCircle',[6,6,5],'Color','white');
dilatorMat=dilatorMatRGB(:,:,1)>0.3;
%get the vastSubs structure
if loadAll
dsObjPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\dsObj.mat';
obiPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\obI.mat';
tisPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\tis.mat';
fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\fvLibrary\';
%dsObjPath='G:\Data\MATLAB\0827_analysis\Volumes\0827\Merge\dsObj.mat';
load(dsObjPath);
load(obiPath);
load(tisPath);
end 
%load(dsObjPath);
vRes = obI.em.dsRes;

%Get the cids of all the bpcs
bpcCids=tis.cells.cids(tis.cells.type.cellTypeMat(:,7)==1);
bpcVxls=getCidVox(bpcCids,1,dsObj,tis);

if testFigs
    figure();
    hold on
    for bpcIt=1:length(bpcVxls)
        curVxls=bpcVxls{bpcIt};
        scatter3(curVxls(:,1),curVxls(:,2),curVxls(:,3),1,'.');
        drawnow
    end
end

allHistDat=getHisto(tis,fvDir);

%make an array to hold the object ids
bpcVxls=cell(1,length(bpcCids));
bpcVxlsDil=cell(1,length(bpcCids));
bpcBBs=cell(1,length(bpcCids));
targ = find(obI.cell.name == bpcCids(i));
isCell = obI.cell.obIDs{targ};
parfor curBpcIt=1:length(bpcCids)
    bpcCids(curBpcIt)
    curBpcVxls=[];
    for curObIt=1:length(obI.nameProps.cids)
        %if sum(obI.nameProps.cids{curObIt}==bpcCids(curBpcIt))>0
        if ismember(bpcCids(curBpcIt),obI.nameProps.cids{curObIt})
            curBpcVxls=[curBpcVxls;vastSubs{curObIt}];
        end
    end
    minX=min(curBpcVxls(:,1));
    minY=min(curBpcVxls(:,2));
    minZ=min(curBpcVxls(:,3));
    maxX=max(curBpcVxls(:,1));
    maxY=max(curBpcVxls(:,2));
    maxZ=max(curBpcVxls(:,3));
    curBpcBB=[minX minY minZ maxX maxY maxZ];
    anchor=[minX-100,minY-100,minZ-100];
    curBpcVxlsRel=curBpcVxls-anchor;
    curBpcVxlsRelDil=[0 0 0];
    for curSlice=100:maxZ-minZ+100
        %curSlice
        curSlicePxls=curBpcVxlsRel(curBpcVxlsRel(:,3)==curSlice,:);
        dilMat=zeros(maxX-minX+200,maxY-minY+200);
        inds = sub2ind([maxX-minX+200,maxY-minY+200],curSlicePxls(:,1),curSlicePxls(:,2));
        dilMat(inds) = 1;
        outputMat=imdilate(dilMat,dilatorMat);
        [row,col]=ind2sub(size(outputMat),find(outputMat>0));
        outputResults=horzcat(row,col,repmat(curSlice,[length(row) 1]));
        curBpcVxlsRelDil=[curBpcVxlsRelDil;outputResults];
    end
    curBpcVxlsDil=curBpcVxlsRelDil(2:end,:)+anchor;
    bpcBBs{curBpcIt}=curBpcBB;
    bpcVxls{curBpcIt}=curBpcVxls;
    bpcVxlsDil{curBpcIt}=curBpcVxlsDil;
end

%%For reference - fast distance matric
if 0
ds = 10;
vox = bpcVxls{40};
voxS = round(vox/ds);
maxVox= max(voxS,[],1);
voxInd = sub2ind(maxVox,voxS(:,1),voxS(:,2),voxS(:,3));
voxInd = unique(voxInd);
[Y X Z] = ind2sub(maxVox,voxInd);
vox = [Y X Z];

vox1 = vox;
vox2 = vox;

dif1 = vox1(:,1) - vox2(:,1)';
dif2 = vox1(:,2) - vox2(:,2)';
dif3 = vox1(:,3) - vox2(:,3)';
dist = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
end

bpcVoxZMode=[];
for i=1:length(bpcCids)
    bpcVoxZMode = [bpcVoxZMode;mode(bpcVxls{i}(:,3))];
end

bpcVoxZMean=[];
for i=1:length(bpcCids)
    bpcVoxZMean = [bpcVoxZMean;mean(bpcVxls{i}(:,3))];
end
%get the number of mip4 voxels in each of the cids
bpcVoxNums=[];
for i=1:length(bpcCids)
    vxlNum=length(bpcVxls{i});
    bpcVoxNums=[bpcVoxNums;vxlNum];
end
    
%get xy centers
centroids = {};
for curBpcIt=1:length(bpcCids)
    if ~isempty(bpcVxls{curBpcIt})
        centroids{curBpcIt}=mean(bpcVxls{curBpcIt},1);
        
    end
end

%get the combinations of bpc cids
comparisonList = nchoosek(1:length(bpcCids),2);
comparisonListCids = bpcCids(comparisonList);

%get the distance between centroids for each pair
centDistList=zeros(length(comparisonList),1);
for curCompIt=1:length(comparisonList)
    bpcA=comparisonList(curCompIt,1);
    bpcB=comparisonList(curCompIt,2);
    centA=centroids{bpcA};
    centB=centroids{bpcB};
    dist=sqrt((centA(1)-centB(1))^2+(centA(2)-centB(2))^2);
    centDistList(curCompIt)=dist;
end

%compute the overlap for each comparison pair
comparisonOverlapVxl = zeros(length(comparisonList),1);
parfor curCompIt=1:length(comparisonList)
    comparisonList(curCompIt,:)
    bpcIdxA=comparisonList(curCompIt,1);
    bpcIdxB=comparisonList(curCompIt,2);
    bpcVxlsA=bpcVxlsDil{bpcIdxA};
    bpcVxlsB=bpcVxlsDil{bpcIdxB};
    overlap=intersect(bpcVxlsA,bpcVxlsB,'rows');
    %totalOverlap=0;
    %make an empty slice list
    %sliceArray=zeros(1050,1);
    %get the slices where CHANGE TO BOTH KARL of the two cells has pixels
    %unique(bpcVxlsA(:,3)
    %sliceArray = intersect(unique(bpcVxlsA(:,3)),unique(bpcVxlsB(:,3)));
        
    if 0
        for curSliceIt=1:length(sliceArray)
        curSlice=sliceArray(curSliceIt);
        curSlice
        %overlapArray=zeros(4000,4000,1050);
        sliceVxlsA=bpcVxlsA(bpcVxlsA(:,3)==curSlice,:);
        sliceVxlsB=bpcVxlsB(bpcVxlsB(:,3)==curSlice,:);
        %This can be sped up a lot if I stay in voxel lists
        sizeX = 4000;
        sizeY = 4000;
        imA=zeros(sizeY,sizeX);
        imB=zeros(sizeY,sizeX);
        %JOSH, the syntax here could be faster
        
        indA = sub2ind([sizeY sizeX],sliceVxlsA(:,1),sliceVxlsA(:,2));
        imA(indA) = 1;
        indB = sub2ind([sizeY sizeX],sliceVxlsB(:,1),sliceVxlsB(:,2));
        imB(indB) = 1;
        %         for pxl=1:length(sliceVxlsA(:,1))
        %             imA(sliceVxlsA(pxl,1),sliceVxlsA(pxl,2))=1;
        %         end
        %         for pxl=1:length(sliceVxlsB(:,1))
        %             imB(sliceVxlsB(pxl,1),sliceVxlsB(pxl,2))=1;
        %         end
        %create the dilation matrix. KARL - Make this a circle.
        tic
        SE = strel('disk',15);
        dilator=ones(dilX,dilY);
        tic
        imDilA=imdilate(imA,dilator);
        toc
        
        imDilB=imdilate(imB,dilator);
        tic
        %imDilA = conv2(imA,dilator);
        toc
        imOverlap=imDilA&imDilB;
        %overlapArray(:,:,curSlice)=imOverlap;
        slicOverlap=sum(imOverlap(:));
        totalOverlap=totalOverlap+slicOverlap;
        end
    end
    comparisonOverlapVxl(curCompIt)=length(overlap);
end