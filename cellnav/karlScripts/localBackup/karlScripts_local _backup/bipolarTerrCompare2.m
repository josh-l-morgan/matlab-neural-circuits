
bpcSizeCutoff=5000;

%bipolar overlap calculation
while 0
testFigs=0;
loadAll=1;
IPLhistStepSz=0.01;
%size of dilation in pixels (at mip4)
dilatorMat=zeros(11,11);
dilatorMatRGB=insertShape(dilatorMat,'FilledCircle',[6,6,5],'Color','white');
dilatorMat=dilatorMatRGB(:,:,1)>0.3;
%get the vastSubs structure
if loadAll
    dsObjPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\dsObj.mat';
    obiPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\obI.mat';
    tisPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat';
    fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
    %dsObjPath='G:\Data\MATLAB\0827_analysis\Volumes\0827\Merge\dsObj.mat';
    load(dsObjPath);
    load(obiPath);
    load(tisPath);
    %allHistDat=getHisto(tis,fvDir);
end
%load(dsObjPath);
vRes = obI.em.dsRes;

%Get the cids of all the bpcs
allTypes=cid2type(tis.cells.cids,tis);
%bpcCids=tis.cells.cids(tis.cells.type.cellTypeMat(:,7)==1);
bpcCids=tis.cells.cids(allTypes{1}==7);
%fix that 1008 has no voxels
bpcCids(2)=[];


bpcVxls=getCidVox(bpcCids,1,dsObj,tis);
for i=1:length(bpcVxls)
    curVxls=bpcVxls{i};
    curVxls=curVxls(curVxls(:,1)>1,:);
    bpcVxls{i}=curVxls;
end
bpcVxlDepths=cell(size(bpcVxls));
bpcHistDat=zeros(size(bpcVxls,2),1/IPLhistStepSz);
bpcMeanDepths=zeros(size(bpcVxls,2),1);
bpcSizes=zeros(size(bpcVxls,2),1);
for i=1:length(bpcVxlDepths)
    curVxls=bpcVxls{i};
    bpcSizes(i)=size(curVxls,1);
    dsVxls=double(curVxls)/10;
    [t,y,curDepths]=getIPLdepth(dsVxls(:,3),dsVxls(:,1),dsVxls(:,2),[],[]);
    bpcVxlDepths{i}=curDepths;
    curHistDat=histcounts(curDepths,[0:IPLhistStepSz:1]);
    bpcHistDat(i,:)=curHistDat;
    trmHistDat=curHistDat;
    trmHistDat(trmHistDat<1000)=0;
    trmDepthDat=curDepths(curDepths>0.2&curDepths<0.8);
    meanDepth=mean(trmDepthDat);
    if isnan(meanDepth)
        meanDepth=0;
    end
    bpcMeanDepths(i)=meanDepth;
    
end

if testFigs
    midList=[];
    f1=figure();
    hold on
    %f2=figure();
    %hold on
    for bpcIt=1:length(bpcVxls)
        curVxls=bpcVxls{bpcIt};
        curVxlDepths=bpcVxlDepths{bpcIt}*1000;
        if size(curVxls,1)>bpcSizeCutoff
            figure(f1);
            curVxls=horzcat(curVxls(:,1:2),curVxlDepths);
            curVxls=curVxls(curVxls(:,3)>50,:);
            h=scatter3(curVxls(:,1),curVxls(:,2),curVxls(:,3),1,'.');
            alpha = 0.25;
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
            %zlim([0 0.7]);
            drawnow
            curHist=bpcHistDat(bpcIt,:);
            %figure(f2);
            %plot(curHist);
            %xline(bpcMeanDepths(bpcIt)*100);
            if bpcMeanDepths(bpcIt)>0.4&&bpcMeanDepths(bpcIt)<0.5
                midList=[midList;bpcCids(bpcIt)];
            end
        end
    end
    xlim([800 2000])
    ylim([700 1900])
    zlim([0 900])
    
end

if spinnit
    fileDir='D:\work\gif\';
    stepAng=3;
    fileName=[fileDir 'test'];
    view(i-1*stepAng,0);
        drawnow
        set(gca,'visible','off')
    testFrame=getframe(f1);
    testFrame=testFrame.cdata;
    d4Fig=zeros(size(testFrame,1),size(testFrame,2),size(testFrame,3),(360/stepAng));
    for i=1:360/stepAng
        fileName=[fileDir 'im' num2str(i,'%03.f') '.jpg'];
        figure(f1);
        view(i-1*stepAng,0);
        drawnow
        set(gca,'visible','off')
        fim = getframe(f1);
        %testIm = fim.cdata;
        testIm = frame2im(fim);
        %(:,:,1);
        %saveas(testIm,fileName,'jpg');
        imwrite(testIm,fileName);
        d4Fig(:,:,:,i)=testIm;
        %d4Fig(:,:,1,i)=testIm(:,:,1);
        %d4Fig(:,:,2,i)=testIm(:,:,2);
        %d4Fig(:,:,3,i)=testIm(:,:,3);
%         if i==1
%         imwrite(testIm,fileName,'gif','LoopCount',Inf,'DelayTime',0.04);
%         else
%         imwrite(testIm,fileName,'gif','WriteMode','append','DelayTime',0.04);
%         end
    end
    ft=figure();
    hold on
    %sp1=subplot(3,1,1);
    %sp2=subplot(3,2,1);
    %sp3=subplot(3,3,1);
    if 0
    for i=1:30 %360/stepAng
        curImDat=d4Fig(:,:,:,i);
        imshow(curImDat);
        drawnow
        pause(0.01)
%         rim=d4Fig(:,:,1,i);
%         gim=d4Fig(:,:,2,i);
%         bim=d4Fig(:,:,3,i);
%         %figure(sp1)
%         subplot(1,3,1);
%         imshow(rim);
%         hold on
%         drawnow
%         %figure(sp2);
%         %pause(.04)
%         subplot(2,3,1);
%         imshow(gim);
%         hold on
%         drawnow
%         %pause(.04)
%         %figure(sp3);
%         subplot(3,3,1);
%         imshow(bim);
%         hold on
%         drawnow
        %pause(.04)
    end   
    end
    %imwrite(d4Fig,fileName,'gif','LoopCount',Inf,'DelayTime',0.04);
end


while 0
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
end

%get the dilated images
bpcDilIms=zeros(2000,2000,length(bpcVxls));
bpcBoundPolys={};
emptyIm=zeros(2000,2000);

for i=1:length(bpcVxls)
    curVxls=bpcVxls{i};
    curProj=unique(curVxls(:,1:2),'rows');
    %curHull=convhull(double(curProj));
    curHull=boundary(double(curProj));
    curPoly=polyshape(curProj(curHull(1:end-1),1),curProj(curHull(1:end-1),2));
    bpcBoundPolys{i}=curPoly;
    curIm=emptyIm;
    %curIm(sub2ind(size(emptyIm),curProj(:,2),curProj(:,1)))=1;
    %dilIm=conv2(curIm,dilatorMat,'same');
    curMask=poly2mask(curPoly.Vertices(:,1),curPoly.Vertices(:,2), ...
        size(curIm,1),size(curIm,2));
    %finIm=dilIm&curMask;
    bpcDilIms(:,:,i)=curMask;
    dilIm=conv2(curMask,dilatorMat,'same');
    bpcDilIms(:,:,i)=dilIm;
end

if 0
    figure();
    hold on
    xlim([700 2000]);
    ylim([700 2000]);
    for i=1:length(bpcVxls)
        if bpcSizes(i)>bpcSizeCutoff
            curPoly=bpcBoundPolys{i};
            curVerts=curPoly.Vertices;
            p0=patch(curVerts(:,1),curVerts(:,2),rand([1 3]));
            p0.FaceAlpha=0.1;
            %plot(curPoly);
            hold on
        end
    end
end



bpcVoxZMode=[];
for i=1:length(bpcVxls)
    bpcVoxZMode = [bpcVoxZMode;mode(bpcVxls{i}(:,3))];
end

bpcVoxZMean=[];
for i=1:length(bpcVxls)
    bpcVoxZMean = [bpcVoxZMean;mean(bpcVxls{i}(:,3))];
end
%get the number of mip4 voxels in each of the cids
bpcVoxNums=[];
for i=1:length(bpcVxls)
    vxlNum=length(bpcVxls{i});
    bpcVoxNums=[bpcVoxNums;vxlNum];
end

%get xy centers
centroids = {};
for curBpcIt=1:length(bpcVxls)
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
comparisonOverlapVxl = zeros(length(comparisonList),5);
for curCompIt=1:length(comparisonList)
    %curCompIt
    %comparisonList(curCompIt,:)
    bpcIdxA=comparisonList(curCompIt,1);
    bpcIdxB=comparisonList(curCompIt,2);
    %[bpcCids(bpcIdxA) bpcCids(bpcIdxB)]
    overlap=sum(bpcDilIms(:,:,bpcIdxA).*bpcDilIms(:,:,bpcIdxB),'all');
    %bpcVxlsA=bpcVxlsDil{bpcIdxA};
    %bpcVxlsB=bpcVxlsDil{bpcIdxB};
    %overlap=intersect(bpcVxlsA,bpcVxlsB,'rows');
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
    comparisonOverlapVxl(curCompIt,1)=overlap;
    comparisonOverlapVxl(curCompIt,2)=bpcMeanDepths(bpcIdxA);
    comparisonOverlapVxl(curCompIt,3)=bpcMeanDepths(bpcIdxB);
    comparisonOverlapVxl(curCompIt,4:5)=[bpcCids(bpcIdxA) bpcCids(bpcIdxB)];
end

overlapMat=zeros(length(bpcCids),length(bpcCids));
for i=1:length(comparisonList)
    bpcIdxA=comparisonList(i,1);
    bpcIdxB=comparisonList(i,2);
    overlapMat(bpcIdxA,bpcIdxB)=comparisonOverlapVxl(i,1);
end

distMat=zeros(length(bpcCids),length(bpcCids));
for i=1:length(comparisonList)
    bpcIdxA=comparisonList(i,1);
    bpcIdxB=comparisonList(i,2);
    distMat(bpcIdxA,bpcIdxB)=centDistList(i);
end

%figure out which ones are close but do not overlap
compDat=horzcat(comparisonOverlapVxl(:,1),centDistList,comparisonOverlapVxl(:,2:5));
[comparisonSrtd,srtIdx]=sortrows(compDat,1,'descend');

% figure();
% tiledlayout('flow');
% image(distMat,'CDataMapping','scaled');
% hold on
% image(overlapMat,'CDataMapping','scaled');

%% OFF and ON layers separate
compDat=horzcat(comparisonOverlapVxl(:,1),centDistList,comparisonOverlapVxl(:,2:5));
offDat=compDat(compDat(:,3)<0.42&compDat(:,4)<0.42,:);
onDat=compDat(compDat(:,3)>0.42&compDat(:,4)>0.42,:);
testFigs2=0;
if testFigs2
    tf2=figure();
    OFFbpcs=offDat(:,5:6);
    OFFbpcs=OFFbpcs(:);
    OFFbpcs=unique(OFFbpcs);
    testTypeDat=cid2type(OFFbpcs,tis);
    hold on
    for i=1:length(OFFbpcs)
        curCid=OFFbpcs(i);
        vxls=getCidVox(curCid,1,dsObj,tis);
        vxls=vxls{1};
        scatter3(vxls(:,1),vxls(:,2),vxls(:,3),'.');
    end
end


%% groups the cells by overlap and distance
logLap=log(offDat(:,1));
logLap(logLap<0.01)=0.01;
rootDist=sqrt(offDat(:,2));
dist=offDat(:,2);
[v,rankLap]=sort(offDat(:,1),'descend');
[v,rankDist]=sort(offDat(:,2),'ascend');

offCompCids=comparisonListCids(compDat(:,3)<0.42&compDat(:,4)<0.42,:);
labels=cell(length(offDat(:,1)),1);
for i=1:length(offDat(:,1))
    curPair=offCompCids(i,:);
    labels{i}=num2str(curPair);
end

figure();
hold on
scatter(dist,logLap);
text(dist+0.25, logLap,labels,'FontSize',6) 
ylabel('log(voxelOverlap)');
xlabel('centroid distance (um)');

%offProLaps=offDat(:,1)./offDat(:,2);
%% try to creep through the bpcs and make a matrix
% gg=[1006 1033 2008 2033 3637 1084 1065 1074];
% gg=[gg 1179 3116 1061 1137 1083]; %added 033122
% bg=[1028 1128 1111 2110 3632 1195 3640 1138 6027 ];
%gg=[];
%bg=[];
%goodComps=[[1033 1074];[1033 1065];[1006 1033];[1006 3637];[1006 2008];[2008 3637];[1006 2033];[1065 2033];[1084 2033];[1084 2008]];
%badComps=[[1006 1128];[1128 1083];[2013 2014];[2013 2015]];


%cidComps=nchoosek(gg,2);
cidComps=nchoosek(OFFbpcs,2);
cidComps=sortrows(cidComps,1,'ascend');

goodCompIdxs=[];
for i=1:length(goodComps)
    %curComp=cidComps(i,:);
    curComp=goodComps(i,:);
    [a,idx]=ismember(curComp,offCompCids,'rows');
    if ~isempty(idx)
        goodCompIdxs=[goodCompIdxs;idx];
    end
end
goodCompIdxs=goodCompIdxs(goodCompIdxs~=0);


badCompIdxs=[];
for i=1:length(badComps)
    %curComp=cidComps(i,:);
    curComp=goodComps(i,:);
    [a,idx]=ismember(curComp,offCompCids,'rows');
    if ~isempty(idx)
        badCompIdxs=[badCompIdxs;idx];
    end
end
badCompIdxs=badCompIdxs(badCompIdxs~=0);

hold on

scatter(dist(goodCompIdxs),logLap(goodCompIdxs),50,'go','filled');
scatter(dist(badCompIdxs),logLap(badCompIdxs),50,'ro','filled');


%%
seed=1065;
seed=1121;
%seed=1074;
%seed=1208;
%seed=3116;
seedList=[1065, 1074, 1208, 3116];
seedList=[1033,1105,2033,3632,3640,1208,1179,1074,1137,1138,3310,3334,6012];
compIdxs=find(offDat(:,5)==seed|offDat(:,6)==seed);
curDat=offDat(compIdxs,:);
[x,srtIdx]=sort(curDat(:,2),'ascend');
curDatSrtd=curDat(srtIdx,:);
testGrp=curDatSrtd(curDatSrtd(:,2)<300,5:6);
testGrp=testGrp(:);
testGrp=testGrp(testGrp~=seed);
testGrp=testGrp(~ismember(testGrp,gg));
testGrp=testGrp(~ismember(testGrp,bg));

truths=[];
lies=[];
validComps=[];
badComps=[];

%% testing

plotTestGrp=1;
if plotTestGrp
    figure();
    hold on
    for w=1:length(seedList)
        curSeed=seedList(w);
        seedVxls=bpcVxls{bpcCids==curSeed};
        seedPoly=bpcBoundPolys{bpcCids==curSeed};
        
        compIdxs=find(offDat(:,5)==curSeed|offDat(:,6)==curSeed);
        curDat=offDat(compIdxs,:);
        [x,srtIdx]=sort(curDat(:,2),'ascend');
        curDatSrtd=curDat(srtIdx,:);
        testGrp=curDatSrtd(curDatSrtd(:,2)<300,5:6);
        testGrp=testGrp(:);
        testGrp=testGrp(testGrp~=curSeed);
        testGrp=testGrp(~ismember(testGrp,gg));
        testGrp=testGrp(~ismember(testGrp,bg));
        y=1;
        r='';
        while true
            if y>length(testGrp)
                break
            end
            curPart=testGrp(y);
            compDatIdx=find(sum(ismember(offDat(:,5:6),curSeed),2)>0& ...
                sum(ismember(offDat(:,5:6),curPart),2)>0);
            CompOverlap=offDat(compDatIdx,1);
            if CompOverlap>0
                [curSeed curPart]
                %offDat(compIdxs,1)
                
                clf
                scatter3(seedVxls(:,1),seedVxls(:,2),seedVxls(:,3),1,'m.');
                hold on
                seedTextLoc=centroids{bpcCids==curSeed};
                text(seedTextLoc(1),seedTextLoc(2),seedTextLoc(3)+50,num2str(curSeed));
                partVxls=bpcVxls{bpcCids==curPart};
                scatter3(partVxls(:,1),partVxls(:,2),partVxls(:,3),1,'c.');
                partTextLoc=centroids{bpcCids==curPart};
                text(partTextLoc(1),partTextLoc(2),partTextLoc(3)+50,num2str(curPart));
                partPoly=bpcBoundPolys{bpcCids==curPart};
                curVerts=partPoly.Vertices;
                p2=patch(curVerts(:,1),curVerts(:,2),rand([1 3]));
                p2.FaceColor=[1 1 .8];
                
                curVerts=seedPoly.Vertices;
                p1=patch(curVerts(:,1),curVerts(:,2),rand([1 3]));
                %p1.FaceColor=(rand([1 3])+.5)./1.5;
                p1.FaceColor=[.6 1 1];
                
                interPoly=intersect(seedPoly,partPoly);
                %interBound=boundary(interPoly);
                curVerts=interPoly.Vertices;
                goodVerts=sum(isnan(curVerts),2)==0;
                if ~isempty(goodVerts)
                    upatch=patch(curVerts(goodVerts,1),curVerts(goodVerts,2),rand([1 3]));
                    upatch.FaceColor=[1 0 1];
                end
                zlim([0 400]);
                centCent=mean([seedTextLoc;partTextLoc]);
                xlim([centCent(1)-200 centCent(1)+200]);
                ylim([centCent(2)-200 centCent(2)+200])
                view(0,90);
                drawnow
                q1="g=good; b=bad; v=valid; n=no touchie; c=go back; x=exit";
                r=input(q1,"s");
                if r=='x'
                    break
                elseif r=='g'
                    gg=[gg; curPart];
                    truths=[truths;[curSeed, curPart]];
                elseif r=='b'
                    bg=[bg; curPart];
                    lies=[lies;[curSeed, curPart]];
                elseif r=='n'
                    badComps=[badComps;[curSeed, curPart]];
                elseif r=='c'
                    y=y-2;
                end
            end
            y=y+1;
        end
    end
end

%% Maximize the goodness of the mosaics
%truths=[];
%lies=[];


putCidList=[truths(:);lies(:)];
putCidList=unique(putCidList);

end

%% another Test
badCompOverlap=zeros(length(badComps),1);
for i=1:length(badComps)
    curComp=badComps(i,:);
    compDatIdx=find(sum(ismember(offDat(:,5:6),curComp(1)),2)>0& ...
        sum(ismember(offDat(:,5:6),curComp(2)),2)>0);
    badCompOverlap(i)=offDat(compDatIdx,1);
end

goodCompOverlap=zeros(length(truths),1);
for i=1:length(truths)
    curComp=goodComps(i,:);
    compDatIdx=find(sum(ismember(offDat(:,5:6),curComp(1)),2)>0& ...
        sum(ismember(offDat(:,5:6),curComp(2)),2)>0);
    goodCompOverlap(i)=offDat(compDatIdx,1);
end

%%



outputTest=cell(10,10);
outputHists=zeros(10,10,40);
% figure();
% tiledlayout('flow');
% hold on
for n=10:5:30

    for w=2:2:20
    
testCombs=nchoosek(putCidList,n);
testCombs=sort(testCombs,2);
testCombs=unique(testCombs,'rows');

TS=sort(truths,2);
LS=sort(lies,2);

combRating=zeros(length(testCombs),1);
%figure();
%hold on
for combIt=1:length(testCombs)
    curComb=testCombs(combIt,:);
    curCombList=nchoosek(curComb,2);
    CCL=sort(curCombList,2);
    good=sum(ismember(TS,CCL,'rows'));
    bad=sum(ismember(LS,CCL,'rows'));
    %scatter(good,bad);
    combRating(combIt)=good-(bad*w);
    
end
% nexttile
% histogram(combRating);
curHist=histcounts(combRating,[-20:20]);
% drawnow
best=find(combRating==max(combRating));
outputHists(n,w,:)=curHist;
outputTest{n,w}=testCombs(best,:);
    end
end

figure();
hold on;
for i=1:10
    curHistDat=outputHists(i,10,:);
    curHistDat=curHistDat(:);
    plot(curHistDat);
end

