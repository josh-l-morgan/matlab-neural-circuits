
bpcVoxCutoff=15000;
HCcolmap=load('HCcolmap.mat');
HCcolmap=HCcolmap.HCcolmap;
HCcolmap=HCcolmap./255;
jetCmap=colormap(jet);
jetTest=1;
if jetTest
    jt1=figure();
    x=5:5:255;
    scatter(x,repmat([1],[length(x),1]),100,jetCmap(x,:),'filled');
    xticks(x);
end
%bipolar overlap calculation
while true
    testFigs=0;
    loadAll=1;
    IPLhistStepSz=0.01;
    %size of dilation in pixels (at mip4)
    dilatorMat=zeros(11,11);
    dilatorMatRGB=insertShape(dilatorMat,'FilledCircle',[6,6,5],'Color','white');
    dilatorMat=dilatorMatRGB(:,:,1)>0.3;
    %get the vastSubs structure
    if loadAll
        srcDir='C:\work\ml\fakeLib\';
        vastSubsPath=[srcDir 'Volumes\Final\Merge\vastSubs.mat'];
        dsObjPath=[srcDir 'Volumes\Final\Merge\dsObj.mat'];
        obiPath=[srcDir 'Volumes\Final\Merge\obI.mat'];
        tisPath=[srcDir 'Volumes\Final\Analysis\fvLibrary\tis.mat'];
        fvDir=[srcDir 'Volumes\Final\Analysis\fvLibrary\'];
        %dsObjPath='G:\Data\MATLAB\0827_analysis\Volumes\0827\Merge\dsObj.mat';
        load(vastSubsPath);
        load(dsObjPath);
        load(obiPath);
        load(tisPath);
        %allHistDat=getHisto(tis,fvDir);
        clear dsObjPath obiPath tisPath vastSubsPath
        curTis=tis;
    end
    %load(dsObjPath);
    vRes = obI.em.dsRes;
    
    %Get the cids of all the bpcs
    allTypes=cid2type(tis.cells.cids,curTis);
    %bpcCids=tis.cells.cids(tis.cells.type.cellTypeMat(:,7)==1);
    bpcCids=tis.cells.cids(allTypes{1}==7);
    %bpcCids(bpcCids==2108)=[];
    bpcCids=bpcCids(bpcCids<4000|bpcCids>5000);
    bpcTypeDat=cid2type(bpcCids,curTis);
    bpcSubtypes=bpcTypeDat{3};
    bpcSubtypes(bpcSubtypes==0)=find(contains(tis.cells.type.subTypeNames{7},{' '}));
    subTypeIDList=unique(bpcSubtypes);
    subTypeNameList=tis.cells.type.subTypeNames{7}(subTypeIDList);
    subTypeNames=tis.cells.type.subTypeNames{7}(bpcSubtypes)';
    [ind, bpcColID]=ismember(bpcSubtypes,subTypeIDList);
    
    bpcVxls=getCidVox(bpcCids,1,dsObj,curTis);
    for i=1:length(bpcVxls)
        curVxls=bpcVxls{i};
        curVxls=double(curVxls);
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
    
    bpcDilIms=zeros(2000,2000,length(bpcVxls));
    bpcBoundPolys={};
    emptyIm=zeros(2000,2000);
    rangmult=3;
    for i=1:length(bpcVxls)
        curVxls=bpcVxls{i};
        curProj=unique(curVxls(:,1:2),'rows');
        %HERE IS WHERE TO REMOVE OUTLIER VOXELS
        xmean=mean(curProj(:,1));
        ymean=mean(curProj(:,2));
        xstdev=std(curProj(:,1));
        ystdev=std(curProj(:,2));
        xrange=[xmean-rangmult*xstdev xmean+rangmult*xstdev];
        yrange=[ymean-rangmult*ystdev ymean+rangmult*ystdev];
        goodIndsX=find(curProj(:,1)>xrange(1)&curProj(:,1)<xrange(2));
        goodIndsY=find(curProj(:,2)>yrange(1)&curProj(:,2)<yrange(2));
        goodInds=unique([goodIndsX;goodIndsY]);
        curProj=curProj(goodInds,:);
        %figure();
        %scatter(curProj(:,1),curProj(:,2));
        %title(num2str(bpcCids(i)));
        %curHull=convhull(double(curProj));
        curHull=boundary(double(curProj));
        curPoly=polyshape(curProj(curHull(1:end-1),1),curProj(curHull(1:end-1),2));
        if sum(isnan(curPoly.Vertices(:)))==0
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
    end
    
    %% how to draw the perimeters
    while 0
        curPoly=bpcBoundPolys{i};
        curVerts=curPoly.Vertices;
        %curFV=
        %p0=patch(curVerts(:,1),curVerts(:,2),rand([1 3]));
        p0=patch(curVerts(:,1),curVerts(:,2),BCcolmap(bpcColID(i),:));
        p0.FaceAlpha=0.1;
        %plot(curPoly);
        hold on
    end
    
    %% Get the other needed matrices
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
    
    bpcBigEnough=bpcVoxNums>bpcVoxCutoff;
    
    
    
    
    %get xy centers
    centroids = {};
    for curBpcIt=1:length(bpcVxls)
        if ~isempty(bpcVxls{curBpcIt})
            centroids{curBpcIt}=mean(bpcVxls{curBpcIt},1);
            
        end
    end
    
    %get the combinations of bpc cids
    filterBool=1;
    if filterBool
        if 0
            turmap=colormap(turbo);
            figure();
            hold on
            for g=1:length(subTypeIDList)
                curSubType=subTypeIDList(g);
                scatter(bpcMeanDepths(bpcTypeDat{3}==curSubType),rand(1),10,turmap(g*10,:));
            end
        end
        onBpcCids=bpcMeanDepths>0.4;
        offBpcCids=bpcMeanDepths<0.4;
        onBpcCidListFinal=find(bpcBigEnough&onBpcCids);
        offBpcCidListFinal=find(bpcBigEnough&offBpcCids);
        combList={onBpcCidListFinal,offBpcCidListFinal};
        
    end
%     onBpcCidListFinal=type2cid({'bpc','bpc','bpc','bpc'},{'bc5i','bc5o','bc5t','bcon'},curTis);
%     offBpcCidListFinal=type2cid({'bpc','bpc','bpc','bpc'},{'bc3a','bc3b','bc4','bcoff'},curTis);
%     onBpcCidListFinal=[onBpcCidListFinal{1};onBpcCidListFinal{2};onBpcCidListFinal{3};onBpcCidListFinal{4}];
%     offBpcCidListFinal=[offBpcCidListFinal{1};offBpcCidListFinal{2};offBpcCidListFinal{3};offBpcCidListFinal{4}];
%     bigEnoughCids=bpcCids(bpcSizes>bpcVoxCutoff);
%     onBpcCidListFinal=intersect(onBpcCidListFinal,bigEnoughCids);
%     offBpcCidListFinal=intersect(offBpcCidListFinal,bigEnoughCids);
%     for r=1:length(onBpcCidListFinal)
%         findCid=onBpcCidListFinal(r);
%         onBpcCidListFinal(r)=find(bpcCids==findCid);
%     end
%     for r=1:length(offBpcCidListFinal)
%         findCid=offBpcCidListFinal(r);
%         offBpcCidListFinal(r)=find(bpcCids==findCid);
%     end
    %% doing this for on and off
    results=struct();
    for i=1:2
        curList=combList{i};
        comparisonList = nchoosek(1:length(curList),2);
        comparisonListCids = bpcCids(curList(comparisonList));
        results(i).comparisonCids=comparisonListCids;
        %comparisonList=comparisonListCids; %I don't like it, but it keeps me from having to change the variable for all the other loops
        
        
        %get the distance between centroids for each pair
        centDistList=zeros(length(comparisonListCids),1);
        for curCompIt=1:length(comparisonListCids)
            bpcA=comparisonListCids(curCompIt,1);
            bpcB=comparisonListCids(curCompIt,2);
            centA=centroids{find(bpcCids==bpcA)};
            centB=centroids{find(bpcCids==bpcB)};
            if ~isempty(centA) & ~isempty(centB)
                dist=sqrt((centA(1)-centB(1))^2+(centA(2)-centB(2))^2);
                centDistList(curCompIt)=dist;
            end
        end
        results(i).centDistList=centDistList;
        
        %compute the overlap for each comparison pair
        comparisonOverlapVxl = zeros(length(comparisonListCids),5);
        for curCompIt=1:length(comparisonListCids)
            %curCompIt
            %comparisonList(curCompIt,:)
            bpcIdxA=comparisonListCids(curCompIt,1);
            bpcIdxB=comparisonListCids(curCompIt,2);
            bpcCidIdxA=find(bpcCids==bpcIdxA);
            bpcCidIdxB=find(bpcCids==bpcIdxB);
            
            %[bpcCids(bpcIdxA) bpcCids(bpcIdxB)]
            overlap=sum(bpcDilIms(:,:,bpcCidIdxA).*bpcDilIms(:,:,bpcCidIdxB),'all');
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
            comparisonOverlapVxl(curCompIt,2)=bpcMeanDepths(bpcCidIdxA);
            comparisonOverlapVxl(curCompIt,3)=bpcMeanDepths(bpcCidIdxB);
            comparisonOverlapVxl(curCompIt,4:5)=[bpcCids(bpcCidIdxA) bpcCids(bpcCidIdxB)];
            if mod(curCompIt,100)==0
                curCompIt
            end
        end
        results(i).comparisonOverlapVxl=comparisonOverlapVxl;
        
        overlapMat=zeros(length(curList),length(curList));
        for q=1:length(comparisonList)
            bpcIdxA=comparisonList(q,1);
            bpcIdxB=comparisonList(q,2);
            overlapMat(bpcIdxA,bpcIdxB)=comparisonOverlapVxl(q,1);
        end
        results(i).overlapMat=overlapMat;
        
        distMat=zeros(length(curList),length(curList));
        for s=1:length(comparisonList)
            bpcIdxA=comparisonList(s,1);
            bpcIdxB=comparisonList(s,2);
            distMat(bpcIdxA,bpcIdxB)=centDistList(s);
        end
        results(i).distMat=distMat;
        
        %figure out which ones are close but do not overlap
        compDat=horzcat(comparisonOverlapVxl(:,1),centDistList,comparisonOverlapVxl(:,2:5));
        [comparisonSrtd,srtIdx]=sortrows(compDat,1,'descend');
        results(i).comparisonSrtd=comparisonSrtd;
        results(i).srtIdx=srtIdx;
    end
    
    % Sort the things from the previous section.
    
    
    %% Parse and visualize the previously generated data
    f11=figure();
    hold on
    tl10=tiledlayout(f11,2,1);
    for j=1:2
        logLap=log(results(j).comparisonOverlapVxl(:,1));
        logLap(logLap<0.01)=0.01;
        results(j).logLap=logLap;
        rootDist=sqrt(results(j).centDistList);
        dist=results(j).centDistList;
        [v,rankLap]=sort(results(j).comparisonOverlapVxl(:,1),'descend');
        [v,rankDist]=sort(results(j).centDistList,'ascend');
        compCids=results(j).comparisonCids;
        %offCompCids=comparisonListCids(compDat(:,3)<0.42&compDat(:,4)<0.42,:);
        labels=cell(length(results(j).centDistList),1);
        for k=1:length(results(j).centDistList)
            curPair=compCids(k,:);
            labels{k}=num2str(curPair);
        end
        noLap=logLap;
        noLap(noLap==0.01)=inf;
        overDist=sqrt(((dist/10).*(dist/10))+(noLap.*noLap));
        results(j).overDist=overDist;
        nexttile(j);
        scatter(dist/10,logLap);
        results(j).distLapRatio=(dist./logLap);
        text(dist/10+0.025, logLap,labels,'FontSize',6)
        ylabel('log(voxelOverlap)');
        xlabel('centroid distance (um)');
        xlim([0 50])
        ylim([1 20])
    end
    
    
    
    %% Visual inspection of results
    manualComp=cell(2,1);
    for m=2:-1:1
        %only looking at the 3a
        curResults=results(m);
        manualComp{m}=repmat(9,[length(curResults.distLapRatio),1]);
        %get the order of the overlap ratios
        [srtd,indx]=sort(curResults.centDistList);
        comparFig=figure();
        
        for n=1:length(curResults.centDistList)%curResults.overDist(curResults.overDist<50))
            curCompCids=curResults.comparisonCids(indx(n),:);
            curCompTypes=cid2type(curCompCids,curTis);
            if ismember(3,curCompTypes{3})
            compareMorph(comparFig,curCompCids,fvDir);
            
            title('0=diff, 1=same, 2=not sure');
            inp=input(num2str(curCompCids));
            manualComp{m}(indx(n))=inp;
            end
        end
        
    end
    
    %%
    comparFig=figure();
    for i=1:comparisons
        compareMorph(comparFig,comparisons(i,:),fvDir);
        pause()
    end
    
    %% test code
    test=results(2).centDistList(find(results(2).comparisonCids(:,1)==1111 ...
        |results(2).comparisonCids(:,2)==1111));
    
    [dist,inds]=sort(test);
    cids=results(2).comparisonCids(find(results(2).comparisonCids(:,1)==1111 ...
        |results(2).comparisonCids(:,2)==1111),:);
    cidsCln=cids(cids~=1111);
    cidsClnSrtd=cidsCln(inds);
    cidsDistsSrtd=test(inds);
    res=horzcat(cidsClnSrtd,cidsDistsSrtd);
    
    
    
    
    test2=results(2).centDistList(find(results(2).comparisonCids(:,1)==1084 ...
        |results(2).comparisonCids(:,2)==1084));
    
    compFig=figure();
    tileBool=zeros(12,1);
    for y=1:12
        curPartCid=cidsClnSrtd(y);
        compareMorph(compFig,[1111 curPartCid],[]);
        p=input('');
        tileBool(y)=p;
    end
    
    
    
    %% ************************************
    % OLD CODE FROM HERE ON DOWN
    % *************************************
    
    
    
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




%% FV code snippets
curCid=bpcCids(bpcIt);
curFV=load([fvDir num2str(curCid) '.mat']);
curFV=curFV.fv;
curP=patch(curFV);
curP.CData=BCcolmap(bpcColID(bpcIt),:);
curP.LineStyle='none';
curP.EdgeColor='none';
curP.FaceAlpha=0.2;
curP.FaceColor=BCcolmap(bpcColID(bpcIt),:);
%curP.EdgeAlpha=0.2;
%curP.EdgeColor=BCcolmap(bpcColID(bpcIt),:);
%drawnow
curHist=bpcHistDat(bpcIt,:);




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

