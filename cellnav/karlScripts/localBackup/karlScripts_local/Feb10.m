if 0
    figure();
    hold on
    cidList=[2,3,4,5,13,14];
    testVox=getCidVox(cidList,10,dsObj,tis);
    cmap=colormap(turbo);
    for i=1:length(cidList)
        curCol=cmap(round(i/length(cidList)*90),:);
        scatter3(testVox{i}(:,1),testVox{i}(:,2),testVox{i}(:,3),1,curCol)
        
        
    end
    corrPts=ptDat(:,6:8);
    corrPts=corrPts/25;
    corrPts(:,3)=corrPts(:,3)*10;
    scatter3(corrPts(:,2),corrPts(:,1),corrPts(:,3),50,'r','filled');
end

if 0
    analDir='Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\';
    datDir='Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\';
    load([analDir 'fvLibrary\tis.mat']);
end

if 0
    vgcList=[2 3 4 5 13 14];
    targCidList=tis.syn.edges(ismember(tis.syn.edges(:,2),vgcList),:);
    targTypeDat=cid2type(targCidList(:,1),tis);
    rgcTarg=targCidList(targTypeDat{1}==1);
    rgcTargTypeDat=targTypeDat{3}(targTypeDat{1}==1);
    rgcTargSubs=histc(rgcTargTypeDat,0:59);
    
    bpcCidList=tis.syn.edges(ismember(tis.syn.edges(:,1),vgcList),:);
    upTypeDat=cid2type(bpcCidList(:,2),tis);
    bpcUps=bpcCidList(upTypeDat{1}==7);
    bpcTypeDat=upTypeDat{3}(upTypeDat{1}==7);
    
end

if 1
    load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\dsObj.mat')
    load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\tis.mat')
    fullCidList=tis.cells.cids;
    length(fullCidList);
    fullTypeDat=cid2type(fullCidList,tis);
    sum(fullTypeDat{1}==7);
    allBPCcids=fullCidList(fullTypeDat{1}==7);
    allBPCVox=getCidVox(allBPCcids,10,dsObj,tis);
    allBPCVols=zeros(length(allBPCcids),1);
    for i=1:length(allBPCVox)
        allBPCVols(i)=length(allBPCVox{i});
    end
    %figure();
    %histogram(allBPCVols,100);
    %sum(allBPCVols>500);
end



skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
cidList=[2 3 4 5 13 14];

%make the phys images
mergeList=[7:10];
mergedImg=zeros(32,256);
for j=1:length(mergeList)
    curImIt=mergeList(j);
    curImg=physImgs{curImIt};
    mergedImg=mergedImg+curImg;
end

loadSkel=0;
if loadSkel
    allSkels=cell(1,6);
    for i=1:length(cidList)
        
        curCid=cidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end
end

for i=1:length(cidList)
    curCid=cidList(i);
    curSkel=allSkels{i};
    skNodePos=curSkel.sm.arbor.nodes.pos;
    skNodeSize=curSkel.sm.arbor.nodes.rad;
    %curVox=getCidVox(curCid,10,dsObj,tis);
    %curVox=curVox{1};
    curPtDat=ptDat(ptDat(:,3)==curCid,:);
    bottomPtDat=curPtDat(curPtDat(:,1)>1999,:);
    figBool=1;
    figBool2=1;
    
    ROIdat=ROI.ROITable;
    ROIpol=ROI.Polarity;
    ROIpolScld=(-ROIpol+1)*128;
    ROIpolScld(ROIpolScld>256)=256;
    ROIpolScld(ROIpolScld<1)=1;
    coolio=colormap(cool);
    ROIcolor=coolio(round(ROIpolScld),:);
    
    
    if figBool
        figure();
        
        hold on
        %scatter3(curVox(:,1),curVox(:,2),curVox(:,3),1,'.k');
        scatter3(skNodePos(:,1),skNodePos(:,2),skNodePos(:,3),skNodeSize*10,'.k');
        curPts=ptDat(ptDat(:,3)==curCid,6:8);
        curCols=ROIcolor(ptDat(:,3)==curCid,:);
        curPts=curPts./[25 25 2.5];
        scatter3(curPts(:,2),curPts(:,1),curPts(:,3),75,curCols,'o','filled');
        title([string(curCid)+' corrPt=green, OFF input=magenta, ON input=cyan']);
        %get the local BPC inputs
        bpcSynPos=tis.syn.pos(tis.syn.edges(:,1)==curCid,:)*10;
        bpcTypeDat=cid2type(tis.syn.edges(tis.syn.edges(:,1)==curCid,2),tis);
        offSubs=[1:5,15,18,19,24];
        onSubs=[6:12,14,17,20,21];
        offIDs=ismember(bpcTypeDat{3},offSubs);
        onIDs=ismember(bpcTypeDat{3},onSubs);
        scatter3(bpcSynPos(offIDs,1),bpcSynPos(offIDs,2),bpcSynPos(offIDs,3),50,'mv','filled');
        scatter3(bpcSynPos(onIDs,1),bpcSynPos(onIDs,2),bpcSynPos(onIDs,3),50,'cv','filled');
    end
    
    if figBool2
        figure();
        image(mergedImg/prctile(mergedImg(:),90)*256);
        colormap(gray);
        %colormap(turbo);
        hold on
        scatter(bottomPtDat(:,4),bottomPtDat(:,5),25,'ro','filled');
        title(string(curCid));
        
        
    end
    
    
    
    
    
    
end



