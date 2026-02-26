%% load up some data
analDir='Y:\karlsRetina\CellNavLibrary_IxQ\';
volName='AprilMerge';
mergeDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\';
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
obI=load([mergeDir 'obI.mat']);
dsObj=load([mergeDir 'dsObj.mat']);
global tis
%% params
testFig=0;



%% CB palette
cbPalette=[51, 34, 136;  17, 119, 51; 68, 170, 153; 136, 204, 238; ...
    221, 204, 119; 204, 102, 119; 170, 68, 153; 136, 34, 85];
cbPalette=cbPalette/255;
%% VG3 syn location analyses
vgcCidList=[2 3 4 5 13];
figure();
hold on
for vIt=1:length(vgcCidList)
    curVgcCid=vgcCidList(vIt);
    
    %get the vgcVoxels
    vgcIt=find(tis.obI.cell.name==curVgcCid);
    curVgcObList=tis.obI.cell.obIDs{vgcIt};
    curVgcVox=[];
    for curObIt=1:length(obI.obI.nameProps.cids)
        curObCids=obI.obI.nameProps.cids{curObIt};
        if curObCids==curVgcCid
            %curObID=curVgcObList(curObIt);
            curVgcVox=[curVgcVox; dsObj.dsObj(curObIt).subs];
        end
    end
    %test figure
    if testFig==1
        h=scatter3(curVgcVox(1:100:end,2),curVgcVox(1:100:end,1),curVgcVox(1:100:end,3),5,'filled');
        h.MarkerEdgeColor = cbPalette(vIt,:);
        h.MarkerFaceColor = cbPalette(vIt,:);
    end
    inputSynIDs=find(tis.syn.edges(:,1)==curVgcCid);
    outputsynIDs=find(tis.syn.edges(:,2)==curVgcCid);
end

%% load skels (this will take a while)
skelList=[2 3 4 5];
skelDat=struct;
for i=1:length(skelList)
    curCid=skelList(i);
    skelDat(i).SM=struct;
    skelDat(i).SM=getSkelDat(skelDir,curCid);
end
%% of skels and syns

%This is all going to be working from the skelDat, so the data fetching and
%structures wil be kinda different than dealing with tis.dat

figure();
hold on

for vIt=1:length(skelList)
    curVgcCid=skelList(vIt);
    %inputSynIDs=find(tis.syn.edges(:,1)==curVgcCid);
    %outputsynIDs=find(tis.syn.edges(:,2)==curVgcCid);
    curSkelPosDat = skelDat(vIt).SM.sm.nep.pos(:,:);
    %curSkelPosDat = skelDat(vIt).SM.sm.arbor.nodes.pos(:,:);
    %curSkelTrunkPos = curSkelPosDat(curSkelPosDat(:,3)>max(curSkelPosDat(:,3))-2,:);
    %curSkelRootPos = mean(curSkelTrunkPos);
    %rootNodeDistList = abs(sum(curSkelPosDat-curSkelRootPos,2));
    rootNodeID = find(curSkelPosDat(:,3)==max(curSkelPosDat(:,3)));
    rootNodeID = rootNodeID(1);
    %rootNodeID = find(curSkelPosDat==rootNodePos);
    
    synPos = skelDat(vIt).SM.sm.syn.pos(:,:);
    syn2Skel = skelDat(vIt).SM.sm.syn2Skel.skelTopoDist(:,:);
    rootSynDistList = syn2Skel(:,rootNodeID);
    
    %get the colormap set up
    maxSynRootDist = max(rootSynDistList);
    rootSynDistFracList = round(rootSynDistList/maxSynRootDist*255);
    synDistCMap = jet(256);
    cmapMat=synDistCMap(rootSynDistFracList,:);
    
    
    if testFig==1
        h=scatter3(curSkelPosDat(:,2),curSkelPosDat(:,1),curSkelPosDat(:,3),5,'filled');
        h.MarkerEdgeColor = cbPalette(vIt,:);
        h.MarkerFaceColor = cbPalette(vIt,:);
        %j=scatter3(curSkelRootPos(2),curSkelRootPos(1),curSkelRootPos(3),150,'filled');
        %j=scatter3(curSkelPosDat((2),curSkelRootPos(1),curSkelRootPos(3),150,'filled');
        %j.MarkerEdgeColor = cbPalette(vIt,:);
        %j.MarkerFaceColor = cbPalette(vIt,:);
        k=scatter3(synPos(:,2),synPos(:,1),synPos(:,3),50,cmapMat,'filled');
        %k.MarkerEdgeColor = synDistCMap(round(rootSynDistFracList*100),:);
        %k.MarkerFaceColor = synDistCMap(round(rootSynDistFracList*100),:);
    end
    
    
    
end
