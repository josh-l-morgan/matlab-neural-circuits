
%% load up some data
analDir='Y:\karlsRetina\CellNavLibrary_IxQ\';
volName='AprilMerge';
mergeDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\';
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
%skelDir='C:\Users\karl\Documents\MATLAB\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
obI=load([mergeDir 'obI.mat']);
dsObj=load([mergeDir 'dsObj.mat']);


%% params
testFig=0;
makeFig=1;
skelList=[2 3 4 5 13];
vgcCidList=[2 3 4 5 13];
cbPalette=[51, 34, 136;  17, 119, 51; 68, 170, 153; 136, 204, 238; ...
    221, 204, 119; 204, 102, 119; 170, 68, 153; 136, 34, 85];
cbPalette=cbPalette/255;

%%
loadAll=0;
if loadAll==1
    global tis
    %load skels (this will take a while
    skelDat=struct;
    for i=1:length(skelList)
        curCid=skelList(i);
        skelDat(i).SM=struct;
        skelDat(i).SM=getSkelDat(skelDir,curCid);
    end
end
%% VG3 syn location analyses
if testFig==1
    figure();
    hold on
end
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

%% of skels and syns
%This is all going to be working from the skelDat, so the data fetching and
%structures wil be kinda different than dealing with tis.dat
if makeFig==1
    figure();
    hold on
end
for vIt=1:length(skelList)
    curVgcCid=skelList(vIt);
    curSkelPosDat = skelDat(vIt).SM.sm.nep.pos(:,:);
    rootNodeID = find(curSkelPosDat(:,3)==max(curSkelPosDat(:,3)));
    rootNodeID = rootNodeID(1);
    synPos = skelDat(vIt).SM.sm.syn.pos(:,:);
    syn2Skel = skelDat(vIt).SM.sm.syn2Skel.skelTopoDist(:,:);
    rootSynDistList = syn2Skel(:,rootNodeID);
    maxSynRootDist = max(rootSynDistList);
    rootSynDistFracList = round(rootSynDistList/maxSynRootDist*255);
    synDistCMap = jet(256);
    cmapMat=synDistCMap(rootSynDistFracList,:);
    if makeFig==1
        h=scatter3(curSkelPosDat(:,2),curSkelPosDat(:,1),curSkelPosDat(:,3),5,'filled');
        h.MarkerEdgeColor = cbPalette(vIt,:);
        h.MarkerFaceColor = cbPalette(vIt,:);
        k=scatter3(synPos(:,2),synPos(:,1),synPos(:,3),50,cmapMat,'filled');
    end
end
