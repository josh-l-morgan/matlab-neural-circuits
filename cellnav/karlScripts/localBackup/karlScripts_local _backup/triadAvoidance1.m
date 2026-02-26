
%Initial Parameters
vgcCidList=[2 3 4 5 13 14];
loadAll=0; loadSkel=0;

%% loading zone
if loadAll
    %standard loading of stuff
    dsObjPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\dsObj.mat';
    obiPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\obI.mat';
    tisPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat';
    fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
    %fvLib='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
    load('Y:\MATLAB\cellNav\karlScripts\localBackup\karlScripts_local\aliasList.mat');
    HCcolmap=load('HCcolmap.mat');
    HCcolmap=HCcolmap.HCcolmap;
    HCcolmap=HCcolmap./255;
    colTest=0;
    if colTest
        colFig=figure();
        scatter(1:15,repmat(1,[1,15]),50,HCcolmap,'filled');
    end
    load(dsObjPath);
    load(obiPath);
    load(tisPath);
    clear dsObjPath obiPath tisPath
    skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\';
    if ~exist('allSkels')
        if ~exist('D:\work\remove\allSkels.mat')
            load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\allSkels.mat');
        else
            load('D:\work\allSkels.mat');
        end
        if 0
            allSkels=cell(1,length(vgcCidList));
            for i=1:length(vgcCidList)
                curCid=vgcCidList(i);
                skelFileName=['sm_cid' + string(curCid)+'.mat'];
                curSkel=load([skelDir+skelFileName]);
                allSkels{i}=curSkel;
            end
        end
    end
    if ~exist('curTis')
        %global glob tis
        curTis=tis;
    end
end

%% parameters
skelDist=30;
synRad=30;
plotz=0;
plotz2=1;
p=0;

%% NOT WORKING loop
if 0
    for i=1:length(allSkels)
        %% load skels, etc.
        curSM=allSkels{i}.sm;
        curCid=curSM.cid;
        curSWC=allSkels{i}.swc;
        
        %get the types of all the synapse partners
        preTypes=cid2type(curSM.syn.edges(:,2),curTis);
        postTypes=cid2type(curSM.syn.edges(:,1),curTis);
        %get the fv library object
        %curCidFv=load([fvDir num2str(curCid) '.mat']);
        curCidFv=curSM.nep.fv;
        % get depth info for nodes/syns/etc
        [inl,gcl,synDepths]=getIPLdepth(curSM.syn.pos(:,3),curSM.syn.pos(:,1),curSM.syn.pos(:,2),[],[]);
        [inl,gcl,nepDepths]=getIPLdepth(curSM.nep.pos(:,3),curSM.nep.pos(:,1),curSM.nep.pos(:,2),[],[]);
        
        %find bpc inputs and get the partner cids
        bpcInIDs2=find(curSM.syn.preClass==7);
        bpcInIDs=find(preTypes{1}==7);
        %find the IDs of those synapses
        bpcInObIDs=curSM.syn.obID(bpcInIDs);
        bpcInsynIDs=curSM.syn.synID(bpcInIDs);
        
        %loop through those
        for j=1:length(bpcInObIDs)
            curobID=bpcInObIDs(j);
            curSynID=bpcInsynIDs(j);
            %get all the other parts of the ribbon synapse
            ribIDs=find(curTis.syn.edges(:,3)==curobID);
            if ~isempty(ribIDs)
                bpcCid=curTis.syn.edges(ribIDs(1),2);
                typeCheck=cid2type(bpcCid,curTis);
                if typeCheck{1}==7
                    partnercids=curTis.syn.edges(ribIDs,1);
                    notSelfCids=partnercids(partnercids~=curCid);
                    if sum(notSelfCids>0)>0
                        for k=1:length(notSelfCids)
                            curPartCid=notSelfCids(k);
                            curPartTypeDat=cid2type(curPartCid,curTis);
                            curPartType=curPartTypeDat{1};
                            curPartSubtype=curPartTypeDat{3};
                            %is it an alpha??
                            if curPartSubtype==23
                                curPartCid
                                p=p+1;
                            end
                        end
                    end
                end
            end
        end
    end
end
%% That didn't work

allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
vgcCidList=[2 3 4 5 10 13 14 20];
bpcPre=find(allPreTypes{1}==7);
vgcPost=find(ismember(curTis.syn.edges(:,1),vgcCidList)');
bpcPreIDs=curTis.syn.edges(bpcPre,3);
vgcPostIDs=curTis.syn.edges(vgcPost,3);
b2vIDs=intersect(bpcPreIDs,vgcPostIDs);
% find the RGCs that share ribbons
b2vEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),b2vIDs),:);
b2vEdgeTypes=cid2type(b2vEdges(:,1),curTis);
% set the target type to t-off-a
targType=[1 23];

% if length(targType)==2
%     [curTis.cells.type.typeNames{targType(1)} ...
% curTis.cells.type.subTypeNames{targType(1)}(targType(2))]
% else
% curTis.cells.type.typeNames{targType(1)}
% end

if length(targType)==1
    targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1),3);
else
    targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1)&b2vEdgeTypes{3}==targType(2),3);
end
targRibNames=curTis.obI.nameProps.names(targRibIds)';
targRibEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),targRibIds),:);
targRibLocs=curTis.obI.colStruc.anchors(targRibIds,:);

%There are 9 x17x shared ribbon inputs with alphas?!? That's insane
%three are with cid2

%% Go through each of the shared ribbons and get an idea of the inputs that are near.
synRad=30;
fAlpha=0.2;
vgcCidList=[2 3 4 5 13 14];
for n=1:length(targRibIds)
    curRibID=targRibIds(n);
    bpcCid=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,2);
    bpcCid=bpcCid(1);
    targets=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,1);
    vgcCid=intersect(targets,vgcCidList);
    rgcCid=targets(targets~=vgcCid);
    curSM=allSkels{find(vgcCidList==vgcCid)}.sm;
    curCid=curSM.cid;
    curSWC=allSkels{find(vgcCidList==vgcCid)}.swc;
    curSynID=find(curSM.syn.obID==curRibID);
    if ~isempty(curSynID)
        curSynPos=curSM.syn.pos(curSynID,:);
        curCloseNode=curSM.syn2Skel.closest(curSynID);
        dist2syns=curSM.syn2Skel.syn2SynDist(:,curSynID);
        %figure(); histogram(dist2syns,[0:100]);
        
        %Here is the spot to figure out whether the euc or lin dist is bigger
        pos = curSM.syn.pos;
        numSyn = size(pos,1);
        dif1 = pos(:,1)-pos(:,1)';
        dif2 = pos(:,2)-pos(:,2)';
        dif3 = pos(:,3)-pos(:,3)';
        dif = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
        dist2synsEuc=dif(:,curSynID);
        %figure(); scatter(dist2syns,dist2synsEuc);
        dist2synsComb=max(horzcat(dist2syns,dist2synsEuc),[],2);
        curPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
        curPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
        curInputs=find(curSM.syn.edges(:,1)==curCid);
        curOutputs=find(curSM.syn.edges(:,2)==curCid);
        closeSyns=find(dist2synsComb<synRad);
        nearInputs=intersect(closeSyns,curInputs);
        nearOutputs=intersect(closeSyns,curOutputs);
        
        debugTriad=1;
        if debugTriad
            winRad=20;
            df1=figure();
            hold on
            
            plottables=struct();
            curVGCfv=load([fvDir num2str(curCid) '.mat']);
            curRGCfv=load([fvDir num2str(rgcCid) '.mat']);
            curBPCfv=load([fvDir num2str(bpcCid) '.mat']);
            
            curVGCfv.fv.vertices=curVGCfv.fv.vertices(:,[2 3 1]);
            curRGCfv.fv.vertices=curRGCfv.fv.vertices(:,[2 3 1]);
            curBPCfv.fv.vertices=curBPCfv.fv.vertices(:,[2 3 1]);
            
            dp1=patch(curVGCfv.fv);
            dp1.FaceColor=[.25 0 .25];
            dp1.FaceAlpha=fAlpha;
            dp1.LineStyle='none';
            
            dp2=patch(curRGCfv.fv);
            dp2.FaceColor=[0 .25 .25];
            dp2.FaceAlpha=fAlpha;
            dp2.LineStyle='none';
            
            dp3=patch(curBPCfv.fv);
            dp3.FaceColor=[.25 .25 0];
            dp3.FaceAlpha=fAlpha;
            dp3.LineStyle='none';
            
            xlim([curSynPos(1)-winRad curSynPos(1)+winRad]);
            ylim([curSynPos(2)-winRad curSynPos(2)+winRad]);
            zlim([curSynPos(3)-winRad curSynPos(3)+winRad]);
            
            scatter3(curSynPos(1), curSynPos(2), curSynPos(3),50,'m+');
            %go through the nearby synapses and see if they are of interest
        end
        
        for j=1:length(nearInputs)
            curInput=nearInputs(j);
            curPreCid=curSM.syn.edges(curInput,2);
            curPreType=curPreTypes{1};
            curInputPos=curSM.syn.pos(curInput,:);
            curPreType(curInput)
            if curPreType(curInput)==7
                scatter3(curInputPos(1), curInputPos(2), curInputPos(3),50,'c+');
                if curPreCid==bpcCid
                    scatter3(curInputPos(1), curInputPos(2), curInputPos(3),50,'co','filled');
                end
            end
        end
        
        for j=1:length(nearInputs)
            curInput=nearInputs(j);
            curPreCid=curSM.syn.edges(curInput,2);
            curPostCid=curSM.syn.edges(curInput,1);
            curPreType=curPreTypes{1};
            curInputPos=curSM.syn.pos(curInput,:);
            if curPostCid==rgcCid
                scatter3(curInputPos(1), curInputPos(2), curInputPos(3),50,'ro');
            end
        end
    end
end


%% new triad avoidance code
%need to run top part first to get the ribbon locations
targRibNames;
targRibEdges;
targRibLocs;
for n=1:length(targRibIds)
    curRibID=targRibIds(n);
    bpcCid=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,2);
    bpcCid=bpcCid(1);
    targets=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,1);
    vgcCid=intersect(targets,vgcCidList);
    rgcCid=targets(targets~=vgcCid);
    curSM=allSkels{find(vgcCidList==vgcCid)}.sm;
    curCid=curSM.cid;
    curSWC=allSkels{find(vgcCidList==vgcCid)}.swc;
    curSynID=find(curSM.syn.obID==curRibID);
    if ~isempty(curSynID)
        
        %code here
    else
        'synapse not found on skeleton'
        n
        curRibID
        curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,:)
        %targRibEdges(n,:)
        %targRibLocs(n,:)
    end

end
