
%Initial Parameters
vgcCidList=[2 3 4 5 13 14];
loadAll=0; loadSkel=0;
lc=20;

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
    %HCcolmap=newCols./255;
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
        if  exist('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\allSkels.mat');
            load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\allSkels.mat');
        elseif exist('C:\Users\karlf\Documents\Lab\data\allSkels.mat')
            load('C:\Users\karlf\Documents\Lab\data\allSkels.mat');
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
            for i=1:length(allSkels)
                curSkel=allSkels{i};
                curSWC=nep2swc(curSkel.sm.nep);
                allSkels{i}.swc=curSWC;
                
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
debugPlot=0;
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
vgcCidList=[2 3 4 5 13 14];
bpcPre=find(allPreTypes{1}==7);
amcPre=find(allPreTypes{1}==0|allPreTypes{1}==8);
vgcPost=find(ismember(curTis.syn.edges(:,1),vgcCidList)');
bpcPreIDs=curTis.syn.edges(bpcPre,3);
amcPreIDs=curTis.syn.edges(amcPre,3);
vgcPostIDs=curTis.syn.edges(vgcPost,3);
b2vIDs=intersect(bpcPreIDs,vgcPostIDs);
a2vIDs=intersect(amcPreIDs,vgcPostIDs);
% find the RGCs that share ribbons
b2vEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),b2vIDs),:);
b2vEdges=b2vEdges(~ismember(b2vEdges(:,1),vgcCidList),:);
b2vEdgeTypes=cid2type(b2vEdges(:,1),curTis);

testList=curTis.syn.edges(find(allPreTypes{1}'==7&ismember(curTis.syn.edges(:,1),vgcCidList)),3);

% set the target type to t-off-a
% if length(targType)==2
%     [curTis.cells.type.typeNames{targType(1)} ...
% curTis.cells.type.subTypeNames{targType(1)}(targType(2))]
% else
% curTis.cells.type.typeNames{targType(1)}
% end
targRibLocMat=struct();
% get the ribbons shared between rgcs and vgc
targType=[1 23];
if length(targType)==1
    targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1),3);
else
    targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1)&b2vEdgeTypes{3}==targType(2),3);
end
targRibNames=curTis.obI.nameProps.names(targRibIds)';
targRibEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),targRibIds),:);
targRibLocs=curTis.obI.colStruc.anchors(targRibIds,:);
targRibLocMat(1).name='all rgc';
targRibLocMat(1).targIDs=targRibIds;
targRibLocMat(1).targNames=targRibNames;
targRibLocMat(1).targEdges=targRibEdges;
targRibLocMat(1).targLocs=targRibLocs;

%get the ones for amcs
targType=[0];
if length(targType)==1
    targRibIds_a=b2vEdges(b2vEdgeTypes{1}==targType(1),3);
else
    targRibIds_a=b2vEdges(b2vEdgeTypes{1}==targType(1)&b2vEdgeTypes{3}==targType(2),3);
end
targRibNames_a=curTis.obI.nameProps.names(targRibIds_a)';
targRibEdges_a=curTis.syn.edges(ismember(curTis.syn.edges(:,3),targRibIds_a),:);
targRibLocs_a=curTis.obI.colStruc.anchors(targRibIds_a,:);

targType=[8];
if length(targType)==1
    targRibIds_b=b2vEdges(b2vEdgeTypes{1}==targType(1),3);
else
    targRibIds_b=b2vEdges(b2vEdgeTypes{1}==targType(1)&b2vEdgeTypes{3}==targType(2),3);
end
targRibNames_b=curTis.obI.nameProps.names(targRibIds_b)';
targRibEdges_b=curTis.syn.edges(ismember(curTis.syn.edges(:,3),targRibIds_b),:);
targRibLocs_b=curTis.obI.colStruc.anchors(targRibIds_b,:);

targRibIds=vertcat(targRibIds_a,targRibIds_b);
targRibNames=vertcat(targRibNames_a,targRibNames_b);
targRibEdges=vertcat(targRibEdges_a,targRibEdges_b);
targRibLocs=vertcat(targRibLocs_a,targRibLocs_b);

targRibLocMat(2).name='all amc';
targRibLocMat(2).targIDs=targRibIds;
targRibLocMat(2).targNames=targRibNames;
targRibLocMat(2).targEdges=targRibEdges;
targRibLocMat(2).targLocs=targRibLocs;

if 0
    figure();
    hold on
    colOrd=[9 13];
    for i=1:length(targRibLocMat)
        locList=targRibLocMat(i).targLocs;
        scatter3(locList(:,1),locList(:,2),locList(:,3),50,HCcolmap(colOrd(i),:),'filled')
    end
end
% get the ones for that specific RGC
% targType=[1];
% if length(targType)==1
%     targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1),3);
% else
%     targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1)&b2vEdgeTypes{3}==targType(2),3);
% end
% targRibNames=curTis.obI.nameProps.names(targRibIds)';
% targRibEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),targRibIds),:);
% targRibLocs=curTis.obI.colStruc.anchors(targRibIds,:);
% targRibLocMat(1).name='all rgc';
% targRibLocMat(1).targIDs=targRibIds;
% targRibLocMat(1).targNames=targRibNames;
% targRibLocMat(1).targEdges=targRibEdges;
% targRibLocMat(1).targLocs=targRibLocs;

%% Go through each of the shared ribbons and get an idea of the inputs that are near.
while 0
fAlpha=0.2;
lc=20;
binWid=0.1;
vgcCidList=[2 3 4 5 13 14];
figure(); scatter(1:15,repmat(1,[1 15]),100,HCcolmap,'o','filled');
%synTypColMat=ones(2,10);
synTypColMat=[[4 15 1 1 1 1 1 14 4 1 13 12 11 11]; ...
    [5 9  1 1 1 1 1 3  5 1 13 12 11 11]];
synTypPosMat=[[1.1 1 1 1 1 1 1 1.2 1.1 1 1 1 1 1]; ...
    [.9 .8  1 1 1 1 1 .7  .9 1 1 1 1 1]];
plotSub=1;

f0=figure();
hold on
% f1=figure();
% hold on

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
    
    %figureout the influence, etc.
    curDistMat=curSM.syn2Skel.syn2SynDist;
    preTypes=cid2type(curSM.syn.edges(:,2),curTis);
    postTypes=cid2type(curSM.syn.edges(:,1),curTis);
    infMat=exp(-curDistMat./lc);
    %bpcSynIDs=find(preTypes{1}==7);
    %amcSynIDs=find(curSM.syn.edges(:,1)==curCid&(preTypes{1}==8|preTypes{1}==0)');
    
    %     vMat=zeros(size(curSM.syn2Skel.syn2SynDist));
    %     vMat(bpcSynIDs,:)=1;
    %     vMat(amcSynIDs,:)=-1;
    %     infMat=vMat.*exp(-curDistMat./lc);
    %     EIbalMat=sum(infMat,1);
    
    
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
        infMatComb=exp(-dist2synsComb./lc);
        curPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
        curPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
        curInputs=find(curSM.syn.edges(:,1)==curCid);
        curOutputs=find(curSM.syn.edges(:,2)==curCid);
        closeSyns=find(dist2synsComb<synRad);
        nearInputs=intersect(closeSyns,curInputs);
        nearOutputs=intersect(closeSyns,curOutputs);
        
        %%
        %         curDistMat=curSM.syn2Skel.syn2SynDist;
        %         vMat=zeros(size(curSM.syn2Skel.syn2SynDist));
        %         vMat(bpcSynIDs,:)=1;
        %         vMat(amcSynIDs,:)=-1;
        %         infMat=vMat.*exp(-curDistMat./lc);
        %         EIbalMat=sum(infMat,1);
        plotSynDat=1;
        if plotSynDat
            %             figure();
            %             hold on
            %             title(num2str(curSynID))
            for o=1:length(closeSyns)
                curCloseSyn=closeSyns(o);
                preType=preTypes{1}(o);
                preSub=preTypes{3}(o);
                postType=postTypes{1}(o);
                postSub=postTypes{3}(o);
                if curSM.syn.edges(o,1)==curCid
                    dir=1; %input
                    col=preType+1;
                    figure(f0)
                    scatter(dist2synsComb,synTypPosMat(dir,col)+rand()/10,25, ...
                        HCcolmap(synTypColMat(dir,col),:),'v','filled');
                    %figure(f1)
                    %scatter(dist2synsComb,synTypPosMat(dir,col)+rand()/10,25, ...
                    %    HCcolmap(synTypColMat(dir,col),:),'v','filled');
                    %drawnow
                else
                    dir=2; %output
                    col=postType+1;
                    figure(f0)
                    scatter(dist2synsComb,synTypPosMat(dir,col)+rand()/10,25, ...
                        HCcolmap(synTypColMat(dir,col),:),'^','filled');
                    %                     scatter((1-infMatComb(o)),synTypPosMat(dir,col)+rand()/10,25, ...
                    %                         HCcolmap(synTypColMat(dir,col),:),'^','filled');
                    %figure(f1)
                    %scatter(dist2synsComb,synTypPosMat(dir,col)+rand()/10,25, ...
                    %    HCcolmap(synTypColMat(dir,col),:),'^','filled');
                    %drawnow
                end
                %scatter((10-10*infMatComb(o)),0.9,25,HCcolmap(synTypColMat(dir,col),:),'^','filled')
            end
            drawnow
            ylim([0 2])
        end
        %%
        debugTriad=0;
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
end

end
%% new triad avoidance code
while 0
%need to run top part first to get the ribbon locations
targRibNames;
targRibEdges;
targRibLocs;

%set up the density bins
binEdges=[0.25:0.25:50];
distBin=1;
g2g = zeros(length(targRibIds),4,length(binEdges));
nearLength = g2g;
countG = zeros(4,1);
for n=1:length(targRibIds)
    curRibID=targRibIds(n);
    bpcCid=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,2);
    bpcCid=bpcCid(1);
    targets=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,1);
    vgcCid=intersect(targets,vgcCidList);
    if ~isempty(vgcCid)
        vgcCid=vgcCid(1);
        rgcCid=targets(~ismember(targets,vgcCid));
        curSM=allSkels{find(vgcCidList==vgcCid)}.sm;
        curCid=curSM.cid;
        curSWC=allSkels{find(vgcCidList==vgcCid)}.swc;
        curSynID=find(curSM.syn.obID==curRibID);
        if ~isempty(curSynID)
            
            %         This is Josh's density code
            %         for b = 1:(length(bRange))
            %             [ y x] = find((synDG2>=(bRange(b)-distBin/2)) & (synDG2<(bRange(b)+distBin/2)));
            %             g2g(g,hT,b) =  g2g(g,hT,b) + length(x);
            %         end
            %code here
            curSynPos=curSM.syn.pos(curSynID,:);
            curCloseNode=curSM.syn2Skel.closest(curSynID);
            dist2syns=curSM.syn2Skel.syn2SynDist(curSynID,:);
            dist2syns(curSynID) = inf;
            dist2skels=curSM.syn2Skel.syn2SkelDist(curSynID,:);
            nodeLengths=curSM.nep.props.nodeLength;
            pos = curSM.syn.pos;
            numSyn = size(pos,1);
            dif1 = pos(:,1)-pos(:,1)';
            dif2 = pos(:,2)-pos(:,2)';
            dif3 = pos(:,3)-pos(:,3)';
            dif = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
            dist2synsEuc=dif(:,curSynID);
            %figure(); scatter(dist2syns,dist2synsEuc);
            dist2synsComb=max(horzcat(dist2syns',dist2synsEuc),[],2);
            infMatComb=exp(-dist2synsComb./lc);
            curPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
            curPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
            curInputs=find(curSM.syn.edges(:,1)==curCid);
            curOutputs=find(curSM.syn.edges(:,2)==curCid);
            closeSyns=find(dist2synsComb<synRad);
            nearInputs=intersect(closeSyns,curInputs);
            nearOutputs=intersect(closeSyns,curOutputs);
            
            %go through the close synapses and make the density of the
            %different types
            bpcInIDs=intersect(nearInputs,find(curPreTypes{1}==7));
            amcInIDs=intersect(nearInputs,find(curPreTypes{1}==0|curPreTypes{1}==8));
             
            rgcOutSynIDs=intersect(nearOutputs,find(curPostTypes{1}==1));
            amcOutIDs=intersect(nearOutputs,find(curPreTypes{1}==0|curPostTypes{1}==8));
            
            synIDcell={bpcInIDs,amcInIDs,rgcOutSynIDs,amcOutIDs};
                        
            for hT = 1:length(synIDcell) %check against all types
                %get the distances from those synapses to all the other synapses
                countG(hT)=countG(hT)+1;
                synDG2 = dist2synsComb(synIDcell{hT});
                for b = 1:(length(binEdges))
                    %get the indices and items for everything within a binwidth of that range
                    [ y x] = find((synDG2>=(binEdges(b)-distBin/2)) & (synDG2<(binEdges(b)+distBin/2)));
                    %add the number of things at that range to the data
                    g2g(n,hT,b) =  g2g(n,hT,b) + length(x);
                end
            end
            for b = 1:(length(binEdges)) %get length for density
                %find all the skeleton nodes at that range from the point
                [ y x] = find((dist2skels>=(binEdges(b)-distBin/2)) & (dist2skels<(binEdges(b)+distBin/2)));
                %that is the denominator for the density!
                %how many synapses at that range / sum of the nodes at that range
                nearLength(n,:,b) = nearLength(n,:,b) + sum(nodeLengths(x));
            end
            
            
            
        else
            'synapse not found on skeleton'
            n
            curRibID
            curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,:)
            %targRibEdges(n,:)
            %targRibLocs(n,:)
        end
    end
end

g2g = g2g ./ permute(repmat(countG,[1 size(g2g,1) size(g2g,3)]),[2 1 3]);
nearLength = nearLength ./ permute(repmat(countG,[1 size(g2g,1) size(g2g,3)]),[2 1 3]);
g2g(isnan(g2g))=0;
g2g=g2g.*100;
nearLength(isnan(nearLength))=0;
g2gN = g2g ./ nearLength;
g2gN(isnan(g2gN))=0;


test=sum(g2gN,1);
test=squeeze(test);
test=test./countG;
means=mean(test(:,round(0.5*length(binEdges)):round(0.9*length(binEdges))),2);
figure();
hold on;
for i=1:4
    plotz(i)=plot(binEdges,test(i,:)./100,'LineWidth',2);
    %yline(means(i)/100);
end

for i=1:4
    %plot(binEdges,test(i,:)./100,'LineWidth',2);
    yline(means(i)/100,':');
end

legend([plotz(:)],{'bpcIn','amcIn','rgcOut','amcOut'})

%% diff
testRat=test_rgcPart-test_amcPart;
figure();
hold on;
for i=1:4
    plot(binEdges,testRat(i,:)./100,'LineWidth',2);
end
legend({'bpcIn','amcIn','rgcOut','amcOut'})
end
%% same but for the same rgc
%% new triad avoidance code
%need to run top part first to get the ribbon locations
targRibNames;
targRibEdges;
targRibLocs;

groupNum=5;
descr='allRGC';
%set up the density bins
binEdges=[0.5:0.25:30];
distBin=2;
g2g = zeros(length(targRibIds),groupNum,length(binEdges));
g2g2 = g2g;
nearLength = g2g;
countG = zeros(groupNum,1);
graphDenom=zeros(groupNum,1);
ribSynPos = zeros(length(targRibIds),3);
whichVG3 = zeros(groupNum,1);
tots=zeros(1,groupNum);
if debugPlot
    fdp=figure();
end
for n=1:length(targRibIds)
    curRibID=targRibIds(n);
    bpcCid=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,2);
    bpcCid=bpcCid(1);
    targets=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,1);
    vgcCid=intersect(targets,vgcCidList);
    if ~isempty(vgcCid)
        vgcCid=vgcCid(1);
        rgcCid=targets(~ismember(targets,vgcCid));
        curSM=allSkels{find(vgcCidList==vgcCid)}.sm;
        curCid=curSM.cid;
        curSWC=allSkels{find(vgcCidList==vgcCid)}.swc;
        curSynID=find(curSM.syn.obID==curRibID);
        if ~isempty(curSynID)
            
            %         This is Josh's density code
            %         for b = 1:(length(bRange))
            %             [ y x] = find((synDG2>=(bRange(b)-distBin/2)) & (synDG2<(bRange(b)+distBin/2)));
            %             g2g(g,hT,b) =  g2g(g,hT,b) + length(x);
            %         end
            %code here
            curSynPos=curSM.syn.pos(curSynID,:);
            ribSynPos(n,:)=curSynPos(1,:);
            whichVG3(n)=curCid;
            curCloseNode=curSM.syn2Skel.closest(curSynID);
            dist2syns=curSM.syn2Skel.syn2SynDist(curSynID,:);
            dist2syns(curSynID) = inf;
            %randomSkelNode=randi([1 max(curSM.syn2Skel.closest)])
            %dist2syns=curSM.syn2Skel.syn2SkelDist(:,randomSkelNode)';
            dist2skels=curSM.syn2Skel.syn2SkelDist(curSynID,:);
            %dist2skels=curSM.skel2skel.linDist(randomSkelNode,:);

            nodeLengths=curSM.nep.props.nodeLength;
            pos = curSM.syn.pos;
            numSyn = size(pos,1);
            dif1 = pos(:,1)-pos(:,1)';
            dif2 = pos(:,2)-pos(:,2)';
            dif3 = pos(:,3)-pos(:,3)';
            dif = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
            dist2synsEuc=dif(:,curSynID);
            %figure(); scatter(dist2syns,dist2synsEuc);
            dist2synsComb=max(horzcat(dist2syns',dist2synsEuc),[],2);
            infMatComb=exp(-dist2synsComb./lc);
            curPostCids=curSM.syn.edges(:,1);
            curPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
            curPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
            curInputs=find(curSM.syn.edges(:,1)==curCid);
            curOutputs=find(curSM.syn.edges(:,2)==curCid);
            closeSyns=find(dist2synsComb<synRad);
            nearInputs=intersect(closeSyns,curInputs);
            nearOutputs=intersect(closeSyns,curOutputs);
            
            %go through the close synapses and make the density of the
            %different types
            bpcInIDs=intersect(nearInputs,find(curPreTypes{1}==7));
            amcOutIDs=intersect(nearOutputs,find(curPostTypes{1}==0|curPostTypes{1}==8));
            rgcOutSynIDs=intersect(nearOutputs,find(curPostTypes{1}==1));
            rgcFFSynIDs=intersect(rgcOutSynIDs,find(ismember(curPostCids,rgcCid)));
            rgcOutSynIDs=rgcOutSynIDs(~ismember(rgcOutSynIDs,rgcFFSynIDs));
            amcInIDs=intersect(nearInputs,find(curPreTypes{1}==0|curPreTypes{1}==8));
            synIDcell={bpcInIDs,amcInIDs,rgcOutSynIDs,rgcFFSynIDs,amcOutIDs};
            tots=tots+[length(bpcInIDs),length(amcInIDs),length(rgcOutSynIDs),length(rgcFFSynIDs),length(amcOutIDs)];
            for hT = 1:length(synIDcell) %check against all types
                %get the distances from those synapses to all the other synapses
                countG(hT)=countG(hT)+1;
                graphDenom(n,hT)=length(synIDcell{hT});
                synDG2 = dist2synsComb(synIDcell{hT});
                if debugPlot
                    figure(fdp);
                    scatter(synDG2,repmat(hT,size(synDG2)),20,'filled');
                    hold on
                end
                for b = 1:(length(binEdges))
                    %get the indices and items for everything within a binwidth of that range
                    [ y x ] = find((synDG2>=(binEdges(b)-distBin/2)) & (synDG2<(binEdges(b)+distBin/2)));
                    %add the number of things at that range to the data
                    g2g(n,hT,b) =  g2g(n,hT,b) + length(x);
                    g2g2(n,hT,b) = g2g2(n,hT,b) + length(x); %/length(synIDcell{hT}) ;
                end
            end
            for b = 1:(length(binEdges)) %get length for density
                %find all the skeleton nodes at that range from the point
                [ y x] = find((dist2skels>=(binEdges(b)-distBin/2)) & (dist2skels<(binEdges(b)+distBin/2)));
                %that is the denominator for the density!
                %how many synapses at that range / sum of the nodes at that range
                nearLength(n,:,b) = nearLength(n,:,b) + sum(nodeLengths(x));
            end
            
            
            
        else
            'synapse not found on skeleton'
            n
            curRibID
            curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,:)
            %targRibEdges(n,:)
            %targRibLocs(n,:)
        end
    end
    if debugPlot
                            yticks([1:5]);
                    yticklabels({'bpcIn','amcIn','rgcOut','rgcFeedForward','amcOut'});
    end
end
countMat=g2g;
denomTot=sum(countMat,[3]);
g2g2(isnan(g2g2))=0;
g2g22=g2g2./sum(g2g2,3);
g2g22(isnan(g2g22))=0;
g2g = g2g ./ permute(repmat(countG,[1 size(g2g,1) size(g2g,3)]),[2 1 3]);
nearLength = nearLength ./ permute(repmat(countG,[1 size(g2g,1) size(g2g,3)]),[2 1 3]);
g2g(isnan(g2g))=0;
g2g=g2g.*100;
nearLength(isnan(nearLength))=0;
g2gN = g2g ./ nearLength;
g2gN(isnan(g2gN))=0;

if debugPlot
    f9=figure(); hold on; for i=1:size(g2gN,1); scatter(i,squeeze(g2gN(i,1,1))./100,25); drawnow; end;
    title('bpcIn Density at 0um+bin/2 for each shared ribbon');
end

test=sum(g2gN,1);
test=squeeze(test);
%test=test./countG;
totalDens=sum(test,2);
means=mean(test(:,round(0.5*length(binEdges)):round(0.9*length(binEdges))),2);
f10=figure();
sgtitle(['percent of synapse density near shared bipolar input with ' descr ' (n=' num2str(countG(1)) ')']);
hold on;
subplot(1,3,1:2)
hold on
for i=1:5%[1 2 3 5]% 
    plotz(i)=plot(binEdges,test(i,:).*100./totalDens(i),'LineWidth',2);
    %yline(means(i)/100);
end

% for i=1:5%[1 2 3 5]
%     %plot(binEdges,test(i,:)./100,'LineWidth',2);
%     yline(means(i)/100,':');
% end
legend({'bpcIn','amcIn','rgcOut','rgcFeedForward','amcOut'})
%legend({'bpcIn','amcIn','rgcOut','amcOut'})
xlabel('distance (um)');
ylabel('% of total synapse density');
xlim([0 20])

subplot(1,3,3)
hold on
for i=1:5%[1 2 3 5]%
    plotz(i)=plot(binEdges,test(i,:).*100./totalDens(i),'LineWidth',2);
    %yline(means(i)/100);
end

% for i=1:5%[1 2 3 5]
%     %plot(binEdges,test(i,:)./100,'LineWidth',2);
%     yline(means(i)/100,':');
% end
%legend({'bpcIn','amcIn','rgcOut','rgcFeedForward','amcOut'})
xlim([0 5])
xlabel('distance (um)');
title('0-5 microns')
