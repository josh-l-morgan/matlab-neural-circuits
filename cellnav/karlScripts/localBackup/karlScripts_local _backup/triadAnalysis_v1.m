
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

%% find the synapse types
% get the types of all the pre and post synaptic cells
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
% list of vgc that I have good skels for
vgcCidList=[2 3 4 5 13 14];
% get the IDs for the different kinds of synapses
bpcPre=find(allPreTypes{1}==7);
amcPre=find(allPreTypes{1}==0|allPreTypes{1}==8);
vgcPost=find(ismember(curTis.syn.edges(:,1),vgcCidList)');
bpcPreIDs=curTis.syn.edges(bpcPre,3);
amcPreIDs=curTis.syn.edges(amcPre,3);
vgcPostIDs=curTis.syn.edges(vgcPost,3);
% combine lists (like bpc-pre, vgc-post)
b2vIDs=intersect(bpcPreIDs,vgcPostIDs);
a2vIDs=intersect(amcPreIDs,vgcPostIDs);
%get edge data for everything that is a bipolar to vgc
b2vEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),b2vIDs),:);
% remove the vgcs from the edge data so that we can search for AMC partners
b2vEdges=b2vEdges(~ismember(b2vEdges(:,1),vgcCidList),:);
%get the type data so that we can search through for IDs
b2vEdgeTypes=cid2type(b2vEdges(:,1),curTis);

%% make the target ID lists (similar to "s" in neighborhood script)
targRibLocMat=struct();
% get the ribbons shared between rgcs and vgc

% 1. RGC partners
targType=[1];
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

% 2. AMC partners 
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

% figure to test the positions of the different shared synapses.
if 0
    figure();
    hold on
    colOrd=[9 13];
    for i=1:length(targRibLocMat)
        locList=targRibLocMat(i).targLocs;
        scatter3(locList(:,1),locList(:,2),locList(:,3),50,HCcolmap(colOrd(i),:),'filled')
    end
end


%% Density recovery

% I could multiply the 
%calc
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

%% same but for the same rgc
%% new triad avoidance code
%need to run top part first to get the ribbon locations
targRibNames;
targRibEdges;
targRibLocs;

groupNum=5;

%set up the density bins
binEdges=[0.5:0.5:30];
distBin=2;
g2g = zeros(length(targRibIds),groupNum,length(binEdges));
nearLength = g2g;
countG = zeros(groupNum,1);
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
            amcInIDs=intersect(nearInputs,find(curPreTypes{1}==0|curPreTypes{1}==8));
            synIDcell={bpcInIDs,amcInIDs,rgcOutSynIDs,rgcFFSynIDs,amcOutIDs};
            
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
f10=figure();
sgtitle('Density Recovery from VG3/AMC ribbons')
hold on;
subplot(1,3,1:2)
hold on
for i=1:5%[1 2 3 5]
    plotz(i)=plot(binEdges,test(i,:)./100,'LineWidth',2);
    %yline(means(i)/100);
end

for i=1:5%[1 2 3 5]
    %plot(binEdges,test(i,:)./100,'LineWidth',2);
    yline(means(i)/100,':');
end
legend({'bpcIn','amcIn','rgcOut','rgcFeedForward','amcOut'})
%legend({'bpcIn','amcIn','rgcOut','amcOut'})
xlabel('distance (um)');
ylabel('density (#syn/um)');
xlim([0 20])



subplot(1,3,3)
hold on
for i=1:5%[1 2 3 5]
    plotz(i)=plot(binEdges,test(i,:)./100,'LineWidth',2);
    %yline(means(i)/100);
end

for i=1:5%[1 2 3 5]
    %plot(binEdges,test(i,:)./100,'LineWidth',2);
    yline(means(i)/100,':');
end
%legend({'bpcIn','amcIn','rgcOut','rgcFeedForward','amcOut'})
xlim([0 5])
xlabel('distance (um)');
title('0-5 microns')
