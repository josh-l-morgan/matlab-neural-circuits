dsObjPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\dsObj.mat';
obiPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\obI.mat';
tisPath='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat';
fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
%fvLib='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
HCcolmap=load('HCcolmap.mat');
HCcolmap=HCcolmap.HCcolmap;
HCcolmap=HCcolmap./255;
colTest=0;
if colTest
    colFig=figure();
    scatter(1:15,repmat(1,[1,15]),50,HCcolmap,'filled');
end
%dsObjPath='G:\Data\MATLAB\0827_analysis\Volumes\0827\Merge\dsObj.mat';
load(dsObjPath);
load(obiPath);
load(tisPath);
%allHistDat=getHisto(tis,fvDir);
clear dsObjPath obiPath tisPath

skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\';
vgcCidList=[2 3 4 5 13 14];
loadSkel=0;
if ~exist('allSkels')
    if ~exist('D:\work\allSkels.mat')
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


%% parameters
skelDist=30;
synRad=30;
plotz=0;
plotz2=1;
%% loop
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
    %get the tips and forks and such from the edgeList
    curPred=curSWC.pred(curSWC.arbor2swcID)+1;
    curPredInv=zeros(size(curPred));
    curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
    curSM.pred=curPred;
    allEdges=curSM.arbor.edges;
    uniqueNodes=unique(allEdges(:));
    nodeCounts=histcounts(allEdges,uniqueNodes);
    tipIDs=find(nodeCounts==1);
    forkIDs=find(nodeCounts>2);
    distFromRoot=curSM.skel2skel.linDist(tipIDs,curSM.nep.seedNode);
    [srtd srtIdx] = sort(distFromRoot,'descend');
    tipIDsSrtd=tipIDs(srtIdx');
    %go through the tips and get the distances to the other synapses
    if plotz2
        f2=figure();
        hold on
    end
    for j=1:length(tipIDsSrtd)
        if plotz
            f1=figure();
            hold on
        end
        curTipID=tipIDsSrtd(j);
        dists2syns=curSM.syn2Skel.syn2SkelDist(:,curTipID);
        %figure(); histogram(dists2syns,[0:150]);
        tipPos=curSM.nep.pos(curTipID,:);
        tipZ=nepDepths(curTipID,:);
        if plotz
            if plotz
                figure(f1);
                scatter3(tipPos(1),tipPos(2),tipPos(3),250,'yp','filled');
            end
            curNode=curTipID;
            nearNode=find(curSM.skel2skel.linDist(curNode,:)<skelDist);
            nearTipNode=intersect(nearNode,tipIDs);
            nearNodePos=curSM.nep.pos(nearNode,:);
            nearTipNodePos=curSM.nep.pos(nearTipNode,:);
            if plotz
                scatter3(nearTipNodePos(:,1),nearTipNodePos(:,2),nearTipNodePos(:,3),25,'m+');
                [bboxa,bboxb]=bounds([nearNodePos;tipPos]);
                xlim([bboxa(1)-5 bboxb(1)+5]);
                ylim([bboxa(2)-5 bboxb(2)+5]);
                zlim([bboxa(3)-5 bboxb(3)+5]);
                %scatter3(nearNodePos(:,1),nearNodePos(:,2),nearNodePos(:,3),5,'k.');
                fvPatch=patch(curCidFv);
                fvPatch.EdgeColor='none';
                fvPatch.FaceColor=[.2 .2 .2];
                fvPatch.FaceAlpha=.1;
            end
        end
        %find all nearby synapses
        closeSyns=find(dists2syns<synRad);
        %get positions of close syns
        closeSynPos=curSM.syn.pos(closeSyns,:);
        if 1
            %find the kinds of synapses so we can plot them correctly
            closeIn=intersect(curSM.syn.isIn,closeSyns);
            closeOut=intersect(curSM.syn.isOut,closeSyns);
            closeOutTypes=postTypes{1}(closeOut);
            closeInTypes=preTypes{1}(closeIn);
            
            
            amcIn=closeIn(closeInTypes==8|closeInTypes==0);
            bpcIn=closeIn(closeInTypes==7);
            amcOut=closeOut(closeOutTypes==8|closeOutTypes==0);
            rgcOut=closeOut(closeOutTypes==1);
            if plotz2
                figure(f2);
                if ~isempty(amcIn)
                    scatter(dists2syns(amcIn),repmat(j,[length(dists2syns(amcIn)) 1]), ...
                        20,'mv','filled');
                end
                if ~isempty(bpcIn)
                    scatter(dists2syns(bpcIn),repmat(j,[length(dists2syns(bpcIn)) 1]), ...
                        20,'cv','filled');
                end
                if ~isempty(amcOut)
                    scatter(dists2syns(amcOut),repmat(j,[length(dists2syns(amcOut)) 1]), ...
                        20,'m^','filled');
                end
                if ~isempty(rgcOut)
                    scatter(dists2syns(rgcOut),repmat(j,[length(dists2syns(rgcOut)) 1]), ...
                        20,'g^','filled');
                end
                
            end
            if plotz
                
                synScale=75;
                scaleSyns=0;
                synSizes=[synScale synScale synScale synScale];
                if scaleSyns
                    synSizes(1)=synScale./dists2syns(amcIn);
                    synSizes(2)=synScale./dists2syns(bpcIn);
                    synSizes(3)=synScale./dists2syns(amcOut);
                    synSizes(4)=synScale./dists2syns(rgcOut);
                end
                
                closeSynCols=repmat([0 0 0],[length(closeSyns) 1]);
                closeSynCols(ismember(closeSyns,closeIn),:)=repmat(HCcolmap(4,:),[length(closeIn) 1]);
                closeSynCols(ismember(closeSyns,closeOut),:)=repmat(HCcolmap(9,:),[length(closeOut) 1]);
                %plot'em
                %closeSynScat=scatter3(closeSynPos(:,1),closeSynPos(:,2),closeSynPos(:,3),100./dists2syns(closeSyns),closeSynCols,'o','filled');
                closeSynScat=scatter3(curSM.syn.pos(amcIn,1),curSM.syn.pos(amcIn,2),curSM.syn.pos(amcIn,3),synSizes(1),'mv','filled');
                closeSynScat=scatter3(curSM.syn.pos(bpcIn,1),curSM.syn.pos(bpcIn,2),curSM.syn.pos(bpcIn,3),synSizes(2),'cv','filled');
                closeSynScat=scatter3(curSM.syn.pos(amcOut,1),curSM.syn.pos(amcOut,2),curSM.syn.pos(amcOut,3),synSizes(3),'m^','filled');
                closeSynScat=scatter3(curSM.syn.pos(rgcOut,1),curSM.syn.pos(rgcOut,2),curSM.syn.pos(rgcOut,3),synSizes(4),'g^','filled');
                %closeSynScat.Marker=cell(size(closeSynPos,1),1);
            %           closeSynScat.CData=
                %             for l=1:size(closeSynPos,1)
                %                 curSynPos=closeSynPos(l,:);
                %                 p=plot3([curSynPos(3) tipPos(3)],[curSynPos(1) tipPos(1)], ...
                %                     [curSynPos(2) tipPos(2)],'--','LineWidth',25/dists2syns(closeSyns(l)));
                %                 p.LineWidth
                %             end
            end
        end
        %%
        drawnow
    end
    title(num2str(curCid))
end

%% Do the branches go back and forth, or do they stay put?

branchLengthCutoff=10;
plotSyns=1;

if plotz3
    curFig=figure();
    hold on;
end

for i=1:length(allSkels)
    subplot(2,3,i);
    hold on
    curSM=allSkels{i}.sm;
    curCid=curSM.cid;
    curSWC=allSkels{i}.swc;
    % get depth info for nodes/syns/etc
    [inl,gcl,synDepths]=getIPLdepth(curSM.syn.pos(:,3),curSM.syn.pos(:,1),curSM.syn.pos(:,2),[],[]);
    [inl,gcl,nepDepths]=getIPLdepth(curSM.nep.pos(:,3),curSM.nep.pos(:,1),curSM.nep.pos(:,2),[],[]);
    %get the tips and forks and such from the edgeList
    curPred=curSWC.pred(curSWC.arbor2swcID)+1;
    curPredInv=zeros(size(curPred));
    curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
    curSM.pred=curPred;
    allEdges=curSM.arbor.edges;
    uniqueNodes=unique(allEdges(:));
    nodeCounts=histcounts(allEdges,uniqueNodes);
    tipIDs=find(nodeCounts==1);
    forkIDs=find(nodeCounts>2);
    distFromRoot=curSM.skel2skel.linDist(tipIDs,curSM.nep.seedNode);
    [srtd srtIdx] = sort(distFromRoot,'descend');
    tipIDsSrtd=tipIDs(srtIdx');
    plotz3=1;
    
    preTypes=cid2type(curSM.syn.edges(:,2),curTis);
    postTypes=cid2type(curSM.syn.edges(:,1),curTis);
    vgcPre=curSM.syn.edges(:,2)==curCid;
    
    usedNodes=[];
    briter=1;
    liter=0;
    for j=1:length(tipIDsSrtd)
        curTipID=tipIDsSrtd(j);
        tipPos=curSM.nep.pos(curTipID,:);
        tipZ=nepDepths(curTipID,:);
        curNode=curTipID;
        curTipDist=0;
        curBranchList=curNode;
        abandon=0;
        lengthsList=0;
        cumLength=0;
        while ~abandon
            nextNode=curPredInv(curNode);
            usedNodes=[usedNodes;curNode];
            if nextNode==0 | ismember(nextNode,usedNodes)
                abandon=1;
                break
            end
            %if ismember(curNode,forkIDs)
            %    abandon=1;
            %    break
            %end
            curNodeRad=curSM.arbor.nodes.rad(curNode);
            curEdgeLength=curSM.nep.props.nodeLength(curNode);
            curTipDist=curTipDist+curEdgeLength;
            cumLength=cumLength+curEdgeLength;
            lengthsList=[lengthsList;cumLength];
            curNode=nextNode;
            curBranchList=[curBranchList;curNode];
            %plot out the thing
        end
        
        
        if curTipDist>branchLengthCutoff
            %length( curBranchList )
            %curTipDist
            %j
            %for k=1:length(curBranchList)
            %    depNode=curBranchList(k);
            curDepths=nepDepths(curBranchList);
            curEdges=lengthsList;
            
            
            plotSyns=1;
            synSize=25;
            if plotSyns
                for nodeIt=1:length(curBranchList)
                    curNode=curBranchList(nodeIt);
                    if ismember(curNode,curSM.syn2Skel.closest)
                        [ind, m]=find(curSM.syn2Skel.closest==curNode);
                        preType=preTypes{1}(ind);
                        postType=postTypes{1}(ind);
                        if vgcPre(ind)
                            %output
                            if postType==1
                                scatter(liter+curEdges(nodeIt),curDepths(nodeIt),synSize,'g^','filled');
                            else
                                scatter(liter+curEdges(nodeIt),curDepths(nodeIt),synSize,'m^','filled');
                            end
                            
                        else
                            %input
                            if preType==7
                                scatter(liter+curEdges(nodeIt),curDepths(nodeIt),synSize,'cv','filled');
                            else
                                scatter(liter+curEdges(nodeIt),curDepths(nodeIt),synSize,'mv','filled');
                            end
                        end
                    end
                end
            end
            plot(liter+curEdges,curDepths);
            drawnow
            liter=curTipDist+liter+10;
            %briter=briter+1;
        end
        
        
    end
    title(num2str(curCid));
    ylabel('IPL depth %');
    xlabel('Arbor Branch Segments');
    xticks([]);
    yline(0.47)
    xlim([0 liter]);
    ylim([0 0.7]);
end
%% Notes
% How to get the density?
%
%
%
%
%
%
%
%
%
%
%
%

