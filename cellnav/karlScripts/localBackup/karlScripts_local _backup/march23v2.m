skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
vgcCidList=[2 3 4 5 13 14];
loadSkel=0;
if ~exist('allSkels')
    allSkels=cell(1,6);
    for i=1:length(vgcCidList)
        curCid=vgcCidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end
end
if ~exist('curTis')
    global glob tis
    curTis=tis;
end
%% test new func
loadAll=0;
if loadAll
    for i=1:length(allSkels)
        curSkel=allSkels{i};
        curSWC=nep2swc(curSkel.sm.nep);
        allSkels{i}.swc=curSWC;
    end
end
%% STRUCTS

%% PARAMS
skipPlot=1;
randomizeTipOrder=1;
stopFork=0;
reuseNodes=0;
draw3dBranches=0;
drawBranches=0;
drawDensity=0;
throwSmall=0;
%% Data processing loop
tipDistCutoff=100;
tipDatDist={};
distRef={"bpc","amc","syn","rgc"}; %"rand","rand","rand","rand","rand"};
 randRep="rand";
 distRef2=repmat(randRep,[10000 1]);
 distRef2=cellstr(distRef2);
 %distRef={distRef2{:} distRef{:}};
distRef=distRef2;

if skipPlot~=1
f5=figure();
f6=figure();
end

uberDat=[];
uberRawDat={};
for i=1:length(distRef)
    distRefs=[];
    for j=1:length(allSkels)
        tipDatCid={};
        %get the skeleton data set up for going from the tips up the
        %branches to the root.
        curSkel=allSkels{j};
        curSM=curSkel.sm;
        curSWC=curSkel.swc;
        curCid=curSM.cid;
        curPred=curSWC.pred(curSWC.arbor2swcID)+1;
        curPredInv=zeros(size(curPred));
        curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
        curSM.pred=curPred;
        allEdges=curSM.arbor.edges;
        uniqueNodes=unique(allEdges(:));
        nodeCounts=histcounts(allEdges,uniqueNodes);
        tipIDs=find(nodeCounts==1);
        forkIDs=find(nodeCounts>2);
        closestNodes=curSM.syn2Skel.closest;
        synDistFromRoot=curSM.skel2skel.linDist(closestNodes,1);
        distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
        [srtd srtIdx] = sort(distFromRoot,'descend');
        tipIDsSrtd=tipIDs(srtIdx');
        % get all the synapse types and their closest nodes
        
        % get the distances from those nodes to the tips
        %synTipDists=curSM.syn2Skel.syn2SkelDist(:,tipIDs);
        synTipDists=curSM.skel2skel.linDist(closestNodes,tipIDs);
        randomSkelNodes=randperm(length(curSM.skel2skel.linDist(1,:)),length(curSM.syn2Skel.closest));
        %randomSkelNodes=randperm(length(curSM.skel2skel.linDist(1,:)),500);
        %randomSkelNodes=randomSkelNodes(find(~ismember(randomSkelNodes,closestNodes)));
        randomSkelNodes=randomSkelNodes(~ismember(randomSkelNodes,tipIDs));
        randomSkelDists=curSM.skel2skel.linDist(randomSkelNodes,tipIDs);
        randDistFromRoot=curSM.skel2skel.linDist(randomSkelNodes,1);
        if 0
            f3=figure();
            swarmchart(repmat(1,length(closestNodes)),closestNodes)
            hold on
            swarmchart(repmat(2,length(randomSkelNodes)),randomSkelNodes)
            swarmchart(repmat(3,length(tipIDs)),tipIDs)
        end
        %         figure(f2);
        %         histogram(closestNodes);
        %         histogram(randomSkelNodes);
        
        synPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
        synPreCid=curSM.syn.edges(:,2);
        synPostCid=curSM.syn.edges(:,1);
        synPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
        
        curRef=distRef{i};
        
        if curRef=="bpc"
            synTipDists=synTipDists(synPreTypes{1}==7,:);
        elseif curRef=="amc"
            synTipDists=synTipDists(synPostCid==curCid&synPreTypes{1}'~=7,:);
        elseif curRef=="syn"
            1;
        elseif curRef=="rgc"
            synTipDists=synTipDists(synPreCid==curCid&synPostTypes{1}'==1,:);
        elseif curRef=="rand"
            %synTipDists=randomSkelDists(randperm(length(randomSkelDists),sum(synPreTypes{1}==7)),:);
            synTipDists=randomSkelDists;
        end
        
        %     figure();
        %     tiledlayout('flow')
        
        bins=[0:.25:10];
        
        
        dists=[];
        for k=1:length(synTipDists(:,1))
            closestTipID=find(synTipDists(k,:)==min(synTipDists(k,:)));
            distToNearestTip=min(synTipDists(k,:));
            closestTipNodeID=tipIDs(closestTipID);
            curSM.syn.pos(k,:);
            curSM.arbor.nodes.pos(closestTipNodeID,:);
            dists=[dists;distToNearestTip];
            
        end
        %size(dists)
        distRefs=[distRefs;dists];
        %figure();
        %nexttile
        
        
    end
    curHistCts=histcounts(distRefs,bins);
    uberRawDat{i}=distRefs;
    if skipPlot~=1
    plotHistCts=curHistCts/sum(curHistCts);
    plotHistCts=smooth(plotHistCts);
    figure(f5);
    if curRef=="bpc"||curRef=="amc"||curRef=="rgc"||curRef=="syn"
        plot(bins(2:end),plotHistCts,'LineWidth',3);
    else
        plot(bins(2:end),plotHistCts,':','LineWidth',1);
    end
    hold on
    
    figure(f6);
    %cp=cdfplot(distRefs);
    %run after
    %cp=cdfplot(distRefs);
    if curRef=="bpc"||curRef=="amc"||curRef=="rgc"||curRef=="syn"
        cp=cdfplot(distRefs);
        cp.LineWidth=3;
    else
        cp=cdfplot(distRefs);
        cp.LineWidth=1;
        cp.LineStyle=':';
    end
    hold on
    end
uberDat(i,:)=plotHistCts;   
%i
end

if skipPlot~=1
figure(f5);
%lgd=legend(distRef);
%lgd.String(3)={'all syns'};
%lgd.String(1:1000)={''};

title('distance from synapses and random skeleton nodes to the nearest branch tip');
xlabel('dist (um)');
ylabel('freq');

figure(f6);
%lgd2=legend(distRef);
%lgd2.String(3)={'all syns'};
end

stdevs=std(uberDat,0,2);
stdevs=std(uberDat,0,1);
means=mean(uberDat,1);
lows=means-(2*stdevs);
highs=means+(2*stdevs);

figure(f5);
hold on
plot(bins(2:end),lows,'LineWidth',1);
plot(bins(2:end),highs,'LineWidth',1);
%lgd.String(5)={'95% low'};
%lgd.String(6)={'95% high'};

figure(f6);
hold on;
cp=cdfplot(lows);
cp.LineWidth=1;
cp=cdfplot(highs);
cp.LineWidth=1;


while 0 
    figure();
    tiledlayout('flow','Padding','none');
    hold on
    for i=1:40
        nexttile
        histogram(uberDat(:,i),100);
        
        
    end
    
end

while 0
    %%
    while 0
        synTipDists=curSM.syn2Skel.syn2SkelDist(:,tipIDs);
        randomSkelNodes=(randperm(length(curSM.arbor.nodes.rad),length(find(synPreTypes{1}==7))));
        randomSkelDists=curSM.skel2skel.linDist(randomSkelNodes,tipIDs);
        
        synPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
        synPostCid=curSM.syn.edges(:,1);
        synPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
        
        synTipDists=synTipDists(synPreTypes{1}==7,:);
        %synTipDists=synTipDists(synPostCid==curCid&synPreTypes{1}'~=7,:);
        synTipDists=randomSkelDists;
        
        figure();
        tiledlayout('flow')
        
        bins=[0:.25:10];
        
        dists=[];
        for k=1:length(synTipDists(:,1))
            closestTipID=find(synTipDists(k,:)==min(synTipDists(k,:)));
            distToNearestTip=min(synTipDists(k,:));
            closestTipNodeID=tipIDs(closestTipID);
            curSM.syn.pos(k,:)
            curSM.arbor.nodes.pos(closestTipNodeID,:)
            dists=[dists;distToNearestTip];
            
        end
        %figure();
        nexttile
        histogram(dists,bins)
    end
    
    %% do it againf ro random
    %start with random nodes and see what it looks like
    
    
    %get random skel nodes at similar zdepths
    
    
    
    synTipDists=curSM.syn2Skel.syn2SkelDist(:,tipIDs);
    
    synPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
    synPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
    
    dists=[];
    for k=1:length(synTipDists(:,1))
        closestTipID=find(synTipDists(k,:)==min(synTipDists(k,:)));
        distToNearestTip=min(synTipDists(k,:));
        closestTipNodeID=tipIDs(closestTipID);
        curSM.syn.pos(k,:)
        curSM.arbor.nodes.pos(closestTipNodeID,:)
        dists=[dists;distToNearestTip];
        
    end
    %figure();
    nexttile
    histogram(dists)
    
    
    
    
    
    % find the nearest tip
    % graph out vs random node selection or vs amc inputs (randomly select
    % a similar number to the bpc total synapse count)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    while 0
        if randomizeTipOrder
            graphTipsSrtd=tipIDsSrtd(randperm(length(tipIDsSrtd)));
        else
            graphTipsSrtd=tipIDsSrtd;
        end
        curCidDat={};
        usedNodeList=[];
        %Actually go through the tips and build out the structure
        for k=1:length(graphTipsSrtd)
            curTipHistDat=[];
            previousEdgeID=0;
            curTipDist=0;
            curTip=graphTipsSrtd(k);
            if curTip~=1
                curNode=curTip;
                if curNode>length(curSM.nep.props.nodeLength);
                    'PROBLEM'
                end
                curEdgeLength=curSM.nep.props.nodeLength(curNode);
                curTipHistDat=[curNode,9,curEdgeLength,curTipDist];
                while curNode~=0 && curTipDist<tipDistCutoff
                    curEdgeLength=curSM.nep.props.nodeLength(curNode);
                    nextNode=curPredInv(curNode);
                    
                    
                    curTipDist=curTipDist+curEdgeLength;
                    closeSyns=find(curSM.syn2Skel.closest==curNode);
                    if ~isempty(closeSyns)
                        for m=1:length(closeSyns)
                            curSyn=closeSyns(m);
                            curSynType=getSMsynType(curSyn,curCid,curSM,curTis);
                            curTipHistDat=[curTipHistDat;[curNode,curSynType,curEdgeLength,curTipDist]];
                        end
                    end
                    %stop if next node has been used already
                    if ~reuseNodes
                        if ismember(nextNode,usedNodeList)
                            break
                        end
                    end
                    %stop if the current node is a fork
                    
                    if ismember(curNode,forkIDs)
                        curTipHistDat=[curTipHistDat;[curNode,6,curEdgeLength,curTipDist]];
                        if stopFork
                            break
                        end
                    else
                        curTipHistDat=[curTipHistDat;[curNode,0,curEdgeLength,curTipDist]];
                        
                    end
                    curNode=nextNode;
                    usedNodeList=[usedNodeList curNode];
                end
                if throwSmall
                    if curTipDist>tipDistCutoff
                        tipDatCid{k}=curTipHistDat;
                    else
                        tipDatCid{k}=[];
                    end
                else
                    tipDatCid{k}=curTipHistDat;
                end
            end
        end
        tipDatDist{j}=tipDatCid;
    end
    tipDat{i}=tipDatDist;
    
    
    %% Setup figures
    if draw3dBranches
        %     branch3d=figure();
        %     branch3dtl=tiledlayout('flow');
        %     title(branch3dtl,'3d arbors');
        %     hold on
        for i=1:length(tipArray)
            branch3d=figure();
            branch3dtl=tiledlayout('flow','TileSpacing','none','Padding','compact');
            hold on
            title(branch3dtl,'3d arbors');
            hold on
            for j=1:length(allSkels)
                curSkel=allSkels{j};
                curSM=curSkel.sm;
                CANP=curSM.arbor.nodes.pos;
                nexttile([1 1]);
                hold on
                view(0,90)
                for k=1:length(tipDat{i}{j})
                    curBranchDat=tipDat{i}{j}{k};
                    if ~isempty(curBranchDat)
                        scatter3(CANP(curBranchDat(:,1),1), ...
                            CANP(curBranchDat(:,1),2), ...
                            CANP(curBranchDat(:,1),3), ...
                            5,'filled');
                    end
                end
                title(num2str(allSkels{j}.sm.cid));
            end
        end
    end
    
    %%
    if drawBranches
        branch2d=figure();
        branch2dtl=tiledlayout('flow');
        title(branch2dtl,'branches');
    end
    
    %%
    if drawDensity
        %make the density plot structure
        branch2d=figure();
        branch2dtl=tiledlayout('flow');
        title(branch2dtl,'branches');
        
    end
    %% TEST BLOCK
    tipTally=[0 0 0 0 0];
    TTD=tipDat{2};
    for curCidIt=1:length(TTD)
        curCidDat=TTD{curCidIt};
        for curTipIt=1:length(curCidDat)
            curTipDat=curCidDat{curTipIt};
            if ~isempty(curTipDat)
                firstForkTipID=find(curTipDat(:,2)==6);
                if curTipDat(firstForkTipID,4)>2
                    tippyTipDat=curTipDat(1:firstForkTipID-1,:);
                    tippyTipOrder=tippyTipDat(tippyTipDat(:,2)>0&tippyTipDat(:,2)~=9,:);
                    tippyTipHist=histcounts(tippyTipOrder(:,2),[0.5:1:5.5]);
                    tipTally=vertcat(tipTally,tippyTipHist);
                    
                    %legend([p1 p2 p3 p4 p5],{'bpc in','amc in','rgc out','amc out','unk'});
                end
            end
        end
    end
    
    %make pie graphs of tips with and without bpcs on them.
    bpcTips=tipTally(tipTally(:,1)>0,:);
    noBpcTips=tipTally(tipTally(:,1)==0,:);
    
    bpcPieDat=sum(bpcTips(:,2:5),1);
    noBpcPieDat=sum(noBpcTips(:,2:5),1);
    
    labs={'amc in','rgc out','amc out','unk'};
    pieFig=figure();
    
    pieTL=tiledlayout(1,2,'TileSpacing','tight','Padding','none');
    
    nexttile
    pie(bpcPieDat);
    sta=title('tips With BPC');
    %posA=get(sta,'position');
    %set(sta,'position',-posA);
    %sta.HorizontalAlignment='right';
    
    nexttile
    pie(noBpcPieDat);
    stb=title('tips Without BPC');
    %posB=get(stb,'position');
    %set(stb,'position',-posB);
    %stb.HorizontalAlignment='right';
    
    lgd = legend(labs);
    lgd.Layout.Tile = 'east';
    
    %suptitle('What other synapses are on tips?');
    bigtitle=title(pieTL,'What other synapses are on tips?');
    bigtitle.VerticalAlignment='bottom';
    
    
end