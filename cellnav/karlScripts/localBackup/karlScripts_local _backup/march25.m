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
cumulativeDat={};
uberDat=[];
uberRawDat={};
tipDatDist={};
distRef={"bpc","amc","rgc"}; %"rand","rand","rand","rand","rand"};
randRep="rand";
distRef2=repmat(randRep,[10000 1]);
distRef2=cellstr(distRef2);
%distRef={distRef2{:} distRef{:}};
distRef=distRef2;

if skipPlot~=1
    f5=figure();
    f6=figure();
end


for i=1:length(distRef)
    distRefs=[];
    i
    for j=1:length(allSkels)
        tipDatCid={};
        %get the skeleton data set up for going from the tips up the
        %branches to the root.
        curSkel=allSkels{j};
        curSM=curSkel.sm;
        %curSWC=curSkel.swc;
        curCid=curSM.cid;
        %curPred=curSWC.pred(curSWC.arbor2swcID)+1;
        %curPredInv=zeros(size(curPred));
        %curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
        %curSM.pred=curPred;
        allEdges=curSM.arbor.edges;
        uniqueNodes=unique(allEdges(:));
        nodeCounts=histcounts(allEdges,uniqueNodes);
        tipIDs=find(nodeCounts==1);
        forkIDs=find(nodeCounts>2);
        closestNodes=curSM.syn2Skel.closest;
        %synDistFromRoot=curSM.skel2skel.linDist(closestNodes,1);
        %distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
        %[srtd srtIdx] = sort(distFromRoot,'descend');
        %tipIDsSrtd=tipIDs(srtIdx');
        % get all the synapse types and their closest nodes
        
        % get the distances from those nodes to the tips
        %synTipDists=curSM.syn2Skel.syn2SkelDist(:,tipIDs);
        synTipDists=curSM.skel2skel.linDist(closestNodes,tipIDs);
        
        randomSkelNodes=randperm(length(curSM.skel2skel.linDist(1,:)),length(curSM.syn2Skel.closest));
        %randomSkelNodes=randperm(length(curSM.skel2skel.linDist(1,:)),500);
        %randomSkelNodes=randomSkelNodes(find(~ismember(randomSkelNodes,closestNodes)));
        randomSkelNodesB=randomSkelNodes(~ismember(randomSkelNodes,closestNodes));
        randomSkelNodes=randomSkelNodes(~ismember(randomSkelNodes,tipIDs));
        randomSkelDists=curSM.skel2skel.linDist(randomSkelNodes,tipIDs);
        %randDistFromRoot=curSM.skel2skel.linDist(randomSkelNodes,1);
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

if skipPlot~=1
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
end

%% work on the cumulative graphs for the illustrator file.

%sort the distances in the random block
superDist=[];
for i=1:1000
    curDists=uberRawDat{i};
    superDist=[superDist;curDists];
end

%SDsrtd=sort(superDist);
plotPts=zeros(100,2);
len=length(SDsrtd);
stepSz=0.25; %calc step in percent
plotPts=zeros(100/stepSz,2);
for i=1:floor(100/stepSz)
    plotPts(i,:)=([i*stepSz/100,SDsrtd(round(i*stepSz/100*len))]);

    
end
f0=figure();
a1=axes();
p1=plot(plotPts(:,2),plotPts(:,1));
p1.LineWidth=2;
hold on
dsX=5;
for synIt=1:3
    curDists=uberRawDat{synIt};
    plotPts=zeros(length(curDists),2);
    CDsrtd=sort(curDists);
    for j=1:length(CDsrtd)
        plotPts(j,:)=([CDsrtd(j),j/length(CDsrtd)]);
    end
    plotPts=plotPts(1:dsX:end,:);
    p2=plot(plotPts(:,1),plotPts(:,2));
    p2.LineWidth=2;
end
legend({'rand','bpc','amc','rgc'});

%%
while 0
    figure();
    tiledlayout('flow','Padding','none');
    hold on
    for i=1:40
        nexttile
        histogram(uberDat(:,i),100);
        
        
    end
    
end