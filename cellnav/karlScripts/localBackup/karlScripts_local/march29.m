
% Make sure that you have cellNav open and on the AprilMerge volume so that
% the tis.mat is loaded correctly. Alternatively, the tis.mat could be
% loaded separately.

%number of random iterations
randLength=100;

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

%% PARAMS
skipPlot=1;
randomizeTipOrder=1;
stopFork=0;
reuseNodes=0;
draw3dBranches=0;
drawBranches=0;
drawDensity=0;
throwSmall=0;

%% Depths
nodeColMap=turbo(100);
depthCheck=0;
if depthCheck
    figure();
    hold on
    for i=1:length(allSkels)
        curSkel=allSkels{i};
        curSM=curSkel.sm;
        nodePos=curSM.arbor.nodes.pos/10;
        %doublecheck that these nodes are in the same space as the tis
        [n,m,nodeDepth]=getIPLdepth(nodePos(:,3),nodePos(:,1),nodePos(:,2),[],[]);
        %testSyns=curTis.syn.pos(curTis.syn.edges(:,1)==curSM.cid,:);
        allSkels{i}.sm.arbor.nodes.depth=nodeDepth;
        nodeDepthRnd=(nodeDepth+.05)*.95;
        nodeDepthRnd(nodeDepthRnd<0.01)=0.01;
        nodeDepthRnd(nodeDepthRnd>.99)=.99;
        nodeCols=nodeColMap(round(nodeDepthRnd*100),:);
        %synDepth=getIPLdepth(
        
        scatter3(nodePos(1:10:end,1),nodePos(1:10:end,2),nodePos(1:10:end,3),3,nodeCols(1:10:end,:));
    end
    %scatter3(testSyns(:,1),testSyns(:,2),testSyns(:,3),'m.');
end

%% Data processing loop
uberDat=[];
uberRawDat={};
tipDatDist={};
distRef=["bpc","amcIn","amcOut","rgc"]; %"rand","rand","rand","rand","rand"};
distRef=cellstr(distRef');
randRep="rand";
distRef2=repmat(randRep,[randLength 1]);
distRef2=cellstr(distRef2);
distRef={distRef2{:} distRef{:}};
%distRef=distRef2;

if skipPlot~=1
    f5=figure();
    f6=figure();
end

%these are the numbers of synapses in the smallest class (bpc / rgc)
randNodeNums=[140,75,74,24,76,52];

cumulativeDat=cell(length(distRef),length(allSkels));
for i=1:length(distRef)
    distRefs=[];
    i
    for j=1:length(allSkels)
        tipDatCid={};
        curSkel=allSkels{j};
        curSM=curSkel.sm;
        curCid=curSM.cid;
        allEdges=curSM.arbor.edges;
        [o,p, nodeDepths]=getIPLdepth(curSM.arbor.nodes.pos(:,3)/10, ...
            curSM.arbor.nodes.pos(:,1)/10, ...
            curSM.arbor.nodes.pos(:,2)/10, [], []);
        uniqueNodes=unique(allEdges(:));
        nodeCounts=histcounts(allEdges,uniqueNodes);
        tipIDs=find(nodeCounts==1);
        forkIDs=find(nodeCounts>2);
        closestNodes=curSM.syn2Skel.closest;
        closestNodeDepths=nodeDepths(closestNodes);
        synTipDists=curSM.skel2skel.linDist(closestNodes,tipIDs);
        
        %potentially normalize the random nodes by nodeLength
        
        randNodes=zeros([length(closestNodes),1]);
        for t=1:length(closestNodeDepths)
            curNode=closestNodes(t);
            curNodeDepth=closestNodeDepths(t);
            similarNodes=find(nodeDepths>(curNodeDepth-0.05)&nodeDepths<(curNodeDepth+0.05));
            % Try it also with nodeLength
            randNodes(t)=randsample(similarNodes,1,true,curSM.nep.props.nodeLength(similarNodes));
        end
        %h=histcounts(nodeDepths(randNodes),[0:0.025:0.9]);
        %h=histcounts(nodeDepths(closestNodes),[0:0.025:0.9]);
        %plot(h);
        shortRand=randsample(randNodes,randNodeNums(j));
        randomSkelDists=curSM.skel2skel.linDist(shortRand,tipIDs);
        %randomSkelDists=randsample(randomSkelDists,randNodeNums(j));
        %randDistFromRoot=curSM.skel2skel.linDist(randomSkelNodes,1);
        if 0
            f3=figure();
            swarmchart(repmat(1,length(closestNodes)),closestNodes)
            hold on
            swarmchart(repmat(2,length(randNodes)),randNodes)
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
        %curRef
        if curRef=="bpc"
            synTipDists=synTipDists(synPreTypes{1}==7,:);
        elseif curRef=="amcIn"
            synTipDists=synTipDists(synPostCid==curCid&synPreTypes{1}'~=7,:);
        elseif curRef=="amcOut"
            synTipDists=synTipDists(synPreCid==curCid&synPostTypes{1}'==8|synPostTypes{1}'==0,:);    
        elseif curRef=="syn"
            1;
        elseif curRef=="rgc"
            synTipDists=synTipDists(synPreCid==curCid&synPostTypes{1}'==1,:);
        elseif curRef=="rand"
            %synTipDists=randomSkelDists(randperm(length(randomSkelDists),sum(synPreTypes{1}==7)),:);
            %03/28 Make the random node list the length of the shortest syn
            %type
            synTipDists=randomSkelDists; %(randNodeNums);
        end
        
        %     figure();
        %     tiledlayout('flow')
        
        bins=[0:.25:10];
        
        
        %dists=[];
        %         for k=1:length(synTipDists(:,1))
        %             %closestTipID=find(synTipDists(k,:)==min(synTipDists(k,:)));
        %             distToNearestTip=min(synTipDists(k,:));
        %             %closestTipNodeID=tipIDs(closestTipID);
        %             %curSM.syn.pos(k,:);
        %             %curSM.arbor.nodes.pos(closestTipNodeID,:);
        %             dists=[dists;distToNearestTip];
        %
        %         end
        dists=min(synTipDists,[],2);
        %size(dists)
        distRefs=[distRefs;dists];
        %figure();
        %nexttile
        cumulativeDat{i,j}=dists;
        
        
    end
    %end
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

%% test plots
if 0
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
end
%% work on the cumulative graphs for the illustrator file.
figure();

stepSz=0.05;
SLL=4;
dat=zeros(size(cumulativeDat(end-SLL+1:end,:),1),10/stepSz);

for n=length(distRef)-SLL+1:length(distRef)
    curDists=[];
    for m=1:length(allSkels)
        curDistsSm=cumulativeDat{n,m};
        curDists=[curDists;curDistsSm];
    end
    k=1;
    for i=0:stepSz:10-stepSz
        dat(mod(n,length(distRef)-SLL+1)+1,k)=sum(curDists<i)/length(curDists);
        k=k+1;
    end
    
end

cumulativeDatBackup=cumulativeDat(1:end-SLL,:);
uberRawDatBackup=uberRawDat(1:end-SLL);

cumulativeDat=cumulativeDat(end-SLL+1:end,:);
uberRawDat=uberRawDat(end-SLL+1:end);


%dat is the structure that has all the confidence intervals on it. The long
%axis is the number of random samples, and the shorter axis is the search
%range (in this case 0:0.05:9.95 ; 200 bins)

stepSz=0.05;
datRand=zeros(size(cumulativeDatBackup,1),10/stepSz);

for n=1:length(cumulativeDatBackup)
    curDists=[];
    for m=1:length(allSkels)
        curDistsSm=cumulativeDatBackup{n,m};
        curDists=[curDists;curDistsSm];
    end
    k=1;
    for i=0:stepSz:10-stepSz
        datRand(n,k)=sum(curDists<i)/length(curDists);
        k=k+1;
    end
end

%%
lows=mean(datRand,1)-3*std(datRand,0,1);
highs=mean(datRand,1)+3*std(datRand,0,1);
% %5/4=1.25? 


clf
p1=plot(0:stepSz:10-stepSz,dat(1,:),'LineWidth',3,'Color',[0 0.8 0]);
hold on
p2=plot(0:stepSz:10-stepSz,dat(2,:),'LineWidth',3,'Color',[.8 0 0]);
p3=plot(0:stepSz:10-stepSz,dat(3,:),'LineWidth',3,'Color',[.8 0 .8]);
p4=plot(0:stepSz:10-stepSz,dat(4,:),'LineWidth',3,'Color',[0 0 .8]);
patch1=patch([0:stepSz:10-stepSz fliplr(0:stepSz:10-stepSz)], ...
    [lows fliplr(highs)],[0 0 0], 'FaceAlpha', 0.1);

%legend({'bpc','amcIn','amcOut','rgc','rand'});
ylabel('cumulative %');
xlabel('distance um');
xlim([0 10]);
ylim([0 1]);

if 0
    %fDir = uigetdir;
    fDir='Y:\PUBLICATIONS\VG3\Figure\Draft1\Sources\';
    filename = [fDir 'test']
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])

end

%% Extra code
if 0
    while 0
        
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
    end
end