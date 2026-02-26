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

%% check the orientation of the tip ends.
sprayPatternAnalysis=0;
if sprayPatternAnalysis
    f2=figure();
    tl2=tiledlayout('flow');
    title(tl2,"Orientation of neurite tips");
    for i=1:length(allSkels)
        %pause();
        quickHistDat=[];
        curCidHistDat={};
        %nexttile([1 1]);
        nexttile
        hold on
        colCodes=[[1 0 0];[0 1 0];[0 0 1];[1 0 1];[0 1 1]];
        
        curSkel=allSkels{i};
        curCid=curSkel.sm.cid;
        title(num2str(curCid));
        curSM=curSkel.sm;
        curSWC=curSkel.swc;
        curPred=curSWC.pred(curSWC.arbor2swcID)+1;
        curPredInv=zeros(size(curPred));
        curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
        curSM.pred=curPred;
        allEdges=curSM.arbor.edges;
        uniqueNodes=unique(allEdges(:));
        nodeCounts=histcounts(allEdges,uniqueNodes);
        tipIDs=find(nodeCounts==1);
        forkIDs=find(nodeCounts>2);
        usedNodeList=[];
        sumVect=[0 0 0];
        %need to sort tipIDs
        distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
        [srtd srtIdx] = sort(distFromRoot,'descend');
        tipIDsSrtd=tipIDs(srtIdx');
        for k=1:length(tipIDsSrtd)
            %curTipHistDat={};
            previousEdgeID=0;
            curTip=tipIDsSrtd(k);
            if curTip~=1
                %drawnow
                %xline(k);
                curTipHistDat=[];
                curNode=curTip;
                beginningPos=curSM.arbor.nodes.pos(curNode,:);
                curTipDist=0;
                while curTipDist<20
                    nextNode=curPredInv(curNode);
                    if nextNode==0 %| ismember(curNode,usedNodeList)
                        break
                    end
                    if ismember(curNode,forkIDs)
                        %scatter(k,curTipDist,100,'b_','LineWidth',5);
                        %                     scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                        %                         curSM.arbor.nodes.pos(curNode,2), ...
                        %                         curSM.arbor.nodes.pos(curNode,3),25,'bo','filled');
                        curTipHistDat=[curTipHistDat;[curTip,6,curTipDist]];
                        quickHistDat=[quickHistDat;[curTip,6,curTipDist]];
                        endingPos=curSM.arbor.nodes.pos(curNode,:);
                        break
                    end
                    curNodeRad=curSM.arbor.nodes.rad(curNode);
                    if curTipDist>20
                        graphTipDist=20;
                    else
                        graphTipDist=curTipDist;
                    end
                    %scatter(k,curTipDist,curNodeRad*50,cool8(round(graphTipDist*12.5)+1,:));
                    %                 scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                    %                     curSM.arbor.nodes.pos(curNode,2), ...
                    %                     curSM.arbor.nodes.pos(curNode,3), ...
                    %                     curNodeRad*10,cool8(round(graphTipDist*12.5)+1,:));
                    curEdgeLength=curSM.skel2skel.linDist(nextNode,curNode);
                    curTipDist=curTipDist+curEdgeLength;
                    %drawnow
                    closeSyns=find(curSM.syn2Skel.closest==curNode);
                    if ~isempty(closeSyns)
                        getSyns=0;
                        if getSyns
                            for m=1:length(closeSyns)
                                curSyn=closeSyns(m);
                                curSynType=getSMsynType(curSyn,curCid,curSM,curTis);
                                curTipHistDat=[curTipHistDat;[curTip,curSynType,curTipDist]];
                                quickHistDat=[quickHistDat;[curTip,curSynType,curTipDist]];
                                %scatter(k,curTipDist,5,colCodes(curSynType,:),'filled');
                                scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                                    curSM.arbor.nodes.pos(curNode,2), ...
                                    curSM.arbor.nodes.pos(curNode,3), ...
                                    curNodeRad*12,cool8(round(graphTipDist*12.5)+1,:),'filled');
                                
                            end
                        end
                        
                    end
                    usedNodeList=[usedNodeList curNode];
                    curNode=nextNode;
                    %pause();
                end
                
                if curTipDist<20
                    deltPos=endingPos-beginningPos;
                    plot3([0 deltPos(1)],[0 deltPos(2)],[0 deltPos(3)]);
                    if curTipDist>3
                        normVect=deltPos/sqrt(sum(deltPos.^2));
                        sumVect=vertcat(sumVect,normVect);
                    end
                end
            end
            curCidHistDat{k}=curTipHistDat;
        end
        curCid
        %sumVect
        
        %         title(num2str(curSM.cid));
        %         p1=scatter(0,0,1,colCodes(1,:),'filled');
        %         p2=scatter(0,0,1,colCodes(2,:),'filled');
        %         p3=scatter(0,0,1,colCodes(3,:),'filled');
        %         p4=scatter(0,0,1,colCodes(4,:),'filled');
        %         p5=scatter(0,0,1,colCodes(5,:),'filled');
        %legend([p1 p2 p3 p4 p5],{'bpc in','amc in','rgc out','amc out','unk'});
        %ylim([0,10]);
        winSize=20;
        xlim([-winSize,winSize]);
        ylim([-winSize,winSize]);
        zlim([-winSize,winSize]);
        avgVect=sum(sumVect,1)/size(sumVect,1);
        
        plot3([0 avgVect(1)*50],[0 avgVect(2)*50],[0 avgVect(3)*50],':','Color','b','LineWidth',5);
        view(90,0);
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        allTipHistDat{i}=curCidHistDat;
        allQuickHistDat{i}=quickHistDat;
        
    end
    
end

%% PARAMS
circSize=0;
graphAll=1;
randomTipOrder=0;
tipDistCutoff=5;
noPlot=1;
reProc=1;
plot1=0;
%% classify the order of things in the tips
uberUsedNodes={};
uberDists={};
tipArray=[1:20];
for tipIt=1:length(tipArray)
    tipDistCutoff=tipArray(tipIt);
if reProc
showTipType=1;
if showTipType
    allTipHistDat={};
    if plot1
        f3=figure();
        fl2=tiledlayout('flow');
    end
    %title(tl2,"Order of Operations at tips");
    allDists=zeros(length(allSkels),1);
    for i=1:length(allSkels)
        %pause();
        quickHistDat=[];
        curCidHistDat={};
        curCid=curSkel.sm.cid;
        if plot1
            nexttile([4 1]);
            hold on
            title(num2str(curCid));
        end
        %nexttile
        colCodes=[[1 0 0];[0 1 0];[0 0 1];[1 0 1];[0 1 1]];
        curSkel=allSkels{i};
        curSM=curSkel.sm;
        curSWC=curSkel.swc;
        curPred=curSWC.pred(curSWC.arbor2swcID)+1;
        curPredInv=zeros(size(curPred));
        curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
        curSM.pred=curPred;
        allEdges=curSM.arbor.edges;
        uniqueNodes=unique(allEdges(:));
        nodeCounts=histcounts(allEdges,uniqueNodes);
        tipIDs=find(nodeCounts==1);
        forkIDs=find(nodeCounts>2);
        usedNodeList=[];
        sumVect=[0 0 0];
        %need to sort tipIDs
        distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
        [srtd srtIdx] = sort(distFromRoot,'descend');
        tipIDsSrtd=tipIDs(srtIdx');
        if graphAll
            graphTipsSrtd=tipIDsSrtd;
        else
            goodTips=TOI{i};
            goodTips=find(goodTips==1);
            graphTipsSrtd=tipIDsSrtd(goodTips);
        end
        if randomTipOrder
            graphTipsSrtd=graphTipsSrtd(randperm(length(graphTipsSrtd)));
        end
        totalDistCovered=0;
        for k=1:length(graphTipsSrtd)
            %curTipHistDat={};
            previousEdgeID=0;
            curTip=graphTipsSrtd(k);
            cidsTipsOrd={};
            if curTip~=1
                if ~noPlot
                    drawnow
                end
                %xline(k);
                curNode=curTip;
                curTipHistDat=[curNode,9,0];
                
                beginningPos=curSM.arbor.nodes.pos(curNode,:);
                curTipDist=0;
                abandon=0;
                while ~abandon & curTipDist<tipDistCutoff
                    nextNode=curPredInv(curNode);
                    if nextNode==0 | ismember(nextNode,usedNodeList)
                        abandon=1;
                        break
                    end
                    if ismember(curNode,forkIDs)
                        if ~noPlot
                            scatter(k,curTipDist,25,'b_','LineWidth',1);
                        end
                        %                     scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                        %                         curSM.arbor.nodes.pos(curNode,2), ...
                        %                         curSM.arbor.nodes.pos(curNode,3),25,'bo','filled');
                        %                     quickHistDat=[quickHistDat;[curTip,6,curTipDist]];
                        endingPos=curSM.arbor.nodes.pos(curNode,:);
                        curTipHistDat=[curTipHistDat;[curNode,6,curTipDist]];
                        quickHistDat=[quickHistDat;[curTip,6,curTipDist]];
                        %break
                    end
                    curNodeRad=curSM.arbor.nodes.rad(curNode);
                    if curTipDist>20
                        graphTipDist=20;
                    else
                        graphTipDist=curTipDist;
                    end
                    %scatter(k,curTipDist,curNodeRad*50,cool8(round(graphTipDist*12.5)+1,:));
                    %                 scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                    %                     curSM.arbor.nodes.pos(curNode,2), ...
                    %                     curSM.arbor.nodes.pos(curNode,3), ...
                    %                     curNodeRad*10,cool8(round(graphTipDist*12.5)+1,:));
                    curEdgeLength=curSM.skel2skel.linDist(nextNode,curNode);
                    curTipDist=curTipDist+curEdgeLength;
                    %drawnow
                    closeSyns=find(curSM.syn2Skel.closest==curNode);
                    if ~isempty(closeSyns)
                        getSyns=1;
                        if getSyns
                            for m=1:length(closeSyns)
                                curSyn=closeSyns(m);
                                curSynType=getSMsynType(curSyn,curCid,curSM,curTis);
                                curTipHistDat=[curTipHistDat;[curNode,curSynType,curTipDist]];
                                quickHistDat=[quickHistDat;[curNode,curSynType,curTipDist]];
                                if ~noPlot
                                    scatter(k,curTipDist,15,colCodes(curSynType,:),'filled');
                                end
                                %scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                                %    curSM.arbor.nodes.pos(curNode,2), ...
                                %    curSM.arbor.nodes.pos(curNode,3), ...
                                %    curNodeRad*12,cool8(round(graphTipDist*12.5)+1,:),'filled');
                                
                            end
                        end
                        
                    end
                    usedNodeList=[usedNodeList curNode];
                    curNode=nextNode;
                    %pause();
                end
                
                if curTipDist<tipDistCutoff | abandon
                    curTipHistDat=[];
                else
                    totalDistCovered=totalDistCovered+curTipDist;
                end
                
            end
            curCidHistDat{k}=curTipHistDat;
        end
        curCid
        %sumVect
        
        if plot1
        title(num2str(curSM.cid));
        p1=scatter(0,0,1,colCodes(1,:),'filled');
        p2=scatter(0,0,1,colCodes(2,:),'filled');
        p3=scatter(0,0,1,colCodes(3,:),'filled');
        p4=scatter(0,0,1,colCodes(4,:),'filled');
        p5=scatter(0,0,1,colCodes(5,:),'filled');
        %legend([p1 p2 p3 p4 p5],{'bpc in','amc in','rgc out','amc out','unk'});
        ylim([0,11]);
        end
        winSize=20;
        %xlim([-winSize,winSize]);
        %ylim([-winSize,winSize]);
        %zlim([-winSize,winSize]);
        avgVect=sum(sumVect,1)/size(sumVect,1);
        
        %plot3([0 avgVect(1)*50],[0 avgVect(2)*50],[0 avgVect(3)*50],':','Color','b','LineWidth',5);
        %view(90,0);
        %xlabel('X');
        %ylabel('Y');
        %zlabel('Z');
        allTipHistDat{i}=curCidHistDat;
        allQuickHistDat{i}=quickHistDat;
        allDists(i)=totalDistCovered;
        if plot1
        a=gca;
        a.XGrid='on';
        a.XTick=[0:length(allTipHistDat{i})];
        end
    end
    
end
end


if 0
    for i=1:length(allTipHistDat)
        figure();
        axes();
        hold on
        curSkel=allSkels{i};
        curSM=curSkel.sm;
        curSWC=curSkel.swc;
        curPred=curSWC.pred(curSWC.arbor2swcID)+1;
        curPredInv=zeros(size(curPred));
        curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
        curCidTipDat=allTipHistDat{i};
        for j=1:length(curCidTipDat)
            drawnow
            curCol=rand(1,3);
            curTipDat=curCidTipDat{j};
            if ~isempty(curTipDat)
                curNode=curTipDat(1:1);
                while curNode~=curTipDat(end,1)
                    nextNode=curPredInv(curNode);
                    if nextNode==0
                        break
                    end
                    
                    scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                        curSM.arbor.nodes.pos(curNode,2), ...
                        curSM.arbor.nodes.pos(curNode,3), ...
                        5,curCol);
                    curNode=nextNode;
                end
            end
        end
        title(num2str(curSM.cid))
    end
end


%% Make new histograms

if 1
allQuickHistDat={};
for i=1:length(allTipHistDat)
    curCid=allTipHistDat{i};
    curCidHistDat=[];
    for j=1:length(curCid)
    curTip=curCid{j};
    curCidHistDat=vertcat(curCidHistDat,curTip);
    end
    allQuickHistDat{i}=curCidHistDat;
end


a=allQuickHistDat{1};
b=allQuickHistDat{2};
c=allQuickHistDat{3};
d=allQuickHistDat{4};
e=allQuickHistDat{5};
f=allQuickHistDat{6};

totHistDat=vertcat(a,b,c,d,e,f); %think about skipping cid5
hista=histcounts(totHistDat(totHistDat(:,2)==1,3),[0:1:25]);
histb=histcounts(totHistDat(totHistDat(:,2)==2,3),[0:1:25]);
histc=histcounts(totHistDat(totHistDat(:,2)==3,3),[0:1:25]);
histe=histcounts(totHistDat(totHistDat(:,2)==4,3),[0:1:25]);
histf=histcounts(totHistDat(totHistDat(:,2)==5,3),[0:1:25]);
histg=histcounts(totHistDat(totHistDat(:,2)==6,3),[0:1:25]);


shista=smooth(hista);
shistb=smooth(histb);
shistc=smooth(histc);
shiste=smooth(histe);
shistf=smooth(histf);
shistg=smooth(histg);

figure();
tl1=tiledlayout('flow');
nexttile
hold on
plot(1:length(shista),hista,'LineWidth',4);
plot(1:length(shista),histb,'LineWidth',4);
plot(1:length(shista),histc,'LineWidth',4);
plot(1:length(shista),histe,'LineWidth',4);
plot(1:length(shista),histf,'LineWidth',4);
plot(1:length(shista),histg,'LineWidth',2);
%legend({'bpc in','amc in','rgc out','amc out','unk','fork'});
if tipDistCutoff<10
xlim([0 11])
else
    xlim([0 tipDistCutoff+1]);
end
title(['raw histogram for tip length ' num2str(tipDistCutoff) 'um']);
%legend({'bpc in','amc in','rgc out','amc out','unk','fork'});

nexttile
%figure();
hold on
plot(1:length(shista),shista,'LineWidth',4);
plot(1:length(shista),shistb,'LineWidth',4);
plot(1:length(shista),shistc,'LineWidth',4);
plot(1:length(shista),shiste,'LineWidth',4);
plot(1:length(shista),shistf,'LineWidth',4);
plot(1:length(shista),shistg,'LineWidth',2);
title(['smooth histogram for tip length ' num2str(tipDistCutoff) 'um']);
if tipDistCutoff<10
xlim([0 11])
else
    xlim([0 tipDistCutoff+1]);
end
%legend({'bpc in','amc in','rgc out','amc out','unk','fork'});


end
uberDists{tipDistCutoff}=allDists;


end


%% show the amounts of arbor covered by the different tip lengths
if 1
    figure();
    title('how much branch is included?');
    hold on
    graphDat=zeros(length(uberDists),length(uberDists{10}));
    for i=1:length(uberDists)
        if ~isempty(uberDists{i})
            graphDat(i,:)=uberDists{i}';
            
        end
    end
    
    for i=1:length(uberDists{1})
        plot(1:length(graphDat(:,i)),graphDat(:,i));
        %hold on
    end
    
    
    
end


%% show the orders

if circSize
    tipCutOff=5;
    k=0;
    TOI={};
    for i=1:length(allTipHistDat)
        curCidHistDat=allTipHistDat{i};
        curCidTOI=zeros(length(curCidHistDat),1);
        for j=1:length(curCidHistDat)
            curTipDat=curCidHistDat{j};
            if sum(curTipDat(find(curTipDat(:,2)==6),3)>tipCutOff)>0
                curCidTOI(j)=1;
            end
        end
        TOI{i}=curCidTOI;
        
    end
end