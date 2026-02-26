%load skeletons if necessary
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
%% main
coolColMap=cool(101);
cool8=cool(256);
if 1
    allTipHistDat={};
    allQuickHistDat={};
    %tiledlayout(2,3);
    
    figure();
    tl=tiledlayout('flow');
    %xlabel('tip# (starting most distal to root)')
    %ylabel('distance from tip')
    
    hold on
    %title(num2str(curCid));
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    
    for i=5%:length(allSkels)
        %pause();
        quickHistDat=[];
        curCidHistDat={};
        %nexttile([1 1]);
        nexttile
        hold on
        colCodes=[[1 0 0];[0 1 0];[0 0 1];[1 0 1];[0 1 1]];
        
        curSkel=allSkels{i};
        curCid=curSkel.sm.cid;
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
        %need to sort tipIDs
        distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
        [srtd srtIdx] = sort(distFromRoot,'descend');
        tipIDsSrtd=tipIDs(srtIdx');
        for k=1:length(tipIDsSrtd)
            %curTipHistDat={};
            previousEdgeID=0;
            curTip=tipIDsSrtd(k);
            if curTip~=1
                drawnow
                %xline(k);
                curTipHistDat=[];
                curNode=curTip;
                curTipDist=0;
                while curTipDist<200
                    nextNode=curPredInv(curNode);
                    if nextNode==0 | ismember(nextNode,usedNodeList)
                        break
                    end
                    if ismember(curNode,forkIDs)
                        %scatter(k,curTipDist,100,'b_','LineWidth',5);
                        scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                            curSM.arbor.nodes.pos(curNode,2), ...
                            curSM.arbor.nodes.pos(curNode,3),25,'bo','filled');
                        quickHistDat=[quickHistDat;[curTip,6,curTipDist]];
                    end
                    curNodeRad=curSM.arbor.nodes.rad(curNode);
                    if curTipDist>20
                        graphTipDist=20;
                    else
                        graphTipDist=curTipDist;
                    end
                    %scatter(k,curTipDist,curNodeRad*50,cool8(round(graphTipDist*12.5)+1,:));
                    scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                        curSM.arbor.nodes.pos(curNode,2), ...
                        curSM.arbor.nodes.pos(curNode,3), ...
                        curNodeRad*10,cool8(round(graphTipDist*12.5)+1,:));
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
                
            end
            curCidHistDat{k}=curTipHistDat;
        end
        
        %         title(num2str(curSM.cid));
        %         p1=scatter(0,0,1,colCodes(1,:),'filled');
        %         p2=scatter(0,0,1,colCodes(2,:),'filled');
        %         p3=scatter(0,0,1,colCodes(3,:),'filled');
        %         p4=scatter(0,0,1,colCodes(4,:),'filled');
        %         p5=scatter(0,0,1,colCodes(5,:),'filled');
        %legend([p1 p2 p3 p4 p5],{'bpc in','amc in','rgc out','amc out','unk'});
        %ylim([0,10]);
        
        allTipHistDat{i}=curCidHistDat;
        allQuickHistDat{i}=quickHistDat;
        
    end
end

a=allQuickHistDat{1};
b=allQuickHistDat{2};
c=allQuickHistDat{3};
e=allQuickHistDat{5};
f=allQuickHistDat{6};

totHistDat=vertcat(a,b,c,e,f);
hista=histcounts(totHistDat(totHistDat(:,2)==1,3),[0:0.25:10]);
histb=histcounts(totHistDat(totHistDat(:,2)==2,3),[0:0.25:10]);
histc=histcounts(totHistDat(totHistDat(:,2)==3,3),[0:0.25:10]);
histe=histcounts(totHistDat(totHistDat(:,2)==4,3),[0:0.25:10]);
histf=histcounts(totHistDat(totHistDat(:,2)==5,3),[0:0.25:10]);
histg=histcounts(totHistDat(totHistDat(:,2)==6,3),[0:0.25:10]);


shista=smooth(hista);
shistb=smooth(histb);
shistc=smooth(histc);
shiste=smooth(histe);
shistf=smooth(histf);
shistg=smooth(histg);

figure();
hold on
plot(1:40,hista,'LineWidth',4);
plot(1:40,histb,'LineWidth',4);
plot(1:40,histc,'LineWidth',4);
plot(1:40,histe,'LineWidth',4);
plot(1:40,histf,'LineWidth',4);
plot(1:40,histg,'LineWidth',2);
legend({'bpc in','amc in','rgc out','amc out','unk','fork'});

figure();
hold on
plot(1:40,shista,'LineWidth',4);
plot(1:40,shistb,'LineWidth',4);
plot(1:40,shistc,'LineWidth',4);
plot(1:40,shiste,'LineWidth',4);
plot(1:40,shistf,'LineWidth',4);
plot(1:40,shistg,'LineWidth',2);
legend({'bpc in','amc in','rgc out','amc out','unk','fork'});


for q=0
    if 0
        tf1=figure();
        tiledlayout(2,3);
        for i=1:length(allSkels)
            %pause();
            curSkel=allSkels{i};
            curCid=curSkel.sm.cid;
            curSM=curSkel.sm;
            colMat=zeros(length(curSM.nep.nodeTipDist),3);
            distMat=curSM.nep.nodeTipDist;
            %distMat(distMat>10)=10;
            colMat=coolColMap(102-(round(distMat*10)+1),:);
            curSWC=curSkel.swc;
            curPred=curSWC.pred(curSWC.arbor2swcID)+1;
            curPredInv=zeros(size(curPred));
            curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
            curSM.pred=curPred;
            growFig=1;
            if growFig
                %figure();
                %axes();
                nexttile
                hold on
                curNodeIDList=[0];
                
                while ~isempty(curNodeIDList)
                    newNodeID=find(ismember(curPredInv,curNodeIDList));
                    scatter3(curSM.nep.pos(newNodeID,1),curSM.nep.pos(newNodeID,2), ...\
                        curSM.nep.pos(newNodeID,3),5,colMat(newNodeID,:));
                    %newNodeID=find(ismember(curPred,curNodeIDList));
                    curNodeIDList=find(ismember(curPredInv,newNodeID));
                    pause(0.01);
                    
                end
            end
            title(num2str(curSM.cid));
            
        end
    end
    
    %%
    if 0
        tf2=figure();
        tiledlayout(2,3);
        for i=1:length(allSkels)
            if 1
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
                
                %firstNodeID=find(curSkel.sm.arbor.nodes.pos(:,3)==min(curSkel.sm.arbor.nodes.pos(:,3)));
                %rootNode=firstNodeID(1);
                
                tipStruct=struct();
                nodeDistTip=zeros(size(curSM.nep.nodes));
                nodeDistTip(:)=10;
                for k=1:length(tipIDs)
                    tipStruct(k).tipID=tipIDs(k);
                    tipStruct(k).nodeList=[];
                    tipStruct(k).lengthList=[];
                    tipStruct(k).totalLength=0;
                    previousEdgeID=0;
                    curTip=tipIDs(k);
                    if curTip~=1
                        curNode=curTip;
                        curTipDist=0;
                        while ~ismember(curNode,forkIDs)
                            curEdgeID=find(any(allEdges == curNode, 2));
                            curEdgeID=curEdgeID(curEdgeID~=previousEdgeID);
                            curEdgeNodes=allEdges(curEdgeID,:);
                            if curEdgeID>length(curSM.nep.edgeLength)
                                break
                            end
                            curEdgeLength=curSM.nep.edgeLength(curEdgeID);
                            nextNode=curEdgeNodes(curEdgeNodes~=curNode);
                            tipStruct(k).nodeList=[tipStruct(k).nodeList; curNode];
                            tipStruct(k).lengthList=[tipStruct(k).lengthList; curEdgeLength];
                            nodeDistTip(curNode)=curTipDist;
                            curTipDist=curTipDist+curEdgeLength;
                            curNode=nextNode;
                            previousEdgeID=curEdgeID;
                        end
                        tipStruct(k).root=curNode;
                        tipStruct(k).totalLength=curTipDist;
                    end
                end
            end
            allSkels{i}.tipDat=tipStruct;
            allSkels{i}.sm.nep.nodeTipDist=nodeDistTip;
            curSM.cid
            sum(nodeDistTip>10)
            %figure();
            nexttile
            histogram(nodeDistTip,[0:0.25:10]);
            hold on;
            title(num2str(curSM.cid));
            
            
        end
    end
    
    %% do VG3 branches terminate with synapses?
    % If so, are they output / input / both?
    if 0
        tf2=figure();
        tiledlayout(2,3);
        tipHood=struct();
        for i=1:length(allSkels)
            curSkel=allSkels{i};
            curCid=curSkel.sm.cid;
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
            allTipSyns={};
            synTypeNames={'bpc in','amc in','rgc out','amc out','unk'};
            if graphAll
                figure();
                title(num2str(curCid));
                hold on
                colCodes=[[1 0 0];[0 1 0];[0 0 1];[1 0 1];[0 1 1]];
                p1=scatter(0,0,1,colCodes(1,:),'filled');
                p2=scatter(0,0,1,colCodes(2,:),'filled');
                p3=scatter(0,0,1,colCodes(3,:),'filled');
                p4=scatter(0,0,1,colCodes(4,:),'filled');
                p5=scatter(0,0,1,colCodes(5,:),'filled');
                %L(m).Annotation=synTypeNames{m};
                %end
                %legend({'bpc in','amc in','rgc out','amc out','unk'});
            end
            for j=1:length(tipIDs)
                curTipSyns=[];
                curTipID=tipIDs(j);
                curTipSynDist=curSM.syn2Skel.syn2SkelDist(:,curTipID);
                previousEdgeID=0;
                if curTip~=1
                    curTipDist=0;
                    curNode=curTipID;
                    while ~ismember(curNode,forkIDs)
                        curEdgeID=find(any(allEdges == curNode, 2));
                        curEdgeID=curEdgeID(curEdgeID~=previousEdgeID);
                        curEdgeNodes=allEdges(curEdgeID,:);
                        if curEdgeID>length(curSM.nep.edgeLength)
                            break
                        end
                        curEdgeLength=curSM.nep.edgeLength(curEdgeID);
                        nextNode=curEdgeNodes(curEdgeNodes~=curNode);
                        %find the synapses that are closeby
                        curNodeSyns=find(curSM.syn2Skel.closest==curNode);
                        if ~isempty(curNodeSyns)
                            curTipSyns=vertcat(curTipSyns,[repmat(curNode,length(curNodeSyns),1),curNodeSyns,repmat(curTipDist,length(curNodeSyns),1)]);
                            if graphAll
                                for k=1:length(curNodeSyns)
                                    curSyn=curNodeSyns(k);
                                    curEdges=curSM.syn.edges(curSyn,:);
                                    edgeTypes=cid2type(curEdges(1:2),curTis);
                                    synType=5;
                                    if curEdges(1)==curCid %input
                                        if edgeTypes{1}(2)==7
                                            synType=1; %BPC IN
                                        elseif curEdges(2)==0 | edgeTypes{1}(2)==8
                                            synType=2; %AMC IN
                                        end
                                    elseif curEdges(2)==curCid %output
                                        if edgeTypes{1}(1)==1
                                            synType=3; %RGC OUT
                                        elseif curEdges(1)==0 | edgeTypes{1}(1)==8
                                            synType=4; %AMC OUT
                                        end
                                    end
                                    scatter(j,curTipDist,15,colCodes(synType,:),'filled');
                                end
                            end
                        end
                        nodeDistTip(curNode)=curTipDist;
                        curTipDist=curTipDist+curEdgeLength;
                        curNode=nextNode;
                        previousEdgeID=curEdgeID;
                    end
                end
                allTipSyns{j}=curTipSyns;
                %legend({'bpc in','amc in','rgc out','amc out','unk'});
                legend([p1 p2 p3 p4 p5],{'bpc in','amc in','rgc out','amc out','unk'});
                %ylim([0,8]);
            end
            
        end
    end
    
end
% for each tip, get the closest synapses. Plot them out colorcoded with a
% log distance from the tip. See if anything jumps out.

function typeInfo=getSMsynType(synID,curCid,curSM,tisDat)
synType=5;
curEdges=curSM.syn.edges(synID,:);
edgeTypes=cid2type(curEdges(1:2),tisDat);
if curEdges(1)==curCid %input
    if edgeTypes{1}(2)==7
        synType=1; %BPC IN
    elseif curEdges(2)==0 | edgeTypes{1}(2)==8
        synType=2; %AMC IN
    end
elseif curEdges(2)==curCid %output
    if edgeTypes{1}(1)==1
        synType=3; %RGC OUT
    elseif curEdges(1)==0 | edgeTypes{1}(1)==8
        synType=4; %AMC OUT
    end
end
typeInfo=synType;
end

