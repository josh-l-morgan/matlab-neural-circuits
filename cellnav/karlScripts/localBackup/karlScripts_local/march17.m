%before this, run tipsy and justthetip. Them data structures are needed.

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

goodTipFig=0;
if goodTipFig
    tf4=figure();
    tl4=tiledlayout('Flow','Padding', 'none', 'TileSpacing', 'compact');
end
    for i=1:length(allSkels)
        if goodTipFig
        nexttile
        scatter3(1000,1000,300,1,'b.');
        hold on
        view(0,90);
        title(num2str(curCid));
        end
        curSkel=allSkels{i};
        curSM=curSkel.sm;
        curCid=curSM.cid;
        curSWC=curSkel.swc;
        
        curPred=curSWC.pred(curSWC.arbor2swcID)+1;
        curPredInv=zeros(size(curPred));
        curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
        allEdges=curSM.arbor.edges;
        uniqueNodes=unique(allEdges(:));
        nodeCounts=histcounts(allEdges,uniqueNodes);
        tipIDs=find(nodeCounts==1);
        forkIDs=find(nodeCounts>2);
        usedNodeList=[];
        distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
        [srtd srtIdx] = sort(distFromRoot,'descend');
        tipIDsSrtd=tipIDs(srtIdx');
        longTips=TOI{i};
        for k=1:length(tipIDsSrtd)
            drawnow
            if longTips(k)==1
                tipCol=[1 0 1];
            else
                tipCol=[0 1 1];
            end
            previousEdgeID=0;
            curTip=tipIDsSrtd(k);
            if curTip~=1
                %drawnow
                curNode=curTip;
                iterator=0;
                while iterator<500
                    nextNode=curPredInv(curNode);
                    if nextNode==0 | ismember(curNode,usedNodeList)
                        break
                    end
                    %drawnow
                    if ismember(nextNode,forkIDs)
                        tipCol=[0 1 1];
                    end
                    %scatter(k,curTipDist,curNodeRad*50,cool8(round(graphTipDist*12.5)+1,:));
                    if mod(iterator,5)==0
                    if goodTipFig
                        scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                        curSM.arbor.nodes.pos(curNode,2), ...
                        curSM.arbor.nodes.pos(curNode,3), ...
                        10,tipCol);
                    end
                    end
                    if ~isempty(closeSyns)
                        getSyns=0;
                        if getSyns
                            for m=1:length(closeSyns)
                                curSyn=closeSyns(m);
                                curSynType=getSMsynType(curSyn,curCid,curSM,curTis);
                                curTipHistDat=[curTipHistDat;[curTip,curSynType,curTipDist]];
                                quickHistDat=[quickHistDat;[curTip,curSynType,curTipDist]];
                                %scatter(k,curTipDist,5,colCodes(curSynType,:),'filled');
                                if goodTipFig
                                    scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                                    curSM.arbor.nodes.pos(curNode,2), ...
                                    curSM.arbor.nodes.pos(curNode,3), ...
                                    25,colCodes(curSynType,:),'filled');
                                end
                            end
                        end
                    end
                    iterator=iterator+1;
                    usedNodeList=[usedNodeList curNode];
                    curNode=nextNode;
                end
                
            end
        end
        curCid
    end