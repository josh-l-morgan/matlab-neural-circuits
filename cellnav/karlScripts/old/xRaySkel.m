%% xRaySkel
% This will do some diagnostics on the skeletons to make sure that
% everything is as it should be for proceeding to the analyses. Get it?

%% load skeletons if necessary
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
cidList=[2 3 4 5 13 14];
loadSkel=0;
if loadSkel
    allSkels=cell(1,6);
    for i=1:length(cidList)
        
        curCid=cidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end
end

%% Go through skels
for i=1:length(allSkels)
    curSM=allSkels{i}.sm;
    curSyns=curSM.syn;
    targetTypeDat=cid2type(curSyns.post,curTis);
    
    
    
        figure();
        hold on
        scatter3(curSM.arbor.nodes.pos(:,1),curSM.arbor.nodes.pos(:,2), ...
            curSM.arbor.nodes.pos(:,3),curSM.arbor.nodes.rad*5,'.k');
        scatter3(curSyns.pos(:,1)*10,curSyns.pos(:,2)*10,curSyns.pos(:,3)*10, ...
            25,synCol,'o','filled');
        title('Average ratio of topological to euclidean distance of synapse to all other synapses?');
    end
    
end




%% skip skip skip
if 0
    
    
for i=1:length(allSkels)
    curSM=allSkels{i}.sm;
    curSyns=curSM.syn;
    dumbFig=0;
    if dumbFig=1
        curSynTopo=curSM.syn2Skel.syn2SynDist;
        curSynEuc=pdist2(curSyns.pos,curSyns.pos);
        distRatio=curSynTopo./curSynEuc;
        distRatio(distRatio==Inf)=0;
        synDistRatioMean=mean(distRatio,2);
        curMap=turbo(256);
        synCol=curMap(round(synDistRatioMean*50),:);
        badSyns=find(distRatio>10);
        figure();
        hold on
        scatter3(curSM.arbor.nodes.pos(:,1),curSM.arbor.nodes.pos(:,2), ...
            curSM.arbor.nodes.pos(:,3),curSM.arbor.nodes.rad*5,'.k');
        scatter3(curSyns.pos(:,1)*10,curSyns.pos(:,2)*10,curSyns.pos(:,3)*10, ...
            25,synCol,'o','filled');
        title('Average ratio of topological to euclidean distance of synapse to all other synapses?');
    end
    
end

end