while 0
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
    end
end
%% get the list of the synapses that I want to graph
bpcLists=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc'},{'bc3a','bc3b','bc4','bc5o','bc5i','bc5t'},curTis);

rgcLists=type2cid({'rgc','rgc','rgc','rgc','rgc'},{'4ow','5ti','37','6sw','63'},curTis);

cidCellArrayA=cell(2,1);
cidCellArrayA{1}=vertcat(bpcLists{1},bpcLists{2},bpcLists{3});
cidCellArrayA{2}=vertcat(bpcLists{4},bpcLists{5},bpcLists{6});
axisNamesA={'OFF bpc','ON bpc'};

cidCellArrayB=cell(3,1);
cidCellArrayB{1}=rgcLists{1};
cidCellArrayB{2}=rgcLists{2};
cidCellArrayB{3}=rgcLists{3};
cidCellArrayB{4}=rgcLists{4};
cidCellArrayB{5}=rgcLists{5};
axisNamesB={'4ow','5ti','37','6sw','63'};

useSkels=[1 2 3 5 6];
for i=1:length(useSkels)
    curSkel=useSkels(i);
    SMcellArray{i}=allSkels{curSkel}.sm;
end

[outputDat stats]=getDRG(cidCellArrayA,cidCellArrayB,SMcellArray,[0:1:30]);

%% go through the outputs and plot the input neighborhoods

plotz=1;
synRad=10;
for i=1:length(SMcellArray)
    curSM=SMcellArray{i};
    curCidFv=curSM.nep.fv;
    for m=1:size(outputDat(1).locsB,2)
        outputIDs=outputDat(i).idsB{m};
        f1=figure();
        for j=1:length(outputIDs)
            hold off
            hold on
            curTipID=outputIDs(j);
            tipPos=outputDat(i).locsB{m}(j,:);
            dists2syns=curSM.syn2Skel.syn2SynDist(:,curTipID);
            %dists2synsA=dists2syns(inputIDs);
            %closeAids=find(dists2synsA<synRad);
            %figure(); histogram(dists2syns,[0:150]);
            %tipPos=curSM.nep.pos(curTipID,:);
            %tipZ=nepDepths(curTipID,:);
            figure(f1);
            scatter3(tipPos(1),tipPos(2),tipPos(3),250,'cp','filled');
            plotRad=10;
            xlim([tipPos(1)-plotRad tipPos(1)+plotRad ])
            ylim([tipPos(2)-plotRad tipPos(2)+plotRad ])
            zlim([tipPos(3)-plotRad tipPos(3)+plotRad ])
            axis vis3d
            fvPatch=patch(curCidFv);
            fvPatch.EdgeColor='none';
            fvPatch.FaceColor=[.2 .2 .2];
            fvPatch.FaceAlpha=.1;
            %find all nearby synapses
            %closeSyns=find(dists2syns<synRad);
            %get positions of close syns
            % closeSynPos=curSM.syn.pos(closeSyns,:);
            %scatter3(inputLocs(closeAids,1),inputLocs(closeAids,2),inputLocs(closeAids,3),100,'mo','filled');
            acols=[[1 0 1];[1 .5 0]];
            for n=1:size(outputDat(1).locsA,2)
                inputIDs=outputDat(i).idsA{n};
                inputLocs=outputDat(i).locsA{n};
                dists2synsA=dists2syns(inputIDs);
                closeAids=find(dists2synsA<synRad);
                scatter3(inputLocs(closeAids,1),inputLocs(closeAids,2),inputLocs(closeAids,3),100,acols(n,:),'filled');
            end
            while 0
                if 0
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
            end
            pause();
        end
    end
    %title(num2str(curCid))
end