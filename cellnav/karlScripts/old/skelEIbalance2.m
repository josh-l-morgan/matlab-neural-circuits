%this will graph the skeleton of a cell and show the EI balance in a few
%different ways.

skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
cidList=[2 3 4 5 13 14];

%% get the badSynList
getBadSyns=1;
if getBadSyns
    global tis
    curTis=tis;
    badSyn=find(curTis.syn.pos(:,1)<1);
    badSynPos=curTis.syn.pos(badSyn,:);
    badSynObID=curTis.syn.obID(badSyn);
    badSynName=curTis.obI.nameProps.oldNames(badSynObID);
    allSynPos=curTis.syn.pos;
    distMat=pdist2(allSynPos,allSynPos);
    noDistSyns=find(distMat<0.1);
    self=sub2ind(size(distMat),[1:size(distMat,1)],[1:size(distMat,2)]);
    noDistSynsClean=noDistSyns(~ismember(noDistSyns,self));
    [badx,bady]=ind2sub(size(distMat),noDistSynsClean);
    badInds=horzcat(badx,bady);
    realBadInds=[];
    badIndsCln=badInds(~ismember(badInds(:,2),badSyn),:);
    for i=1:size(badIndsCln,1)
        curCompSynIDs=badIndsCln(i,:);
        badEdges=curTis.syn.edges(curCompSynIDs,:);
        if sum(badEdges(1,1:2)==badEdges(2,1:2))==2
            realBadInds=[realBadInds;curCompSynIDs(1);curCompSynIDs(2)];
            %a=curTis.syn.pos(curComSynIDs(1),[2 1 3])
            %l=curTis.syn.edges(curCompSynIDs(1),:)
            %a=curTis.syn.pos(curCompSynIDs(1),[2 1 3])
            %b=a.*[250 250 25];
            %clipboard('copy',b);
            %pause;
        end
    end
    noPosSyns=badSyn;
    realBadInds=[realBadInds;noPosSyns];
    removeInds=unique(realBadInds);
end


loadSkels=0;
if loadSkels
    allSkels=cell(length(cidList),1);
    for i=1:length(cidList)
        curCid=cidList(i);
        skelFN=['sm_cid' + string(curCid)+'.mat'];
        allSkels{i}=load([skelDir + skelFN]);
    end
end

l=10;
labels={'inh','bal','exc'};
labelTicks=[0,0.5,1];
fig=figure();

cmap=colormap(cool);
hold on
%colorbar(ax);

for l=0.1:0.1:20
for i=1:length(cidList)
    sgtitle(['length const= '+string(l)]);
    %colorbar;
    curSM=allSkels{i}.sm;
    %I want the syn2Skel topo distance
    curDistMat=curSM.syn2Skel.skelTopoDist;
    %get the input synapses
    curCellSynIDs=curSM.syn.synID;
    goodSMsynIDs=~ismember(curSM.syn.synID,removeInds);
    %curCellSynIDs=curCellSynIDs(~ismember(curCellSynIDs,removeInds));
    inputIDs=curSM.syn.edges(:,1)==cidList(i);
    inputIDs=inputIDs';
    upstreamCids=curSM.syn.edges(:,2);
    upstreamTypes=cid2type(upstreamCids,tis);
    bpcInput=upstreamTypes{1}==7;
    bpcInput2=bpcInput&inputIDs;
    amcInput=inputIDs&~bpcInput;
    %divide into bpc and amc
    
    %% make the skelDat for E/I
    
    skelNodeCols=zeros(length(curSM.arbor.nodes.pos(:,1)),3);
    %v(x)=v(max)*e^(-x/l))
    vMat=zeros(size(curDistMat));
    vMat(bpcInput,:)=1;
    vMat(amcInput,:)=-1;
    
    
    %% length constant loop
    %figure();
    curPlot=subplot(2,3,i);
    %hold on
    %axis vis3d
    %colorbar
    %for l=0.1
        infMat=vMat.*exp(-curDistMat./l);
        EIbalMat=sum(infMat,1);
        %figure(); histogram(EIbalMat);
        EIcolNum=(EIbalMat+10)/20*256;
        EIcolNum(EIcolNum<1)=1;
        EIcolNum(EIcolNum>256)=256;
        skelNodeCols=colmap(round(EIcolNum),:);
        %% make a fig
        skNodePos=curSM.arbor.nodes.pos;
        skNodeSize=curSM.arbor.nodes.rad;
        cla
        scatter3(skNodePos(:,1),skNodePos(:,2),skNodePos(:,3),skNodeSize*20,skelNodeCols,'.');
        branchPos=curSM.arbor.nodes.pos(branchNodeList,:);
        hold on
        bpcInputLocs=curSM.syn.pos(bpcInput,:)*10;
        amcInputLocs=curSM.syn.pos(amcInput,:)*10;
        scatter3(bpcInputLocs(:,1),bpcInputLocs(:,2),bpcInputLocs(:,3),40,'mo','filled','MarkerEdgeColor','k');
        scatter3(amcInputLocs(:,1),amcInputLocs(:,2),amcInputLocs(:,3),40,'co','filled','MarkerEdgeColor','k');
        scatter3(branchPos(:,1),branchPos(:,2),branchPos(:,3),30,'ok','filled');
        
        %hold on
        %title(['length const= '+string(l)]);
%         curPlot.XAxisLine = 'off';
%         curPlot.XLabel = '';
%         curPlot.YLabel = '';
%         curPlot.Arrow = 'off';
%         curPlot.Origin = [-Inf -Inf 0];
        xticks([]);
        yticks([]);
        zticks([]);
        box off;
        view([0 0]);
        set(curPlot,'color',[.94 .94 .94]);
        %pause(0.1);
        %figure();
        %trisurf(curSM.show.voxFV.faces,curSM.show.voxFV.vertices(:,1),curSM.show.voxFV.vertices(:,2),curSM.show.voxFV.vertices(:,3))
    %end
    %set(gca,'visible','off')
    title(string(cidList(i)));
    %colorbar(ax);
end
ax=axes(fig,'visible','off');
cbar=colorbar(ax);
cbar.Ticks=labelTicks;
cbar.TickLabels=labels;
cbar.Position=[0.93,0.1,0.025,0.8];
if l==0.1
    pause;
end
filename=['eibal_',num2str(l*10,'%04.f'),'.png'];
saveas(fig,filename)
end
