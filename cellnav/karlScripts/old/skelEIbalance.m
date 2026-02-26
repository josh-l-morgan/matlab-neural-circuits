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

for i=1:length(cidList)
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
    cmap=colormap(cool);
    skelNodeCols=zeros(length(curSM.arbor.nodes.pos(:,1)),3);
    %v(x)=v(max)*e^(-x/l))
    vMat=zeros(size(curDistMat));
    vMat(bpcInput,:)=1;
    vMat(amcInput,:)=-1;
    infMat=vMat.*exp(-curDistMat./l);
    EIbalMat=sum(infMat,1);
    figure(); histogram(EIbalMat);
    EIcolNum=(EIbalMat+30)/60*256;
    skelNodeCols=colmap(round(EIcolNum),:);
    
    % for each of the nodes, get the distances to all syns and divide the
    % input syns into amc and bpc. Set those as -1 and +1 and then do the
    % modelling with the length constant. Add in that you can do a cylinder
    % with a radius the average of the radii of the two adjacent nodes.
    % Look up how that influences how the signal travels.
    
    
    
    % For a first shot, though, just do the varying length constants with
    % the E/I influence without any radius data and see what it looks like.
    % Make a movie of the thing spinning and see if it looks awesome
    
    
    
    %% make a fig
    figure();
    skNodePos=curSM.arbor.nodes.pos;
    skNodeSize=curSM.arbor.nodes.rad;
    scatter3(skNodePos(:,1),skNodePos(:,2),skNodePos(:,3),skNodeSize*20,skelNodeCols,'.');
    hold on
    axis vis3d
    %figure();
    %trisurf(curSM.show.voxFV.faces,curSM.show.voxFV.vertices(:,1),curSM.show.voxFV.vertices(:,2),curSM.show.voxFV.vertices(:,3))
        
    
end