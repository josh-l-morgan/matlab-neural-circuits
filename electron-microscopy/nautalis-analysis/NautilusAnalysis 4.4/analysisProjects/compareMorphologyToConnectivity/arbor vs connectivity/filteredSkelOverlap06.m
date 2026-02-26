%%Make predictions of how many synapses each pair of cells should form
%%based of the distribution of their skeletons. -
%-Skeletons are filtered



%% Plot axon to tcr
clear all
MPN = GetMyDir;
load([MPN 'obI.mat'])
shouldPause = 10;
anchorScale = [-.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];

skelOverlapPred.predictionType = 'skelOverlap';

%% variables

firstPredThresh = 0;
binWidth = 1;


seedList = [108 201]
crossoverAxons = [2032	2033	2034	2035]
noSkel = [2014 1026]

muchTraced = [106 107 108 109 111 112 117 120 123 129 133 134 148 156 159 162 163 169 ...
    170 201 203 205 206 207 210 212 213 215 216 218];
skelOverlapPred.muchTraced = muchTraced;

%% Get cells

useList = obI2cellList_seedInput(obI,seedList);
axList = useList.preList;
cellList = useList.postList;
synMat = useList.con;

skelOverlapPred.useList = useList;

synapses = obI.nameProps.edges;
edges = synapses(:,1:2);


axList = setdiff(axList,noSkel);
disp(sprintf('Results calculated without %d',noSkel));

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end

skelOverlapPred.con = con;

%%  Get Skeletons

for i = 1 : length(axList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,axList(i));
    load(fileName)
    axSkel(i).skel = cellStruct.skel;
    axSkel(i).name = axList(i);
end


for i = 1 : length(cellList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,cellList(i));
    load(fileName)
    cellSkel(i).skel = cellStruct.skel;
    cellSkel(i).name = cellList(i);
end

%% Create filtered node list
minLength = 3; %minimum branch length to be counted 
clear axNodes cellNodes

for i = 1:length(axSkel)
   
    skel = axSkel(i).skel;
    nodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keepNodes = [];
    for b = 1:length(skel.bones);
        bone = skel.bones(b);
        L = 0;
        for e = 1:size(bone.edges,1);
            edgeNodes = nodeSubs(bone.edges(e,:),:);
            edgeDif = diff(edgeNodes,1);
            edgeLength = sqrt(sum(edgeDif.^2));
            L = L + edgeLength;
        end
        disp(L)
        boneLengths(b) = L;
        if L>minLength
            keepNodes = cat(2,keepNodes,bone.nodes);
        end
    end
    keepNodes = unique(keepNodes);
    
    
    allNodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keptNodeSubs = allNodeSubs(keepNodes,:);
    scatter(allNodeSubs(:,3),allNodeSubs(:,2),'.','k')
    hold on
    scatter(keptNodeSubs(:,3),keptNodeSubs(:,2),'.','r')
    hold off
    pause(.01)
    
    axNodes{i} = keptNodeSubs;
end

for i = 1:length(cellSkel)
    
    
    
     skel = cellSkel(i).skel;
    nodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keepNodes = [];
    for b = 1:length(skel.bones);
        bone = skel.bones(b);
        L = 0;
        for e = 1:size(bone.edges,1);
            edgeNodes = nodeSubs(bone.edges(e,:),:);
            edgeDif = diff(edgeNodes,1);
            edgeLength = sqrt(sum(edgeDif.^2));
            L = L + edgeLength;
        end
        disp(L)
        boneLengths(b) = L;
        if L>minLength
            keepNodes = cat(2,keepNodes,bone.nodes);
        end
    end
    keepNodes = unique(keepNodes);
    
    
    allNodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keptNodeSubs = allNodeSubs(keepNodes,:);
    scatter(allNodeSubs(:,3),allNodeSubs(:,2),'.','k')
    hold on
    scatter(keptNodeSubs(:,3),keptNodeSubs(:,2),'.','r')
    hold off
    pause(.01)
    
    cellNodes{i} = keptNodeSubs;
    
    
end




%% Map existing synapses
rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors = scaleSubs(double(rawSynAnchors),anchorScale);

goodSyn = zeros(size(synAnchors,1),1)>0;

for i = 1:size(con,1);
    for p = 1:size(con,2)
        conDists{i,p} = [];
    end
end
goodCon = con * 0;

for i = 1:size(synAnchors,1)
    if sum(cellList == edges(i,1)) & sum(axList == edges(i,2))
%         preSkel = axSkel(find(axList == edges(i,2))).skel;
%         postSkel = cellSkel(find(cellList == edges(i,1))).skel;
%         
%         preNodes = scaleSubs(double(preSkel.node2subs),voxelScale);
%         postNodes = scaleSubs(double(postSkel.node2subs),voxelScale);
        
        %     preNodes = preSkel.node2subs;
        %     preNodes  = [preNodes(:,1) * voxelScale(1) preNodes(:,2) * voxelScale(2) preNodes(:,3) * voxelScale(3)];
        %     postNodes = postSkel.node2subs;
        %     postNodes  = [postNodes(:,1) * voxelScale(1) postNodes(:,2) * voxelScale(2) postNodes(:,3) * voxelScale(3)];
        %
        preTarg = find(axList == edges(i,2));
        preNodes = axNodes{preTarg};
        postTarg = find(cellList == edges(i,1));
        postNodes = cellNodes{postTarg};
        postDists = sqrt((postNodes(:,1)-synAnchors(i,1)).^2 + (postNodes(:,2)-synAnchors(i,2)).^2 + ...
            (postNodes(:,3)-synAnchors(i,3)).^2);
        
        preDists = sqrt((preNodes(:,1)-synAnchors(i,1)).^2 + (preNodes(:,2)-synAnchors(i,2)).^2 + ...
            (preNodes(:,3)-synAnchors(i,3)).^2);
        
        if sum([min(preDists) min(postDists)])<10
            goodSyn(i) = 1;
        
        synNodeId = [find(postDists == min(postDists),1) find(preDists == min(preDists),1)];
        synNodeSub = [postNodes(synNodeId(1),:); preNodes(synNodeId(2),:)];
        synDif = diff(synNodeSub,1);
        synDist = sqrt(sum(synDif.^2));
        
        conDists{preTarg,postTarg} = cat(1,conDists{preTarg,postTarg},synDist);
        goodCon(preTarg,postTarg) = goodCon(preTarg,postTarg) + 1;
        end
        
    end
end

skelOverlapPred.conDists = conDists;
skelOverlapPred.goodCon = goodCon;

%% Calculate overlap
%%. 1) node cloud, 2)continous, 3) threshold vs gaussian



%%Get distribution of distances for relevant cells
allDists  = zeros(10000,1);
lastDist = 0;
histBin = [0:binWidth:100];
for i = 1:length(axList)
    disp(sprintf('checking axon %d of %d',i,length(axList)))
    preSkel = axSkel(i).skel;
    preNodes = scaleSubs(preSkel.node2subs,voxelScale);
    %preNodes  = [preNodes(:,1) * voxelScale(1) preNodes(:,2) * voxelScale(2) preNodes(:,3) * voxelScale(3)];
    tic
    parfor p = 1: length(cellList)
        postSkel = cellSkel(p).skel;
        postNodes = scaleSubs(postSkel.node2subs,voxelScale);
        %postNodes  = [postNodes(:,1) * voxelScale(1) postNodes(:,2) * voxelScale(2) postNodes(:,3) * voxelScale(3)];
        
        pairDists = zeros(size(preNodes,1),size(postNodes,1));
        for n = 1:size(preNodes,1)
            pairDists(n,:) = sqrt((postNodes(:,1)-preNodes(n,1)).^2 + (postNodes(:,2)-preNodes(n,2)).^2 + ...
                (postNodes(:,3)-preNodes(n,3)).^2);
        end
        histDists = histc(pairDists(:),histBin);
        axCellHists{i,p} = histDists;
        
        %
        %         nextDist = lastDist + size(pairDists,1) * size(pairDists,2);
        %         if nextDist > length(allDists)
        %             allDists(length(allDists) * 2) = 0;
        %         end
        %         %nodeDists{i,p} = pairDists;
        %         allDists(lastDist+1:nextDist) = pairDists(:);
        %         lastDist = nextDist;
    end
    toc
end

%% Get one synapse prediction filter

sumHists = zeros(1,length(histBin));
synDists = [];
for i = 1:length(axList)
    for p = 1:length(cellList)
        if sum(find(muchTraced== cellList(p)));  % only use traced post synaptic cells
            synDists = cat(1,synDists,  conDists{i,p});
            sumHists = sumHists + axCellHists{i,p}';
        end
    end
end

synHistDist = hist(synDists,histBin);
nodeFrac = synHistDist./sumHists;


firstPredSyn = con * 0;
for i = 1:length(axList)
    for p = 1:length(cellList)
        pairCount= axCellHists{i,p} ;
        firstPredSyn(i,p) = sum(pairCount'.*nodeFrac);
    end
end


%%Repeat using filter
useOverlap = firstPredSyn>=firstPredThresh;
sumHists = zeros(1,length(histBin));
synDists = [];
for i = 1:length(axList)
    for p = 1:length(cellList)
        if sum(find(muchTraced== cellList(p))) & useOverlap(i,p);  % only use traced post synaptic cells
            synDists = cat(1,synDists,  conDists{i,p});
            sumHists = sumHists + axCellHists{i,p}';
        end
    end
end

synHistDist = hist(synDists,histBin);
nodeFrac = synHistDist./sumHists;


predSyn = con * 0;
for i = 1:length(axList)
    for p = 1:length(cellList)
        pairCount= axCellHists{i,p} ;
        predSyn(i,p) = sum(pairCount'.*nodeFrac);
    end
end


skelOverlapPred.predSyn = predSyn;
skelOverlapPred.firstPredSyn = firstPredSyn;
skelOverlapPred.histBin = histBin;
skelOverlapPred.axCellHists = axCellHists;
skelOverlapPred.nodeFrac = nodeFrac;
skelOverlapPred.usedOverlap = useOverlap;
skelOverlapPred.firstPredThresh = firstPredThresh;
%}


%% Save data
binWidth
fileName = sprintf('%sskelOverlapPred_%d.mat',MPN,binWidth);
save(fileName,'skelOverlapPred')



%% display
clf
colList = ['rgbcmyk']
r = 1
    %   scatCol = [r/length(allPred) 1-r/length(allPred) .3]
    %   scatter(predSyn(:),con(:),15,scatCol,'filled','o');
    difPred{r} = con - predSyn;
    histDif{r} = hist(difPred{r}(:),[-10:10])
    
    subplot(2,1,1)
    bar(histDif{r},colList(r))
    image(con*20)
    %   text([10 10000],num2str(binRes(r)),'EdgeColor',colList(r))
    subplot(2,1,2)
    conPref = (con-predSyn)./(predSyn+1)
    conPref(isnan(conPref)) = 0;
    image(conPref*100+100)
    
    
    
    checkPred = predSyn(predSyn ~=0);
    
    
hold off




%% Run predictPreferencesFrom....mat

%{

%% Predict from similarity
% 1) could predict based on the number of alternate paths.  2) could
% predict based on similarity  3) run it from the seed cell perspective
clf
seedCells = [108 201];
predSyn = allPred{r};
conPref = (con-predSyn)./(predSyn + 1);
preSeedPrefMat = con*0;
postSeedPrefMat = con*0;
preSeesPostPref = zeros(size(con,1), size(con,2), 2); %connection to seed count of each post with each pre removed
preSeesCountPref =  zeros(size(con,1), size(con,2), 2);
for i = 1: length(axList)
    
    usePre = setdiff(1:length(axList),i);
    sampCon = con(usePre,:);
    %sampPref =
    for s = 1:length(seedCells)
        targ = find(cellList==seedCells(s));
        simMat = sqrt(sampCon .* repmat(sampCon(:,targ),[1 size(con,2)]));
        image(simMat*20);
        sumMat(s,:) = sum(simMat,1);
        preSeed(s) = con(i,targ);
    end
    %preSeed = preSeed/sum(preSeed);
    postSeedPref = sumMat(1,:)./sum(sumMat,1) % tcr pref for seed minus the RGC of interest
    if sum(isnan(postSeedPref))
        postSeedPref(isnan(postSeedPref)) = .5;
        'must deal with singly innervated cells'
    end
    preSeesPostPref(i,:,:) = sumMat';
    preSeesCountPref(i,:,1) = preSeed(1);
    preSeesCountPref(i,:,2) = preSeed(2);


    preSeedPref = preSeed(1)/sum(preSeed);
    bar(sumMat');
    preSeedPrefMat(i,:) = preSeedPref;
    postSeedPrefMat(i,:) = postSeedPref;
end

%%
  
for i = 1:size(con,1)
    for p = 1:size(con,2)
        if con(i,p)
            % scatter(preSeesPostPref(i,p,1),preSeesPostPref(i,p,2)),con(i,p),[preSeedPrefMat
        end
    end
end




%%
clf
hold on
for i = 1:size(con,1)
    for p = 1:size(con,2);
        if (con(i,p))
        scatter(postSeedPrefMat(i,p),preSeedPrefMat(i,p),con(i,p)*500 + rand,'b')
        end
        
        if conPref(i,p)>0
             scatter(postSeedPrefMat(i,p),preSeedPrefMat(i,p),conPref(i,p)*500 + rand,'g')
        end
        if conPref(i,p)<0
                scatter(postSeedPrefMat(i,p),preSeedPrefMat(i,p),abs(conPref(i,p))*500 + rand,'r')

        end
                
    end
end
ylim([-1 2])
xlim([-1 2])
       
%%
    conPref = (con-predSyn);

clf
hold on
scatterShapes = 'ods*.';
for i = 1:size(con,1)
    for p = 1:size(con,2);
        if (con(i,p))
        %scatter(postSeedPrefMat(i,p),con(i,p),50,colList(round(preSeedPrefMat(i,p))+1))
        end
        
        scatter(postSeedPrefMat(i,p)+rand/10,conPref(i,p),50,colList(round(preSeedPrefMat(i,p))+1))

                
    end
end
ylim([-20 20])
xlim([-.1 1.1])

%}


%%

%% show cells
for i = 1:length(axList)
    disp(sprintf('checking axon %d of %d',i,length(axList)))
    preSkel = axSkel(i).skel;
    preNodes = scaleSubs(preSkel.node2subs,voxelScale);
    %preNodes  = [preNodes(:,1) * voxelScale(1) preNodes(:,2) * voxelScale(2) preNodes(:,3) * voxelScale(3)];
    tic
    for p = 1: length(cellList)
        
        postSkel = cellSkel(p).skel;
        postNodes = scaleSubs(postSkel.node2subs,voxelScale);
        %postNodes  = [postNodes(:,1) * voxelScale(1) postNodes(:,2) * voxelScale(2) postNodes(:,3) * voxelScale(3)];
        
        scatter(preNodes(:,3),preNodes(:,2),'.','k')
        hold on
        scatter(postNodes(:,3),postNodes(:,2),'.','r')
        hold off
        disp(sprintf('ax %d and cell %d',axList(i),cellList(p)))
        ylim([0 400])
        xlim([0 300])
        
        pause(.01)
        
        
    end
end

% subplot(1,1,1)
% colSyn = cat(3,predSyn*100,con*100,con*0);
% image(uint8(colSyn))
% scatter(predSyn(:),con(:))

%}


%% Analysis by preference



