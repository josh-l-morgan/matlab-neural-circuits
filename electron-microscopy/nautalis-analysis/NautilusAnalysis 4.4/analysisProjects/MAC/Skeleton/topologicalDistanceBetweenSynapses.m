

%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    MPN = 'Z:\joshm\LGNs1\Exports\MAC_export\MAC_merge_mat\';
    WPN = MPN;
    DPN = [MPN 'data\'];
    if ~exist(DPN,'dir'),mkdir(DPN);end
    
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903 125];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
    % useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end
binaryMat = 0;

postList = allEdges(:,2);
postList = unique(postList(postList>0));
%load([springRes 'turkey_plusOne125.mat'])

cells = obI.cell.name(obI.cell.isCell>0);
cellList = intersect(postList,cells);
%% Get skeletons
skelList = cellList;
minTip = 1;
minSpace = 1;

for i = 1 : length(skelList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,skelList(i));
    load(fileName)
    %%
    skel = cellStruct.arbor;
    node2subs = skel.nodes.node2subs;
    bones = skel.branches;
    
    edges = cat(1,bones.edges, skel.bridges);
    nodes = 1:size(node2subs,1); % will go wrong if gaps
    
    nodePos = node2subs(nodes,:);
    voxelScale = obI.em.dsRes; %adjusted for skeleton down samp
    nodePos = scaleNodes(nodePos,voxelScale);
    nep.nodePos = nodePos;
    
    nep.nodes = nodes;
    nep.edges = edges;
    
    
    
    %% condition skeleton
    
    showNep(nep)
    
    nep = uniteSkel(nep);%% fix breaks in skeleton
%     nep = groupNepEdges(nep); %% group edges into bones
%     nep = edgeGroup2bones(nep);
%     nep = bonesVsSpurs(nep,minTip); %% sort bones and spurs
%     
    %%redo without spurs
    %nep.edges = cat(1,nep.bones.edges);
    %nep = groupNepEdges(nep);
    %nep = edgeGroup2bones(nep);
    %showNepBones(nep,2);
    
    %%Simplify skeleton
%     nep = simpleNep(nep,minSpace);
    % showNepBones(nep)
    %nep = cleanBones(nep);
    %showNepBones(nep);
    
    skels(i).id = skelList(i);
    skels(i).nep = nep;
    
end

%% Load synapses

nepSyn = obISyn2nep(obI);
%nepSyn = linkSyn2Skel(nepSyn,skels);



%% One at a time


syn = nepSyn.cellEdges
cellList = [skels.id];

for i = 1: length(cellList)
    disp(sprintf('running %d of %d',i,length(cellList)));
    cellID = cellList(i);
    skel = skels(i).nep;
    sPos = skel.nodePos;
    
    isSyn = find(syn(:,2) == cellID);
    
    if length(isSyn)>1
        synPos = nepSyn.nodePos(isSyn,:);
        syns = 1:length(isSyn);
        
        pairs = nchoosek(syns,2);
        clear pLength eucDist keepPaths keepLengths
                
        for p = 1:size(pairs,1)
            
            pos1 = synPos(pairs(p,1),:);
            pos2 = synPos(pairs(p,2),:);
            
            eucDist(p) = sqrt((pos1(1)-pos2(1)).^2 + (pos1(2)-pos2(2)).^2 + ...
                (pos1(3)-pos2(3)).^2 );
            
            
            syn2Skel = sqrt((sPos(:,1)-pos1(1)).^2 + (sPos(:,2)-pos1(2)).^2 + ...
                (sPos(:,3)-pos1(3)).^2 );
            link1 = (find(syn2Skel==min(syn2Skel),1));
            
            syn2Skel = sqrt((sPos(:,1)-pos2(1)).^2 + (sPos(:,2)-pos2(2)).^2 + ...
                (sPos(:,3)-pos2(3)).^2 );
            link2 = (find(syn2Skel==min(syn2Skel),1));
            
            %% find path
            if link1 ~= link2
                edges = skel.edges;
                nodePos = sPos;
                startNode = link1;
                stopNode = link2;
                
                
                doubleCheck = 2; % factor of spreads past first hitting target
                reps = size(edges,1);
                nodeNum = size(nodePos,1);
                pred = zeros(nodeNum,1);
                dists = inf(nodeNum,1);
                dists(startNode) = 0;
                lengths = sqrt((nodePos(edges(:,1),1)-(nodePos(edges(:,2),1))).^2 + ...
                    (nodePos(edges(:,1),2)-(nodePos(edges(:,2),2))).^2 + ...
                    (nodePos(edges(:,1),3)-(nodePos(edges(:,2),3))).^2);
                
                firstHit = 0;
                for r = 1:reps*10
                    
                    
                    oldDists = [dists(edges(:,1)) dists(edges(:,2))]; %current distance of nodes in edges
                    newDists = oldDists(:,[2 1]) + repmat(lengths,[1 2]); %Distance of nodes if connected through edge
                    closer = newDists<oldDists; %nodes who would benefit from given edge
                    
                    [e a] = find(closer);
                    b = 3-a;
                    n = edges(sub2ind(size(edges),e,a));
                    n2 = edges(sub2ind(size(edges),e,b));
                    
                    for ce = 1:length(n)
                        if newDists(e(ce),a(ce)) < dists(n(ce));
                            pred(n(ce)) = n2(ce);
                            dists(n(ce)) = newDists(e(ce),a(ce));
%                             scatter3(sPos(n2(ce),1),sPos(n2(ce),2),sPos(n2(ce),3),'y','filled')
%                             scatter3(sPos(n(ce),1),sPos(n(ce),2),sPos(n(ce),3),'y','filled')
%                             
%                             pause(.01)
                        end
                    end
                    
                    
                    %% Check for finish
                    if sum(n == stopNode) & ~firstHit;
                        firstHit = r;
                    end
                    
                    if firstHit & (r > (firstHit * doubleCheck))
                        break
                    end
                    
                    
                    % Show
                    %             predNodes = pred(pred>0);
                    %             scatter3(sPos(predNodes,1),sPos(predNodes,2),sPos(predNodes,3),'y','filled')
                    %             pause(.01)
                    
                    
                end
                
                
                %% read path
                
                lastNode = stopNode;
                clear path
                for pa = 1:reps
                    %             scatter3(sPos(lastNode,1),sPos(lastNode,2),sPos(lastNode,3),'b','filled')
                    %             pause(.1)
                    predNode = pred(lastNode);
                    path(pa,:) = [predNode lastNode];
                    lastNode = predNode;
                    if lastNode == startNode
                        break
                    end
                end
                path = flipud(path);
                
                
                
                scatter3(sPos(:,1),sPos(:,2),sPos(:,3),2,'k','filled')
                hold on
                
                %% path length
                clear pLengths
                for pa = 1:size(path,1)
                    pLengths(pa) = sqrt((sPos(path(pa,1),1) - sPos(path(pa,2),1)).^2 + ...
                        (sPos(path(pa,1),2) - sPos(path(pa,2),2)).^2 +...
                        (sPos(path(pa,1),3) - sPos(path(pa,2),3)).^2);
                    plot3(sPos(path(pa,:),1),sPos(path(pa,:),2),...
                        sPos(path(pa,:),3),'r','linewidth',4)
                end
                
                pLength(p) = sum(pLengths);
                keepPaths{p} = path;
                keepLengths{p} = pLengths;
                
                %% show path
                
                scatter3(pos1(1),pos1(2),pos1(3),'r','filled')
                scatter3(pos2(1),pos2(2),pos2(3),'g','filled')
                scatter3(sPos(link1,1),sPos(link1,2),sPos(link1,3),'m','filled')
                scatter3(sPos(link2,1),sPos(link2,2),sPos(link2,3),'c','filled')
                
                scatter3(sPos(path(:),1),sPos(path(:),2),sPos(path(:),3),'b','filled')
                sum(pLengths)
                eucDist(p);
                pathRat = sum(pLengths)/eucDist(p);
                pause(.1)
                %pause
                hold off
                
                
            else % same node on arbor
                pLength(p) = 0;
                 keepPaths{p} = [];
                keepLengths{p} = [];
            end
            
        end %run all pairs
        
    else % no pairs on cell
        pLength = [];
        eucDist = [];
        keepPaths = [];
        keepLengths = [];
    end
    
    cellDists(i).pLength = pLength;
    cellDists(i).eucDist = eucDist;
    cellDists(i).cellName = skelList(i);
    cellDists(i).keepPaths = keepPaths;
    cellDists(i).keepLengths = keepLengths;
    cellDists(i).sPos = sPos;
end

topoPairDists = cellDists;
save([DPN 'topoPairDists4.mat'],'topoPairDists')


%% Plot 

c = 0;
clear eucDists pLengths
for t = 1:length(topoPairDists)
    pLength = topoPairDists(t).pLength;
    eucDist = topoPairDists(t).eucDist;
    for p = 1:length(pLength)
        c = c+1;
        eucDists(c) = eucDist(p);
        pLengths(c) = pLength(p);
    end
end

scatter(eucDists,pLengths,'k','o','filled')
hold on
plot([0 30],[0 30],'r')
hold off




if 0
    %% test edges
    edges = skel.edges;
    startNode = 1
    
    scatter(sPos(:,1),sPos(:,2),'k','filled')
    hold on
    
    lastNode = startNode;
    clear path
    use = edges * 0 + 1;
    
    for pa = 1:size(edges,1)
        scatter3(sPos(lastNode,1),sPos(lastNode,2),sPos(lastNode,3),'b','filled')
        pause(.01)
        predNode = [];
        for n = 1:length(lastNode)
            [e a] = find((edges==lastNode(n)) & ( use));
            b = 3-a;
            predNode = cat(1,predNode,edges(sub2ind(size(edges),e,b)));
            %path(pa,:) = [predNode lastNode];
            use(e,:) = 0;
        end
        lastNode = unique(predNode);
        
        
    end
    
    
    hold off
    
    
    
    
end



