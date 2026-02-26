



%% Get data
loadData = 1;
if loadData
    clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end
binaryMat = 0;


%load([springRes 'turkey_plusOne125.mat'])


%% Get skeletons
skelList = 1009;
minSpace = 30;

for i = 1 : length(skelList)
    fileName = sprintf('%sskel\\mat\\%d.mat',MPN,skelList(i));
    load(fileName)
    %%
    skel = cellStruct.arbor;
    node2subs = skel.nodes.node2subs;
    bones = skel.branches;
    
    edges = cat(1,bones.edges, skel.bridges);
    nodes = unique(edges(:));
    nodePos = node2subs(nodes,:);
    
    nep.nodes = nodes;
    nep.edges = edges;
    nep.nodePos = nodePos;
    
    %% fix breaks in skeleton
    
    nep2 = uniteSkel(nep)
%     
%     
%     subplot(1,1,1),clf
%     hold on
%     showEdges =  edges;%edges(groupEdges == i,:);
%     sub1 = node2subs(showEdges(:,1),:);
%     sub2 = node2subs(showEdges(:,2),:);
%     for e = 1:size(showEdges,1)
%         plot([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)],'color','k')
%     end
%     pause(.01)
% 
%     addedEdges = [];
%     while 1
%         %%group nodes
%         g = 0;
%         checkEdge = ones(size(edges,1),1);
%         checkNodes = ones(size(nodes));
%         groupNodes = zeros(size(nodes));
%         
%         while sum(checkNodes)
%             newNodes = nodes(find(checkNodes,1));
%             g = g+1;
%             while ~isempty(newNodes)
%                 oldNodes = newNodes;
%                 newNodes = [];
%                 for n = 1:length(oldNodes);
%                     checkNodes(nodes == oldNodes(n)) = 0;
%                     groupNodes(nodes == oldNodes(n)) = g;
%                     [foundEdges p] = find((edges==oldNodes(n)) & [checkEdge checkEdge]);
%                     checkEdge(foundEdges) = 0;
%                     newNodes = [newNodes; setdiff(edges(foundEdges,:),oldNodes(n))];
%                 end
%             end
%         end
%         sprintf('found %d disconnected groups',max(groupNodes))
%         if max(groupNodes) == 1, break, end
%         
%         %%find tips and link disconnected objects by tips
%         uE = unique(edges(:));
%         hE = hist(edges(:),uE);
%         tips = uE(hE==1);
%         minDist = zeros(size(tips));
%         minNode = zeros(size(tips));
%         for t = 1:length(tips)
%             nodeTarg = find(nodes == tips(t));
%             nodeGroup = groupNodes(nodeTarg);
%             tPos = nodePos(nodeTarg,:);
%             dists = sqrt((nodePos(:,1)-tPos(1)).^2 + (nodePos(:,2)-tPos(2)).^2 + ...
%                 (nodePos(:,3)-tPos(3)).^2 );
%             dists(groupNodes == nodeGroup) = Inf;
%             minDist(t) = min(dists);
%             minNode(t) = nodes(find(dists==minDist(t),1));
%         end
%         
%         bestDist = min(minDist)
%         bestTarg = find(minDist==bestDist,1);
%         
%         newEdge = [tips(bestTarg) minNode(bestTarg)]
%         edges = cat(1,edges,newEdge);
%         addedEdges = cat(1,addedEdges,newEdge);
%         
%     end    
%         
%     showEdges =  newEdge;%edges(groupEdges == i,:);
%     sub1 = node2subs(showEdges(:,1),:);
%     sub2 = node2subs(showEdges(:,2),:);
%     for e = 1:size(showEdges,1)
%         plot([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)],'color','r','linewidth',3)
%     end
%     pause(.01)
%         
%     hold off
    
    %% group edges into bones

    
    %freeze tips and branch points
    hE = hist(edges(:),nodes);
    fixedNode = hE~=2;
    
    checkNodes = ~fixedNode;
    checkEdge = ones(size(edges,1),1);
    groupNodes = zeros(size(nodes));
    groupEdges = zeros(size(edges,1),1);
    listFixed = nodes(fixedNode>0);
    g = 0;
    while sum(checkNodes)
        newNodes = nodes(find(checkNodes,1));
        g = g+1;
        while ~isempty(newNodes)
            oldNodes = newNodes;
            newNodes = [];
            for n = 1:length(oldNodes);
                checkNodes(nodes == oldNodes(n)) = 0;
                groupNodes(nodes == oldNodes(n)) = g;
                
                [foundEdges p] = find((edges==oldNodes(n)) & [checkEdge checkEdge]);
                groupEdges(foundEdges) = g;
                checkEdge(foundEdges) = 0;
                newNodes = [newNodes; setdiff(edges(foundEdges,:),oldNodes(n))];
            end
            newNodes = setdiff(newNodes,listFixed);
        end
    end
  
    %%label zeros
    zeroEdge = find(groupEdges ==0);
    zeroNames = [1:length(zeroEdge)] + max(groupEdges);
    groupEdges(zeroEdge) = zeroNames;
    
    
    subplot(1,1,1),clf
    hold on
    colmap = hsv(100);
    for i = 0:max(groupEdges);
        showEdges =  edges(groupEdges == i,:);
        sub1 = node2subs(showEdges(:,1),:);
        sub2 = node2subs(showEdges(:,2),:);
        rCol = colmap(ceil(rand*100),:)/2;
        for e = 1:size(showEdges,1)
            
            plot([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)],'color',rCol,'linewidth',1)
            
        end
    end
    pause(.01)

    hold off
    
    
    %% min tip
    minTip = 5;
    edgeGroups = unique(groupEdges);    
    hE = hist(edges(:),nodes);

    for i = 1:length(edgeGroups)
        
        gEdges = edges(groupEdges == edgeGroups(i),:);
        hCE = hist(gEdges(:),nodes);
        ends = nodes(hCE==1);
        isTip = 0; 
        if isempty(ends)
            isTip = 1;
            startN = gEdges(1,1);
        elseif hE(find(nodes == ends(1))) == 1
            startN = ends(1);
            isTip = 1;
        elseif hE(find(nodes == ends(2))) == 1
            startN = ends(2);
            isTip = 1;
        else
            isTip = 0;
            startN = ends(1);
        end
       
        newNode = startN;
        c = 0;
        checkEdges = ones(size(gEdges));
        newEdges = zeros(size(checkEdges));
        L = zeros(size(checkEdges,1),1);
        tog = [2 1];
        while 1
            c = c+1;
            oldNode = newNode;
            [y x] = find((gEdges == oldNode) & checkEdges,1);
            if ~isempty(y)
                
                checkEdges(y,:) = 0;
                newNode = gEdges(y,tog(x));
                sub1 = node2subs(nodes==oldNode,:);
                sub2 = node2subs(nodes==newNode,:);
                %             scatter(sub1(e,1), sub1(e,2),'b')
                %             scatter(sub2(e,1),sub2(e,2),'r','^')
                
                L(c) = sqrt((sub1(1)-sub2(1)).^2 +  (sub1(2)-sub2(2)).^2 + ...
                    (sub1(3)-sub2(3)).^2);
                newEdges(c,:) = [oldNode newNode];
            else
                break
            end
        end
          
         bRun(i).edges = newEdges;
    bRun(i).lengths = L;
    bRun(i).isTip = isTip;
        
    end
    
   useRun = ones(length(bRun),1); 
   for b = 1:length(bRun)
       if bRun(b).isTip
            if sum(bRun(b).lengths)<minTip
                useRun(b) = 0;
            end
       end
   end
    
  spurs = bRun(~useRun>0);
  runs = bRun(useRun>0);
  
  
  %% Regroup
  
  
  
  
  
   
   
    
  return

    
    
    %%
    
    
    
    
    
    
    
    
    
    simpleEdges = skel.bridges;
    
    bNodes = cat(1,bones.edges);
    hNodes =  hist(bNodes(:),unique(bNodes(:)));
    hHist = hist(hNodes, unique(hNodes));
    
    parentSubs = node2subs(parents,:);
    scatter(parentSubs(:,1),parentSubs(:,2))
    
    clf
    for b = 1:length(bones)
        
        %bRun = [ bones(b).edges; bones(b).base bones(b).parent] ;
        bRun = [bones(b).edges] ;
        sub1 = node2subs(bRun(:,1),:);
        sub2 = node2subs(bRun(:,2),:);
        %        scatter(sub1(:,1),sub1(:,2),'b','+')
        %        hold on
        %        scatter(sub2(:,1),sub2(:,2),'r','.')
        
        clear newEdges eL
        
        newEdges(1,1) = bRun(1,1);
        edgeCount = 1;
        L = 0;
        
        for e = 1:size(bRun,1)-1
            %             scatter(sub1(e,1), sub1(e,2),'b')
            %             scatter(sub2(e,1),sub2(e,2),'r','^')
            
            L = L + sqrt((sub1(e,1)-sub2(e,1)).^2 +  (sub1(e,2)-sub2(e,2)).^2 + ...
                (sub1(e,3)-sub2(e,3)).^2);
            if (L >= minSpace) || fixedNode(bRun(e,2))>0
                %if  fixedNode(bRun(e,2))>0
                
                eL(edgeCount) = L;
                L = 0;
                newEdges(edgeCount,2) = bRun(e,2);
                edgeCount = edgeCount + 1;
                newEdges(edgeCount,1) = bRun(e,2);
                %               scatter(sub2(e,1),sub2(e,2),900,'.','g')
            end
            
        end
        newEdges(end,2) = bRun(end,2);
        
        %         %%Get new length between nodes
        %         sub1 = node2subs(newEdges(:,1),:);
        %         sub2 = node2subs(newEdges(:,2),:);
        %          hold on
        %         Ls = zeros(size(newEdges,1),1);
        %         for e = 1:size(newEdges,1)
        %
        %             plot([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)])
        %             Ls(e) = sqrt((sub1(e,1)-sub2(e,1)).^2 +  (sub1(e,2)-sub2(e,2)).^2 + ...
        %                 (sub1(e,3)-sub2(e,3)).^2);
        %
        %         end
        %           pause(.4)
        
        simpleEdges = cat(1,simpleEdges,newEdges);
        
        %scatter(cellNodePos(:,1),cellNodePos(:,2),'.'),pause(.01)
    end %bones
    
    %%Get new length between nodes
    sub1 = node2subs(simpleEdges(:,1),:);
    sub2 = node2subs(simpleEdges(:,2),:);
    hold on
    Ls = zeros(size(simpleEdges,1),1);
    for e = 1:size(simpleEdges,1)
        
        plot([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)],'k')
        Ls(e) = sqrt((sub1(e,1)-sub2(e,1)).^2 +  (sub1(e,2)-sub2(e,2)).^2 + ...
            (sub1(e,3)-sub2(e,3)).^2);
        
    end
    pause(.4)
    
    %%
    
end


%% Filter nodes
nodes = useList.nodes;
seedFilt = 108;
seedPos = find(nodes==seedFilt);
con = useList.con;
symCon = con + con';
[cellGroup] = segmentCon(symCon);
seedGroup = cellGroup(seedPos);
nodes = nodes(cellGroup == seedGroup);

nodeType = nodes*0;
for i = 1:length(nodes)
    nodeType(i) = useList.nodeType(useList.nodes==nodes(i));
end
nodeIDs = nodes;




nodeNum = length(nodeIDs);
lookUpID(nodeIDs+1) = 1:length(nodeIDs);


allEdges(:,3) = 1;
allEdges(allEdges(:,1) == plusOne,3) = 1;
allEdges(allEdges(:,2) == plusOne,3) = 1;

%% Set color
nodeCol = zeros(nodeNum,3);

if 1
    
    [nodeCol, useCol] = getList_seedColor(seedList,nodeIDs) ;
    nodeCol = nodeCol * .8;
    nodeCol(sum(nodeCol,2) == 0,:) = .3;
    
    seedCol = [1 0 0; 0 1 0; .5 0 1 ; 0 .5 1];
    for s = 1:length(seedList)
        targ = find(nodeIDs == seedList(s));
        nodeType(targ) = nodeType(targ)+2;
        nodeCol(targ,:) = seedCol(s,:);
    end
    
    showCol = permute(useCol.col,[1 3 2])
    
    image(uint8(showCol*256))
    useCol.comb{10}
    
    targ = find(nodeIDs == plusOne);
    nodeCol(targ,:) = [10 10 10]
    
end


%%Set color according to attribute
if 0
    %[colorList cellCol] = getAttributes(obI);
    
    load([MPN 'cb2d.mat'])
    colorList = cb2d.IDs;
    nodeCol = nodeCol + .2;
    
    colPropRaw = cb2d.areaUM;
    
    nodeProp = [];nodePropRef = [];
    for i = 1:length(nodeIDs)
        targ = find(colorList == nodeIDs(i));
        if ~isempty(targ)
            if length(targ)>1
                'too many targets'
                colorList(targ)
            end
            nodeProp= [nodeProp colPropRaw(targ(1))];
            nodePropRef = [nodePropRef i];
            
        end
    end
    
    colProp = nodeProp-min(nodeProp)+1;
    colProp = ceil(colProp/max(colProp) * 100);
    
    colTable = jet(100);
    cellCol = colTable(colProp,:);
    
    nodeCol(nodePropRef,:) = cellCol;
    nodeType(nodePropRef) = nodeType(nodePropRef) + 2;
    
    cat(2,nodeIDs', nodeCol);
end

%%Label specific cells
if 0
    crossoverAxons = [2032	2033	2034	2035]
    gotList = getList_giantBoutons(MPN);
    labelCells = gotList;
    for i = 1:length(labelCells)
        nodeCol(find(nodeIDs==labelCells(i)),:) = nodeCol(find(nodeIDs==labelCells(i)),:) + .5;
    end
end


nodeCol(nodeCol>1) = 1;
nodeCol(nodeCol<0) = 0;

%% Create springDat (need nodeCol, nodeIDs,

springIn.nodeIDs = nodeIDs;
springIn.allEdges = allEdges;
springIn.allWeights = allEdges ;
springIn.nodeCol = nodeCol;
springIn.nodeType = nodeType;
springIn.seedList = seedList;

springDat = skelSpringParameters_01(springIn);

if binaryMat
    springDat.edges.ew = springDat.edges.ew>0;
end

%% Run springs
movDir  = [springDir 'seedMovie_?\'];
if ~exist(movDir,'dir'), mkdir(movDir), end

for rerun = 1: 1
    allResults{rerun} = runSkelSprings(springDat);
end

%%  print eps
if 0
    springDir = 'D:\LGNs1\Analysis\springDat\skelSpring\';
    if ~exist(springDir,'dir'), mkdir(springDir), end
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'net125_white';
    
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
end



%% Save result
if 0
    
    springRes = 'D:\LGNs1\Analysis\springDat\results\';
    if ~exist(springRes,'dir'), mkdir(springRes), end
    
    result = allResults{end};
    save([springRes 'turkey_plusOne125.mat'],'result')
    
    
end






