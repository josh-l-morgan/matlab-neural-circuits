

clear all
springRes = 'D:\LGNs1\Analysis\springDat\results\';
%springRes = 'D:\LGNs1\Analysis\springDat\results\randNoSeedsFilt1';


%load([springRes 'res_snapLong01.mat'])
%load([springRes 'res_snapLongWithSeeds.mat'])
%load([springRes 'res_snapLongWithSeedsRandomEdgeRemoval.mat'])
%load([springRes 'res_snapLongNoSeedsAllSynRemoval.mat'])
%load([springRes 'res_snapLongNoSeedsLowDispersion.mat'])
%load([springRes 'res_snapLongNoSeedsLowNoise.mat'])
% load([springRes 'res_snapLongNoSeedsHighDispersion.mat'])
% load([springRes 'res_snapLongNoSeedsVeryLosDispersion'])
% load([springRes 'res_snapLongNoSeeds4'])
% load([springRes 'res_snapLongNoSeedsNice'])
%load([springRes 'res_snapLongNoSeedsRand1'])
%load([springRes 'res_realNoSeedFilt1Wait10c'])

% load('D:\LGNs1\Analysis\springDat\results\res_snapLongNoSeedsAll.mat')


%%Standard run
load('D:\LGNs1\Analysis\springDat\results\res_snapLongNoSeedsAll_2.mat')
%% Rerun with 4 seeds
load('D:\LGNs1\Analysis\springDat\results\res_snapLongNo4Seeds.mat')
load('D:\LGNs1\Analysis\springDat\results\IB_binSyn.mat')


stopOnes = 0;

results = result;%allResults{1};

snapTime = results.snapTime;
nodeIDs = results.nodeIDs;
cellGroups = results.cellGroups;

cladeCounts = max(cellGroups,[],2);

num = length(nodeIDs);
steps = length(snapTime);
layers = 1:steps;

cladeIDmat = cellGroups + repmat(layers' * num*2,[1 num]);
maxG = max(cellGroups(:));
cladeIDs = unique(cladeIDmat(:));
newGroupIDs = 1:length(cladeIDs);
lookupID(cladeIDs) = newGroupIDs;
cladeIDmat = lookupID(cladeIDmat);

cladeIDs = newGroupIDs;
cladeNum = length(cladeIDs);



%%
clear cladedNodes e2 height
for i = 1:cladeNum
    [y x] = find(cladeIDmat==i);
    g = cellGroups(y(1),x(1));
    xPos(i) = (g * (num/cladeCounts(y(1))))-(num/cladeCounts(y(1))/2);
    cladedNodes{i} = nodeIDs(x);
    yPos(i) = y(1);
    
    
    if y(1)>1
        upvec = cladeIDmat(y(1)-1,x);
        e2(i) = unique(upvec);
        recY(i) = y(1);
    end
end


%e2(1) = 1;
allEdges = [cladeIDs(find(e2>0,1):end)' e2(find(e2>0,1):end)'];





%% Get node order
cladeOrder{1} = unique(cladeIDmat(1,:));
for h = 1:max(yPos)-1
    
    isUp = cladeOrder{h};
    %isUp = find(yPos == h);
    %downNum = cladeCounts(h+1);
    
    nextOrder = [];
    for u = 1:length(isUp)
        isDown = find(e2 == isUp(u));
        nextOrder = cat(2,nextOrder,isDown);
    end
    cladeOrder{h+1} = nextOrder;
    L = length(nextOrder);
    xPos(nextOrder) = (((1:L)/L)-(1/L/2))*num;
    
    %    if L~= sum(recY(h))
    %        'bad'
    %        return
    %    end
end

%%Correct node order
lastGroupOrder = cladeOrder{end};
nodeOrder = [cladedNodes{lastGroupOrder}]';
lookupNodes(nodeOrder) = 1:length(nodeOrder);

%%Reorder nodes in clades
for i = 1:length(cladeIDs)
    
    theseNodes = cladedNodes{i};
    theseOrders = lookupNodes(theseNodes);
    [sortNodes idx] = sort(theseOrders,'ascend');
    cladedNodes{i} = theseNodes(idx);
    xPos(i) = mean(theseOrders);
end


%% groupMat

nodeNum = size(cellGroups,2);
cutNum = size(cellGroups,1);
relMat = zeros(nodeNum);
for i = 1:nodeNum
    for p = 1:nodeNum
        for h = cutNum:-1:1;
            hit = (cellGroups(h,i) == cellGroups(h,p));
            if hit
                %i,p,hit,pause
                relMat(i,p) = h;
                break
            end
        end
    end
end
relMat = relMat/cutNum;

for i = 1:length(nodeIDs)
    transNode(i) = find(nodeOrder == nodeIDs(i));
    transNode(i) = find(nodeIDs == nodeOrder(i));
    
end

rawRelMat = relMat;
relMat = relMat(transNode,:);
relMat = relMat(:,transNode);

image(relMat*100)
colormap(bluered(cutNum))
colormap(gray(100))

cladeResults.nodeOrder = nodeOrder;
cladeResults.relMat = relMat;
cladeResults.rawRelMat = rawRelMat;
%% get nodeColor
nodeCol = zeros(length(nodeOrder),3);

% 
% seedList = [ 108 201 109 907 903];
% useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
% synEdges = obI.nameProps.edges(:,[2 1]);
% 
% seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .0 1; 0 0 1];
% for s = 1:length(seedList)
%     useCol = seedCol(s,:);
%     for n = 1:length(nodeOrder)
%         nodeCol(n,:) = nodeCol(n,:) + useCol *  ...
%             double(sum((synEdges(:,1)==nodeOrder(n)) & (synEdges(:,2)) == seedList(s))>0);
%     end
% end

%% get nodeColor


checkIDs = result.useList.preList;
checkProp = checkIDs;



%% Reshape color

if 0
colProp = checkProp-min(checkProp);
colProp = ceil(colProp/max(colProp) * 99)+1;
colTable = bluered(100);
cellCol = colTable(colProp,:);

else
   colMap = hsv(256);
    colMap = cat(1,[0 0 0],colMap);
    checkNum = length(checkIDs);
    cellCol = colMap(ceil((1:checkNum)*256/checkNum),:);
    %cellCol = cellCol(randperm(checkNum),:);
    
end
    
nodeCol = zeros(length(nodeOrder),3);
useCol = zeros(length(nodeOrder),1);
    for n = 1:length(nodeOrder)
       targ = find(checkIDs == nodeOrder(n),1);
       if isempty(targ)
           nodeCol(n,:) = [0 0 0];
           useCol(n) = 0;
       else
          nodeCol(n,:) = cellCol(targ,:) ;
          useCol(n) = 1;
       end
    end


%% get clade color

cladeCol = zeros(length(cladeIDs),3);

if 1
    
    for i = 1:length(cladeIDs)
        theseNodes = lookupNodes(cladedNodes{i});
        theseColors = nodeCol(theseNodes,:);
        useThese = find(useCol(theseNodes));
        
        if isempty(useThese)
            newColor = [.1 .1 .1];
        else
            newColor = mean(theseColors(useThese,:),1);
        end
        cladeCol(i,:) = newColor;
    end
      
else
    
    for i = 1:length(cladeIDs)
        theseNodes = lookupNodes(cladedNodes{i});
        theseColors = nodeCol(theseNodes,:);
        sumColors = sum(theseColors,1);
        maxColor = max(sumColors);
        if maxColor
            newColor = sumColors/maxColor;
        else
            newColor = [.5 .5 .5];
        end
        cladeCol(i,:) = newColor;
    end
end

if stopOnes
    
    terminals = []; terminalNum = 0;
    terminals = zeros(length(nodeOrder),2);
    if 1 %% Blackout singles
        for i = find(e2>0,1):length(cladeIDs)
            parentClade = cladedNodes{e2(i)};
            if length(parentClade) <2
                cladeCol(i,:) = cladeCol(e2(i),:) * 0;
                targ = lookupNodes(parentClade);
                if ~sum(terminals(targ,:))
                    terminals(targ,2) = xPos(e2(i));
                    terminals(targ,1) = yPos(e2(i));
                    terminalColor(targ,:) = cladeCol(e2(i),:);
                end
            end
        end
    end
end
%% Show data
%scatter(xPos,yPos,'.')

clf



for e = 1:size(allEdges,1)
    plot([xPos(allEdges(e,1)) xPos(allEdges(e,2))],...
        [yPos(allEdges(e,1)) yPos(allEdges(e,2))],'color',cladeCol(allEdges(e,1),:),...
        'LineWidth',2)
    hold on
end

if stopOnes
    sg = scatter(terminals(:,2),terminals(:,1),8,terminalColor,'filled');
    %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
    %set(sg,'MarkerEdgeColor','w');
end
xlim([-1 num+1])


axis on
set(gcf,'color','w')
set(gca,'color','k')
%
%         for p = 1:length(pl)
%             set(pl(p),'color',edgeCol(p,:),'LineWidth',edges.width(p));
%         end
hold off


%% save results
cladeResults.results = results;
cladeResults.cladeIDs = cladeIDs;
cladeResults.e2 = e2;
cladeResults.allEdges = allEdges;
cladeResults.cladedNodes = cladedNodes;
cladeResults.yPos = yPos;
cladeResults.xPos = xPos

%{
save([springRes 'res_snapLong01_randClade1.mat'],'cladeResults')

%}



if 0
    %%
    
    itag  = 'IB_clade';
    springDir = 'D:\LGNs1\Analysis\springDat\figDraft1b\';
    
    figWin = [10 10 1900 400];
    outerWindow  = get(gca,'OuterPosition');
    axisWin  = [0.02 0.02 0.96 0.96];
    
    set(gcf,'Position',figWin)
    set(gca, 'Position', axisWin)
    set(gcf, 'InvertHardCopy', 'off');
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 1000 400])
    
    
    imageName = sprintf('%s%s.png',springDir,itag);
    %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    
    epsName = sprintf('%s%s.eps',springDir,itag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
    
end












