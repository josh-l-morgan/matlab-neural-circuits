
%% Plot axon to tcr
clear all
MPN = GetMyDir;
load([MPN 'obI.mat'])
shouldPause = 1;


plotLimits = [-400 -200;0 400];
set(gcf,'Position',[50 50 1800 900])
%% variables
dimOrder = [1 3 2];
seedList = [108 201]

crossoverAxons = [2032	2033	2034	2035]

%% Get cells
tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
linList = obI.nameProps.cellNum(obI.nameProps.lin);
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

preList = [];
for i = 1:length(seedList)
    isPost = edges(:,1) == seedList(i);
    preList = [preList ;edges(isPost,2)];
end
axList = intersect(unique(preList),rgcList);
axList = axList((axList>=1000) & (axList<5000));

% axList = setdiff(axList, crossoverAxons)
% axList = axList((axList>=1000) & (axList<5000));
% 
% axList = axList(axList<2000);
% axList = [axList crossoverAxons]



postList = [];
for i = 1:length(axList)
    isPre = edges(:,2) == axList(i);
    postList = [postList; edges(isPre,:)];
end
cellList = intersect(unique(postList),tcrList);
cellList = cellList((cellList>0) & (cellList<1000));

%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end

%% Pick colors
cmap = [0 0 0; 1 0 0 ; 0 1 0; 0 0 1];

preColIdx = size(cmap,1)+1:size(cmap,1)+length(axList);
cellColIdx = max(preColIdx)+1:max(preColIdx)+length(cellList);

%%Axon color
PreConToSeed = zeros(size(con,1),3);
for i = 1:length(seedList)
    PreConToSeed(:,i) = con(:,find(cellList == seedList(i)));
end

PreCol(:,1) = PreConToSeed(:,1)/mean(PreConToSeed(PreConToSeed(:,1)>0,1))*.25;
PreCol(:,2) = PreConToSeed(:,2)/mean(PreConToSeed(PreConToSeed(:,2)>0,2))*.25;
% PreCol(:,3) = 1-(PreCol(:,1)-PreCol(:,2));
% PreCol(PreCol<0) = 0;
PreCol(:,3) = .5;

cmap(preColIdx,:) = PreCol;

PostConToSeed = zeros(size(con,2),3);
    for p = 1:length(cellList)
        for i = 1:length(axList)
        PostConToSeed(p,:) = PostConToSeed(p,:) +  sqrt(PreConToSeed(i,:)*con(i,p));
        %PostConToSeed(p,:) = PostConToSeed(p,:) +  min(PreConToSeed(i,:)*con(i,p));
    end
end

PostCol(:,1) = PostConToSeed(:,1)/mean(PostConToSeed(PostConToSeed(:,1)>0,1))*100;
PostCol(:,2) = PostConToSeed(:,2)/mean(PostConToSeed(PostConToSeed(:,2)>0,2))*100;
PostCol(:,3) = 0;
cmap(cellColIdx,:) = PostCol;

cmap(cmap>1) = 1;
colormap(cmap)


showA = [axList' PreCol];
showC = [cellList' PostConToSeed];



%% get anchors
hold on

dSamp = [.016 0.016 0.032];
rawAnchors = zeros(length(cellList),3);
for i = 1:length(cellList)
    rawAnchors(i,:) = obI.cell.anchors(find(obI.cell.names == cellList(i)),:);
end
anchors(:,dimOrder(1)) = rawAnchors(:,1).*dSamp(1)*-1;
anchors(:,dimOrder(2)) = rawAnchors(:,2).*dSamp(2);
anchors(:,dimOrder(3)) = rawAnchors(:,3).*dSamp(3);


%anchors = round(anchors);
%anchors(anchors<1) = 1;

%synapses = cat(1,[ 0 0 0], synapses);
rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors(:,dimOrder(1)) = rawSynAnchors(:,1).*dSamp(1)*-1;
synAnchors(:,dimOrder(2)) = rawSynAnchors(:,2).*dSamp(2);
synAnchors(:,dimOrder(3)) = rawSynAnchors(:,3).*dSamp(3);



%% axon anchors
%%Get axon information
hold on
seedColList = [1 0 0; 0 1 0];

for i = 1:length(axList)
    allAnchors =  synAnchors(edges(:,2) == axList(i),:);
    axAnchor(i,:) = mean(allAnchors,1);
end

%%Axon anchor Circle
h = 100;
meanAnchor = mean(axAnchor,1);
for i = 1:length(axList)
    oldAnch = axAnchor(i,:);
    slope = oldAnch - meanAnchor;
    x = slope(2);
    y = slope(1);
    c = sqrt((x^2 +y^2)/(h^2));
    circAxAnchor(i,:) = ...
        [oldAnch(1) - y/c oldAnch(2) - x/c oldAnch(3)];
    scatter(y,x,'r')
    scatter(y/c,x/c,'g')
    sqrt((y/c)^2+(x/c)^2)
end
%circAxAnchor = axAnchor;





%%

clf
subplot(1,1,1)

hold on
cellCount = zeros(length(cellList),3);
for i = 1:length(axList);
    for p = 1:length(cellList)
        
        axCol = PreCol(i,:);
        conWeight = con(i,p);
        if conWeight>0
            
            axCol(axCol>1) = 1;
            plotY = [circAxAnchor(i,1) anchors(p,1)];
            plotX = [circAxAnchor(i,2) anchors(p,2)];
            lineWidth = max(1,sqrt(conWeight));
            plot(plotX, plotY,'LineWidth',lineWidth,'color',axCol);
           % plot(plotX, plotY,'LineWidth',1,'color',axCol);
            
        end
        
    end
end
% cellColor(:,1) = cellColor(:,1)/mean(cellColor(:,1))/4;
% cellColor(:,2) = cellColor(:,2)/mean(cellColor(:,2))/4;
% cellColor(:,3) = .5;
%


ylim(plotLimits(1,:));
xlim(plotLimits(2,:));
myVectorFile = [MPN 'conGraph\draftLines.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')


%%

if shouldPause
pause
clf
end

%colormap(PreCol)
%scatter(circAxAnchor(:,2),circAxAnchor(:,1),[1:length(axList)]',40, 'filled')
% scatter(circAxAnchor(:,2),circAxAnchor(:,1),50,cmap(preColIdx,:),'filled','d',...
%     'MarkerEdgeColor','w','LineWidth',1)
hold on
for i = 1:length(axList)
    plot([circAxAnchor(i,2) circAxAnchor(i,2)],[circAxAnchor(i,1) circAxAnchor(i,1)],'color',cmap(preColIdx(i),:),'Marker','d',...
        'MarkerFaceColor',cmap(preColIdx(i),:),'MarkerEdgeColor','w','MarkerSize',8)
end

ylim(plotLimits(1,:));
xlim(plotLimits(2,:)); 


if shouldPause
pause
clf
end

hold on
%%f



%scatter(cbPos(:,2),cbPos(:,1),30,'filled','b')
% scatter(anchors(:,2),anchors(:,1),100,cellColIdx+1,'filled','o',...
%     'MarkerEdgeColor','w','LineWidth',1)


 scatterBar(PostConToSeed(:,1:2),anchors)
hold on


% 
% for i = 1:length(cellList)
%     text(anchors(i,2),anchors(i,1),num2str(cellList(i)),...
%         'FontSize',6,'HorizontalAlignment','center','FontWeight','bold','color','w')
% end

highlightCell = seedList

cbPos = zeros(length(highlightCell),3);
seedCol = zeros(length(highlightCell),1);
seedCol = {'r','g','b'};
for i = 1:length(highlightCell)
    targ = find(cellList == highlightCell(i));
    cbPos(i,:) = anchors(targ,:);
    seedColIdx(i) = cellColIdx(targ);
    scatter(cbPos(i,2),cbPos(i,1),250,seedCol{i},'filled','o',...
    'MarkerEdgeColor','w','LineWidth',1)
end

for i = 1:length(highlightCell)
    text(cbPos(i,2)+1,cbPos(i,1),num2str(highlightCell(i)),...
        'FontSize',6,'HorizontalAlignment','center','FontWeight','bold')
end

ylim(plotLimits(1,:));
xlim(plotLimits(2,:)); 
myVectorFile = [MPN 'conGraph\draftcirc.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')


%% synDistribution
 showSynPostList = [201]
 preToShowCell = unique(edges(edges(:,1)== showSynPostList,2))
 [preToShowCell ia ib] = intersect(axList,preToShowCell);
 

clf
PreConToConTo = zeros(length(axList),3);
for i = 1:size(con,1)
    for p = 1:size(con,2)
       PreConToConTo(i,:) =  PreConToConTo(i,:) + sqrt(con(i,p) * PostConToSeed(p,:));
    end
end

barHeight = 15;
 scatterBar(PreConToConTo(ia,1:2),circAxAnchor(ia,:),15)
 
 
 ylim(plotLimits(1,:));
xlim(plotLimits(2,:)); 

 
 %%
 %clf
 
 

isPost = find(edges(:,1) == showSynPostList)
synPos = synAnchors(isPost,:);
scatter(synPos(:,2),synPos(:,1),30,'k','filled')
 
 hold on
 seedShapes = {'o','s'};
 for p = 1:length(showSynPostList)
     for i = 1:length(axList)
         isAx = edges(:,2) == axList(i);
         isSeed = edges(:,1) == showSynPostList(p);
         synPos = synAnchors(isAx & isSeed,:);
         synCol = PreConToConTo(i,:) ;
         synCol(synCol>1) = 1;
         scatter(synPos(:,2),synPos(:,1),20,synCol,'filled','Marker',seedShapes{p})
        
         
     end
 end
 
 
 
ylim(plotLimits(1,:));
xlim(plotLimits(2,:)); 


%%
%clf
 
  for p = 1:length(showSynPostList)
     for i = 1:length(axList)
         isAx = edges(:,2) == axList(i);
         isSeed = edges(:,1) == showSynPostList(p);
         synPos = synAnchors(isAx & isSeed,:);
         synCol = PreConToConTo(i,:) ;
         synCol(synCol>1) = 1;
        % scatter(synPos(:,2),synPos(:,1),50,synCol,'filled','Marker',seedShapes{p})
         hold on
         for s = 1:size(synPos,1)
             axCol = synCol;
             plotY = [circAxAnchor(i,1) synPos(s,1)];
             plotX = [circAxAnchor(i,2) synPos(s,2)];
             lineWidth = 1;
             plot(plotX, plotY,'LineWidth',lineWidth,'color',axCol);
             
         end
         
     end
 end
 
 
 
 ylim(plotLimits(1,:));
xlim(plotLimits(2,:)); 




