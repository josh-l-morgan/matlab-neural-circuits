
%% Plot axon to tcr
clear all
MPN = GetMyDir;
load([MPN 'obI.mat'])

%% variables
DimOrder = [1 2 3];
seedList = [108 201]

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
cmap = [0 0 0];

preColIdx = 2:length(axList)+1;
cellColIdx = max(preColIdx)+1:max(preColIdx)+length(cellList);

%%Axon color
PreConToSeed = zeros(size(con,1),3);
for i = 1:length(seedList)
    PreConToSeed(:,i) = con(:,find(cellList == seedList(i)));
end

PreCol(:,1) = PreConToSeed(:,1)/mean(PreConToSeed(PreConToSeed(:,1)>0,1))*.5;
PreCol(:,2) = PreConToSeed(:,2)/mean(PreConToSeed(PreConToSeed(:,2)>0,2))*.5;
PreCol(:,3) = 0;
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

dSamp = [.016 0.016 0.030];
anchors = zeros(length(cellList),3);
for i = 1:length(cellList)
    anchors(i,:) = obI.cell.anchors(find(obI.cell.names == cellList(i)),:);
end
anchors(:,1) = anchors(:,1).*dSamp(1)*-1;
anchors(:,2) = anchors(:,2).*dSamp(2);
anchors(:,3) = anchors(:,3).*dSamp(3);


%anchors = round(anchors);
%anchors(anchors<1) = 1;

%synapses = cat(1,[ 0 0 0], synapses);
synAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors(:,1) = synAnchors(:,1).*dSamp(1)*-1;
synAnchors(:,2) = synAnchors(:,2).*dSamp(2);
synAnchors(:,3) = synAnchors(:,3).*dSamp(3);



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
    % scatter(y,x)
    % scatter(y/c,x/c)
    %sqrt((y/c)^2+(x/c)^2)
end
circAxAnchor = axAnchor;





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
            
            axCol(3) = 1;
            axCol(axCol>1) = 1;
            plotY = [circAxAnchor(i,1) anchors(p,1)];
            plotX = [circAxAnchor(i,2) anchors(p,2)];
            plot(plotX, plotY,'LineWidth',ceil(sqrt(conWeight)),'color',axCol);
            plot(plotX, plotY,'LineWidth',1,'color',axCol);
            
        end
        
    end
    pause(.01)
end
% cellColor(:,1) = cellColor(:,1)/mean(cellColor(:,1))/4;
% cellColor(:,2) = cellColor(:,2)/mean(cellColor(:,2))/4;
% cellColor(:,3) = .5;
%


ylim([-400 -200]);
xlim([50 450]);
myVectorFile = [MPN 'conGraph\draftLines.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')

%%

clf

%colormap(PreCol)
%scatter(circAxAnchor(:,2),circAxAnchor(:,1),[1:length(axList)]',40, 'filled')
scatter(circAxAnchor(:,2),circAxAnchor(:,1),50,preColIdx,'filled','d',...
    'MarkerEdgeColor','w','LineWidth',1)


%%
hold on
%%f


%scatter(cbPos(:,2),cbPos(:,1),30,'filled','b')
scatter(anchors(:,2),anchors(:,1),100,cellColIdx+1,'filled','o',...
    'MarkerEdgeColor','w','LineWidth',1)
for i = 1:length(cellList)
    text(anchors(i,2)+1,anchors(i,1)+2,num2str(cellList(i)),...
        'FontSize',6,'HorizontalAlignment','center','FontWeight','bold')
end
% 
% highlightCell = [108 201]
% cbPos = zeros(length(highlightCell),3);
% seedCol = zeros(length(highlightCell),1);
% for i = 1:length(highlightCell)
%     targ = find(cellList == highlightCell(i));
%     cbPos(i,:) = anchors(targ,:);
%     seedColIdx(i) = cellColIdx(targ);
% end
% scatter(cbPos(:,2),cbPos(:,1),250,seedColIdx,'filled','o',...
%     'MarkerEdgeColor','w','LineWidth',1)
% for i = 1:length(highlightCell)
%     text(cbPos(i,2)+1,cbPos(i,1),num2str(highlightCell(i)),...
%         'FontSize',6,'HorizontalAlignment','center','FontWeight','bold')
% end

ylim([-400 -200]);
xlim([50 450]);

myVectorFile = [MPN 'conGraph\draftcirc.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')




