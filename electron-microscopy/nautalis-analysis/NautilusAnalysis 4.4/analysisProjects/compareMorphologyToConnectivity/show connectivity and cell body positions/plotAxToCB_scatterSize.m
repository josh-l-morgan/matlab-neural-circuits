
%% Display connections between axons and thalamocortical cells using the cell body positions of thalamocortical cells


clear all
load('MPN.mat')
load([MPN 'obI.mat'])
shouldPause = 0;

conGraphDir = [MPN 'conGraph\'];
if ~exist(conGraphDir,'dir'),mkdir(conGraphDir),end

h = 100;
dimOrder = [2 3 1];
dSamp = [-.0184 0.016 0.030];

plotLimits = [-500 -200;-50 450; 0 400];
plotLimits = [-600 -0;0 600; 0 600];

set(gcf,'Position',[50 50 1000 1000])
%% variables
seedList = [108 201]
seedChannel = [1 2 3];

crossoverAxons = [2032	2033	2034	2035]

%% Get cells
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

useList = obI2cellList_seedInput_RGC_TCR(obI,seedList)
axList = useList.preList;
cellList = useList.postList;
con = useList.con;

cellPref = seedPreferences(seedList,useList)


%% Pick colors
cmap = [0 0 0; 1 0 0 ; 0 1 0; 0 0 1];

preColIdx = size(cmap,1)+1:size(cmap,1)+length(axList);
cellColIdx = max(preColIdx)+1:max(preColIdx)+length(cellList);

%%Axon color
PreConToSeed = zeros(size(con,1),3);
for i = 1:length(seedList)
    PreConToSeed(:,seedChannel(i)) = con(:,find(cellList == seedList(i)));
end

PreCol(:,1) = PreConToSeed(:,1)/mean(PreConToSeed(PreConToSeed(:,1)>0,1))*.25;
PreCol(:,2) = PreConToSeed(:,2)/mean(PreConToSeed(PreConToSeed(:,2)>0,2))*.25;
% PreCol(:,3) = 1-(PreCol(:,1)-PreCol(:,2));
% PreCol(PreCol<0) = 0;
PreCol(:,3) = .5;
PreCol(isnan(PreCol)) = 0;

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

rawAnchors = zeros(length(cellList),3);
for i = 1:length(cellList)
    rawAnchors(i,:) = obI.cell.anchors(find(obI.cell.name == cellList(i)),:);
end
anchors(:,1) = rawAnchors(:,dimOrder(1)).*dSamp(dimOrder(1));
anchors(:,2) = rawAnchors(:,dimOrder(2)).*dSamp(dimOrder(2));
anchors(:,3) = rawAnchors(:,dimOrder(3)).*dSamp(dimOrder(3));


%anchors = round(anchors);
%anchors(anchors<1) = 1;

%synapses = cat(1,[ 0 0 0], synapses);
rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors(:,1) = rawSynAnchors(:,dimOrder(1)).*dSamp(dimOrder(1));
synAnchors(:,2) = rawSynAnchors(:,dimOrder(2)).*dSamp(dimOrder(2));
synAnchors(:,3) = rawSynAnchors(:,dimOrder(3)).*dSamp(dimOrder(3));


%% axon anchors
%%Get axon information
hold on
seedColList = [1 0 0; 0 1 0];

for i = 1:length(axList)
    allAnchors =  synAnchors(edges(:,2) == axList(i),:);
    axAnchor(i,:) = mean(allAnchors,1);
end

%%Axon anchor Circle

meanAnchor = mean(anchors,1);
for i = 1:length(axList)
    oldAnch = axAnchor(i,:);
    slope = meanAnchor - oldAnch ;
    x = slope(2);
    y = slope(1);
    c = sqrt((x^2 +y^2)/(h^2));
    circAxAnchor(i,:) = ...
        [meanAnchor(1) - y/c meanAnchor(2) - x/c meanAnchor(3)];
    scatter(y,x,'r')
    scatter(y/c,x/c,'g')
    sqrt((y/c)^2+(x/c)^2);
end

%%Evenly spaced anchors
targetAngles = [1:length(axList)]/(length(axList)) * 2 * pi;
actualAngles = atan2(axAnchor(:,1)-meanAnchor(1),axAnchor(:,2)-meanAnchor(2));
[sortAngles aix] = sort(actualAngles,'ascend');
newAngles(aix) = targetAngles-pi;
circAxAnchor= [sin(newAngles)' * h + meanAnchor(1) ...
    cos(newAngles)' * h + meanAnchor(2) ...
    axAnchor(:,3)];
scatter(circAxAnchor(:,2),circAxAnchor(:,1))



circAxAnchor = axAnchor;





%% Scatter pie Sizes

if shouldPause
    pause
    
end
clf

%
%
% scatter(cbPos(:,2),cbPos(:,1),30,'filled','b')
% scatter(anchors(:,2),anchors(:,1),100,cellColIdx+1,'filled','o',...
%     'MarkerEdgeColor','w','LineWidth',1)


seedLink = cellPref.sharedSyn;
maxNum = 1;
sumLinks = sum(seedLink,1);
[sortLinks idx] = sort(sumLinks,'descend');

for s = 1:size(anchors,1)
    
    i = idx(s);
    clear p
    vals = seedLink(1:2,i);
    totVal = sum(vals);
    vals = vals/sum(vals);
%     if sum(vals)>= maxNum
%         vals = [vals * maxNum/sum(vals); 0];
%     else
%         vals = [vals; maxNum - sum(vals)];
%     end
    p.vals = vals;
    %p.vals = [10 3 1];
    p.color = [1 0 0; 0 .2 1; 0 0 0];
    p.center = anchors(i,1:2);
    p.radius = sqrt(totVal/pi);
    myPie(p)
    hold on
    
end
hold off
set(gca,'color',[0 0 0])
set(gcf,'color',[1 1 1])
xlim([0 600])
ylim([0 600])


myFileName = 'D:\LGNs1\Analysis\connectivityGraphs\scatterPies\pieSize1.eps'
% saveas(gcf,'myFileName','epsc2');
% print('-r150','-depsc2', myFileName);
% plot2svg(myFileName,gcf,'png')
print(gcf, myFileName, '-depsc2','-painters')


'return'
return



%% Scatter pie

if shouldPause
    pause
    
end
clf

%
%
% scatter(cbPos(:,2),cbPos(:,1),30,'filled','b')
% scatter(anchors(:,2),anchors(:,1),100,cellColIdx+1,'filled','o',...
%     'MarkerEdgeColor','w','LineWidth',1)


seedLink = cellPref.sharedSynNorm;
maxNum = 1;
for i = 1:size(anchors,1)
    
    clear p
    vals = seedLink(1:2,i);
    if sum(vals)>= maxNum
        vals = [vals * maxNum/sum(vals); 0];
    else
        vals = [vals; maxNum - sum(vals)];
    end
    p.vals = vals;
    %p.vals = [10 3 1];
    p.color = [1 0 0; 0 .5 1; 0 0 0];
    p.center = anchors(i,1:2);
    p.radius = 5;
    myPie(p)
    hold on
    
end
hold off
set(gca,'color',[.1 .1 .1])
set(gcf,'color',[1 1 1])

myFileName = 'D:\LGNs1\Analysis\connectivityGraphs\scatterPies\pie4.eps'
% saveas(gcf,'myFileName','epsc2');
% print('-r150','-depsc2', myFileName);
% plot2svg(myFileName,gcf,'png')
print(gcf, myFileName, '-depsc2','-painters')


'return'
return



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


ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal

myVectorFile = [conGraphDir 'draftLines.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')


%%

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
    scatter(cbPos(i,2),cbPos(i,1),250,seedCol{seedChannel(i)},'filled','o',...
        'MarkerEdgeColor','w','LineWidth',1)
end

for i = 1:length(highlightCell)
    text(cbPos(i,2)+1,cbPos(i,1),num2str(highlightCell(i)),...
        'FontSize',6,'HorizontalAlignment','center','FontWeight','bold')
end

ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal

myVectorFile = [MPN 'conGraph\draftcirc.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')

if shouldPause
    pause
    clf
end
%%

if shouldPause
    pause
    clf
    
end
set(gca,'color','k')
%colormap(PreCol)
%scatter(circAxAnchor(:,2),circAxAnchor(:,1),[1:length(axList)]',40, 'filled')
% scatter(circAxAnchor(:,2),circAxAnchor(:,1),50,cmap(preColIdx,:),'filled','d',...
%     'MarkerEdgeColor','w','LineWidth',1)
hold on
for i = 1:length(axList)
    %     plot([circAxAnchor(i,2) circAxAnchor(i,2)],[circAxAnchor(i,1) circAxAnchor(i,1)],'color',cmap(preColIdx(i),:),'Marker','d',...
    %         'MarkerFaceColor',cmap(preColIdx(i),:),'MarkerEdgeColor','w','MarkerSize',8)
    %
    scatter(circAxAnchor(i,2),circAxAnchor(i,1),50,cmap(preColIdx(i),:),'filled','d',...
        'MarkerEdgeColor','w','LineWidth',1)
    
end

ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal


if shouldPause
    pause
    clf
end

hold on
%%clf



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
    scatter(cbPos(i,2),cbPos(i,1),250,seedCol{seedChannel(i)},'filled','o',...
        'MarkerEdgeColor','w','LineWidth',1)
end

for i = 1:length(highlightCell)
    text(cbPos(i,2)+1,cbPos(i,1),num2str(highlightCell(i)),...
        'FontSize',6,'HorizontalAlignment','center','FontWeight','bold')
end

ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal

myVectorFile = [MPN 'conGraph\draftcirc.eps'];
print(gcf, myVectorFile, '-depsc2','-painters')

if shouldPause
    pause
    clf
end

%% synDistribution
showSynPostList = [seedList(1)]
%plotLimits = [-450 -150;-50 450; 0 400];

preToShowCell = unique(edges(edges(:,1)== showSynPostList,2))
[preToShowCell ia ib] = intersect(axList,preToShowCell);


clf
PreConToConTo = zeros(length(axList),3);
for i = 1:size(con,1)
    for p = 1:size(con,2)
        PreConToConTo(i,:) =  PreConToConTo(i,:) + sqrt(con(i,p) * PostConToSeed(p,:));
    end
end


meanAnchor = anchors(find(cellList == showSynPostList),:);
subAxAnchors = axAnchor(ia,:);
targetAngles = [1:length(ia)]/(length(ia)) * 2 * pi;
actualAngles = atan2(subAxAnchors(:,1)-meanAnchor(1),subAxAnchors(:,2)-meanAnchor(2));
[sortAngles aix] = sort(actualAngles,'ascend');
newAngles = [];
newAngles(aix) = targetAngles-pi;
circSubAxAnchor= [sin(newAngles)' * h + meanAnchor(1) ...
    cos(newAngles)' * h + meanAnchor(2) ...
    subAxAnchors(:,3)];
scatter(circSubAxAnchor(:,2),circSubAxAnchor(:,1))


barHeight = 30;
scatterBar(PreConToConTo(ia,1:2),circSubAxAnchor,barHeight)


ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal


%%
clf



% isPost = find(edges(:,1) == showSynPostList)
% synPos = synAnchors(isPost,:);
% scatter(synPos(:,2),synPos(:,1),30,'k','filled')
%
hold on
seedShapes = {'o','s'};
for p = 1:length(showSynPostList)
    for i = 1:length(axList)
        isAx = edges(:,2) == axList(i);
        isSeed = edges(:,1) == showSynPostList(p);
        synPos = synAnchors(isAx & isSeed,:);
        synCol = PreConToConTo(i,:) ;
        synCol = synCol/sum(synCol)*10;
        synCol(3) = 0
        synCol(synCol>1) = 1;
        %scatter(synPos(:,2),synPos(:,1),10,'k','filled','Marker',seedShapes{p},'MarkerEdge',synCol)
        scatter(synPos(:,2),synPos(:,1),10,synCol,'filled','Marker',seedShapes{p})
        
        
    end
end



ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal


%%
clf

for p = 1:length(showSynPostList)
    for i = 1:length(axList)
        isAx = edges(:,2) == axList(i);
        isSeed = edges(:,1) == showSynPostList(p);
        synPos = synAnchors(isAx & isSeed,:);
        synCol = PreConToConTo(i,:);
        synCol = synCol/sum(synCol)*10;
        synCol(3) = 0;
        synCol(synCol>1) = 1;
        % scatter(synPos(:,2),synPos(:,1),50,synCol,'filled','Marker',seedShapes{p})
        hold on
        for s = 1:size(synPos,1)
            axCol = synCol;
            targ = find(preToShowCell ==axList(i))
            plotY = [circSubAxAnchor(targ,1) synPos(s,1)];
            plotX = [circSubAxAnchor(targ,2) synPos(s,2)];
            lineWidth = 1;
            plot(plotX, plotY,'LineWidth',lineWidth,'color',axCol);
            
        end
        
    end
end



ylim(plotLimits(dimOrder(1),:));
xlim(plotLimits(dimOrder(2),:));
axis equal



