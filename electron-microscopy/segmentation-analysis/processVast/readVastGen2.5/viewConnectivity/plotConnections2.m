

load([MPN 'obI.mat'])


DimOrder = [1 2 3];

clf
subplot(1,1,1)
%% plot edges
hold on

dSamp = [.016 0.016 0.030];

anchors = obI.colStruc.anchors;
anchors(:,1) = anchors(:,1).*dSamp(1);
anchors(:,2) = anchors(:,2).*dSamp(2);
anchors(:,3) = anchors(:,3).*dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;

synapses = obI.nameProps.edges;
%synapses = cat(1,[ 0 0 0], synapses);
synAnchors = obI.colStruc.anchors(synapses(:,3),:);
edges = synapses(:,1:2);

[con sortPre sortPost] = edge2con(edges);
sharedPre = zeros(length(sortPost));
for i = 1:size(con,1)
    for p = 1:size(con,1)
        if i~=p
            sharedPre(i,p) = sum(min(con(:,i),con(:,p)));
        end
    end
end

cellList = [conTo(1).tcrList conTo(2).tcrList]
seedList = [108 201]
for i = 1:length(seedList);
    i
    for p = 1:length(cellList)
%         
%         con1 = sortPost(i);
%         con2 = sortPost(p);
%         
%     con1 = find(obI.cell.names==synapses(i,1));
%     con2 = find(obI.cell.names==synapses(i,2));
%     

    w1 = find(sortPost == seedList(i));
    w2 = find(sortPost == cellList(p));
    conWeight = sharedPre(w1,w2);
    if conWeight>0
    
    con1 = find(obI.cell.names==seedList(i));
    con2 = find(obI.cell.names==cellList(p));

    plotY = [anchors(con1,1) anchors(con2,1)];
    plotX = [anchors(con1,2) anchors(con2,2)];
  
    plot(plotX, plotY,'LineWidth',ceil(sqrt(conWeight)));
    
    end
    
    end
     pause(.01)
end
hold off


%%
synapses;

hold on
%%f
cbPos = zeros(length(cellList),3);
for i = 1:length(cellList)
    targ = find(obI.cell.name == cellList(i));
    cbPos(i,:) = anchors(targ,:);
end
scatter(cbPos(:,2),cbPos(:,1),30,'filled')


cellList = [201 108 ];
cbPos = zeros(length(cellList),3);
for i = 1:length(cellList)
    targ = find(obI.cell.name == cellList(i));
    cbPos(i,:) = anchors(targ,:);
end
scatter(cbPos(:,2),cbPos(:,1), 100, 'filled','k')
hold off

