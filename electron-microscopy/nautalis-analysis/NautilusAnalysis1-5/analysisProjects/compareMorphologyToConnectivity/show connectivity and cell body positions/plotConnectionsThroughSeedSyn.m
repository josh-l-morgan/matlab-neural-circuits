

load([MPN 'obI.mat'])


DimOrder = [1 2 3];

cellList = [conTo(1).tcrList conTo(2).tcrList]
seedList = [108 201]
axRad = 100;

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
synAnchors(:,1) = synAnchors(:,1).*dSamp(1);
synAnchors(:,2) = synAnchors(:,2).*dSamp(2);
synAnchors(:,3) = synAnchors(:,3).*dSamp(3);
edges = synapses(:,1:2);


seedColList = [1 0 0; 0 1 0];
for i = 1:length(seedList)

    isSeed = edges(:,1) == seedList(i);
    targ = find(obI.cell.names == seedList(i));
    synAnchor = synAnchors(isSeed,1:2)
        anchor = anchors(targ,1:2);
        anchor = mean(synAnchor,1);
recSeedAnchors(i,:) = anchor;

synNum = size(synAnchor,1);
    plotY = [synAnchor(:,1) repmat(anchor(1),[synNum 1])];
    plotX = [synAnchor(:,2) repmat(anchor(2),[synNum 1])];
    plot(plotX',plotY','Color',seedColList(i,:),'LineWidth',1)


end
centerPoint = mean(recSeedAnchors,1);
hold off

%% axon anchors

hold on
seedColList = [1 0 0; 0 1 0];
allPre = unique(edges(:,2));
PreConToSeed = zeros(length(allPre),length(seedList));
for i = 1:length(allPre)
    clear conToSeed
   axAnchors = [];
    for s = 1:length(seedList)
        isSeed = (edges(:,1) == seedList(s)) & (edges(:,2) == allPre(i));
        axAnchors = cat(1,axAnchors,synAnchors(isSeed,:));
        conToSeed(s) = sum(isSeed);
    end
    PreConToSeed(i,:) = conToSeed;
    if sum(conToSeed)
   allAnchors =  synAnchors(edges(:,2) == allPre(i),:);
   axAnchor(i,:) = mean(allAnchors,1);
   synNum = size(axAnchors,1);
    plotY = [axAnchors(:,1) repmat(axAnchor(i,1),[synNum 1])];
    plotX = [axAnchors(:,2) repmat(axAnchor(i,2),[synNum 1])];
    axCol = [0 0 1 ]
    axCol(1:length(seedList)) = conToSeed/sum(conToSeed);
    plot(plotX',plotY','Color',axCol,'LineWidth',1)

end
end




%%
[con sortPre sortPost] = edge2con(edges);
sharedPre = zeros(length(sortPost));
for i = 1:size(con,1)
    for p = 1:size(con,1)
        if i~=p
            sharedPre(i,p) = sum(min(con(:,i),con(:,p)));
        end
    end
end


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


cbPos = zeros(length(seedList),3);
for i = 1:length(seedList)
    targ = find(obI.cell.name == seedList(i));
    cbPos(i,:) = anchors(targ,:);
end
scatter(cbPos(:,2),cbPos(:,1), 100, 'filled','k')
hold off

