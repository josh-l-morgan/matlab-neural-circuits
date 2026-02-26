function[results] = runSkelSprings(springDat,results,movDir)

subplot(1,1,1)
showSide = 1;

group = springDat.group;
groupNum = springDat.groupNum;
param = springDat.param;
nodes = springDat.nodes;
edges = springDat.edges;
seed = springDat.seed;

seedList = seed.list;
try synPref = seed.seedPref; end;
try isPref = seed.isPref; end;
try snap = param.snap;end


figWin = [700 10 1200 1000];
outerWindow  = get(gca,'OuterPosition');
axisWin  = [0.02 0.02 0.96 0.96];

set(gcf,'Position',figWin)
set(gca, 'Position', axisWin)


%% set defaults

defaults = {'k' '1'; 'dampen' '.5'; 'noise' '[10 1 .1 0]' ; 'disperse' '5'; 'reps' '2000';...
    'fsize' '201'; 'speedMax' '1'; 'centerSpring' '2'; 'imageFreq' '100';...
    'startRepulse' 'reps/2';  'repulse' '5'; 'minDist' '5'; 'zeroSeeds' '1';'doSnap'  '0';'adaptDisperse' '1';...
    'bkGround','''k''';'edgeColor','''w''';'bgVec','[0 0 0]';'nodeDim','[2 1]'};

for i = 1:size(defaults,1)
    if isfield(param,defaults{i,1})
        eval(sprintf('%s = param.%s;',defaults{i,1},defaults{i,1}));
    else
        eval(sprintf('%s = %s;',defaults{i,1},defaults{i,2}));
    end
end

minSpring  = centerSpring;



%% Show 


            pl = plot([nodeX(e2(useEdge)) nodeX(e1(useEdge))]',[nodeY(e2(useEdge)) nodeY(e1(useEdge))]');
            hold on
            
            for p = 1:length(pl)
                set(pl(p),'color',edgeCol(useEdge(p),:),'LineWidth',edges.width(useEdge(p)));
            end
        end
        
        
        useEdge = find(~zeroEdge & edges.display);
        pl = plot([nodeX(e2(useEdge)) nodeX(e1(useEdge))]',[nodeY(e2(useEdge)) nodeY(e1(useEdge))]');
        hold on
        
        axis square
        axis off
        set(gcf,'color',bkGround)
        
        
        for p = 1:length(pl)
            set(pl(p),'color',edgeCol(useEdge(p),:),'LineWidth',edges.width(useEdge(p)));
        end
        
        
        










return


%% set fixed
wind = [0 fsize; 0 fsize];
center = [fsize/2 fsize/2];

%% configure nodes

%%Get node data
%     postIDs = useList.postList;
%     preIDs = useList.preList;
nodeIDs = nodes.labelIDs;
nodeNum = length(nodeIDs);
nodeType = nodes.type;
nodeCol = nodes.color;
nodeMarker = nodes.marker;
pushes = nodes.pushes;
isPushed = nodes.isPushed;
mass = nodes.mass;

if isfield(nodes,'nodeEdgeColor')
    nodeEdgeColor = nodes.nodeEdgeColor;
else
    nodeEdgeColor = nodeCol*0+1;
end


usePrev = 0;
if exist('results','var')
    if ~isempty(results)
        usePrev = 1;
    end
end




if usePrev
    move = 0;
    for i = 1:length(nodeIDs)
        
        targ = find(results.nodeIDs == nodeIDs(i),1);
        if ~isempty(targ)
            nodeX(i,1) = results.nodeX(targ);
            nodeY(i,1) = results.nodeY(targ);
        else
            nodeX(i,1) = fsize/2;
            nodeY(i,1) = fsize/2;
        end
        
    end
    
elseif param.usePos
    buff = 5;
    move = 1;
    %nodeDim = [2 1]
    nodeX = zeros(nodeNum,1)+fsize/2;
    nodeY = zeros(nodeNum,1)+fsize/2;
    nodePos = nodes.nodePos(nodes.usePos,:);
    nodePos = nodePos(:,nodeDim);
    nodePos(:,1) = nodePos(:,1) - min(nodePos(:,1));
    nodePos(:,2) = nodePos(:,2) - min(nodePos(:,2));
    maxPos = max(nodePos,[],1);
    
    fieldRat = (fsize-buff*2)./maxPos;
    if max(fieldRat)>0
        fRat = min(fieldRat);
    else
        fRat = max(fieldRat);
    end
    
    
    
    nodeX(nodes.usePos>0) = nodePos(:,1)*fRat+buff;
    nodeY(nodes.usePos>0) = fsize - nodePos(:,2)*fRat-buff;
    
    
    
else
    move = 1;
    %%initialize random positions, zero velocities
    %     nodeX = rand(nodeNum,1)*fsize/2+fsize/4;
    %     nodeY = rand(nodeNum,1)*fsize/2+fsize/4;
    nodeX = zeros(nodeNum,1)+fsize/2;
    nodeY = zeros(nodeNum,1)+fsize/2;
    %     nodeX = rand(nodeNum,1)*fsize/10000+fsize/2;
    %     nodeY = rand(nodeNum,1)*fsize/10000+fsize/2;
    
    
end




nodeXV = zeros(nodeNum,1);
nodeYV = zeros(nodeNum,1);

%%Make connectivity matrix
allEdges = edges.all;
if isfield(edges,'ew')
    
    e1 = edges.e1;
    e2 = edges.e2;
    ew = edges.ew;
    symmetry = edges.symmetry;
    edgeCol = edges.col;
    try edgeWidth = edges.width; end;
    
else
    
    %%Get new edges
    [e1 e2] = find((con.*useCon)>0);
    ei = sub2ind(size(con),e1,e2);
    ew = con(ei);
    edgeCol = nodeCol(e1,:) + nodeCol(e2,:);
    edgeCol(ew==0,:) = edgeCol(ew==0,:)*.2;
    %edgeWidth = ew*0+1;
end


tic

rawCon = zeros(nodeNum,nodeNum);
useCon = rawCon * 0;
for e = 1:length(e1)
    if symmetry(e)
        rawCon(e1(e),e2(e))  = rawCon(e2(e),e1(e)) + ew(e) * .5;
        rawCon(e2(e),e1(e)) = rawCon(e2(e),e1(e)) + ew(e) * .5;
    else
        rawCon(e1(e),e2(e)) = rawCon(e1(e),e2(e)) + ew(e);
    end
end
con  = rawCon;
con = con';
if showSide
    subplot(2,3,1)
    colormap jet(100)
    image(con*20),pause(.01)
end
'con'
toc


    if showSide
        subplot(2,3,4)
        %     showFieldG = (imresize(fieldG,1/fres));
        fieldGsc = fieldG * 256/ max(fieldG(:));
        showFieldG = uint8(cat(3,fieldGsc,fieldDis*1000,fieldGsc));
        image(flipud(showFieldG))
        pause(.1)
    end
    %%Display
    if length(imageFreq) == 1;
        showEdge = ~mod(r-1,imageFreq);
    else
        showEdge = sum(imageFreq==r)>0;
    end
    
    if showEdge | (r == reps) | (r == 1)
        'phys'
        toc
        tic
        hold off
        disp(sprintf('running %d of %d',r,reps))
        if showSide
            subplot(2,3,[2 3 5 6])
        end
        axis square
        axis off
        set(gcf,'color',bkGround)
        
        
      useEdge = find(zeroEdge & edges.display);
        if sum(useEdge)

            pl = plot([nodeX(e2(useEdge)) nodeX(e1(useEdge))]',[nodeY(e2(useEdge)) nodeY(e1(useEdge))]');
            hold on
            
            for p = 1:length(pl)
                set(pl(p),'color',edgeCol(useEdge(p),:),'LineWidth',edges.width(useEdge(p)));
            end
        end
        
        
        useEdge = find(~zeroEdge & edges.display);
        pl = plot([nodeX(e2(useEdge)) nodeX(e1(useEdge))]',[nodeY(e2(useEdge)) nodeY(e1(useEdge))]');
        hold on
        
        axis square
        axis off
        set(gcf,'color',bkGround)
        
        
        for p = 1:length(pl)
            set(pl(p),'color',edgeCol(useEdge(p),:),'LineWidth',edges.width(useEdge(p)));
        end
        
        
        
        
        if edges.text
            for p = 1:length(pl)
                text(mean([nodeX(e2(p))  nodeX(e1(p))]+1),mean([nodeY(e2(p)) nodeY(e1(p))]+1),num2str(ew(p)),...
                    'color',edgeColor,'FontSize',8)
            end
        end
        %colormap(cat(1,[0 0 0],nodeCol))
        colormap(nodeCol)
        for g = 1:groupNum
            if length(group(g).ind) == 3;
                group(g).ind = [group(g).ind group(g).ind];
            end
            
            ind = group(g).ind;
            for s = 1:length(group(g).ind)
                sg = scatter(nodeX(ind(s)),nodeY(ind(s)),group(g).size,nodeCol(ind(s),:),nodeMarker{ind(s)},'filled');
                set(sg,'LineWidth',group(g).lineWidth);
            end
            
            % %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
            %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,group(g).ind,group(g).marker,'filled');

            %%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,group(g).ind,group(g).marker,'filled');
            
            %%set(sg,'MarkerFaceColor',[group(g).ind 3]')
            %%set(sg,'MarkerEdgeColor',edgeColor);
            %set(sg,'LineWidth',group(g).lineWidth);
            %%set(sg,'SizeData',c;
            
            if group(g).text;
                for t = 1:length(group(g).ind)
                    targ = group(g).ind(t);
                    tx = text(nodeX(targ)-0,nodeY(targ)-0,{num2str(nodeIDs(targ))});
                    %tx = text(nodeX,nodeY,num2str(nodeIDs))
                    set(tx,'Color','w')
                end
            end
            
        end
        
        if zeroSeeds
            for s = 1:length(seedNodes);
                targ = seedNodes(s);
                tx = text(nodeX(targ)-4,nodeY(targ)-4,{num2str(nodes.nodeName(targ))});
                set(tx,'Color','w')
            end
        end
        
        ylim([wind(1,1) wind(1,2)])
        xlim([wind(2,1) wind(2,2)])
        set(gca,'color','k')
        hold off
        pause(.1)
        
        if exist('movDir','var')
            saveSpring(movDir)
        end
        'disp'
        toc
        tic
    end
    
    if ~move,break,end
end

results.nodeX = nodeX;
results.nodeY = nodeY;
results.nodeIDs = nodeIDs;
results.snapTime = snapTime;
results.cellGroups = cellGroups;
if exist('movDir','var')
    saveSpring(movDir)
end
