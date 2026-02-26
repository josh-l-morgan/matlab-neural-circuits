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

%%Zero out seeds
isSeed = zeros(length(nodeIDs),1);
zeroEdge = zeros(length(e2),1)>0;
seedNodes = [];
if zeroSeeds
    for s = 1:length(seedList);
        targ = find((nodes.nodeName == seedList(s)) & (nodes.type == 1));
        if ~isempty(targ)
            isSeed(targ) = 1;
            seedNodes(length(seedNodes)+1) = targ;
            
            prefDifs(s,:) = con(targ,:)+con(:,targ)';
            
            con(targ,:) = 0;
            con(:,targ) = 0;
            isE1 = find(e1==targ);
            isE2 = find(e2 == targ);
            edgeCol(isE1,:) = edgeCol(isE1,:) * .20 + repmat(bgVec,[length( isE1) 1]) * .75;
            edgeCol(isE2,:) = edgeCol(isE2,:) * .20 + repmat(bgVec,[length( isE2) 1]) * .75;
            ew([isE1(:); isE2(:)]) = 0.001;
            zeroEdge(e2 == targ) = 1;
            
        end
        
    end
    
    
    zeroEdge = zeroEdge>0;
    
end

%%Bad use
useCon = con*0;
[uy ux] = find(useCon+1);
ui = sub2ind(size(useCon),uy,ux);
%useCon(ui(ux>uy)) = 1;
useCon(ui(~(ux==uy))) = 1;

%con = con/max(con(:));

%%Condition variables
edgeCol(edgeCol<0) = 0;
edgeCol(edgeCol>1) = 1;
nodeCol(nodeCol<0) = 0;
nodeCol(nodeCol>1) = 1;

%%%%%
symCon = con + con';
[cellGroup] = segmentCon(symCon);
cellGroups = cellGroup;
snapTime = 0;


[usePre usePost useW] = find(con>0);
usePre = e1';
usePost = e2';
useW = ew;

fres = 1;
fieldDis = zeros(fsize*fres);

tic
%kern = gaus3d([100 100 1],minDist);
%kern = kern/max(kern(:));

kern1 = tent((fsize*fres)*2+1)*(1/sum(pushes))* disperse;
kern2 = tent(minDist*2*fres)* repulse;
kern3 = addKern(kern1,kern2);
%image(uint8(repmat(kern,[1 1 3]))*1000)

for r = 1:reps
    %%find all relative positions
    %     preXmat = repmat(nodeX,[1,nodeNum]);
    %     postXmat = preXmat';%repmat(nodeX',[nodeNum,1]);
    %     preYmat = repmat(nodeY,[1,nodeNum]);
    %     postYmat = preYmat';%repmat(nodeY',[nodeNum,1]);
    %
    %%make push field
    nodeX(nodeX<1) = 1;
    nodeY(nodeY<1) = 1;
    nodeX(nodeX>fsize) = fsize;
    nodeY(nodeY>fsize) = fsize;
    nodeInds = sub2ind([fsize*fres fsize*fres],round(nodeY*fres),round(nodeX*fres));
    pushInds = nodeInds(pushes>0);
    uInds = unique(pushInds);
    hInds = histc(pushInds,uInds);
    fieldDis = fieldDis*0;
    fieldDis(uInds) = hInds;
    
    if r>=startRepulse
        %fieldG = fastCon(fieldDis,kern1);
        kern = kern2;
    else
        %fieldG = fastCon(fieldDis,kern1);
        kern = kern1;
    end
    fieldG = convn(fieldDis,kern,'same');
    
    
    %[Gx,Gy] = imgradientxy(fieldG,'intermediate');
    [Gx,Gy] = imgradientxy(fieldG,'sobel');
    
    nodeGx = Gx(nodeInds);%.* (fieldG(nodeInds)-max(kern(:)));
    nodeGy = Gy(nodeInds);%.* (fieldG(nodeInds)-max(kern(:)));
    
    sumPushX = nodeGx .* isPushed *  -1 ;%
    sumPushY = nodeGy .* isPushed *  -1;
    
    %%Calculate forces
    
    preXmat = nodeX(usePre);
    preYmat = nodeY(usePre);
    postXmat = nodeX(usePost);
    postYmat = nodeY(usePost);
    
    difX = postXmat - preXmat;
    difY = postYmat - preYmat;
    
    dists = (sqrt(difX.^2 + difY.^2));
    
    
    
    rads = atan2(difY,difX);
    rads(isnan(rads)) = 0;
    
    
    %%Calculate spring force
    springForce =  dists .* useW * k;
    
    %%Snap a spring
    if doSnap
        if ~sum(con(:)),break,end
        snapping = (r>snap.start) & ~mod(r,snap.wait);
        if snapping
            if snap.synAll
                con(con>0) = con(con>0) - 1;
                sum(con(:))
            else %if not cut all syn by 1
                
                
                
                isCon = find(con>0);
                conDist = dists(isCon);
                if snap.longest
                    longest = find(conDist == max(conDist));
                else
                    longest = ceil(rand*length(conDist));
                end
                snapit = isCon(longest);
                
                [y x] = ind2sub(size(con),snapit);
                snapEdge = find((e1 == y) & (e2 == x));
                edgeCol(snapEdge,:) = edgeCol(snapEdge,:)*.2;
                con(snapit) = 0;
                
            end %if snapAll
            
            %disp('snap')
            symCon = con + con';
            [cellGroup] = segmentCon(symCon);
            cellGroups = cat(1,cellGroups,cellGroup);
            snapTime = cat(1,snapTime,r);
            groupCount = max(cellGroup(:));
            if isfield(snap,'maxGroup')
                if groupCount >= snap.maxGroup
                    doSnap = 0;
                    disperse = snap.finalDisperse;
                end
            end
        end
    end
    
    
    %springForce(~useCon) = 0;
    springForce = springForce*.1;
    pullX = cos(rads) .* springForce;
    pullY = sin(rads) .* springForce;
    pullX(isnan(pullX)) = 0;
    pullY(isnan(pullY)) = 0;
    
    %     sumPullX = sum(pullX,2) - sum(pullX,1)';
    %     sumPullY = sum(pullY,2) - sum(pullY,1)';
    %
    
    sumPullX = full(sum(sparse(e1,e2,pullX,nodeNum,nodeNum),2)- ...
        sum(sparse(e2,e1,pullX,nodeNum,nodeNum),2)) ;
    sumPullY = full(sum(sparse(e1,e2,pullY,nodeNum,nodeNum),2) - ...
        sum(sparse(e2,e1,pullY,nodeNum,nodeNum),2));
    
    
    %%Calculate Gravity
    centDifX = center(1) - nodeX;
    centDifY = center(2) - nodeY;
    radsCent = atan2(centDifY,centDifX);
    centerDist = sqrt(centDifX.^2 + centDifY.^2);
    if adaptDisperse
        if mean(centerDist)>(fsize * .30),
            disperse = disperse - .01;
        elseif mean(centerDist)<(fsize * .2)
            disperse = disperse + .01;
        end
    elseif 0
        if mean(centerDist)>(fsize * .30),
            centerSpring = centerSpring + .001;
        elseif mean(centerDist)<(fsize * .2)
            centerSpring = centerSpring - .001;
        end
        centerSpring = max(minSpring,centerSpring);
        
    end
    centerForce = centerDist/fsize .* centerSpring;
    centerForce(centerDist>(fsize/2)) =  centerForce(centerDist>(fsize/2))*100000000;
    gravX = cos(radsCent) .* centerForce;
    gravY = sin(radsCent) .* centerForce;
    
    
    %%Add forces
    
    nodeXV = (nodeXV + sumPullX + sumPushX + gravX); %% add forces
    nodeYV = (nodeYV + sumPullY+ sumPushY + gravY);
    
    
    %%Regulate speed
    nodeXV = nodeXV * dampen;
    nodeYV = nodeYV * dampen;
    speed = sqrt(nodeXV.^2 + nodeYV.^2);
    speedLim = speed;
    speedLim(speed>speedMax) = speedMax;
    speedRat  = speedLim./speed;
    speedRat(isnan(speedRat)) = 0;
    nodeXV = nodeXV.* speedRat;
    nodeYV = nodeYV.* speedRat;
    
    %%Add noise
    useNoise = noise(ceil(r/reps*length(noise))); %% select noise epoch
    nodeXV = nodeXV + rand(nodeNum,1) * useNoise - useNoise/2;  %% add noise
    nodeYV = nodeYV  + rand(nodeNum,1) * useNoise - useNoise/2;
    
    if move
        %%Change position
        nodeX = nodeX + nodeXV .* nodes.nodeMove;
        nodeY = nodeY + nodeYV .* nodes.nodeMove;
        
        %%Position Seeds
        if zeroSeeds
            for s = 1:length(seedList);
                targ = find(nodes.nodeName == seedList(s));
                if ~isempty(targ)
                    prefDif = prefDifs(s,:);
                    meanX = sum(nodeX.*prefDif')/sum(prefDif);
                    meanY = sum(nodeY.*prefDif')/sum(prefDif);
                    rad  = atan2(meanY-center(1), meanX- center(2));
                    newY = sin(rad) * (fsize/2)+center(1);
                    newX = cos(rad) * (fsize/2) + center(2);
                    
                    nodeX(targ) = newX;
                    nodeY(targ) = newY;
                end
                
                
                
            end
        end
    end
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
