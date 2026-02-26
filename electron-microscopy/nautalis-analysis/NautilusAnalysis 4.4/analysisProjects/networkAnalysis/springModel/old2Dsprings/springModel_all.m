%function[] = springModel(uselist,modPar);



% %%
% F = -kX
% k = stiffness
% X = proportional to distance;

%% Get data
loadData = 1;
if loadData
    clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [108 201];
    useList = obI2cellList_seedInput(obI,seedList);
    %useList = obI2cellList_all(obI);

    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
end

springDir = 'D:\LGNs1\Analysis\springDat3\'
if ~exist(springDir,'dir'), mkdir(springDir), end



%% shape attributes

    
    postSynPref = seedPref.sharedSynNorm(1,:)./sum(seedPref.sharedSynNorm,1);
    preSynPref = seedPref.ax2seed(1,:)./sum(seedPref.ax2seed,1);
    synPref = [preSynPref postSynPref];
    %synPref = preSynPref;
    grouping = [preSynPref * 0 + 1   postSynPref * 0 + 2];
    isPref = ~isnan(synPref);
    
    nodeNum = length(isPref);
    nodeCol = zeros(nodeNum,3);
    nodeCol(isPref,1) = synPref(isPref);
    nodeCol(isPref,3) = 1-synPref(isPref);
    nodeCol(~isPref,2) = 1;
    %set(0,'DefaultAxesColorOrder',nodeCol)
    
  	seed.list = seedList;
    seed.seedPref = synPref;
    seed.isPref = isPref;

    edges.all = allEdges;
    
    nodes.labelIDs = [useList.preList useList.postList];
    nodes.grouping1 = zeros(nodeNum,1);
%     nodes.grouping1(nodeIDs<1000) = 1;
%     nodes.grouping1((nodeIDs>=1000) & (nodeIDs < 10000)) = 2;
%     nodes.grouping1(nodeIDs>=10000) = 3;
%     
    nodes.grouping1 = grouping;
    nodes.color = nodeCol;
    
    
    %%set param
    param.zeroSeeds = 1;
    param.k = 1;
    param.noise = [10 1 .1 0];
    param.damp = .5;
    param.disperse = 5;
    param.reps = 2000;
    param.fsize = 200;
    param.speedMax = 1;
    param.centerSpring = 2;
    param.imageFreq = 100;
    
    param.startRepulse = 500;%
    param.repulse = 5; %1
    param.minDist = 5;
   
     %%Make groups
    %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
    
    group(1).ind = nodes.grouping1 == 1;
    group(2).ind = nodes.grouping1 == 2;
    group(3).ind = nodes.grouping1 == 3;
    group(1).marker = '^';
    group(2).marker = 'o';
    group(3).marker = 's';
    group(1).size = 10;
    grooup(2).size = 40;
    group(3).size = 20;
    groupNum = length(group);   
    
    
    springDat.group = group;
    springDat.param = param;
    springDat.nodes = nodes;
    springDat.edges = edges;
    springDat.seed = seed;
    
    for rerun = 1:10 
    
     group = springDat.group;
     groupNum = length(group);
     param = springDat.param;
     nodes = springDat.nodes;
     edges = springDat.edges;
     seed = springDat.seed;
     
     synPref = seed.seedPref;
     isPref = seed.isPref;
        
        
    %% set defaults
    
    defaults = {'k' '1'; 'damp' '.5'; 'noise' '[10 1 .1 0]' ; 'disperse' '5'; 'reps' '2000';...
        'fsize' '200'; 'speedMax' '1'; 'centerSpring' '2'; 'imageFreq' '100';...
        'startRepulse' 'reps/2';  'repulse' '5'; 'minDist' '5'; 'zeroSeeds' '1'} ;
    
    for i = 1:size(defaults,1)
       if isfield(param,defaults{i,1})
           eval(sprintf('%s = param.%s;',defaults{i,1},defaults{i,1}));
       else
          eval(sprintf('%s = %s;',defaults{i,1},defaults{i,2}));
       end        
    end
        
    
    %% set fixed
    wind = [0 fsize; 0 fsize];
    center = [fsize/2 fsize/2];
    
    %% configure nodes
    
    %%Get node data
%     postIDs = useList.postList;
%     preIDs = useList.preList;
    nodeIDs = nodes.labelIDs;
    nodeNum = length(nodeIDs);
    nodeGroup = nodes.grouping1;
    nodeCol = nodes.color;
    
    %%initialize random positions, zero velocities
    nodeX = rand(nodeNum,1)*fsize/2+fsize/4;
    nodeY = rand(nodeNum,1)*fsize/2+fsize/4;
    nodeXV = zeros(nodeNum,1);
    nodeYV = zeros(nodeNum,1);
    
    %%Make connectivity matrix
    allEdges = edges.all;
    rawCon = zeros(nodeNum,nodeNum);
    useCon = rawCon * 0;
    for i = 1:length(nodeIDs)
        for p = 1:length(nodeIDs)
            rawCon(i,p) = rawCon(i,p) + sum((allEdges(:,1) == nodeIDs(i)) & (allEdges(:,2) == nodeIDs(p)));
            rawCon(i,p) = rawCon(i,p) + sum((allEdges(:,1) == nodeIDs(p)) & (allEdges(:,2) == nodeIDs(i)));
            useCon(i,p) = (p>i);
        end
    end
    con  = rawCon;
    
    %%Get new edges
    [e1 e2] = find((con.*useCon)>0);
    
    %%Zero out seeds
    if zeroSeeds
        for s = 1:length(seedList);
            targ = find(nodeIDs == seedList(s))
            con(targ,:) = 0;
            con(:,targ) = 0;
        end
    end
    
    %con = con/max(con(:));
    
    %% %%%
    for r = 1:reps
        
        %%find all relative positions
        preXmat = repmat(nodeX,[1,nodeNum]);
        postXmat = repmat(nodeX',[nodeNum,1]);
        preYmat = repmat(nodeY,[1,nodeNum]);
        postYmat = repmat(nodeY',[nodeNum,1]);
        
        difX = postXmat - preXmat;
        difY = postYmat - preYmat;
        
        dists = (sqrt(difX.^2 + difY.^2));
        
        rads = atan2(difY,difX);
        rads(isnan(rads)) = 0;
        
        
        %%Calculate spring force
        springForce =  dists .* con * k;
        
        springForce(~useCon) = 0;
        pullX = cos(rads) .* springForce;
        pullY = sin(rads) .* springForce;
        pullX(isnan(pullX)) = 0;
        pullY(isnan(pullY)) = 0;
        
        sumPullX = sum(pullX,2) - sum(pullX,1)';
        sumPullY = sum(pullY,2) - sum(pullY,1)';
        
        %%Calculate Gravity
        centDifX = center(1) - nodeX;
        centDifY = center(2) - nodeY;
        radsCent = atan2(centDifY,centDifX);
        centerDist = sqrt(centDifX.^2 + centDifY.^2);
        centerForce = centerDist/fsize .* centerSpring;
        centerForce(centerDist>(fsize/2)) =  centerForce(centerDist>(fsize/2))*100;
        gravX = cos(radsCent) .* centerForce;
        gravY = sin(radsCent) .* centerForce;
        
        
        
        %%Calculate repulsion
        %repulsion = (1./dists .* repulse) + disperse/numNodes;
        %repulsion = (dists<minDist) * repulse + disperse/nodeNum;
        
        if r>=startRepulse
            repulsion = (dists<minDist) * repulse + (fsize./dists) * disperse/nodeNum;
        else
            repulsion = (fsize./dists) * disperse/nodeNum;
        end
        
        %repulsion(repulsion>repulse) = repulse;
        repulsion(isnan(repulsion)) = 0;
        repulsion(~useCon) = 0;
        pushX = cos(rads) .* repulsion;
        pushY = sin(rads) .*repulsion;
        
        sumPushX = -sum(pushX,2) + sum(pushX,1)';
        sumPushY = -sum(pushY,2) + sum(pushY,1)';
        
        %%Add noise
        useNoise = noise(ceil(r/reps*length(noise))); %% select noise epoch
        nodeXV = nodeXV + rand * useNoise - useNoise/2;  %% add noise
        nodeYV = nodeYV  + rand * useNoise - useNoise/2;
        
        %%Add forces
        nodeXV = (nodeXV + sumPullX + sumPushX + gravX); %% add forces
        nodeYV = (nodeYV + sumPullY+ sumPushY + gravY);
                
        %%Regulate speed
        nodeXV = nodeXV * damp;
        nodeYV = nodeYV * damp;
        speed = sqrt(nodeXV.^2 + nodeYV.^2);
        speedLim = speed;
        speedLim(speed>speedMax) = speedMax;
        speedRat  = speedLim./speed;
        nodeXV = nodeXV.* speedRat;
        nodeYV = nodeYV.* speedRat;
        
        %%Change position
        nodeX = nodeX + nodeXV;
        nodeY = nodeY + nodeYV;
        
        %%Position Seeds
        if zeroSeeds
            for s = 1:length(seedList);
                targ = find(nodeIDs == seedList(s));
                if synPref(targ)>.5
                    prefDif = synPref;
                else
                    prefDif = 1- synPref;
                end
                meanX = sum(nodeX.*prefDif')/sum(prefDif);
                meanY = sum(nodeY.*prefDif')/sum(prefDif);
                rad  = atan2(meanY-center(1), meanX- center(2));
                newY = sin(rad) * (fsize/2)+center(1);
                newX = cos(rad) * (fsize/2) + center(2);
                
                nodeX(targ) = newX;
                nodeY(targ) = newY;
                
            end
        end
        
        %%Display
        showEdge = ~mod(r-1,imageFreq);
        if showEdge | (r == reps)
            hold off
            disp(sprintf('running %d of %d',r,reps))
            
            pl = plot([nodeX(e2) nodeX(e1)]',[nodeY(e2) nodeY(e1)]');
            hold on
            
            ei = sub2ind(size(con),e1,e2);
            ev = con(ei);
            edgeCol = nodeCol(e1,:);
            edgeCol(ev==0,:) = edgeCol(ev==0,:)*.2;
            for p = 1:length(pl)
                set(pl(p),'color',edgeCol(p,:));
            end
            
            
            for g = 1:groupNum
                sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
                set(sg,'MarkerEdgeColor','w')
            end
            
            ylim([wind(1,1) wind(1,2)])
            xlim([wind(2,1) wind(2,2)])
            set(gca,'color','k')
            
            hold off
            pause(.1)
        end
        
    end
    
    result.nodeX = nodeX;
    result.nodeY = nodeY;
    
    imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
    saveas(gcf,imageName,'png')
    %print(gcf, imageName, '-depsc2','-painters')
    set(gcf, 'InvertHardCopy', 'off');
    print(gcf, imageName, '-dpng')
    
    
    %datName
    
end
















