function[results] = runSprings(springDat)



     group = springDat.group;
     groupNum = springDat.groupNum;
     param = springDat.param;
     nodes = springDat.nodes;
     edges = springDat.edges;
     seed = springDat.seed;
     
     seedList = seed.list;
	try synPref = seed.seedPref, end;
    try isPref = seed.isPref, end;
        
        
    %% set defaults
    
    defaults = {'k' '1'; 'dampen' '.5'; 'noise' '[10 1 .1 0]' ; 'disperse' '5'; 'reps' '2000';...
        'fsize' '200'; 'speedMax' '1'; 'centerSpring' '2'; 'imageFreq' '100';...
        'startRepulse' 'reps/2';  'repulse' '5'; 'minDist' '5'; 'zeroSeeds' '1';'flatForce' '10'} ;
    
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
    nodeType = nodes.type;
    nodeCol = nodes.color;
    
    %%initialize random positions, zero velocities
    nodeX = rand(nodeNum,1)*fsize/2+fsize/4;
    nodeY = rand(nodeNum,1)*fsize/2+fsize/4;
    
    nodeX = rand(nodeNum,1)*fsize/10000+fsize/2;
    nodeY = rand(nodeNum,1)*fsize/10000+fsize/2;
    
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
        try edgeWidth = edges.width, end;
        
    else
               
        %%Get new edges
        [e1 e2] = find((con.*useCon)>0);
        ei = sub2ind(size(con),e1,e2);
        ew = con(ei);
        edgeCol = nodeCol(e1,:) + nodeCol(e2,:);
        edgeCol(ew==0,:) = edgeCol(ew==0,:)*.2;
        edgeCol(edgeCol<0) = 0;
        edgeCol(edgeCol>1) = 1;
        %edgeWidth = ew*0+1;
        
    end
    
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
        
    %%Zero out seeds
    if zeroSeeds
        for s = 1:length(seedList);
            targ = find(nodeIDs == seedList(s))
            prefDifs(s,:) = con(targ,:)+con(:,targ)';

            con(targ,:) = 0;
            con(:,targ) = 0;
            edgeCol(e1==targ,:) = edgeCol(e1==targ,:)  * .1;
            edgeCol(e2 == targ,:) = edgeCol(e2 == targ,:) * .1;
        
        end
        
    end
    
    %%Bad use
    useCon = con*0;
    [uy ux] = find(useCon+1);
    ui = sub2ind(size(useCon),uy,ux);
    useCon(ui(ux>uy)) = 1;
    
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
        netRad = fsize/2;
        centerForce(centerDist>(netRad)) =  centerForce(centerDist>(netRad))*10000;
        gravX = cos(radsCent) .* centerForce;
        gravY = sin(radsCent) .* centerForce;
        
        %%Calculate FlatForce
        useFlat = flatForce(ceil(r/reps*length(flatForce))); %% select noise epoch
        centDifX = nodeX * 0;
        centDifY = center(2) - nodeY;
        radsCent = atan2(centDifY,centDifX);
        centerDist = sqrt(centDifX.^2 + centDifY.^2);
        centerForce = centerDist/fsize .* useFlat;
        flatX = cos(radsCent) .* centerForce;
        flatY = sin(radsCent) .* centerForce;
        
        %%Calculate repulsion
        %repulsion = (1./dists .* repulse) + disperse/numNodes;
        %repulsion = (dists<minDist) * repulse + disperse/nodeNum;
        
        useDisperse = disperse(ceil(r/reps*length(disperse)));
        if r>=startRepulse
            repulsion = (dists<minDist) * repulse + (fsize./dists) * useDisperse/nodeNum;
        else
            repulsion = (fsize./dists) * useDisperse/nodeNum;
        end
        
        %repulsion(repulsion>repulse) = repulse;
        repulsion(isnan(repulsion)) = 0;
        repulsion(~useCon) = 0;
        pushX = cos(rads) .* repulsion;
        pushY = sin(rads) .*repulsion;
        
        sumPushX = -sum(pushX,2) + sum(pushX,1)';
        sumPushY = -sum(pushY,2) + sum(pushY,1)';
        
       
        
        %%Add forces
        nodeXV = (nodeXV + sumPullX + sumPushX + gravX); %% add forces
        nodeYV = (nodeYV + sumPullY+ sumPushY + gravY + flatY);
                
        %%Regulate speed
        nodeXV = nodeXV * dampen;
        nodeYV = nodeYV * dampen;
        speed = sqrt(nodeXV.^2 + nodeYV.^2);
        speedLim = speed;
        speedLim(speed>speedMax) = speedMax;
        speedRat  = speedLim./speed;
        nodeXV = nodeXV.* speedRat;
        nodeYV = nodeYV.* speedRat;
        
        %%Add noise
        useNoise = noise(ceil(r/reps*length(noise))); %% select noise epoch
        nodeXV = nodeXV + rand(nodeNum,1) * useNoise - useNoise/2;  %% add noise
        nodeYV = nodeYV  + rand(nodeNum,1) * useNoise - useNoise/2;
        
        %%Change position
        nodeX = nodeX + nodeXV;
        nodeY = nodeY + nodeYV;
        
        %%Position Seeds
        if zeroSeeds
            for s = 1:length(seedList);
                targ = find(nodeIDs == seedList(s));
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
        
        %%Display
        showEdge = ~mod(r-1,imageFreq);
        if showEdge | (r == reps)
           subplot(2,1,1)

            
            hold off
            disp(sprintf('running %d of %d',r,reps))
            
            pl = plot([nodeX(e2) nodeX(e1)]',[nodeY(e2) nodeY(e1)]');
            hold on
            
            axis square
            axis off
            set(gcf,'color','k')
            
            for p = 1:length(pl)
                set(pl(p),'color',edgeCol(p,:),'LineWidth',edges.width(p));
            end
            
            %colormap(cat(1,[0 0 0],nodeCol))
            %colormap(nodeCol)
            for g = 1:groupNum
                %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
                sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
                set(sg,'MarkerEdgeColor','w');
                set(sg,'LineWidth',group(g).lineWidth);
                %set(sg,'SizeData',c;
            end
            
            ylim([wind(1,1) wind(1,2)])
            xlim([wind(2,1) wind(2,2)])
            set(gca,'color','k')
            
            hold off
            
            subplot(2,1,2)
            cmap = cat(1,[0 0 0], jet(30));
            colormap(cmap)
            [sort1 idx1 ] = sort(nodeX(group(1).ind));
            order1 = group(1).ind(idx1);
            [sort2 idx2 ] = sort(nodeX(group(2).ind));
            order2 = group(2).ind(idx2);
            
            subCon = con(order1,:);
            subCon = subCon(:,order2);
            
            subConCol = zeros(size(subCon,1),size(subCon,2),3);
            subCol = subCon * 0;
            for c = 1:3
                    subCol(:)= cmap(subCon+1,c);
                    subConCol(:,:,c) = subCol;
            end
            
            colSize = 5;
            preCol = nodeCol(order1,:);
            preCol = permute(preCol,[1 3 2]);
            subConCol = cat(2,repmat(preCol,[1 colSize]),subConCol);
            postCol = nodeCol(order2,:);
            postCol = permute(postCol,[3 1 2]);
            postCol = cat(2,zeros(1,colSize,3),postCol);
            subConCol = cat(1,repmat(postCol,[colSize 1]),subConCol);
            image(uint8(subConCol*256))
            pause(.1)

        end
        
    end
    
    results.nodeX = nodeX;
    results.nodeY = nodeY;
    

