function[col] = showArborChunks3D(arbor)


subs = arbor.nodes.node2subs;
vsubs = arbor.vox.subs;

for i = 1:3
    subs(:,i) = subs(:,i) - arbor.offset(i) + 1;
    vsubs(:,i) = vsubs(:,i) - arbor.offset(i) + 1;
end



%% Render 3D
hold off
clf
if 0
downSamp = 2;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 1;

smallSub = shrinkSub(vsubs,downSamp);
%smallSub = smallSub(:,flipDim);

fv = subVolFV(smallSub,[],renderProps);
%fv.vertices = fv.vertices * downSamp;
fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;

[p] = renderFV(fv,[1 1 1],.1);

end

view([0 0])
axis off
pause(.01)
hold on
       

%%draw edges


downSamp = 1;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;

nCol = hsv(1000);
for b = 1:length(arbor.branches)
    
    edges = arbor.branches(b).edges;
    drawNum = 10;
    steps = [0:drawNum]/drawNum;
    stepE = zeros(length(steps),3);
    
    nodes = [arbor.branches(b).tip; arbor.branches(b).nodes];
    showNode = subs(nodes,:);
    nodeRad = arbor.nodes.nodeRad(nodes);
    scatSize = (nodeRad*10).^2;
    scatter3(showNode(:,1),showNode(:,2),showNode(:,3),'b','.')
    runE = subs(nodes,:);
    plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',1,'color','w')
    
    for n = 1:length(nodes)
        
        nSubs = vsubs(arbor.vox.nearestNode == nodes(n),:);
        
        
        
        nColi = nCol(ceil(rand*size(nCol,1)),:);
        
        if 1
            renderCon(nSubs,[],nColi,.2)
        else
            
            smallSub = shrinkSub(nSubs,downSamp);
            fv = subVolFV(smallSub,[],renderProps);
            fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;
            [p] = renderFV(fv,nColi,.2);
        end
        
    end
    view([0 0])
    pause(.01)
end


%% Draw bridges

if isfield(arbor,'bridges')
    edges = arbor.bridges;
    for b = 1:size(edges,1)
        pList = edges(b,:);
        runE = subs(pList,:);
        plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','y')
    end
    
end

view([0 0])

ax = gca;     
ax.Clipping = 'off';
hold off


