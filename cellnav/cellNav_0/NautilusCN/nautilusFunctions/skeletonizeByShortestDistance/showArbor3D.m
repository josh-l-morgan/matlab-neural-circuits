function[col] = showArbor3D(arbor,nodeCol)


subs = arbor.nodes.pos;
vsubs = arbor.vox.subs;


if ~exist('nodeCol','var')
   nodeCol = zeros(length(arbor.nodes.rad),3);
   nodeCol(:,3) = 1;
end

for i = 1:3
    subs(:,i) = subs(:,i) - arbor.skel.offset(i) + 1;
    vsubs(:,i) = vsubs(:,i) - arbor.skel.offset(i) + 1;
end



%% Render 3D
hold off
clf

if 1
    renderCon(vsubs,arbor.vox.conMat,[1 1 1],.1);
else
    downSamp = 2;
    renderProps.smooth = 0;
    renderProps.resize = 1;
    renderProps.smoothPatch = 1;
    
    smallSub = shrinkSub(vsubs,downSamp);
    %smallSub = smallSub(:,flipDim);
    
    fv = subVolFV(smallSub,[],renderProps);
    %fv.vertices = fv.vertices * downSamp;
    fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;
    
    [p] = renderFV(fv,[1 1 1],.2);
    
end
view([0 0])
axis off
pause(.01)
hold on


%%draw edges

for b = 1:length(arbor.branches)
    
    edges = arbor.branches(b).edges;
    drawNum = 10;
    steps = [0:drawNum]/drawNum;
    stepE = zeros(length(steps),3);
    
    nodes = [arbor.branches(b).tip arbor.branches(b).nodes];
    showNode = subs(nodes,:);
    nodeRad = arbor.nodes.rad(nodes);
    scatSize = (nodeRad*10).^2;
    scatCol = nodeCol(nodes,:);
    scatter3(showNode(:,1),showNode(:,2),showNode(:,3),scatSize,scatCol,'.')
    runE = subs(nodes,:);
    plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',1,'color','w')
    
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


