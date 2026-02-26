function[col] = showArborVox(arbor)


subs = arbor.nodes.pos;
vsubs = arbor.vox.subs;

for i = 1:3
    subs(:,i) = subs(:,i) - arbor.offset(i) + 1;
    vsubs(:,i) = vsubs(:,i) - arbor.offset(i) + 1;
end



%% Render 3D

renderCon(vsubs,arbor.vox.conMat)
hold on

%% draw edges

for b = 1:length(arbor.branches)
    
    %edges = arbor.branches(b).edges;
    %drawNum = 10;
    %steps = [0:drawNum]/drawNum;
    %stepE = zeros(length(steps),3);
    
    nodes = [arbor.branches(b).tip; arbor.branches(b).nodes; arbor.branches(b).base];
    showNode = subs(nodes,:);
    nodeRad = arbor.nodes.nodeRad(nodes);
    scatSize = (nodeRad*10).^2;
    scatSize  = nodeRad*0+1;
    scatter3(showNode(:,1),showNode(:,2),showNode(:,3),'b','.')
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








