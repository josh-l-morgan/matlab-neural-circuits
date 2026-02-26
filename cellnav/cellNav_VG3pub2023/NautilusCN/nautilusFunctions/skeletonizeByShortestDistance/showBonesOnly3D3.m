function[] = showBones3D3(skel)



subs = skel.node2subs;
subs = subs-repmat(skel.offset,[size(subs,1) 1]);

hold on
for tS = 1:length(skel.bones)
    
    nList = [skel.bones(tS).tip skel.bones(tS).nodes skel.bones(tS).base];
    pList = skel.node2surf(nList);
    runE = subs(nList,:);
    scatter3(runE(:,1),runE(:,2),runE(:,3),1,'r','o','filled')
    plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',1,'color','k')
    %scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')
    
    view([45 45])
    pause(.01)
    
end


if isfield(skel,'bridges')
    for b = 1:size(skel.bridges,1)
        
        nList = skel.bridges(b,:);
        pList = skel.node2surf(nList);
        runE = subs(nList,:);
        plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','r')
        
    end
end





