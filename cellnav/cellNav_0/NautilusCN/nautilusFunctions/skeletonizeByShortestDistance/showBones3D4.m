function[] = showBones3D4(subs,skel,nodeCol)


%%

if ~isempty(subs)
    
    if ~exist('projectDim','var')
        projectDim = 1;
    end
    
    dims = [1:3];
    dims = dims(dims ~= projectDim);
    
    offset = min(subs,[],1)-1;
    
    if ~exist('dims','var')
        dims = [1 2];
    end
    
    for i = 1:3
        subs(:,i) = subs(:,i) - offset(i) + 1;
    end
    
    %% Render 3D
    
    if 1
        
        renderCon(subs,[],[1 1 1],.1)
        
        
    else
        downSamp = 4;
        renderProps.smooth = 0;
        renderProps.resize = 1;
        renderProps.smoothPatch = 0;
        
        smallSub = shrinkSub(subs,downSamp);
        %smallSub = smallSub(:,flipDim);
        
        fv = subVolFV(smallSub,[],renderProps);
        %fv.vertices = fv.vertices * downSamp;
        fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;
        
        [p] = renderFV(fv,[1 1 1],.2);
        
    end
    view([0 0])
   
    
else
    
    offset = [0 0 0];
    
end

 axis off
    pause(.01)
    hold on



%%

%
% pred = path.pred;
% if ~exist('sortProp','var') | isempty(sortProp)
%     sortProp = path.lengths;
% end
% [p pIDX] = sort(sortProp,'descend')
%
% if ~exist('tips','var') | isempty(tips)
%     tips = pIDX;
% end
%

if ~exist('prop','var') | isempty(prop)
    prop = ones(size(skel.node2subs,1),1);
end
%
% if ~exist('N','var') | isempty(N)
%     N = length(path.tips);
% end

col = jet(1000);

minProp = min(prop);
maxProp = max(prop);



sSubs = skel.nodePos;
sSubs = sSubs - repmat(offset,[size(sSubs,1) 1]);

for tS = 1:length(skel.bones)
    
    if 0
        clf
        [p] = renderFV(fv,[0 0 1],.4);
        view([0 0])
        axis off
        hold on
    end
    
    
    nList = [skel.bones(tS).tip skel.bones(tS).nodes skel.bones(tS).base ];
    pList = skel.node2surf(nList);
    
    runVal = prop(nList);
    runSize = runVal.^2 * 100;
    runSize(runSize<1) = 1;
    runSize = runSize*0+30;
   
    runCol = nodeCol(nList,:);
    %runE = subs(pList,:);
    runE = sSubs(nList,:);
    scatter3(runE(:,1),runE(:,2),runE(:,3),runSize,runCol,'o','filled')
    plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',1,'color','w')
    %scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')
    
    view([45 45])
    pause(.01)
    
end


if isfield(skel,'bridges')
    for b = 1:size(skel.bridges,1)
        nList = skel.bridges(b,:);
        pList = skel.node2surf(nList);
        runE = subs(pList,:);
        plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','y')
    end
end


ax = gca;
ax.Clipping = 'off';
set(gca,'color',[0 0 0])




