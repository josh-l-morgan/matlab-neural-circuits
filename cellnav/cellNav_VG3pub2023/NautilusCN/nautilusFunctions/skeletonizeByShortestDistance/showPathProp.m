function[] = showPathProp(subs,pathL,useTips,prop)


   % subs = surfVox.subs;
    offset = min(subs,[],1)-1;
    subs = subs - repmat(offset,[size(subs,1) 1]);
    downSamp = 4;
    renderProps.smooth = 0;
    renderProps.resize = 1;
    renderProps.smoothPatch = 0;
    smallSub = shrinkSub(subs,downSamp);
    
    pred = pathL.pred;
    %prop = atMaxAxRad;
      for t= 1:length(useTips)
        clf
        fv = subVolFV(smallSub,[],renderProps);
        fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;
        [p] = renderFV(fv,[0 0 1],.4);
        view([0 0])
        axis off
        hold on
        
   
             startNode = pathL.tips(useTips(t));
        seed = pathL.bases(useTips(t));
         
        pList = startNode;
        for i = 1:size(subs,1)
            prev = pred(pList(end));
            if ~prev,break,end
            pList(i+1) = prev;
            if prev == seed
                break
            end
        end
        
        runVal = prop(pList);
        runSize = runVal.^2 * 100;
        runSize(runSize<1) = 1;
        runE = subs(pList,:);
        scatter3(runE(:,1),runE(:,2),runE(:,3),runSize,'r','o','filled')
        plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','g')
        scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')
        
        view([45 45])
        pause(.1)
        hold  off
    end
    