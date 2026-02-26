function[pathRads] = mapPathRads(d2sMat)

%describe morphology of paths using 26 projections

%%

%pathL = surfVox.pathL;
%pathS = surfVox.path2surf;

d2sMat = d2sMat + .5;
conInd = lookupCon()

meanAll = mean(d2sMat,2);
medAll = median(d2sMat,2);
for i = 1:size(conInd.ax,1)
    
    ax =conInd.ax(i,:);
    mid = conInd.mid{i};
    
    meanAx(:,i) = mean(d2sMat(:,ax),2);
    meanMid(:,i) = mean(d2sMat(:,mid),2);
 
end

midFrac = meanMid./(meanAx + meanMid);
minMidFrac = min(midFrac,[],2)
minMid = min(meanMid,[],2);
maxAx = max(meanAx,[],2);
minAx = min(meanAx,[],2);

[ a b] = find(meanAx == repmat(maxAx,[1 size(meanAx,2)]));
[a uIDX] = unique(a);
atMaxAxDim = b(uIDX);
atMaxAxRad = meanMid(sub2ind(size(meanAx),a,atMaxAxDim));


pathRads.atMaxAxRad = atMaxAxRad;
pathRads.meanAx = meanAx;
pathRads.meanMid = meanMid;
pathRads.meanAll = meanAll;
pathRads.medAll = medAll;


%% Draw from tip

if 0
    
    
    subs = surfVox.subs;
    offset = min(subs,[],1)-1;
    subs = subs - repmat(offset,[size(subs,1) 1]);
    downSamp = 4;
    renderProps.smooth = 0;
    renderProps.resize = 1;
    renderProps.smoothPatch = 0;
    smallSub = shrinkSub(subs,downSamp);
    
    pred = pathL.pred;
    prop = atMaxAxRad;
    [sortProp sortIDX] = sort(pathL.lengths,'descend');
    
    for t= 1:length(pathL.tips)
        startNode = pathL.tips(sortIDX(t));
        seed = pathL.bases(sortIDX(t));
        clf
        fv = subVolFV(smallSub,[],renderProps);
        fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;
        [p] = renderFV(fv,[0 0 1],.4);
        view([0 0])
        axis off
        hold on
        
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
        runE = subs(pList,:);
        scatter3(runE(:,1),runE(:,2),runE(:,3),runSize,'r','o','filled')
        plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','g')
        scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')
        
        view([45 45])
        pause
        hold off
        
    end
    
    
    
end