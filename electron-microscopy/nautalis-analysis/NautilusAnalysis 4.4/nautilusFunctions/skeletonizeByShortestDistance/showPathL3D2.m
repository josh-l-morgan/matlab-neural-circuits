function[col] = showPathL3D2(subs,path,prop,N)


%%
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
hold off
clf
downSamp = 4;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;

smallSub = shrinkSub(subs,downSamp);
%smallSub = smallSub(:,flipDim);

fv = subVolFV(smallSub,[],renderProps);
%fv.vertices = fv.vertices * downSamp;
fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;

[p] = renderFV(fv,[0 0 1],.4);
view([0 0])
axis off
pause(.01)
hold on


%%


pred = path.pred;
if ~exist('sortProp','var') | isempty(sortProp)
    sortProp = path.lengths;
end

if ~exist('prop','var') | isempty(prop)
    prop = path.pathDist;
end

if ~exist('N','var') | isempty(N)
    N = length(path.tips);
end

col = jet(1000);

minProp = min(prop);
maxProp = max(prop);

[sortProp propIDX] = sort(sortProp,'descend');

for tS = 1:N
    
    
    if 0
        clf
        [p] = renderFV(fv,[0 0 1],.4);
        view([0 0])
        axis off
        hold on
    end
    
    
    t = propIDX(tS);
    
    pList = path.tips(t);
    seed = path.bases(t);
    for i = 1:size(subs,1)
        prev = pred(pList(end));
        if ~prev,break,end
        pList(i+1) = prev;
        if prev == seed
            break
        end
    end
    length(pList)
    
    
    runVal = prop(pList);
    runSize = runVal.^2 * 100;
    runSize(runSize<1) = 1;
    runSize = runSize*0+20;
    runColID = round((runVal-minProp)*size(col,1)/maxProp);
    runColID(runColID<1) = 1;
    runColID(runColID>size(col,1)) = size(col,1);
    runCol = col(runColID,:);
    runE = subs(pList,:);
    scatter3(runE(:,1),runE(:,2),runE(:,3),runSize,runCol,'o','filled')
    plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',1,'color','w')
    %scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')
    
    view([45 45])
    pause(.01)
    
end










