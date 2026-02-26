function[col] = showPathL3D(subs,path,prop,N)


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
if ~exist('prop','var')
    lookUpTip(path.tips) = [1:length(path.tips)];
    prop = path.lengths(lookUpTip(path.owner));
end

if ~exist('N','var')
    N = length(path.tips);
end

col = jet(1000);

minProp = min(prop);
maxProp = max(prop);

[sortProp propIDX] = sort(path.lengths,'descend');

for tS = 1:N
    
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
    % plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','g')
    %scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')
    
    view([45 45])
    pause(.01)
    
end










