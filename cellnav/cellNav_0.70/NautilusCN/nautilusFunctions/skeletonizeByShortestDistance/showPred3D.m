function[col] = showPred3D(subs,pred,seed,startNode)

if ~exist('startNode','var')
    startNode = ceil(rand * size(subs,1));
end

offset = min(subs,[],1)-1;

for i = 1:3
    subs(:,i) = subs(:,i) - offset(i) + 1;
end

%% Render 3D
downSamp = 4;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;

smallSub = shrinkSub(subs,downSamp);
%smallSub = smallSub(:,flipDim);

fv = subVolFV(smallSub,[],renderProps);
%fv.vertices = fv.vertices * downSamp;
fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;

[p] = renderFV(fv,[1 1 1],.1);
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


runE = subs(pList,:);
scatter3(runE(:,1),runE(:,2),runE(:,3),300,'r','.')

plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','g')
scatter3(subs(seed,1),subs(seed,2),subs(seed,3),100,'o','w','filled')

view([45 45])
pause(.01)

hold off













