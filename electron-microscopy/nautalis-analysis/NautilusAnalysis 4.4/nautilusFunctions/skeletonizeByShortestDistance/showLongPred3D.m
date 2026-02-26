function[col] = showLongPred3D(subs,path)


seed = path.seed;
dists = path.dist;
pred = path.pred;

if ~exist('startNode','var')
    hasPred = find(pred>0);
    isLong = find(dists(hasPred) == max(dists(hasPred)),1);
    startNode = hasPred(isLong);
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

[p] = renderFV(fv,[0 0 1],.4);
view([0 0])
axis off
hold on


pList = startNode;
if ~isempty(pList)
for i = 1:size(subs,1)
    prev = pred(pList(end));
    if ~prev,break,end
    pList(i+1) = prev;
end


runE = subs(pList,:);
scatter3(runE(:,1),runE(:,2),runE(:,3),300,'r','.')
scatter3(runE(1,1),runE(1,2),runE(1,3),500,'r','.')
scatter3(runE(end,1),runE(end,2),runE(end,3),500,'m','.')


plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','g')
scatter3(subs(seed,1),subs(seed,2),subs(seed,3),1000,'w','.')

view([45 45])
pause(.01)
end
hold off













