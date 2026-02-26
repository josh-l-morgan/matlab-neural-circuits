function[col] = showAllPred3D(subs,pred,startNode)


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

[p] = renderFV(fv,[0 0 1],.4);
view([0 0])
axis off
hold on


noPred = find(pred==0);
pred(noPred) = noPred;
predSubs = subs(pred,:);

subSamp = 1:1:size(subs,1);
scatter3(subs(subSamp,1),subs(subSamp,2),subs(subSamp,3),'.','r');
hold on
scatter3(subs(noPred,1),subs(noPred,2),subs(noPred,3),10,'x','w')
plot3([subs(subSamp,1) subs(pred(subSamp),1)]' ,...
    [subs(subSamp,2) subs(pred(subSamp),2)]',...
    [subs(subSamp,3) subs(pred(subSamp),3)]','linewidth',1,'color','g')

view([45 45])
pause(.01)

hold off













