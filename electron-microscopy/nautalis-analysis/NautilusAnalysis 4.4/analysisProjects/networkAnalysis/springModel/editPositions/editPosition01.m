
springRes = [WPN 'springDat\results\'];
load([springRes 'res_all_randEdge.mat'])
load([springRes 'fourSeeds_edit1.mat'])
load([springRes 'boutWeightDiam01.mat'])
load([springRes 'turkey_plusOne125e.mat'])

load([springRes 'withFourSeeds_edit1.mat'])
load([springRes 'fourSeeds_edit1.mat'])

clf


scatter(result.nodeX,result.nodeY)
%% Rotate
scatter(result.nodeX,result.nodeY)

midX = (max(result.nodeX) - min(result.nodeX))/2 + min(result.nodeX)
midY = (max(result.nodeY) - min(result.nodeY))/2 + min(result.nodeY)

midX = 100;
midY = 100;

difX = result.nodeX-midX;
difY = result.nodeY-midY;

rads = atan2(difX,difY);
dists = sqrt(difX.^2+difY.^2);
rads = rads + .03;

result.nodeY = cos(rads) .* dists + midX;
result.nodeX = sin(rads) .* dists + midY;

scatter(result.nodeX,result.nodeY)
ylim([0 200])
xlim([0 200])



%%
moveCells = [125];
for t = 1:length(moveCells)
targCell = moveCells(t)
targ = find(result.nodeIDs == targCell);

scatter(result.nodeX,result.nodeY)
hold on
scatter(result.nodeX(targ),result.nodeY(targ),'k','filled')
hold off

[x y] = ginput

if isempty(targ)
    result.nodeX(end+1) = x(end);
    result.nodeY(end+1) = y(end);
    result.nodeIDs(end+1) = moveCells(end);
    result.cellGroups(end+1) = 0;
else
    result.nodeX(targ) = x(end);
    result.nodeY(targ) = y(end);
end

scatter(result.nodeX,result.nodeY)
hold on
scatter(result.nodeX(targ),result.nodeY(targ),'k','filled')
hold off
end
%%
if 0
save([springRes 'fourSeeds_edit_125.mat'],'result')
end