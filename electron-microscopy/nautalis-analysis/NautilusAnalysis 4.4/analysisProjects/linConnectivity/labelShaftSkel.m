function[sm] = labelShaftSkel(sm);

proxThresh = 1;

%% Get shaft

MPNshaft = 'D:\LGNs1\Export\export_joshm_dendType_cE\'
load([MPNshaft 'obI.mat'])
load([MPNshaft 'dsObj.mat'])

isShaft = [];
shaftSub = [];
allSub = [];
axSub = [];
targSub = [];
bodySub = [];

for i = 1:length(obI.nameProps.names)
    sub = double(dsObj(i).subs);
    allSub = cat(1,allSub,sub);
    if sum(regexp(obI.nameProps.names{i}, 'Body'))
       bodySub = cat(1,bodySub,sub);
    elseif sum(regexp(obI.nameProps.names{i}, 'shaft'))
        shaftSub = cat(1,shaftSub,sub);
    elseif sum(regexp(obI.nameProps.names{i}, 'axon'))
        axSub = cat(1,axSub,sub);
    else 
        targSub = cat(1,targSub,sub);
    end
end

shaftSub = double(shaftSub) * .1;
axSub = double(axSub) * .1;
targSub = double(targSub) * .1;
allSub = double(allSub) * .1;
bodySub = double(bodySub) * .1;

% scatter3(sm.skelPos(:,1),sm.skelPos(:,2),sm.skelPos(:,3),80,'r','.');

%% Check nodes

isShaft = zeros(size(sm.skelPos,1));
isTarg = isShaft;
isAx = isShaft;

dist = zeros(size(sm.skelPos,1),1);
for i = 1:size(sm.skelPos,1)
    dists = sqrt((shaftSub(:,1) - sm.skelPos(i,1)).^2 + (shaftSub(:,2)-sm.skelPos(i,2)).^2 + ...
        (shaftSub(:,3)-sm.skelPos(i,3)).^2);
    dist(i) = min(dists);
end

sm.dist2shaft = dist;

dist2 = zeros(size(sm.skelPos,1),1);
for i = 1:size(sm.skelPos,1)
    dists = sqrt((axSub(:,1) - sm.skelPos(i,1)).^2 + (axSub(:,2)-sm.skelPos(i,2)).^2 + ...
        (axSub(:,3)-sm.skelPos(i,3)).^2);
    dist2(i) = min(dists);
end


dist3 = zeros(size(sm.skelPos,1),1);
for i = 1:size(sm.skelPos,1)
    dists = sqrt((bodySub(:,1) - sm.skelPos(i,1)).^2 + (bodySub(:,2)-sm.skelPos(i,2)).^2 + ...
        (bodySub(:,3)-sm.skelPos(i,3)).^2);
    dist3(i) = min(dists);
end

sm.isBody = dist3<=proxThresh;
sm.isAx = ~sm.isBody & (dist2<= proxThresh);
sm.isTarg = ~sm.isBody & ~sm.isAx & ( dist > proxThresh);
sm.isShaft = ~sm.isBody & ~sm.isAx & ( dist <= proxThresh);

sm.shaftSub = sm.skelPos(sm.isShaft,:);
sm.bodySub = sm.skelPos(sm.isBody,:);
sm.dist2body = dist3;
sm.axSub = axSub;
sm.dist2ax = dist2;
sm.targSub = sm.skelPos(sm.isTarg,:);


sm.type = double(sm.isShaft);
sm.type(sm.isTarg>0) = 2;
sm.type(sm.isAx>0) = 3;


%% Type lengths

edgeType = sm.type(sm.skelEdges);
axEdge = sum(edgeType == 3,2)==2;
shaftEdge = sum(edgeType == 1,2)==2;
targEdge = sum(edgeType == 2,2)==2;

pos1 = sm.skelPos(sm.skelEdges(:,1),:);
pos2 = sm.skelPos(sm.skelEdges(:,2),:);

L = sqrt(sum((pos2-pos1).^2,2));

sm.L.all = sum(L);
sm.L.ax = sum(L(axEdge));
sm.L.targ = sum(L(targEdge));
sm.L.shaft = sum(L(shaftEdge));




%% show result



if 0
    
    scatter3(shaftSub(:,1),shaftSub(:,2),...
        shaftSub(:,3),180,'b','.');
    hold on
   scatter3(targSub(:,1),targSub(:,2),...
        targSub(:,3),180,'g','.');
   scatter3(axSub(:,1),axSub(:,2),...
        axSub(:,3),180,'r','.');
    hold off
    
end


if 0
    
    scatter3(sm.skelPos(closeShaft,1),sm.skelPos(closeShaft,2),...
        sm.skelPos(closeShaft,3),180,'b','.');
    hold on
    scatter3(sm.skelPos(sm.isAx,1),sm.skelPos(sm.isAx,2),...
        sm.skelPos(sm.isAx,3),180,'r','.');
    scatter3(sm.skelPos(sm.isTarg,1),sm.skelPos(sm.isTarg,2),...
        sm.skelPos(sm.isTarg,3),180,'g','.');
    hold off
    
end

























