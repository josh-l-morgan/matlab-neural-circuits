
if 0
load(['MPN.mat'])

load([MPN '\vastSubs.mat']);
load([MPN '\dsObj.mat']);
load([MPN '\obI.mat']);
end

tag1 = 'cid3097';
nams = obI.colStruc.names;
obIsCid = zeros(length(nams),1);
for i = 1:length(nams)
    nam = nams{i};
    if sum(regexp(nam,tag1));
        obIsCid(i) = 1;
    end
end

synAnchors = double(obI.colStruc.anchors(obIsCid>0,[2 1 3]));
synAnchorsUM = synAnchors .* repmat(obI.em.res,[size(synAnchors,1) 1])/1000;

ax = gca;
cla(ax)
ax.NextPlot = 'add';

if 0

cellVastSub = cat(1,vastSubs{obIsCid>0});

synAnchorsVS = synAnchorsUM  ./ repmat(obI.em.vRes(1,:),[size(synAnchorsUM,1) 1])*1000; 
scatter3(cellVastSub(:,1),cellVastSub(:,2),cellVastSub(:,3),'.','k')
scatter3(synAnchorsVS(:,1),synAnchorsVS(:,2),synAnchorsVS(:,3),'o','r')

else
cellSub = cat(1,dsObj(obIsCid>0).subs);
synAnchorsDS = synAnchorsUM / obI.em.dsRes(1);

scatter3(synAnchorsDS(:,1),synAnchorsDS(:,2),synAnchorsDS(:,3),'o','r')
scatter3(cellSub(:,1),cellSub(:,2),cellSub(:,3),'.','k')
end



% edges = obI.nameProps.edges(:,1:2);
% synIsCid = sum(edges==3097,2);
%synIsCid = obIsCid;
%synAnchors = double(obI.colStruc.anchors(obI.nameProps.edges(synIsCid>0,3),[2 1 3]));
