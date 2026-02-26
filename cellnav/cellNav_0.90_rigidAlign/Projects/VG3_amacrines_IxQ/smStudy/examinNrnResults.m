

try close(examFig); end
examFig = figure;
exAx = gca(examFig)
cid = 2028;
makeMPNcnv
load('MPN.mat')
swcDir = [WPN 'swc\'];

smDir = [WPN 'SMs\'];
swcDir = [WPN 'swc\'];

smxName = sprintf('smx_cid%d.mat',cid);
smx = load([smDir smxName],'syn','nep','syn2Skel');
nrnName = sprintf('nrn_cid%d.mat',cid);
nrn = load([swcDir nrnName]);

swc = smx.nep.swcS;
pred = swc.pred;

%%map synapse to edge
synPos = smx.syn.pos;
pos =  smx.nep.pos;
nrnPos = nrn.pos;
swcPos = swc.pos;


%%Reshape
nrnPos = nrnPos(:,[2 1 3]);
meanPos = mean(pos,1);
%synPos = synPos-repmat(meanPos,[size(synPos,1) 1]);
%pos = pos-repmat(meanPos,[size(pos,1) 1]);
%swcPos = swcPos-repmat(meanPos,[size(swcPos,1) 1]);

%pos = nrnPos(:,[2 1 3]);

W = nrn.maxVolt;

a = smx.nep.swcS.pred-[1:length(smx.nep.swcS.pred)]';

%% test skeleton
pred = swc.pred+1; %adjust index from python to matlab
nodes = swc.nodes' + 1;
last = find(pred<0);
fillin = pred*0>0;
fillin(last) = 1;
while ~isempty(last)
    next = [];
    for n = 1:length(last)
        next = cat(1,next,find(pred==last(n)));
    end
    last = next;
    fillin(last) = 1;
%     scatSkel4 = scatter3(exAx,swcPos(fillin,1),swcPos(fillin,2),...
%         swcPos(fillin,3),50,'.','k');
%     drawnow
end

predNotNode = setdiff(unique(pred),unique(nodes))
predBiggerThanNode = find(pred>=nodes)
extraNodes = length(nodes) - length(unique(nodes))




%% match nodes

nrnMatch = zeros(size(swcPos,1),1);
nrnDist = nrnMatch; 
for i = 1:size(swcPos,1)
    
    dists = sqrt((nrnPos(:,1)-swcPos(i,1)).^2 + ....
        (nrnPos(:,2)-swcPos(i,2)).^2 + ....
        (nrnPos(:,3)-swcPos(i,3)).^2);
    minDist = min(dists);
    nrnDist(i) = minDist;
    nrnMatch = find(dists==minDist,1);    
    
end
lostSwc = find(nrnDist>.1)
lostSwc = setdiff(lostSwc,1)


%% length rad
cla
L = sqrt((swcPos(2:end,1)-swcPos(pred(2:end),1)).^2 + ...
   (swcPos(2:end,2)-swcPos(pred(2:end),2)).^2+ ...
    (swcPos(2:end,3)-swcPos(pred(2:end),3)).^2);

R = swc.rad(2:end);

scatter(exAx,L,R,'k','.')
hold on
scatter(exAx,L(lostSwc-1),R(lostSwc-1),'r','o')
plot([0 1],[0 1])
hold off

%%  show edges


cla
e1 = 2:size(swcPos,1);
e2 = pred(2:end);
plot3(exAx,[swcPos(e1,1) swcPos(e2,1)]', [swcPos(e1,2) swcPos(e2,2)]',...
    [swcPos(e1,3) swcPos(e2,3)]');
set(gca,'Clipping','off')
hold on
    scatSkel5 = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),100,'.','k');
    scatSkel5 = scatter3(exAx,swcPos(lostSwc,1),swcPos(lostSwc,2),swcPos(lostSwc,3),100,'+','k');
hold off






%% show nrnOrder
if 0
    cla
    hold on
    for i = 1:size(nrnPos,1)
        
        scatSkel = scatter3(exAx,nrnPos(i,1),nrnPos(i,2),nrnPos(i,3),50,'r','o','filled','markerfacealph',.3);
        try     scatSkel5 = scatter3(exAx,swcPos(i,1),swcPos(i,2),swcPos(i,3),100,'+','k');
            end
        drawnow
        pause
    end
    
    hold off
end

%% Show weights
cla
hold on
cmap = colorcube(1000);
for i = 1:size(W,1);
    subplot(size(W,1),1,i),cla
    
    plot(W(i,:),'color',cmap(ceil(rand*900),:))
end
hold off

%% Show result
subplot(1,1,1)
exAx = gca
for s = 1:size(W,1)
    
    prop = W(s,:);
    propCol  = prop;
    propCol = propCol-min(propCol)+1;
    propCol = propCol*100/max(propCol);
    
    scatSkel2 = scatter3(exAx,pos(:,1),pos(:,2),pos(:,3),'o','k');
    hold on
    
    scatSkel = scatter3(exAx,nrnPos(:,1),nrnPos(:,2),nrnPos(:,3),50,'o','filled','markerfacealph',.5);
    scatSkel.CData = propCol;
    set(gca,'Clipping','off')
    %scatSkel4 = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),50,'x','k');
    scatSkel5 = scatter3(exAx,swcPos(lostSwc,1),swcPos(lostSwc,2),swcPos(lostSwc,3),100,'+','k');
    
    scatSynAll = scatter3(exAx,synPos(:,1),synPos(:,2),synPos(:,3),150,'o','g','filled','markerfacealph',.3);
    try scatSynS = scatter3(exAx,synPos(s,1),synPos(s,2),synPos(s,3),200,'o','k'); end
    scatSynS = scatter3(exAx,swcPos(1,1),swcPos(1,2),swcPos(1,3),200,'o','r','filled','markerfacealpha',.5);
    drawnow
    hold off
    
    
    pause(1)
    
end






    
    
