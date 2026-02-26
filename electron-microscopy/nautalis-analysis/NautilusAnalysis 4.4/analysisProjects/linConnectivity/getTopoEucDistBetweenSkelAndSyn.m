function[sm] = getTopoEucDistBetweenSkelAndSyn(sm)

% %% Get data
% loadData = 1;
% if loadData
%     %MPN = GetMyDir;
%     load('MPN.mat');
%     WPN = MPN;
%     DPN = [MPN 'data\'];
%     if ~exist(DPN,'dir'),mkdir(DPN);end
%     load([MPN 'obI.mat']);
%     allEdges = obI.nameProps.edges(:,[2 1]);
%     
% end



D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(D, @nUpdateWaitbar);
c = 1;

show = 0;

%% Get skeletons


%%Load skeleton from nep
% load([MPN 'nep\skelNep125.mat'])
% sm.skelNodes = nep.nodes;
% sm.skelEdges = nep.edges;
% sm.skelPos = nep.nodePos
% pos = nep.nodePos;
% edges = nep.edges;

pos = sm.skelPos;
edges = sm.skelEdges;

lengths = sqrt((pos(edges(:,1),1)-(pos(edges(:,2),1))).^2 + ...
    (pos(edges(:,1),2)-(pos(edges(:,2),2))).^2 + ...
    (pos(edges(:,1),3)-(pos(edges(:,2),3))).^2);

pts = sm.pos;


dif = cat(3,pos(:,1)-pts(:,1)',pos(:,2)-pts(:,2)',pos(:,3)-pts(:,3)');
sm.skelEucDist = sqrt(sum(dif.^2,3));

%nepSyn = obISyn2nep(obI);
%nepSyn = linkSyn2Skel(nepSyn,skels);


%% One at a time

syn = [sm.pre sm.post];

cellID = 125;
%skel = nep;
sPos = pos;

isSyn = find((syn(:,2) == cellID) | (syn(:,1) == cellID));



pathDists = sm.skelEucDist;
distThresh = 5;

synPos = sm.pos(isSyn,:);
syns = 1:length(isSyn);

syn2skel.synID = isSyn;
syn2skel.syn = syn(isSyn,:);
syn2skel.synPos = synPos;
syn2skel.pre = sm.pre(isSyn);
syn2skel.post = sm.post(isSyn);
syn2skel.preClass = sm.preClass(isSyn);
syn2skel.postClass = sm.postClass(isSyn);

pairs = nchoosek(syns,2);
clear pLengths eucDist keepPaths keepLengths closestSkel
closestSkel = zeros(size(synPos,1),1);



dif = cat(3,sPos(:,1)-synPos(:,1)',sPos(:,2)-synPos(:,2)',sPos(:,3)-synPos(:,3)');
dists = sqrt(sum(dif.^2,3));
pathDistsTarg = zeros(size(dists,1),size(synPos,1));
parfor p = 1:size(synPos,1)
    if ~mod(p,10)
        disp(sprintf('measuring dist to skel %d of %d.',p,size(synPos,1)))
    end
    pos1 = synPos(p,:);
    
    
%     syn2Skel = sqrt((sPos(:,1)-pos1(1)).^2 + (sPos(:,2)-pos1(2)).^2 + ...
%         (sPos(:,3)-pos1(3)).^2 );
    syn2Skel = dists(:,p);
    link1 = (find(syn2Skel==min(syn2Skel),1));
    closestSkel(p) = link1;
    
    %edges = skel.edges;
    nodePos = sPos;
    startNode = link1;
    
    pp = node2nodeDist(edges,lengths,startNode);
    pathDistsTarg(:,p) = pp.dists;
    
%     if 0
%         scatter3(sPos(:,1),sPos(:,2),sPos(:,3),10,pp.dists,'.')
%         hold on
%         scatter3(pos1(1),pos1(2),pos1(3),100,'r')
%         hold off
%         pause(.01)
%     end
        send(D, p);

end

function nUpdateWaitbar(~)
        waitbar(c/size(synPos,1), h);
        c = c + 1;
end

pathDists(:,isSyn) = pathDistsTarg;



sm.skelDist = max(pathDists,sm.skelEucDist);
syn2skel.closestSkel = closestSkel;
sm.syn2skel = syn2skel;

close(h)
end



