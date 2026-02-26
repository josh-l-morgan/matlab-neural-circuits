function[swcS] = sm2swc(sm)

%%
%%http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
%%https://neuroinformatics.nl/HBP/morphology-viewer/


%%Origin node parent must be -1
%%Parent nodes must occur before children nodes

load('MPN')

nodes = sm.nep.nodes(:);
type = nodes * 0;
pos = sm.nep.pos;
rad = sm.nep.nodeRad(:);
edges = sm.nep.edges;
parent = type*0;
parent(sm.nep.edges(:,1)) = sm.nep.edges(:,2);

offset = min(pos,[],1)+1;
offset = mean(pos,1);

pos = pos-repmat(offset,[size(pos,1) 1]);


datMat = cat(2,nodes,type,pos,rad,parent);



%% Reorder Nodes so that parents are always first


ne = unique(edges(:));
nodeNum = length(ne);
counts = hist(edges(:,1),ne);
uc = unique(counts(:));
cc = hist(counts,uc);

pred = zeros(nodeNum,1); % list parent of each node
free = ones(size(edges,1),1); % make list edges that have not been checked
if isfield(sm.nep,'seedNode')
    seed = sm.nep.seedNode;
else
    seedE = find(edges(:,1) == edges(:,2));
    seed = edges(seedE,1);
    free(seedE) = 0; % clear edge in which seed connects to seed
end
%edges(seed,2) = -1;

scatter3(pos(:,1),pos(:,2),pos(:,3),'.')
hold on
scatter3(pos(seed,1),pos(seed,2),pos(seed,3),'o')
hold off

last = seed; %initialize search
newID = zeros(length(unique(edges(:))),1); % create empty list for new IDs
c = 0; % initialize node counter
for i = 1:size(edges,1) % loop up to number of edges
    newN = []; % initialize matrix for collecting new nodes
    for n = 1:length(last) % Check all discovered and unchecked nodes
        c = c + 1; % intterate counter
        newID(last(n)) = c; %give ID to last node being checked
        next = find((edges(:,2)==last(n)) & free); % find all edges in which checked node is parent
        nextNode = edges(next,1); %grab the child node IDs from edges
        newN = cat(1,newN, nextNode); % collect child nodes of discoverd edges
        free(next) = 0; % clear the discovered edge from being checked again
       if ~isempty(nextNode) % if a new node has been discovered
            pred(nextNode) = c; %set parent to old position space
        end
    end
    last = newN;
    if isempty(last),break,end
end

pred(newID) = pred;
pred(1) = -1;

hold off

swcS.arbor2swcID = newID;
swcS.pred = pred;

clear nodes2 pos2 rad2
nodes2 = 1:length(newID);
pos2(newID,:) = pos;
rad2(newID,1) = rad;
type2(newID,1) = type;


%% Format swc variable
swcStr = [];
line1 = '# sm to swc';

line2 = sprintf('\n# cell ID %d', sm.cid);
line3 = sprintf('\n# position of cell zero point  = %.3f %.3f %.3f ',...
    offset(1), offset(2), offset(3));


swcStr = [swcStr line1 line2 line3 '\n'];
for i = 1 : size(datMat,1)
    newLine = sprintf('\n%d %d %.3f %.3f %.3f %.3f %d',...
        nodes2(i),type2(i),pos2(i,1),pos2(i,2),pos2(i,3),rad2(i),pred(i));
    swcStr = [swcStr newLine];
end

swcS.swcStr = swcStr;
swcS.nodes = nodes2;
swcS.pos = pos2;
swcS.rad = rad2;
swcS.type = type2; 

%%Write SWC

swcDir = [WPN 'swc'];
if ~exist(swcDir,'dir'), mkdir(swcDir); end
fileName = sprintf('%s\\cid%d.swc',swcDir,sm.cid)
fid = fopen(fileName,'w+');
fprintf(fid, swcStr);
fclose(fid);




    




















    
    
    
    