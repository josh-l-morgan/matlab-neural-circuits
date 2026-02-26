function[swcS] = nep2swc(nep)

%%
%%http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
%%https://neuroinformatics.nl/HBP/morphology-viewer/


%%Origin node parent must be -1
%%Parent nodes must occur before children nodes

load('MPN')

nodes = nep.nodes(:);
type = nodes * 0 + 3;
type(1) = 1;
pos = nep.pos;
rad = nep.nodeRad(:);
edges = nep.edges;
parent = type*0;
parent(nep.edges(:,1)) = nep.edges(:,2);

offset = min(pos,[],1)+1;
offset = mean(pos,1);

if 0; %% apply offset
    pos = pos-repmat(offset,[size(pos,1) 1]);
end

datMat = cat(2,nodes,type,pos,rad,parent);

showSteps = 1;

%% Reorder Nodes so that parents are always first

if showSteps
    clf
    plot([pos(edges(:,1),1) pos(edges(:,2),1)]',[pos(edges(:,1),2) pos(edges(:,2),2)]','k')
end

goodEdge = ~((edges(:,1)-edges(:,2))==0);
edges = edges(goodEdge,:);
ne = unique(edges(:));
nodeNum = length(ne);
counts = hist(edges(:,1),ne);
uc = unique(counts(:));
cc = hist(counts,uc);

pred = zeros(nodeNum,1); % list parent of each node
free = ones(size(edges,1),1); % make list edges that have not been checked
if isfield(nep,'seedNode')
    seed = nep.seedNode;
else
    seedE = find(edges(:,1) == edges(:,2));
    seed = edges(seedE,1);
    free(seedE) = 0; % clear edge in which seed connects to seed
end
%edges(seed,2) = -1;
if showSteps
    clf
    scatter3(pos(:,1),pos(:,2),pos(:,3),'.')
    hold on
    scatter3(pos(seed,1),pos(seed,2),pos(seed,3),'o')
    hold off
end


%% Temp trace

last = seed;
useEdge = ones(size(edges,1),1);
trackN = [];
countN = zeros(max(edges(:)),1);
newID = zeros(length(unique(edges(:))),1); % create empty list for new IDs
c = 0;

predRaw = zeros(max(edges(:)),1);
predRaw(last) = -1;
predOld = predRaw;
while 1
    
    newN = [];
    countN(last) = countN(last) + 1;
    
    if showSteps
        clf
        scatter(pos(:,1),pos(:,2),'.','k')
        hold on
        scatter(pos(trackN,1),pos(trackN,2),'.','b')
        drawnow
    end
    
    for n = 1:length(last)
        nid = last(n);
        c = c+1;
        newID(nid) = c;
        isEdges1 = (edges(:,1) == nid) & useEdge;
        isEdges2 = (edges(:,2) == nid) & useEdge;
        useEdge(isEdges1) = 0;
        useEdge(isEdges2) = 0;
        isNode = edges(isEdges1,2);
        isNode = cat(1,isNode,edges(isEdges2,1));
        newN = cat(1,newN,isNode);
        predRaw(isNode) = c;
        predOld(isNode) = nid;
        
        if 0;%showSteps
            scatter(pos(fromN,1),pos(fromN,2),'.','r')
            scatter(pos(isNode,1),pos(isNode,2),'.','g')
            plot([pos(fromN,1) pos(isNode,1)]',[pos(fromN,2) pos(isNode,2)]','r')
            drawnow
            
        end
        fromN = repmat(nid,[length(isNode) 1]);
        
        Ls = sqrt((pos(fromN,1) - pos(isNode,1)).^2 +...
            (pos(fromN,2)  - pos(isNode,2)).^2);
        
        if max(Ls)>1
            Ls
        end
        %         xlim([3 6])
        %         ylim([1 7])
    end
    hold off
    
    if isempty(newN)
        break
    end
    newN = setdiff(newN,last);
    last = unique(newN);
    
    trackN = cat(1,trackN,last);
    
end

if showSteps
    clf
    useP = (pred>0) & (newID>0);
    plot([pos(useP,1) pos(predOld(useP),1)]',[pos(useP,2) pos(predOld(useP),2)]','k')
    pause(.01)
end

if showSteps
    clf
    useP = (pred>0) & (newID>0);
    plot([pos(useP,1) pos(predOld(useP),1)]',[pos(useP,2) pos(predOld(useP),2)]','k')
    pause(.01)
end

pred = predRaw * 0;
pred(newID) = predRaw;
%pred(2:end) = 1; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

hold off

swcS.arbor2swcID = newID;
lookUp = [1:length(newID)]';
swcS.swc2arborID(newID) = lookUp;

clear nodes2 pos2 rad2
nodes2 = 1:length(newID);
pos2(newID,:) = pos;
rad2(newID,1) = rad;% * 0 + .01;
type2(newID,1) = type;
meanNodeRad2(newID,1) = nep.meanNodeRad;% * 0 + .01;


if 1 %% Recondition for testing
   disp('reconditioning for debug')
   %pred(2:end) = 1:(length(pred)-1);
   
   
   %% make spines
   if 0
   L = length(pred);
   pred3 = 1:L;
   rad3 = rad2;
   type3 = zeros(L,1) + 3;
   nodes3 = L+1:L+L;
   pos3 = pos2 + 1;
   
   pred = cat(1,pred,pred3');
   rad2 = cat(1,rad2,rad3);
   type2 = cat(1,type2,type3);
   nodes2 = cat(2,nodes2,nodes3);
   pos2 = cat(1,pos2,pos3);
   end
   
   %% alternate IDs
   type2(2:2:end) = 4;
   
   
   
   
%    pred(:) = pred*0+1;
%    pred(1) = -1;
%    rad2 = rad2*0 + .001; 
%    type2 = type2 * 3;
%    nodes2(:) = 1:length(nodes2);
   %pos2 = round(pos2);
end


if showSteps
    clf
    plot([pos2(2:end,1) pos2(pred(2:end),1)]',[pos2(2:end,2) pos2(pred(2:end),2)]','k')   
end

%% Adjust for python
pred = pred-1;
pred(1) = -1;
nodes2 = nodes2-1;


%% Format swc variable
swcStr = [];
line1 = '# nep to swc';

line2 = sprintf('\n# cell ID %d', 0);
line3 = sprintf('\n# position of cell zero point  = %.3f %.3f %.3f ',...
    offset(1), offset(2), offset(3));


%swcStr = [swcStr line1 line2 line3 '\n']; %% Add/remove header
for i = 1 : length(nodes2)
    
    newLine = sprintf('%d %d %.3f %.3f %.3f %.3f %d\n',...
        nodes2(i),type2(i),pos2(i,1),pos2(i,2),pos2(i,3),rad2(i),pred(i));
    swcStr = [swcStr newLine];
end

swcStr = swcStr(1:end-1);
swcS.pred = pred;
swcS.swcStr = swcStr;
swcS.nodes = nodes2;
swcS.pos = pos2;
swcS.rad = rad2;
swcS.type = type2;
swcS.meanNodeRad = meanNodeRad2;





























