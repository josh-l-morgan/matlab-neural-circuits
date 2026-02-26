function[linkedNodes] = nearNodes(edges,startNodes,ord)

% 
% edges = [1 2; 2 3; 3 4 ; 3 5 ; 3 6]
% startNodes = 1;

%%
checkEdges = find(ones(size(edges,1),1));
oldNodes = startNodes;
linkedNodes = [startNodes];
for i = 1:ord
   
    c = ismember(edges(checkEdges,1),oldNodes);
    newPost = edges(checkEdges(c),2);
    checkEdges(c) = [];
    
    c = ismember(edges(checkEdges,2),oldNodes);
    newPre = edges(checkEdges(c),1);
    checkEdges(c) = [];

    oldNodes = [newPre(:); newPost(:)];
    linkedNodes = unique([linkedNodes(:); newPre(:); newPost(:)]);
    
end
linkedNodes = setdiff(linkedNodes,startNodes);

%%
% checkEdges = find(ones(size(edges,1),1));
% oldNodes = startNodes;
% linkedNodes = [startNodes];
% for i = 1:ord
%    
%     c = ismember(edges(:,1),oldNodes);
%     newPost = setdiff(edges(c,2),linkedNodes);
%     %checkEdges(c) = [];
%     
%     c = ismember(edges(:,2),oldNodes);
%     newPre = setdiff(edges(c,1),linkedNodes);
%     %checkEdges(c) = [];
% 
%     oldNodes = [newPre(:); newPost(:)];
%     linkedNodes = unique([linkedNodes(:); newPre(:); newPost(:)]);
%     
% end
% linkedNodes = unique(setdiff(linkedNodes,startNodes));
% 

