function[useEdge,allEdges2] = uniqueEdges(allEdges,nodes1,nodes2);

% 
% if isempty(allWeights)
%     
%     edges.e1 = [];
%     edges.e2 =[];
%     edges.ew = [];
%     edges.use = [];
%     edges.all = [];
%     
% else
%     
%     allEdges = allWeights(:,1:2);
%     
%     if size(allWeights,2) < 3;
%         synWeight = ones(size(allWeights,1),1);
%     else
%         synWeight = allWeights(:,3);
%     end
%     
%     if ~exist('nodes1','var')
%         nodes1 = unique(allEdges(:));
%     end
%     
%     if ~exist('nodes2','var')
%         nodes2 = nodes1;
%     end
%     
%     
    
    changeIDs = zeros(max(nodes1),1);
    changeIDs(nodes1) = 1:length(nodes1);
    allEdges2 = changeIDs(allEdges);
    if size(allEdges2,2) <2, allEdges2 = allEdges2'; end
    useEdge = sum(allEdges2>0,2)==2;
    allEdges2 = allEdges2(useEdge>0,:);

%     conSize = max(allEdges2(:));
%     
%     edgeInd = sub2ind([conSize conSize],allEdges2(:,1), allEdges2(:,2));
%     uEdge = unique(edgeInd);
%     ew = histc(edgeInd,uEdge)';
%     [e1 e2] = ind2sub([conSize conSize],uEdge);
%     
%     edges.use = [e1 e2];
%     edges.all = allEdges;
%    
%     
%     edges.e1 =e1;
%     edges.e2 = e2;
%     edges.ew = ew;
    
%     
%     
% end










