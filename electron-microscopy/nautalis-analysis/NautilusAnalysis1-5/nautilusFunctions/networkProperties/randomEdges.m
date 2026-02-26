function[edges] = randomEdges(allEdges,nodes1,nodes2);



if isempty(allEdges)
   

edges.e1 = [];
edges.e2 =[];
edges.ew = [];
edges.use = [];
edges.all = [];

else
changeEdges = allEdges + 1;
maxSubs = max(changeEdges,[],1);

edgeInd = sub2ind(maxSubs,changeEdges(:,1),changeEdges(:,2));

uind = unique(edgeInd);

if length(uind>1)
    hind = hist(edgeInd,uind);
else
    hind = length(edgeInd);
end

[y x] = ind2sub(maxSubs,uind);
y = y-1;
x = x-1;

uedges = [y x];

uweights = hind';

useEdge = zeros(length(y),2);

if exist('nodes2','var')
    for i = 1:length(y)
        targ = find(nodes1==y(i));
        if ~isempty(targ)
            useEdge(i,1) = targ;
        end
        targ = find(nodes2==x(i));
        if ~isempty(targ)
            useEdge(i,2) = targ;
        end
    end
    
    nodeIDs = [nodes1 nodes2];
    preList = nodes1;
    postList = nodes2;
else
    for i = 1:length(y)
        targ = find(nodes1==y(i));
        if ~isempty(targ)
            useEdge(i,1) = targ;
        end
        targ = find(nodes1==x(i));
        if ~isempty(targ)
            useEdge(i,2) = targ;
        end
    end
    nodeIDs = nodes1;
end

goodEdge = sum(useEdge>0,2)==2;

uedges = uedges(goodEdge,:);
uweights = uweights(goodEdge);





edges.e1 = useEdge(goodEdge,1);
edges.e2 = useEdge(goodEdge,2);
edges.ew = uweights;
edges.use = uedges;
edges.all = allEdges;
end










