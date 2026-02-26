function[edges] = uniqueEdges(allWeights,nodes1,nodes2);


if isempty(allWeights)
    
    edges.e1 = [];
    edges.e2 =[];
    edges.ew = [];
    edges.use = [];
    edges.all = [];
    
else
    
    allEdges = allWeights(:,1:2);
    
    if size(allWeights,2) < 3;
        synWeight = ones(size(allWeights,1),1);
    else
        synWeight = allWeights(:,3);
    end
    
    if ~exist('nodes1','var')
        nodes1 = unique(allEdges(:));
    end
    
    if ~exist('nodes2','var')
        nodes2 = nodes1;
    end
    
    useEdge = zeros(size(allEdges,1),1);
    c = 0;
    for y = 1:length(nodes1)
        for x = 1:length(nodes2)
            getEdges = (allEdges(:,1) == nodes1(y)) & ...
                (allEdges(:,2) == nodes2(x));
            eW = sum(synWeight(getEdges));
            if eW>0
                c = c+1;
                e1(c) = y;
                e2(c) = x;
                ew(c) = eW;
                useEdge = useEdge + getEdges;
                
            end
        end
    end
    
    edges.e1 = e1;
edges.e2 = e2;
edges.ew = ew;
edges.use = [e1 e2];
edges.all = allEdges;
    
    
end










