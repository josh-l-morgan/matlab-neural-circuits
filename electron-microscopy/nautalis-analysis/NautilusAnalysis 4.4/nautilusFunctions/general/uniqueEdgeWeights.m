function[edges] = uniqueEdges(allEdges,nodes1,nodes2);


if isempty(allEdges)
   
edges.e1 = [];
edges.e2 =[];
edges.ew = [];
edges.use = [];
edges.all = [];

else
    
  if ~exist('nodes2','var')
      nodes2 = nodes1;
  end
    
    useEdge = zeros(size(allEdges,1),1);
    c = 0;
   for y = 1:length(nodes1)
       for x = 1:length(nodes2)
           getEdges = (allEdges(:,1) == nodes1(y)) & ...
               (allEdges(:,2) == nodes2(x));
           eW = sum(allEdges(getEdges,3));
           if eW>0
               c = c+1;
               e1(c) = nodes1(y);
               e2(c) = nodes2(x);
               ew(c) = eW;
               useEdge = useEdge + getEdges;
               
           end
       end
   end
   
end
  
edges.e1 = e1;
edges.e2 = e2;
edges.ew = e2;
edges.use = uedges;
edges.all = allEdges;
end










