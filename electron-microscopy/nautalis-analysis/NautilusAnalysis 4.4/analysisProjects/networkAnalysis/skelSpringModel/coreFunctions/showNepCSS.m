 function[] = showNep(nep,col,wdth);
%%
 showEdges = nep.edges;
 nodePos = nep.nodePos;
 useNodes = nep.useNode;
 
 %% Show edges
  clf, hold on
    %%edgeTypeCode: {'pre'  'post'  'skel'  'cell2skel'}
 col = [0 0 1; .6 .6 0; .6 .1 .1; 0 0 0];
 edgeOrder = [1 2 3 4];
 wdth = [1 1 1 1];
 
 for i = 1:length(edgeOrder)
     t = edgeOrder(i);
    showEdges = nep.edges(nep.edgeType == t,:);
     sub1 = nodePos(showEdges(:,1),:);
    sub2 = nodePos(showEdges(:,2),:);
     for e = 1:size(showEdges,1)
        plot([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)],'color',col(i,:),'linewidth',wdth(i))
        
    end
 end
 
 %%Show nodes
 %%nodeTypeCode: {'cell'  'syn'  'skel'}
 col = [0 .7 0; 1 1 0; 1 0 0; 0 0 0];
 nodeOrder = [1 2 3 4];
 wdth = [30 20 8 10];
 
 for i = 1:4
    sub1 = nep.nodePos((nep.nodeType == nodeOrder(i))& useNodes,:);
    scatter(sub1(:,1), sub1(:,2),wdth(i),'MarkerFaceColor',col(i,:),'MarkerEdgeColor','k');     
 end
 
 
 
 
 
 