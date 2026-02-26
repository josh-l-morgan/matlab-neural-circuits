 function[] = showNep(nep,col,wdth);

 showEdges = nep.edges;
 nodePos = nep.pos;
 

 
 
 if ~exist('col','var')
     col = 'k';
 end
 if ~exist('wdth','var')
     wdth = 1;
 end
 
 sub1 = nodePos(showEdges(:,1),:);
 sub2 = nodePos(showEdges(:,2),:);
 for e = 1:size(showEdges,1)
     plot3([sub1(e,1) sub2(e,1)], [sub1(e,2) sub2(e,2)], [sub1(e,3) sub2(e,3)],'color',col,'linewidth',wdth)
     hold on
 end
 hold off
 pause(.01)