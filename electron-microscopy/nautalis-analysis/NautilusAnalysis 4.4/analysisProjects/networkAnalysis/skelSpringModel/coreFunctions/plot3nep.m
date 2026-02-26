function[] = plot3nep(nep)

useEdge = find(nep.edgeType == 3);
nodeX = nep.nodePos(:,1);
nodeY = nep.nodePos(:,2);
nodeZ = nep.nodePos(:,3);
e1 = nep.edges(:,1);
e2 = nep.edges(:,2);
edgeCol = nep.edgeCol;
edgeWidth = 3;

pl = plot3([nodeX(e2(useEdge)) nodeX(e1(useEdge))]' ,...
    [nodeY(e2(useEdge)) nodeY(e1(useEdge))]',...
    [nodeZ(e2(useEdge)) nodeZ(e1(useEdge))]');

xlim([0 max(nep.nodePos(:,1))])
ylim([0 max(nep.nodePos(:,2))])
zlim([0 max(nep.nodePos(:,3))])
  axis square
        axis off
        set(gcf,'color','k')
        


for p = 1:length(pl)
    set(pl(p),'color',edgeCol(useEdge(p),:),'LineWidth',edgeWidth);
end
        
        
        