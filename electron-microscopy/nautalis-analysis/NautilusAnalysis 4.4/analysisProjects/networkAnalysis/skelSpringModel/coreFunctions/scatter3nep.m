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
        
            for s = 1:length(group(g).ind)
                sg = scatter(nodeX(ind(s)),nodeY(ind(s)),group(g).size,nodeCol(ind(s),:),nodeMarker{ind(s)},'filled');
                set(sg,'LineWidth',group(g).lineWidth);
            end
            
            % %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,nodeCol(group(g).ind,:),group(g).marker,'filled');
            %sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,group(g).ind,group(g).marker,'filled');

            %%sg = scatter(nodeX(group(g).ind),nodeY(group(g).ind),group(g).size,group(g).ind,group(g).marker,'filled');
            
            %%set(sg,'MarkerFaceColor',[group(g).ind 3]')
            %%set(sg,'MarkerEdgeColor',edgeColor);
            %set(sg,'LineWidth',group(g).lineWidth);
            %%set(sg,'SizeData',c;
            
            if group(g).text;
                for t = 1:length(group(g).ind)
                    targ = group(g).ind(t);
                    tx = text(nodeX(targ)-0,nodeY(targ)-0,{num2str(nodeIDs(targ))});
                    %tx = text(nodeX,nodeY,num2str(nodeIDs))
                    set(tx,'Color','w')
                end
            end