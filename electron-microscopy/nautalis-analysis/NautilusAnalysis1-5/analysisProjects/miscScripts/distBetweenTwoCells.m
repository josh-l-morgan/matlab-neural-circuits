


anchor1 = scaleSubs(obI.cell.anchors(obI.cell.names == 108,:),anchorScale)
anchor2 = scaleSubs(obI.cell.anchors(obI.cell.names == 201,:),anchorScale)

diffs = anchor1-anchor2;
dist = sqrt(sum(diffs.^2));

disp(['distance is ' num2str(dist)])
