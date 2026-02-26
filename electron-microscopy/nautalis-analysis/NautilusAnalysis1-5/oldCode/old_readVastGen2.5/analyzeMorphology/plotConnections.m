

load([MPN 'obI.mat'])


DimOrder = [1 2 3];



anchors = obI.colStruc.anchors;
anchors(:,1) = anchors(:,1)/dSamp(1);
anchors(:,2) = anchors(:,2)/dSamp(2);
anchors(:,3) = anchors(:,3)/dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
synapses;



cellList = [ 108 201]
cbPos = zeros(length(cellList),3);
for i = 1:length(cellList)
    targ = find(obI.cell.name == cellList(i));
    cbPos(i,:) = obI.cell.anchors(targ,:);
end