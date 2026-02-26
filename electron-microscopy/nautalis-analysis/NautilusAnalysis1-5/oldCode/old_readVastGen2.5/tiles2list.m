function[] = tiles2list(TPN)


load([TPN 'sparseTiles.mat']);


pointList.ids = cat(1,sparseTiles.id);
numVox = length(pointList.ids);
pointList.Z = zeros(numVox,1);
pointList.X = zeros(numVox,1);
pointList.Y = zeros(numVox,1);

c = 0;
for i = 1:length(sparseTiles)
    [y x] = ind2sub([sparseTiles(i).Height sparseTiles(i).Width],sparseTiles(i).ind);
    newCount = c + length(y);
    pointList.X(c+1:newCount) = sparseTiles(i).Width * (sparseTiles(i).col-1) + x;
    pointList.Y(c+1:newCount) = sparseTiles(i).Height * (sparseTiles(i).row-1) + y;
    pointList.Z(c+1:newCount) = sparseTiles(i).z;
    c = newCount;
end


[pointList.ids idx] = sort(pointList.ids,'ascend');
pointList.X = pointList.X(idx);
pointList.Y = pointList.Y(idx);
pointList.Z = pointList.Z(idx);


pointList.X = pointList.X - min(pointList.X)+1;
pointList.Y = pointList.Y - min(pointList.Y)+1;
pointList.Z = pointList.Z - min(pointList.Z)+1;


save([TPN 'pointList.mat'],'pointList')