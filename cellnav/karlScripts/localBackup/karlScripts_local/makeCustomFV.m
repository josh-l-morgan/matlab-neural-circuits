function output1=makeCustomFV(objList,dsObj)
objVoxels=vertcat(dsObj(objList).subs);
%figure(); hold on; scatter3(objVoxels(:,1),objVoxels(:,2),objVoxels(:,3),'.');
%gridSize=max(objVoxels)-min(objVoxels);
maxVals=max(objVoxels);
%clipVoxels=objVoxels-min(objVoxels)+1;
V=false(maxVals(1)+1,maxVals(2)+1,maxVals(3)+1);
V(sub2ind(size(V),objVoxels(:,1),objVoxels(:,2),objVoxels(:,3)))=true;
output1=isosurface(V);
if 0
coords=double(objVoxels);
% Create a delaunay triangulation from the coordinates
dt = delaunayTriangulation(coords);
% Extract the unique vertices from the triangulation
vertices = dt.Points;

% Extract the faces from the triangulation
faces = dt.ConnectivityList;

output2.faces=faces;
output2.vertices=vertices;
end

end
