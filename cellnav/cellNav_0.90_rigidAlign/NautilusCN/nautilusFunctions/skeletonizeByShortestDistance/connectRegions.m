function[voxEdges] =  connectRegions(vox, vProp)

%% list of which region connects to which region and how many times. 
%% I want a voxel list that says which region each voxel is closest to

con1 = vox.conMat;
vNum = size(con1,1);
reps = 1;

if size(vProp,2)>size(vProp,1)
    vProp = vProp';
end

%% shift data for easy lookup
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.
ownID = cat(1,0,vProp);
neighborID = ownID; % 
passDist = ownID*0+vNum+10; % keep track of distance to source of value so as to prefer shorter distances
passVox = cat(2,1,2:vNum+1)';
passPred = ownID * 0;

dist1 = vox.conDat.dists;

%% Add a small amount of noise to properties
% propDif = abs(ownID(1:end-1)-ownID(2:end));
% propDif = propDif(propDif>0);
% minDif = min(propDif);
% ownID = ownID + rand(size(ownID))*minDif/2;

%%

for r = 1:reps;
   voxOrder =  randperm(26); %randomize order of voxels to prevent overall drift?
   for iV = 1:26
       i = voxOrder(iV);
       newVox = con2(:,i);
       newDist = passDist(newVox)*0 + dist1(i); %disnance through connection
       newID = ownID(newVox); %property to be dilated
       
       better = (newID ~= ownID) & (newDist<passDist) & (newID>0);
       
       passDist(better) = newDist(better);
       neighborID(better) = newID(better);
       passVox(better) = passVox(newVox(better));
       passPred(better) = newVox(better);
   end
end

%% move passed values out of lookup space
path2max.source = vox.name;
path2max.prop = neighborID(2:end);
path2max.vox = passVox(2:end)-1;
path2max.pred = passPred(2:end) - 1;
path2max.dist = passDist(2:end);



%% connectivity


voxEdges = [vProp path2max.prop];





