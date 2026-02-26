function[path2max] = passMaxAndDist(vox,vProp,reps)

con1 = vox.conMat;
vNum = size(con1,1);

if size(vProp,2)>size(vProp,1)
    vProp = vProp';
end

%% shift data for easy lookup
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.
passProp = cat(1,min(vProp)-1,vProp);
passDist = passProp*0; % keep track of distance to source of value so as to prefer shorter distances
passVox = cat(2,1,2:vNum+1)';
passPred = passProp * 0;

dist1 = vox.conDat.dists;

%% Add a small amount of noise to properties
% propDif = abs(passProp(1:end-1)-passProp(2:end));
% propDif = propDif(propDif>0);
% minDif = min(propDif);
% passProp = passProp + rand(size(passProp))*minDif/2;

%%

for r = 1:reps;
    disp(sprintf('%d of %d',r,reps));
    tic
    pause(.001)
   voxOrder =  randperm(26); %randomize order of voxels to prevent overall drift?
   for iV = 1:26
       i = voxOrder(iV);
       newVox = con2(:,i);
       newDist = passDist(newVox) + dist1(i); %disnance through connection
       newProp = passProp(newVox); %property to be dilizted
       
       bigger = newProp > passProp; %if the connected position + distance to connection is smaller than current distance
       closer = (newProp >= passProp) &  (newDist < passDist);
       better = bigger | closer;
       
       passDist(better) = newDist(better);
       passProp(better) = newProp(better);
       passVox(better) = passVox(newVox(better));
       passPred(better) = newVox(better);
   end
   toc
end

%% move passed values out of lookup space
path2max.source = vox.name;
path2max.prop = passProp(2:end);
path2max.vox = passVox(2:end)-1;
path2max.pred = passPred(2:end) - 1;
path2max.dist = passDist(2:end);


