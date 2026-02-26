function[passVox passProp passPred] = passDist(con1,conDat,vProp,reps)

vNum = size(con1,1);

if size(vProp,2)>size(vProp,1)
    vProp = vProp';
end

%% shift data for easy lookup
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.
passProp = cat(1,min(vProp)-max(vProp),vProp);
passVox = cat(2,1,2:vNum+1)';
passPred = passProp * 0;

dist1 = conDat.dists;

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
       newProp = passProp(newVox); %disnance through connection
       bigger = newProp > passProp; %if the connected position + distance to connection is smaller than current distance
      
       passProp(bigger) = newProp(bigger);
       passVox(bigger) = passVox(newVox(bigger));
       passPred(bigger) = newVox(bigger);
   end
   toc
end

%% move passed values out of lookup space
passProp = passProp(2:end);
passVox = passVox(2:end)-1;
passPred = passPred(2:end) - 1;

