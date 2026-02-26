function[passProp] = quickPassMax(con1,vProp)

vNum = size(con1,1);

if size(vProp,2)>size(vProp,1)
    vProp = vProp';
end

%% shift data for easy lookup
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.
passProp = cat(1,min(vProp)-1,vProp);

for r = 1:vNum;
    
   voxOrder =  randperm(26); %randomize order of voxels to prevent overall drift?
   foundBetter = 0;
   for iV = 1:26
       i = voxOrder(iV);
       newVox = con2(:,i);
       newProp = passProp(newVox); %property to be dilizted
       
       better = find(newProp > passProp); %if the connected position + distance to connection is smaller than current distance
       passProp(better) = newProp(better);
       foundBetter = foundBetter + length(better);
       
   end
   if ~foundBetter,break,end
   
end

passProp = passProp(2:end);
