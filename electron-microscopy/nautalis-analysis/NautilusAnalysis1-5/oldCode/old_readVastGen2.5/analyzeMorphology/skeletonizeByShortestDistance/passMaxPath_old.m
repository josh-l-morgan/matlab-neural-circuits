function[passVox passProp] = passMax(con1,vProp,reps)

vNum = size(con1,1);

if size(vProp,2)>size(vProp,1)
    vProp = vProp';
end

%% shift data for easy lookup
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.
passProp = cat(1,0,vProp);
passVox = cat(2,1,1:vNum)';

%% Add a small amount of noise to properties
propDif = abs(passProp(1:end-1)-passProp(2:end));
propDif = propDif(propDif>0);
minDif = min(propDif);
passProp = passProp + rand(size(passProp))*minDif/2;

%%
% 
% for r = 1:reps;
%     disp(sprintf('%d of %d',r,reps));
%     pause(.0001)
%    for i = 1:26
%        newVox = con2(:,i);
%        newProp = passProp(newVox);
%        
%        bigger = newProp>passProp;
%        
%        passProp(bigger) = newProp(bigger);
%        passVox(bigger) = passVox(newVox(bigger));
%    end
% end
% 

vox26 = repmat(passVox,[1 26]);
prop26 = repmat(passProp,[1 26]);

for r = 1:reps;
    disp(sprintf('%d of %d',r,reps));
    pause(.0001)
    newProp = prop26(con2);
    voxMax = max(newProp,[],2);
    %max26 = repmat(voxMax,1,26);
    %% compare max26, prop26 and newProp
    %bigger =  (max26 > prop26) & (max26 == newProp);
    
    bigger = 
    
    
   for i = 1:26
       newVox = con2(:,i);
       newProp = passProp(newVox);
       
       bigger = newProp>passProp;
       
       passProp(bigger) = newProp(bigger);
       passVox(bigger) = passVox(newVox(bigger));
   end
end



%% move passed values out of lookup space
passProp = passProp(2:end);
passVox = passVox(2:end)-1;

