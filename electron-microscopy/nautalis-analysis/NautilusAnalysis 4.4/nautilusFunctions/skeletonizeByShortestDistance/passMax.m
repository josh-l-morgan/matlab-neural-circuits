function[passProp] = passMax(con1,vProp,reps)

vNum = size(con1,1);

if size(vProp,2)>size(vProp,1)
    vProp = vProp';
end

%% shift data for easy lookup
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.
passProp = cat(1,0,vProp);
passVox = cat(2,1,1:vNum)';

%% 

for r = 1:reps;
    tic
    disp(sprintf('%d of %d',r,reps)); pause(.0001)
    
    prop26 = repmat(passProp,[1 26]);
    newProp = prop26(con2);
    voxMax = max(newProp,[],2);
    bigger = voxMax > passProp;
    passProp(bigger) = newProp(bigger);
    toc
    
end

%% move passed values out of lookup space
passProp = passProp(2:end);
passVox = passVox(2:end)-1;

