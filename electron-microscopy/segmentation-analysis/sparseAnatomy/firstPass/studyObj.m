clear all

SPN = 'C:\Users\joshm\Documents\myWork\LGNs1\jlmHomeSeg\export1Mat\'
TPN = [SPN 'cell1\'];

%% Load object data
load([SPN 'objs.mat'])

%% Pick an object
for i = 1:length(o)
    objSize(i,:) = size(o{i});
end

objNum = find(objSize(:,1) == max(objSize(:,1)));
objSubs = o{objNum};


%% Find center
midObj = median(objSubs,1);

%% turn object into voxel connectivity matrix
tic
disp('Turn object into connectivity matrix')
[conMat conDat] = obj2con(objSubs);
toc


%% Get obj surface
 [objSurf conSurf ] = subs2surf(objSubs, conMat);

 
%% Find Seed
dist2mid = sqrt((objSurf(:,1)-midObj(1)).^2 + (objSurf(:,2)-midObj(2)).^2 + ...
    (objSurf(:,3)-midObj(3)).^2);

minDist = min(dist2mid);
seed = find(dist2mid==minDist,1);

%% find shortest path from seed to surface voxels

tic 
disp('Find shortest path from seed to all surfaces')
[pred dist2seed] = conMat2surfaceShortest(conSurf, conDat,seed);
toc

%% Distance to surface
reps = 2;
vProp  = (sum(conMat(:,1:6)>0,2)>5) * size(conMat,1) * 2; %make surface voxel distance 0;
[passSurfVox passSurfProp passSurfPred] = passDist(conMat,conDat,vProp,reps);

%% Find tips by spreading ditances to seed
reps = 5;
[passTipVox passTipProp passTipPred] = passMaxPath(conSurf,conDat,dist2seed,reps);
tipness = hist(passTipVox,1:size(conSurf,1));
reps = 20;
[passTipVox2 passTipProp2 passTipPred2] = passMaxPath(conSurf,conDat,tipness,reps);
tipness2 = hist(passTipVox2,1:size(conSurf,1));
tips = find((dist2seed>0) & (tipness2>10));




%% Shorten paths

tips = find(tipness2>20);
pathDist = zeros(size(conSurf,1),1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = dist2seed(tip);
    for r = 1:size(conSurf,1)*2
        if tip<1,break,end
       pathDist(tip) = max(pathDist(tip),tipDist); 
       tip =  pred(tip);
    end
end


reps = 20;
[pix2maxDistVox pix2maxDistProp pix2maxDistPred pix2maxDistDist] = passMaxAndDist(conSurf,conDat,pathDist,reps);

bases = pix2maxDistVox(tips); %find new base for each tip
maxPred = pix2maxDistPred; % combine pred lists
maxPred(maxPred <1) = pred(maxPred <1);

tips = find(tipness2>20);
skelTip = zeros(size(conSurf,1),1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = dist2seed(tip);
    for r = 1:size(conSurf,1)*2
        if tip<1,break,end
       skelTip(tip) = tips(t); 
       tip =  maxPred(tip);
    end
end

%%

reps = 20;
[dist2SkelVox dist2SkelProp dist2SkelPred dist2SkelDist] = passMaxAndDist(conSurf,conDat,skelTip,reps);


return

%% Distance between tips
% reps = 60;
% vProp  = (tipness2<10) * size(conMat,1) * 2; %make surface voxel distance 0;
% [tip2tipVox tip2tipProp tip2tipPred] = passDist(conSurf,conDat,vProp,reps);
% 
% vProp  = (tipness2>10) ; %make surface voxel distance 0;
% [tip2tipVox tip2tipProp tip2tipPred] = passMaxPath(conSurf,conDat,tipness,reps);

%% show Path
colormap gray(256)
minSubs = min(objSurf,[],1);
clear moveObj
for i = 1:3
   moveObj(:,i) = objSurf(:,i)-minSubs(i)+1; 
end
maxSubs = max(moveObj,[],1);



%%
objVol = zeros(maxSubs,'uint8');

objVol(sub2ind(maxSubs,moveObj(:,1),moveObj(:,2),moveObj(:,3))) = ...
    mod(skelTip,256);

%objVol(sub2ind(maxSubs,moveObj(tips,1),moveObj(tips,2),moveObj(tips,3))) = 1;

sumObj = sum(objVol,3);
maxObj = max(objVol,[],3);
objIm = sumObj * 10 + ((sumObj>0)*20);
showPath = zeros(maxSubs(1),maxSubs(2));

image(maxObj*1000);
colormap colorcube(256)
image(maxObj)

return
%%
col(:,:,3) = objIm;
%col(:,:,1) = maxTipness;


showPathZ = zeros(maxSubs(1),maxSubs(2));
showPathY = zeros(maxSubs(2),maxSubs(3));
showPathX = zeros(maxSubs(1),maxSubs(3));
for r = 1:length(tips)

tip = tips(r);

for i = 1:20000
    if tip<1
        break
    end
    showPathZ(moveObj(tip,1),moveObj(tip,2)) = showPathZ(moveObj(tip,1),moveObj(tip,2)) + 1;
    showPathY(moveObj(tip,2),moveObj(tip,3)) = showPathY(moveObj(tip,2),moveObj(tip,3)) + 1;
    showPathX(moveObj(tip,1),moveObj(tip,3)) = showPathX(moveObj(tip,1),moveObj(tip,3)) + 1;
    tip = maxPred(tip);
end
col(:,:,2) = showPathZ*100;

image(uint8(col))

pause(.01)
end
return

showPathZ = zeros(maxSubs(1),maxSubs(2));
showPathY = zeros(maxSubs(2),maxSubs(3));
showPathX = zeros(maxSubs(1),maxSubs(3));
for r = 1:length(tips)

tip = bases(r);

for i = 1:20000
    if tip<1
        break
    end
    showPathZ(moveObj(tip,1),moveObj(tip,2)) = showPathZ(moveObj(tip,1),moveObj(tip,2)) + 1;
    showPathY(moveObj(tip,2),moveObj(tip,3)) = showPathY(moveObj(tip,2),moveObj(tip,3)) + 1;
    showPathX(moveObj(tip,1),moveObj(tip,3)) = showPathX(moveObj(tip,1),moveObj(tip,3)) + 1;
    tip = pred(tip);
end
col(:,:,2) = showPathZ*100;

image(uint8(col))

pause(.01)
end


