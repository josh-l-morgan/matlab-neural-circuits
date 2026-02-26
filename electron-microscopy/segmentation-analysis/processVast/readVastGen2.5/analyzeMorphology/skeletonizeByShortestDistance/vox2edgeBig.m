function[pred dist2seed con1] = vox2edgeBig(obj,seed);

%% study

if ~exist('seed','var')
    seed = 1;
end

vNum = size(obj,1);

% %% create region lists
% rsize = 100;
% roundObj = round(obj/rsize)+1;
% regions = cell(max(roundObj));
% roundInd = sub2ind(size(regions),roundObj(:,1),roundObj(:,2),roundObj(:,3));
% numRegions = max(roundInd);
% for i = 1:numRegions
%     regions{i} = find(roundInd == i);
% end

% 
% for i = 1:length(roundInd)
%     if ~mod(i,1000)
%         disp(i)
%     end
%     curInd = roundInd(i);
%    regions{curInd} = [regions{curInd} i]; 
% end
%     
toc


%% group subs by z pos
difZ = obj(2:end,3) - obj(1:end-1,3);
changeZ = cat(1,1,find(difZ>0)+1);

for i = 1:length(changeZ)
   useZ(i,1) = changeZ(max(1,i-1));
   useZ(i,2) = changeZ(min(length(changeZ),i+2))-1;
end
useZ(end-1:end,2) = size(obj,1);

%% turn voxel list into connectivity list
con1 = zeros(vNum,26);
dist1 = con1;
conThresh = 1.9;

% 
% 
% for i = 1:vNum;
%     if sum(changeZ==i);
%         curZ = obj(i,3);
%        nearObj = obj(useZ(curZ,1):useZ(curZ,2),:);
%        disp(sprintf('%02.1f%% done',i/vNum*100))
%        pause(.01)
%     end
%     
%     dists = sqrt((nearObj(:,1)-obj(i,1)).^2 + (nearObj(:,2)-obj(i,2)).^2 ...
%         +(nearObj(:,3) - obj(i,3)).^2);
%     conTo = find((dists >0) & (dists<conThresh));
%     realConTo = conTo + useZ(curZ)-1;
%     con1(i,1:length(conTo)) = realConTo;
%     dist1(i,1:length(conTo)) = dists(conTo);
% end


for z = 1: size(useZ,1)
           nearObj = obj(useZ(z,1):useZ(z,2),:);
           disp(sprintf('%d of %d',z,size(useZ))),pause(.001)
for i = useZ(z,1):useZ(z,2);
    
     dists = sqrt((nearObj(:,1)-obj(i,1)).^2 + (nearObj(:,2)-obj(i,2)).^2 ...
        +(nearObj(:,3) - obj(i,3)).^2);
     conTo = find((dists >0) & (dists<conThresh));
    realConTo = conTo + useZ(curZ)-1;
    con1(i,:) = [realConTo; zeros(26-length(realConTo),1)];
    dist1(i,:) = [dists(conTo); zeros(26-length(realConTo),1)];
end %run i

end %run z subdivisions
    
%%
shiftMat = ones(3,3,3);
shiftMat(2,2,2) = 0;
[y x z] = ind2sub([3 3 3],find(shiftMat))
shift = [y x z];
newSubs = [];
con1 = [];
for z = 1: size(useZ,1)
    disp(sprintf('%d of %d',z,size(useZ,1)));
           nearObj = obj(useZ(z,1):useZ(z,2),:);
           maxNear = max(nearObj,[],1);
           minNear = min(nearObj,[],1);
           nearMat = zeros(maxNear(1) - minNear(1) + 3, maxNear(2) - minNear(2) + 3,3);
           nearObjMove = nearObj;
           nearObjMove(:,1) = nearObjMove(:,1)-minNear(1)+2;
           nearObjMove(:,2) = nearObjMove(:,2)-minNear(2)+2;
           nearObjMove(:,3) = nearObjMove(:,3)-z+2;
           
           nearInd = sub2ind(size(nearMat),nearObjMove(:,1),nearObjMove(:,2),nearObjMove(:,3));
           nearMat(nearInd) = 1;
           
           frameY = size(nearMat,1)-3;
           frameX = size(nearMat,2)-3;
           middleInd = find(nearMat(2:end-1,2:end-1,2));
           conN = zeros(length(middleInd),length(shift));
           for s = 1:size(shift,1)
              cutOut = nearMat(shift(s,1):shift(s,1)+frameY,shift(s,2):shift(s,2)+frameX,shift(s,3));
              conN(:,s) = cutOut(middleInd);
           end
           
           [y x] = ind2sub(size(cutOut),middleInd);
           y = y+ minNear(1)-1;
           x = x+ minNear(2) - 1;
           zList = ones(size(y))* z;
           newSubs = cat(1,newSubs,[y x zList]);
           con1 = cat(1,con1,conN);
           
           
end
%            
%            disp(sprintf('%d of %d',z,size(useZ))),pause(.001)
% for i = useZ(z,1):useZ(z,2);
%     
%      dists = sqrt((nearObj(:,1)-obj(i,1)).^2 + (nearObj(:,2)-obj(i,2)).^2 ...
%         +(nearObj(:,3) - obj(i,3)).^2);
%      conTo = find((dists >0) & (dists<conThresh));
%     realConTo = conTo + useZ(curZ)-1;
%     con1(i,:) = [realConTo; zeros(26-length(realConTo),1)];
%     dist1(i,:) = [dists(conTo); zeros(26-length(realConTo),1)];
% end %run i
% 
% end %run z subdivisions
%     

%% find shortest path assuming no loops
con2 = con1;
dist2seed = ones(1,vNum)*vNum;
dist2seed(seed) = 0;
pred = zeros(vNum,1);

nextCon = seed; %initialize next verticies to try

for vRep = 1:vNum
    lastCon = nextCon;
    nextCon = [];
    for i = 1:length(lastCon) % run all next verticies
        cV = lastCon(i); % grab current vertex
        conTo = con2(cV,:); % find who is connected to current vertex
        con2(cV,:) = 0; % clear connections once vertex is run once
        conTo = conTo(conTo>0); % grab actual connections
        distTo = dist1(cV,conTo>0); %
        betterPath = ( dist2seed(conTo) > (dist2seed(cV)+distTo)); %find connections to which better path has been found
        dist2seed(conTo(betterPath)) = (dist2seed(cV)+1); % Change distances to improve
        pred(conTo(betterPath)) = cV; % mark improved position by predecessor.
        nextCon = cat(2,nextCon,conTo); %  add connected to list for next start verticies.
    end
    
    if ~sum(con2(:,1))
        break
    end
end

