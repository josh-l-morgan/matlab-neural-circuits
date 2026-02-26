function[tipness maxV] = findTips(conMat,dist2seed,seed);
minTipSize = 40;


vNum = size(conMat,1);
maxV = 1:vNum;


startDist = [dist2seed + rand(size(dist2seed)) 0];
startDist = 1:length(startDist);
maxDist = startDist;

con2 = conMat;
con2(con2==0) = length(startDist);


for r = 1:minTipSize
    for i = 1:26
        neighborList = con2(:,i);
        neighborDist = startDist(neighborList);
        newMax = find(((neighborDist - maxDist(1:end-1)) >0) & (maxDist(1:end-1)>0));
        maxDist(newMax) = neighborDist(newMax);
        maxV(newMax) = neighborList(newMax);
    end
    disp(sprintf('%02.1f%% finished',r/minTipSize*100)),pause(.001)
    
end


tipness = hist(maxV,1:vNum);


% 
% 
% %%
% minTipSize = 10;
% vNum = size(conMat,1);
% 
% startDist = dist2seed + rand(size(dist2seed));
% maxDist = startDist;
% maxV = ones(vNum,1)*seed;
% con2 = conMat;
% con2(con2==0) = seed;
% 
% for r = 1:minTipSize
%     for i = 1:26
%         neighborList = con2(:,i);
%         neighborDist = startDist(neighborList);
%         newMax = (neighborDist - maxDist) >0;
%         maxDist(newMax) = neighborDist(newMax);
%         maxV(newMax) = neighborList(newMax);
%     end
%     disp(sprintf('%02.1f%% finished',r/minTipSize*100)),pause(.001)
%     
% end
% 
% %%Surf vox
% % deepVox = sum(conMat(:,1:6)>0,2)>5;
% % maxV(deepVox) = seed;
% tipness = hist(maxV,1:vNum);
% 
