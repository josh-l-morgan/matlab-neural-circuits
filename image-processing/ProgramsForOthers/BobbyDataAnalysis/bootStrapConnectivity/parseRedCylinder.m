


 %load('C:\Users\joshm\Documents\MATLAB\models\bootStrapConnectivity\redCylinder01.mat')

preIn = [];
postIn = [];

%rDat(:,1) = dat(:,1);
%rDat(:,2) = sum(dat(:,2:end),2);
rDat = dat;
numAx = size(dat,1);

for d  = 1:size(rDat,2)
    for a = 1:size(rDat,1) %run through all combinations
        occ = rDat(a,d);
        preIn(length(preIn)+1:length(preIn)+occ) = a;
        postIn(length(postIn)+1:length(postIn)+occ) = d;      
        
    end
end

cons.preIn = preIn;
cons.postIn = postIn;
