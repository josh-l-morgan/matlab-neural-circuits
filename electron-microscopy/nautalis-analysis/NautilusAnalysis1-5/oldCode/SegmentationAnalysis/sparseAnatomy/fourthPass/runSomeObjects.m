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
objectSubs = o{objNum};


connectSubs(objectSubs)