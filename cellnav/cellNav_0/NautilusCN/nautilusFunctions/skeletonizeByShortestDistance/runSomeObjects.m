clear all

SPN = 'D:\LGNs1\segmentation\VAST\S8prelim\export4Mat\'
TPN = [SPN 'cell_1\'];
if ~exist(TPN,'dir'),mkdir(TPN),end

%% Load object data
load([SPN 'obj.mat'])

%% Pick an object
for i = 1:length(obj.subs)
    objSize(i,:) = size(obj.subs{i});
end

objNum = find(objSize(:,1) == max(objSize(:,1)));
objectSubs = obj.subs{objNum};
clear obj

%connectSubs(objectSubs)