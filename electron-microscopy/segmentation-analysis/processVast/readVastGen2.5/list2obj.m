function[] = list2obj(SPN)


TPN = [SPN(1:end-1) '_mat\'];
load([TPN 'pointList.mat']);

uniqueIds = unique(pointList.ids);
countIds = hist(pointList.ids,uniqueIds);
