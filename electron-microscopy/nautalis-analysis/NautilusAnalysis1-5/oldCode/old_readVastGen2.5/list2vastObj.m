function[] = list2vastObj(TPN)

load([TPN 'pointList.mat']);

uniqueIds = unique(pointList.ids);
countIds = hist(pointList.ids,uniqueIds);
idStop = cumsum(countIds);
idStart = [0 idStop(1:end-1)]+1;

pointList = pointList; % lets parfor know pontList is a variable

parfor i = 1:length(uniqueIds);
    
    subs{i} = cat(2,pointList.Y(idStart(i):idStop(i)),...
        pointList.X(idStart(i):idStop(i)), ...
        pointList.Z(idStart(i):idStop(i)));
    
end

vastOb.uniqueIds = uniqueIds;
vastOb.countIds = countIds;
vastOb.size = [max(pointList.Y) max(pointList.X) min(pointList.Z)];    
vastOb.subs = subs;

save([TPN 'vastOb.mat'],'vastOb')