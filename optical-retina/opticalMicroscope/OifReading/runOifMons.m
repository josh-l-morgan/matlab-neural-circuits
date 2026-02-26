

TPN = GetMyDir;
folds = findFolders(TPN);
monList = {}
for i = 1:length(folds)
    dFolds = dir(folds{i});
    for d = 1:length(dFolds)
        mos = regexp(dFolds(d).name,'Mosaic.log');
        if ~isempty(mos)
            monList{length(monList)+1} = folds{i};
            break
        end
    end
end
save([TPN 'monList.mat'],'monList')

for i = 1:length(monList)
    OifMonFun([monList{i} '\']);
end