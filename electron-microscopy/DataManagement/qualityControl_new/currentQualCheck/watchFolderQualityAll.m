clear all
%% get folder
TPN = GetMyDir;
if exist([TPN 'q.mat'])
    load([TPN 'q.mat'])
else
    q.checkedNams = {};
end

%% loop continuously
while 1
    
APN = findFolders(TPN)
    allNams = {}; allTimes = []; allFolders = {};
for f = 1: length(APN)
    [aNams aTimes]= getTiles(APN{f});
    fID = ones(length(aNams),1) * f;
    allFolders = cat(1,allFolders,APN(fID));
    allNams = cat(1,allNams, aNams);
    allTimes = cat(1,allTimes, aTimes);
end


%% check for new names
checkit = ones(length(allNams),1);
for i = 1:length(allNams)
    for c = 1:length(q.checkedNams)
        if strcmp(allNams{i},q.checkedNams{c})
            checkit(i) = 0;
            break
        end
    end
end
newNams = allNams(checkit>0);
newTimes = allTimes(checkit>0,:);
newFolders = allFolders(checkit>0,:);

%% Run quality check
pause(10)
for t = 1:length(newNams)
    L = length(q.checkedNams);
    sprintf('Running image %d of %d.',L+1,length(newNams))
    while 1  % 5 min after file write
        if etime(clock, newTimes(t,:)) >300
            break
        end
    end
    checkFile = [APN{newFolderIDs(t)} '\' newNams{t}];
    q.tile(L+1) = checkFileQual(checkFile);
    q.checkedNams{L+1,1} = newNams{t};
    q.checkedFolders{L+1,1} = newFolders{t};
    save([TPN 'q.mat'],'q');
    
end

%% Organize data
% 


% mos = parseTileNames(checkedNams);

 allQual = [q.tile.quality];
[sortQual ranks] = sort(allQual');%sort(tile.use.quality,'ascend')
worstNams = q.checkedNams(ranks);
worstFolders = q.checkedFolders(ranks);

showWorst = [worstFolders worstNams num2cell(sortQual)]
showTiles = [q.checkedFolders q.checkedNams num2cell(allQual')]
xlswrite([TPN 'worstTiles.xls'],showWorst,'showWorst')
xlswrite([TPN 'worstTiles.xls'],showTiles,'showTiles')




% 
% 
% qualMos = zeros(mos.Dims);
% for i = 1:length(mos.Files)
%     qualMos(mos.Pos(i,1),mos.Pos(i,2)) = tile(mos.Files(i)).quality;
% end
% 
% subplot(1,1,1)
% image(fitH(qualMos)),pause(.01)
% % 
% % qual.mos.qualMos = qualMos;
% % qual.mos.globMos = globMos;
% % qual.mos.satMos = satMos;
% qual.tile = tile;
% qual.nams = checkedNams;
% qual.worst = worstNams;
% qual.mos = mos;
% save([TPN 'qual.mat'],'qual')

'waiting for new file...'
pause(60)
end





