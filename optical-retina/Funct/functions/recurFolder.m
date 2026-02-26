function[listF] = recurFolder(F,listF)
%%CreatesList of  all folders withing a folder. 


%% Add folder to list
if nargin == 0 
    F = GetMyDir;
    F = F(1:length(F)-1);
    listF = {F};
elseif nargin == 1
    listF = {F};
else
    listF{length(listF)+1,1} = F;
end

currentDirectory = F

%% Search for new folders

dirF = dir(F); dirF = dirF(3:length(dirF));

dirFc = struct2cell(dirF);
dirNames = dirFc(1,:);
realDirs = [dirFc{4,:}];
dirNames = dirNames(realDirs);
for i = 1:length(dirNames)
     name = [F '\' dirNames{i}];
     listF = recurFolder(name,listF);
end