clear all
%%
delayCopy = 1;  %Seconds between reading directories and starting file copy

%% get folder
WPN = GetMyDir;  %% Get source directory
TPN = GetMyDir;  %% Get target directory


%% loop continuously
while 1
    
    lastCheckTime = datenum(clock);
    
    %% Read source folder
    APN = findFolders(WPN); % find all folders in source directory
    allNams =  {}; allTimes = []; allFolders = {};  %initialize empty lists
    for f = 1: length(APN)  %run through every found folder
        dAPN = dir(APN{f}); dAPN = dAPN(3:end);
        noDir = ~[dAPN.isdir];
        aFiles = dAPN(noDir);
        
        allNams = [allNams; {aFiles.name}'];
        allTimes = [allTimes; {aFiles.datenum}'];
        fID = ones(length(aFiles),1) * f;
        allFolders = cat(1,allFolders,APN(fID));
    end
    
    
    %% Delay copy
    pause(delayCopy)
    
    %% Write files with times after last check time
    
    
    
    %%
        
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
            
            
        end
        
        
        'waiting for new file...'
        pause(60)
    end
    
    
    
    
    
