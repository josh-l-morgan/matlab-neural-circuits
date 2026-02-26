%%Map disk usage within specific folders and keep track of tags

%% Check folder
SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\'
SPN = 'V:\'

%% ignore folders
skipPN = {"\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\LGN_Adult\LGNs1\"}

%% Save progress in folder
TPN = 'G:\FolderMaps\'; %SPN;
TFN = 'folderMap_C.mat';
disp(sprintf('Mapping directory %s',SPN))

%% Options
tierMaxOne = 10;
maxSubDir = 10; %maximum number of subdirectories a folder can have before it is declared a trunk
maxTB = 100000; % Stop checking path once maxTB is found
saveEveryXMinutes = 10;
shortReportLength = 20;
reportEveryXMinutes = 1/6;
lastReport = datenum(clock);
lastSave = datenum(clock);

% 
% %%Check target folder
% placeHolder = [];
% for i = 1:1000
%     TFN = sprintf('folderMap%04.0f.mat',i);
%     if ~exist([TPN TFN],'file')
%         % save([TPN TFN], 'placeHolder','-V7.3')
%         break
%     end
% end
% disp(sprintf('File name is %s',TFN))
% 


if 0
    disp('Loading previous map')
    load([TPN TFN]);
    nextPath = dat.nextPath;
    paths = dat.paths(1:nextPath,1);
    folders = dat.folders(1:nextPath,1);
    totBytes = dat.totBytes(1:nextPath,1);
    parent = dat.parent(1:nextPath,1);
    tier = dat.tier;
    containsCappedPaths = dat.containsCappedPaths;

    c = dat.c;
    usePaths = dat.usePaths;
    startCheckingPaths = dat.lastPathChecked + 1;

elseif 1 %initialize
    paths = {SPN};
    folders = {SPN};
    totBytes = 0;
    nextPath = 1;
    tier = 1;
    parent = 0;
    c = 0;

    arrayPad = 1000;
    usePaths = [];
    %% find directory
    disp('Finding trunk directories')
    while 1

        c = c + 1;
        if c>size(paths,1), break, end % Stop searching when c is larger than path list
        if isempty(paths{c,1}),break,end % Dont search if path is bland. why?
        isTrunk = 0; %initialize record of trunkness
        if tier(c) == tierMaxOne
            isTrunk = 1;
        end

        %%Get directory information
        d =  dir(paths{c,1});
        newDirID = setdiff(find([d.isdir]),[1 2]); %%Get directories
        newDirs = {d(newDirID).name};
        numSub = length(newDirs);


        %%Collect all bytes from non directory files and add them to directory
        newFileID = find(~[d.isdir]);
        totBytes(c,1) = sum([d(newFileID).bytes]); %record non-directory bytes in file

        %% For each new directory
        if (numSub==0) | (numSub>maxSubDir) %if dead end or mipmap folder
            isTrunk = 1;
        end

        if isTrunk
            usePaths = cat(1,usePaths,c); %record path{c} as trunk
        else % or create space for found sub directories as paths
            for f = 1:length(newDirs);
                dirName = newDirs{f};
                if ~strcmp(dirName,'.') & ~strcmp(dirName,'..')
                    if ~sum(strcmp(dirName,skipPN)) %% ignore some directories
                        nextPath = nextPath+1;

                        %%Padd arrays if necessary
                        if (nextPath + arrayPad) > size(paths,1)
                            paths{nextPath + arrayPad,1} = '';
                            folders{nextPath + arrayPad,1} = '';
                            totBytes(nextPath + arrayPad,1) = 0;
                            tier(nextPath + arrayPad,1) = 0;
                            parent(nextPath + arrayPad,1) = 0;
                        end

                        %%Add new information
                        paths{nextPath,1} = [paths{c,1} newDirs{f} '\'];
                        tier(nextPath,1) = tier(c)+1;
                        parent(nextPath,1) = c;
                        folders{nextPath,1} = newDirs{f};
                        totBytes(nextPath,1) = 0;
                    end
                end
            end
        end

    end

    %%Cut padding
    paths = paths(1:nextPath,1);
    folders = folders(1:nextPath,1);
    tier = tier(1:nextPath,1);
    totBytes = totBytes(1:nextPath,1);
    parent = parent(1:nextPath,1);
    startCheckingPaths = 1; %first index of usePaths to check
    containsCappedPaths = zeros(nextPath,1);
    toc
end

%% Individual trunks
numPaths = length(paths);
for ti = startCheckingPaths:length(usePaths)
    trunkID = usePaths(ti);
    disp(sprintf('mapping %s',paths{trunkID}))
    sub_paths = paths(trunkID);

    c = 0;
    totBytes(trunkID) = 0; % clear previous data for trunk
    arrayPad = 10000;
    nextPath = length(sub_paths);
    while 1

        c = c + 1;
        if c>size(sub_paths,1), break, end
        if totBytes(trunkID,1) > (maxTB * 10^12)
            disp(sprintf('search of path %s capped at %0.2fTB',...
                sub_paths{c},totBytes(trunkID,1)/(10^12)))
            containsCappedPaths(trunkID) = containsCappedPaths(trunkID) + 1;
            break, end % if you have already found the max TB
        if isempty(sub_paths{c,1}),break,end

        %%Get directory information
        d =  dir(sub_paths{c,1});
        newDirID = setdiff(find([d.isdir]),[1 2]); %%Get directories
        newDirs = {d(newDirID).name};

        %%Collect all bytes from non directory files and add them to directory
        newFileID = find(~[d.isdir]);
        totBytes(trunkID,1) = totBytes(trunkID,1) + sum([d(newFileID).bytes]);

        %% For each new directory
        for f = 1:length(newDirs)
            dirName = newDirs{f};
            if ~strcmp(dirName,'.') & ~strcmp(dirName,'..')
                if ~sum(strcmp(dirName,skipPN)) %% ignore some directories

                    nextPath = nextPath+1;
                    %%Padd arrays if necessary
                    if (nextPath + arrayPad) > size(sub_paths,1)
                        sub_paths{nextPath + arrayPad,1} = '';
                    end
                    %%Add new information
                    sub_paths{nextPath,1} = [sub_paths{c,1} newDirs{f} '\'];
                end
            end
        end

        if ((datenum(clock) - lastReport) * 24 * 60) > reportEveryXMinutes
            disp(sprintf('%s',paths{trunkID}))
            disp(sprintf('Path %d of %d. Found %0.6f tb in %s',...
                ti,length(usePaths),totBytes(trunkID,1)/(10^12),paths{trunkID}(length(SPN):end)))
            lastReport = datenum(clock);
        end
    end


    if 1;%((datenum(clock) - lastReport) * 24 * 60) > reportEveryXMinutes
        disp(sprintf('%s',paths{trunkID}))
        disp(sprintf('Path %d of %d. Found %0.6f tb in %s',...
            ti,length(usePaths),totBytes(trunkID,1)/(10^12),paths{trunkID}(length(SPN):end)))

        disp(datetime(clock))
        startSum = datenum(clock);
        dat.SPN = SPN;
        dat.usePaths = usePaths;
        dat.lastPathChecked = ti;
        dat.nextPath = numPaths;
        dat.paths = paths;
        dat.folders = folders;
        dat.totBytes = totBytes;
        dat.parent = parent;
        dat.c = numPaths;
        dat.tier = tier;
        dat.containsCappedPaths = containsCappedPaths;
        report = reportMap(dat);
        listReport = report.listReport;
        shortReport = listReport(1:min(numPaths,shortReportLength),[1 2 3]);
        disp(shortReport)
        pause(.01)
        stopSum = datenum(clock);
        reportTimeMin = (stopSum-startSum) * 24 * 60;
        disp(sprintf('Report took %0.2f minutes',reportTimeMin))
        if reportTimeMin > (saveEveryXMinutes * .2)
            saveEveryXMinutes = saveEveryXMinutes * 2;
            disp(sprintf('Now reporting every %0.2f minutes',saveEveryXMinutes))
        end


        disp('Checking folders...')
    end

    if 1;%((datenum(clock) - lastSave) * 24 * 60) > saveEveryXMinutes
        disp('saving...')
        try
            save([TPN TFN], 'listReport', 'dat','-V7.3')
        catch err
            disp('failed to save')
        end
        lastSave = datenum(clock);
    end

end



%% Final report and save
tic
dat.nextPath = numPaths;
dat.paths = paths;
dat.folders = folders;
dat.totBytes = totBytes;
dat.parent = parent;
dat.c = numPaths;
dat.tier = tier;
report = reportMap(dat);
listReport = report.listReport;
disp(sprintf('checking directory number %d, biggest paths:',c))
shortReport = listReport(1:min(numPaths,30),:);
disp(shortReport)
pause(.01)
toc
tic
disp('saving...')
save([TPN TFN], 'listReport', 'dat','-V7.3')
toc

mediumReport = listReport(1:min(numPaths,100000),:);
disp(sprintf('finished mapping %s',SPN))