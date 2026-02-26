clear all
%%
delayCopySeconds = 60;  %Seconds between reading directories and starting file copy
delayCopyDays = delayCopySeconds/60/60/24;

secondsBetweenFullCheck = 60;
daysBetweenFullCheck = secondsBetweenFullCheck/60/60/24;
lastFullCheck = datenum(clock);
downscaledtargetsize = 2048;
downscaledfileextension ='_downscaled.tif';

%% get folder
WPN = GetMyDir;  %% Get source directory
TPN = GetMyDir;  %% Get target directory
DSPN = GetMyDir;  %% Get downsample directory
startFoldName1 = length(WPN) + 1;
%startFoldName2 = length(WPN)*2 + 1;


%% Read source folder
currentCheckTime =  datenum(clock);%record time at first read

APN = findFolders(WPN); % find all folders in source directory
allNams =  {}; allTimes = []; allFolders = {}; allSizes = [];  %initialize empty lists
for f = 1:length(APN);  %run through every found folder
    dAPN = dir(APN{f}); dAPN = dAPN(3:end);
    noDir = ~[dAPN.isdir];
    aFiles = dAPN(noDir);
     
    allNams = [allNams; {aFiles.name}'];
    allTimes = [allTimes; [aFiles.datenum]'];
    allSizes = [allSizes; [aFiles.bytes]'];
    fID = ones(length(aFiles),1) * f;
    allFolders = cat(1,allFolders,APN(fID));
end

%% Make directories
for i = 1:length(APN)
    source = APN{i};
    dest1 = [TPN source(startFoldName1:end)];
    dest2 = [DSPN source(startFoldName1:end)];
    if ~exist( dest1,'dir')
        mkdir(dest1);
        mkdir(dest2);
    end
end

%% Copy files

for i = 1:length(allNams)
    source = [allFolders{i} '\' allNams{i}];
    dest1 = [TPN source(startFoldName1:end)];
    dest2 = [DSPN source(startFoldName1:end)];
    if ~exist(dest1,'file')
        if  (allTimes(i) - delayCopyDays)>currentCheckTime;
            'waiting for delay time to pass'
            pause(delayCopySeconds)  %wait for delay time to pass
        end
        status = 0;
        while status == 0  %make sure copy succeded
            status = copyfile(source,dest1);
            tiffile =strfind(source, 'tif');
            if tiffile ~=0
            image=imread(source);
            dimage=imresize(image,[downscaledtargetsize downscaledtargetsize]);
            filenamewithpath = source;
            pos=size(filenamewithpath,2); 
                while (pos>0)&&(filenamewithpath(pos)~='\')
                    pos=pos-1;
                end;
            filename=filenamewithpath(pos+1:end);
            targetname=[dest2 downscaledfileextension];
            %targetname=[dest2 filename(1:end-4) downscaledfileextension]
            imwrite(dimage,targetname);
            sprintf('  Writing downscaled image to %s',targetname)
            end
          end
        end
    end


%% loop continuously
while 1
    pause(1)
    %Update after first run through
    lastCheckTime = currentCheckTime; %record time at first read
    currentCheckTime =  datenum(clock);
    noNew = 1;
    %% Read source folder
    APN = findFolders(WPN); % find all folders in source directory
    allNams =  {}; allTimes = []; allFolders = {}; allSizes = [];  %initialize empty lists
    for f = 1: length(APN)  %run through every found folder
        dAPN = dir(APN{f}); dAPN = dAPN(3:end);
        noDir = ~[dAPN.isdir];
        aFiles = dAPN(noDir);
        
        allNams = [allNams; {aFiles.name}'];
        allTimes = [allTimes; [aFiles.datenum]'];
        allSizes = [allSizes; [aFiles.bytes]'];
        fID = ones(length(aFiles),1) * f;
        allFolders = cat(1,allFolders,APN(fID));
    end
    
    %% Make directories
    for i = 1:length(APN)
        source = APN{i};
        dest1 = [TPN source(startFoldName1:end)];
        dest2 = [TPN source(startFoldName1:end)];
        if ~exist( dest1,'dir')
            mkdir(dest1);
        end
        
        if ~exist( dest2,'dir')
            mkdir(dest2);
        end
        noNew = 0;
    end
    
    %% Copy new files
    copyNew = find(allTimes>lastCheckTime);
    for i = 1:length(copyNew)
        source = [allFolders{copyNew(i)} '\' allNams{copyNew(i)}];
        dest1 = [TPN source(startFoldName1:end)];
         dest2 = [DSPN source(startFoldName1:end)];
        if ~exist(dest1,'file')
            if (allTimes(copyNew(i)) + delayCopyDays)>datenum(clock)
                pause((allTimes(copyNew(i)) + delayCopyDays-datenum(clock))/24/60/60);
            end
                status = 0;
            while status == 0  %make sure copy succeded
                status = copyfile(source,dest1);
                sprintf('  Writing image to %s',dest1) %disp(txt);
                tiffile =strfind(source, 'tif');
                if tiffile ~=0
                    image=imread(source);
                    dimage=imresize(image,[downscaledtargetsize downscaledtargetsize]);
                    filenamewithpath = source;
                    pos=size(filenamewithpath,2);
                while (pos>0)&&(filenamewithpath(pos)~='\')
                    pos=pos-1;
                end;
                    filename=filenamewithpath(pos+1:end);
                    targetname=[dest2 filename(1:end-4) downscaledfileextension];
                     imwrite(dimage,targetname);
                     sprintf('  Writing downscaled image to %s',targetname) %disp(txt);
                end;
            end;
        end;
    end;
    
    %% Run full check
    if ( datenum(clock)-lastFullCheck )> daysBetweenFullCheck
        tic
        %'running full check'
        lastFullCheck = datenum(clock);
        
        for i = 1:length(allNams)
            source = [allFolders{i} '\' allNams{i}];
            dest1 = [TPN source(startFoldName1:end)];
            dest2 = [DSPN source(startFoldName1:end)];
            shouldCopy = 0;
            if ~exist(dest1,'file')
                shouldCopy = 1;
            else  %if file exists
                fileInfo = dir(dest1);
                fileInfo = fileInfo(end);
                if fileInfo.datenum<allTimes(i) %replace old file
                    shouldCopy = 1;
                    sprintf('replaced old file %s', allNams{i})
                elseif (fileInfo.datenum == allTimes(i) ) & ...
                        (fileInfo.bytes < allSizes(i)); %replace partial file
                    sprintf('replaced partial file %s', allNams{i})
                end % if overwrite
                
            end %whether dest file exists
            
            if shouldCopy %copy if you should
                status = 0;
                noNew = 1;
                while status == 0  %make sure copy succeded
                    status = copyfile(source,dest1);
                    sprintf('  Writing image to %s',dest1)% disp(txt);
                    tiffile =strfind(source, 'tif');
                    if tiffile ~=0
                        image=imread(source);
                        dimage=imresize(image,[downscaledtargetsize downscaledtargetsize]);
                        filenamewithpath = source;
                        pos=size(filenamewithpath,2);
                        while (pos>0)&&(filenamewithpath(pos)~='\')
                            pos=pos-1;
                        end;
                        filename=filenamewithpath(pos+1:end);
                        targetname=[dest2 filename(1:end-4) downscaledfileextension];
                        imwrite(dimage,targetname);
                        sprintf('  Writing downscaled image to %s',targetname) %disp(txt);
                    end;
                end;
            end;       %run all file names
            %toc
            if noNew
                sprintf('No new files found after full check at %s', datestr(clock))
            else
                sprintf('New files found after full check at %s', datestr(clock))
            end
        end   %if ready for full check
    end



end

