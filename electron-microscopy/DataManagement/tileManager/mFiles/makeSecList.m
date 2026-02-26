
% if matlabpool('size') == 0 % checking to see if my pool is already open
%     matlabpool open 7
% end


load('..\matFiles\chkQ.mat')
load('..\matFiles\logQuals.mat')
load('..\matFiles\noLogQuals.mat')
load('..\matFiles\bottomLine.mat')
load('..\matFiles\tifMap.mat')
load('..\matFiles\utslTiles.mat')


%% Make list of most recent tiles for use

IDs = utslTiles.tiles.ID; % get ID list from Ultrathin section library
recent = cell(length(IDs),1);
missingTiles = {};
for i = 1:length(IDs)
    if mod(i,10000) == 0
        disp(sprintf('Checking tile %d of %d.',i,length(IDs)))
        pause(.01)
    end
    ID = IDs(i);
    ms = []; times = []; tileIndex = [];
    for m = 1:length(tifMap.mon)
        targs = find(tifMap.mon(m).tileIDs == ID);
        if ~isempty(targs)
            ms = cat(1,ms, ones(length(targs),1)*m);
            times = cat(1,times, [tifMap.mon(m).tifTimes(targs)]);
            tileIndex = cat(1,tileIndex,targs);
        end
    end
    if isempty(tileIndex)
        %disp(ID2tile(ID));
        missingTiles{length(missingTiles)+1,1} = ID2tile(ID);
        
    elseif length(tileIndex) ==1
        recent{i,1} = [tifMap.mon(ms).montageDir '\'...
            tifMap.mon(ms).sectionDirs{tileIndex} '\'...
            tifMap.mon(ms).tifNames{tileIndex}];
        if isempty(recent{i,1}),pause,end
    else
        
        targ = find(times == max(times),1); %find most recent
        recent{i,1} = [tifMap.mon(ms(targ)).montageDir '\'...
            tifMap.mon(ms(targ)).sectionDirs{tileIndex(targ)} '\'...
            tifMap.mon(ms(targ)).tifNames{tileIndex(targ)}];
                if isempty(recent{i,1}),pause,end

    end
end


%% Build section list
countSecs =0;
clear us
us.date = datestr(clock);
for w = 1:length(utslTiles.waf)
    for s = 1:length(utslTiles.waf(w).sec)
        countSecs = countSecs + 1;
        us.sec(countSecs).wafNum = w;
        us.sec(countSecs).wafSec = s;
        us.w(w).s(s) = countSecs;
        for r = 1:4
            for c = 1:4
                us.sec(countSecs).IDs(r,c) = sub2ID(w,s,r,c);
            end
        end
    end
end




%% convert IDlist to subs
subs = ID2subs(IDs);

%% Enter recent path into sec list
for i = 1:length(IDs)
    subs = ID2subs(IDs(i));
    secNum = us.w(subs(1)).s(subs(2));
    us.sec(secNum).eval.path2recent{subs(3),subs(4)} = recent{i};
    us.sec(secNum).paths{subs(3),subs(4)} = recent{i};
end



%% collect all directories for each section
allSecFolders = cell(length(us.sec),1);
for m = 1:length(tifMap.mon)
    subs = tifMap.mon(m).subs;
    for i = 1:length(tifMap.mon(m).subs)
        try
            secNum = us.w(subs(i,1)).s(subs(i,2));
        catch err
            secNum = 0;
        end
        if secNum
            paths = allSecFolders{secNum};
            path = [tifMap.mon(m).montageDir '\' ...
                tifMap.mon(m).sectionDirs{i}];
            paths{length(paths)+1} = path;
            allSecFolders{secNum} = paths;
        end
    end
end


for s = 1:length(allSecFolders)
    allFolders = allSecFolders{s}';
    us.sec(s).folders = unique(allFolders);
    
end

%% Enter tifData inot us sec list (including bottom Line)
for m = 1:length(tifMap.mon)
    subs = tifMap.mon(m).subs;
    for i = 1:length(tifMap.mon(m).subs)
        try
            secNum = us.w(subs(i,1)).s(subs(i,2));
        catch err
            secNum = 0;
        end
        if secNum
            paths = us.sec(secNum).folders;
            path = [tifMap.mon(m).montageDir '\' ...
                tifMap.mon(m).sectionDirs{i}];
            pathNum = find(strcmp(paths,path));
            r = tifMap.mon(m).subs(i,3);
            c = tifMap.mon(m).subs(i,4);
            try
                tifCount = length(us.sec(secNum).secDir(pathNum).tifDat(r,c).tif)+1;
            catch err
                tifCount = 1;
            end
            us.sec(secNum).secDir(pathNum).tifDat(r,c).tif(tifCount).name = tifMap.mon(m).tifNames{i};
            us.sec(secNum).secDir(pathNum).tifDat(r,c).tif(tifCount).time = tifMap.mon(m).tifTimes(i);
            us.sec(secNum).secDir(pathNum).tifDat(r,c).tif(tifCount).ID = tifMap.mon(m).tileIDs(i);
            us.sec(secNum).secDir(pathNum).tifDat(r,c).tif(tifCount).bottomLine = ...
                bottomLine.mon(m).lineDat(i,:);
        end
    end
end

%% Set section check status

for s = 1:length(us.sec)
    us.sec(s).needCheck = 0;
end


%% Enter quality values from chkQ into us

for d = 1:length(chkQ.do)
    tile = chkQ.tiles(d);
    for t = 1:length(tile.tif)
        folders{t,1} = [tile.tif(t).montageDir '\' tile.tif(t).sectionDirs];
    end
    
end


%% Enter quality data into us
for d = 1:length(chkQ.do)
    tile = chkQ.tiles(d);
    subs = ID2subs(tile.tileID);
    secNum = us.w(subs(1)).s(subs(2));
    r = subs(3);
    c = subs(4);
    us.sec(secNum).needCheck = 1;
    
    folders = us.sec(secNum).folders;
    for t = 1:length(tile.tif)
        folder = [tile.tif(t).montageDir '\' tile.tif(t).sectionDirs];
        targFold = find(strcmp(folders,folder))
        
        
        tifNames = {us.sec(secNum).secDir(targFold).tifDat(r,c).tif.name};
        tifName = tile.tif(t).tifName;
        targTif = find(strcmp(tifNames,tifName));
        
        us.sec(secNum).secDir(targFold).tifDat(r,c).tif(targTif(1)).quality ...
            = tile.tif(t).qual.quality;
    end
end


%% Enter nolog data into us

for d = 1:length(noLogQuals.tiles)
    tile = noLogQuals.tiles(d);
    subs = ID2subs(tile.tileIDs);
    secNum = us.w(subs(1)).s(subs(2));
    r = subs(3);
    c = subs(4);
    us.sec(secNum).needCheck = 1;
    
    folders = us.sec(secNum).folders;
    for t = 1:length(tile.tif)
        folder = [tile.tif(t).montageDir '\' tile.tif(t).sectionDirs];
        targFold = find(strcmp(folders,folder));
        
        
        tifNames = {us.sec(secNum).secDir(targFold).tifDat(r,c).tif.name};
        tifName = tile.tif(t).tifName;
        targTif = find(strcmp(tifNames,tifName));
        
        us.sec(secNum).secDir(targFold).tifDat(r,c).tif(targTif(1)).quality ...
            = tile.tif(t).qual.quality;
    end
end



%% Record log bood data in Section list us

for s = 1:length(us.sec)
    us.sec(s).eval.bestQuals = zeros(4,4);
    us.sec(s).eval.qualKnown = 0;
end

for w = 1:length(logQuals.waf)
    if logQuals.waf(w).foundBook
        for s = 1:length(logQuals.waf(w).sec)
            secNum = us.w(w).s(s);
            for t = 1:length(logQuals.waf(w).sec(s).tile);
                [ID waf sec r c] = tile2ID(logQuals.waf(w).sec(s).tile(t).tileName);
                us.sec(secNum).eval.bestQuals(r,c) = ...
                    logQuals.waf(w).sec(s).bestQuals(t);
            end
            
        end
    end
end

save('..\matFiles\us.mat','us')


%% Evaluate sections

for s = 1 : length(us.sec)
    if ~mod(s,100)
        disp(sprintf('checking section %d of %d',s,length(us.sec)))
    end
    
    for f = 1:length(us.sec(s).folders)
        useTif = zeros(4,4);
        tifOK = zeros(4,4);
        useQual = zeros(4,4);
        foldName = us.sec(s).folders{f};
        for r = 1:4
            for c = 1:4
                
                try
                    tifs = us.sec(s).secDir(f).tifDat(r,c).tif;
                catch err
                    tifs = [];
                end
                
                if isempty(tifs)
                    useTif(r,c) = 0;
                    tifOK(r,c) = 0;
                    maxTimes(r,c) = 0;
                else
                    
                    %% grab data
                    bL = cat(1,tifs.bottomLine);
                    times = cat(1,tifs.time);
                    try
                        q = [tifs.quality]' ;
                    catch err
                        q = [];
                    end
                    clear pass
                    pass(1:size(bL,1),1) = (bL(:,4)>0);
                    
                    maxTimes(r,c) = max(times); %record latest tif time
                    
                    
                    if isempty(q) %if there is no quality
                        recent = find(times==max(times));
                        if pass(recent(1))
                            useTif(r,c) = recent(1);
                            tifOK(r,c) = 1;
                            useName{r,c} = [foldName '\' tifs(recent(1)).name];
                        else
                            bestLeft = find((times.*pass) == max(times.*pass));
                            useTif(r,c) = bestLeft(1);
                            tifOK(r,c) = 0;
                            useName{r,c} = [foldName '\' tifs(bestLeft(1)).name];
                            us.sec(s).needCheck = 1;
                        end
                    elseif length(q) == length(pass)%quality for every tif
                        best = find(q == max(q));
                        if pass(best)
                            useTif(r,c) = best(1);
                            tifOK(r,c) = 1;
                            useQual(r,c) = q(best(1));
                            useName{r,c} = [foldName '\' tifs(best(1)).name];
                        else
                            bestLeft = find((q.*pass) == max(q.*pass));
                            useTif(r,c) = bestLeft(1);
                            tifOK(r,c) = 0;
                            useQual(r,c) = q(bestLeft(1));
                            useName{r,c} = [foldName '\' tifs(bestLeft(1)).name];
                        end
                        
                    else %something went wrong
                        'Damn, missing data'
                        pause
                    end
                    
                end
                
            end
            
        end
        us.sec(s).secDir(f).tifOK = tifOK;
        us.sec(s).secDir(f).useTif = useTif;
        us.sec(s).secDir(f).useQual = useQual;
        us.sec(s).secDir(f).useName = useName;
        us.sec(s).secDir(f).maxTimes = maxTimes;
        
    end %run each folder
end




%% Choose the best section folder

for s = 1:length(us.sec);
    us.sec(s).eval.checkByEye = 0;
end

for s = 1:length(us.sec)
    if us.sec(s).needCheck;
        foldNum = length(us.sec(s).folders);
        
        if foldNum==1
            bestFold = 1;
            us.sec(s).paths =  us.sec(s).secDir(bestFold).useName;
            us.sec(s).eval.bestQuals = us.sec(s).secDir(bestFold).useQual;
            us.sec(s).eval.tifOK = us.sec(s).secDir(bestFold).tifOK;
        elseif foldNum>1
            us.sec(s).eval.checkByEye = 1;
            
            allQual = zeros(4,4,foldNum);
            startTime = zeros(foldNum,1);
            for f = 1:foldNum
                startTime(f) = max(us.sec(s).secDir(f).maxTimes(:));
                allQual(:,:,f) = us.sec(s).secDir(f).useQual;
            end
            bestFold = find(startTime == max(startTime),1);
            
            if sum(sum(allQual(:,:,1)>0))>7
                wins = []
                for r = 1:4
                    for c = 1:4
                        if sum(allQual(r,c,:))
                            wins = cat(1,wins,...
                                find(allQual(r,c,:) == max(allQual(r,c,:))));
                        end
                    end
                end
                countWins = hist(wins,1:1:max(wins));
                [sortWins order] = sort(countWins,'descend');
                if (sortWins(1)-sortWins(2) )>7
                    bestFold = order(1);  %assign winning folder to bestFold
                end
                
            end
            
             us.sec(s).secDir(bestFold).useName;
            us.sec(s).eval.bestQuals = us.sec(s).secDir(bestFold).useQual;
            us.sec(s).eval.tifOK = us.sec(s).secDir(bestFold).tifOK;
        end
            
    end % in needs check
end

save('..\matFiles\us.mat','us')

%% Review sections






