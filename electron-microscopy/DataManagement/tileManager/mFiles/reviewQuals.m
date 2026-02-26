
if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 7
end


% load('..\matFiles\chkQ.mat')
% load('..\matFiles\logQuals.mat')
% load('..\matFiles\noLogQuals.mat')
% load('..\matFiles\bottomLine.mat')
% load('..\matFiles\tifMap.mat')
% load('..\matFiles\utslTiles.mat')


%% Make list of most recent tiles for use

IDs = utslTiles.tiles.ID; % get ID list from Ultrathin section library
recent = cell(length(IDs),1);
parfor i = 1:length(IDs)
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
        
    elseif length(tileIndex) ==1
       recent{i,1} = [tifMap.mon(ms).montageDir '\'...
           tifMap.mon(ms).sectionDirs{tileIndex} '\'...
           tifMap.mon(ms).tifNames{tileIndex}];
    else 
        
        targ = find(times == max(times),1); %find most recent
              recent{i,1} = [tifMap.mon(ms(targ)).montageDir '\'...
           tifMap.mon(ms(targ)).sectionDirs{tileIndex(targ)} '\'...
           tifMap.mon(ms(targ)).tifNames{tileIndex(targ)}];
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


%% Decide between retakes

bestQual = recent; %start list with most recent files.
useTif.tiles
for d = 1:length(chkQ.do)
    tile = chkQ.tiles(d);
    for t = 1:length(tile.tif)
        folders{t,1} = [tile.tif(t).montageDir '\' tile.tif(t).sectionDirs];
    end   
    uniqueFolders = unique(folders);
    eval.folders = uniqueFolders;
    for uF = 1:length(uniqueFolders)
        sec.folders
        
    end
        
    
    if isempty(tif)
        
    end
        
        
    
end

