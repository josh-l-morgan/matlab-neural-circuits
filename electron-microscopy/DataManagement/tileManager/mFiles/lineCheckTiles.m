%% Enter data locations

load('..\matFiles\tifMap.mat')
load('..\matFiles\bottomLine.mat')

readSize = 25000;
bottomLine.readSize = readSize;

rF = 0; % count read fails
for m = 2:length(tifMap.mon);
    bottomLine.mon(m).date = datestr(clock);
    
    tifNumber = length(tifMap.mon(m).tifNames);
    
    if ~isfield(bottomLine.mon(m),'tileIDs');
        bottomLine.mon(m).lineDat = zeros(tifNumber, 4);
        bottomLine.mon(m).tileID = zeros(tifNumber,1); 
    end
    
    for t = 1:tifNumber
        %if ~tileID(t)
        if ~mod(t,100)
        sprintf('reading tile %d of %d.', t, tifNumber)
        end
        tifPath = [tifMap.mon(m).montageDir '\' tifMap.mon(m).sectionDirs{t} ...
            '\' tifMap.mon(m).tifNames{t}];
        try  Iline = imread(tifPath,'PixelRegion',{[readSize readSize] [1 readSize]});
        catch err
           Iline = 0; 
           rF = rF+1;
           readFail(rF) = t;
        end
            
        bottomLine.mon(m).lineDat(t,:)= ...
            [min(Iline) mean(Iline) std(double(Iline)) max(Iline)];
        bottomLine.mon(m).tileID(t) = tifMap.mon(m).tileIDs(t);
        %end
    end % run tiles
    save('..\matFiles\bottomLine.mat','bottomLine')
    
end