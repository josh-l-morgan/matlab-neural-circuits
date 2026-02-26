%% Enter data locations

load('..\matFiles\tifMap.mat')
load('..\matFiles\bottomLine.mat')

readSize = bottomLine.readSize;

rF = 0; % count read fails

errList = {};
for m = 1:length(tifMap.mon);
    
    
    tifNumber = length(tifMap.mon(m).tifNames);
    
    sumDat = find(~sum(bottomLine.mon(m).lineDat,2))
    
    for sD = 1:length(sumDat)
        %if ~tileID(t)
        t = sumDat(sD);
        
        if ~mod(t,100)
            sprintf('reading tile %d of %d.', sD, length(sumDat))
        end
        
            
            tifPath = [tifMap.mon(m).montageDir '\' tifMap.mon(m).sectionDirs{t} ...
                '\' tifMap.mon(m).tifNames{t}];
            
            for re = 1:3  %try reading 3 times
                readOK = 1;
                try  
                    Iline = imread(tifPath,'PixelRegion',{[readSize readSize] [1 readSize]});
                catch err
                    err
                    readOK = 0;
                    
                    if re ==3
                    errList{length(errList)+1,1} = err;
                    Iline = [0 0 0 0 ];
                    rF = rF+1;
                    readFail(rF) = t;
                    bottomLine.mon(m).readFail(t) = length(errList);
                    end
                end
                if readOK
                    break
                end
            end
            
            bottomLine.mon(m).lineDat(t,:)= ...
                [min(Iline) mean(Iline) std(double(Iline)) max(Iline)];
            bottomLine.mon(m).tileID(t) = tifMap.mon(m).tileIDs(t);
            %end
        end % run tiles
        
end
    
bottomLine.errList = errList;
save('..\matFiles\bottomLine.mat','bottomLine')
        
