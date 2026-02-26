%load('..\matFiles\us.mat')

TPN = GetMyDir
iRange = [13000 13500];
duration = zeros(1:length(us.sec),1);
for s = 1:length(us.sec)
    disp(sprintf('Sampling %d of %d.',s,length(us.sec)))
    paths = {};
    try
        paths = us.sec(s).paths;
        wafNum = us.sec(s).wafNum;
        secNum = us.sec(s).wafSec;
    catch err
        err
    end
    
    for r = 1:4
        for c = 1:4
            newName = sprintf('s%06.0f_w%03.0f_o%03.0f_r%1.0f_c%1.0f.tif',...
                        s,wafNum,secNum,r,c);
            newPath = [TPN newName];
            if ~exist(newPath,'file')
                
                startTime = datenum(clock);
                
                imagePath = paths{r,c};
                if isempty(imagePath)
                    I = zeros(iRange,'uint8');
                elseif ~exist(imagePath,'file') %check to see if source image exists
                    I = zeros(iRange,'uint8')+150;
                else
                    
                    for p = 1:3
                        pass = 1;
                        try
                            tic
                            I = imread(imagePath,'PixelRegion',{iRange, iRange});
                            toc
                        catch err
                            pass = 0;
                            err
                        end
                        if pass
                            break
                        end
                    end
                    
                    if ~pass
                        I = zeros(iRange,'uint8')+300;
                    end
                    
                    for p = 1:3
                        pass = 1;
                        try
                            imwrite(I,newPath);
                            
                        catch err
                            pass = 0;
                            err
                        end
                        if pass
                            break
                        end
                    end % try writing
                    
                end %if there is no image path
                
                stopTime = datenum(clock);
                dur = (stopTime-startTime)*24*60*60
                durations(s) = dur;
            end %if target file already exists
            
        end %columns
        
    end % rows
end % sections
        