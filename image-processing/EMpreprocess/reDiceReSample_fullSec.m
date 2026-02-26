



SPN = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2\'
TPN = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2_reDice4\'

oldRes = 4096;
newRes = 512;
fac = oldRes/newRes;
tile = zeros(newRes);


dSPN = dir(SPN)
secDir = {dSPN(3:end).name};

clear secNum
for i = 1:length(secDir)
    secNum(i) = str2num(secDir{i});
end
[sortNum ia] = sort(secNum,'ascend');
secDir = secDir(ia);


mipDir = {'0'};%, '1', '2', '3', '4', '5', '6', '7', '8'};

for s =1:length(secDir)
    
    disp(sprintf('redicing %d of %d',s,length(secDir)))
    
    oldDir = [SPN secDir{s} '\' '0' '\'];
    newDir = [TPN secDir{s} '\' '0' '\'];
    if ~exist(newDir,'dir')
        mkdir(newDir)
    end
    if exist(oldDir,'dir')
        
        dI = dir([oldDir '*.png']);
        inams = {dI.name};
        rs = zeros(length(inams),1);
        cs = rs;
        I = {};
        tic
        parfor i = 1:length(inams)
            
            nam = inams{i};
            u = regexp(nam,'_');
            d = regexp(nam,'.png');
            rs(i) = str2num(nam(1:u-1));
            cs(i) = str2num(nam(u+1:d-1));
            Ir{i} = imread([oldDir nam]);
        end
        toc
        
     
        
        countI0 = 0;
        I0 = {};
        Irow = [];
        Icol = [];
        for i = 1:length(Ir)
            I = Ir{i};
            [ys xs] = size(I);
            
            cnt = 0;
            yStarts  = zeros(fac^2,1);
            yStops = yStarts;
            xStarts =  yStarts;
            xStops = yStarts;
            newFileNames = {};
            
            for y = 1:fac
                for x = 1:fac
                    yStart = (y-1)* newRes + 1;
                    xStart = (x-1) * newRes + 1;
                    yStop = min(y * newRes,ys);
                    xStop = min(x * newRes, xs);
                    
                    if (yStart<ys) & ( xStart<xs)
                        cnt = cnt+1;
                        yStop = min(yStop,ys);
                        
                        yStarts(cnt) = yStart;
                        yStops(cnt) = yStop;
                        xStarts(cnt) = xStart;
                        xStops(cnt) = xStop;
                        
                        
                        nrs(cnt) = (r)*fac + y-1;
                        ncs(cnt) = (c)*fac + x-1;
                        
                        newFileNames{cnt} = sprintf('%d_%d.png',nrs(cnt),ncs(cnt));
                        
                    end
                    
                end
            end
            
            
            
            if cnt>0
                if ~exist([newDir newFileNames{end}],'file')
                    %I = imread([oldDir nam]);
                    
                    I0t = {};
                    Irowt = zeros(cnt,1);
                    Icolt = Irowt;
                    parfor n = 1:cnt
                        countI0 = countI0 + 1;
                        %                             tile(1:yStops(n)-yStarts(n)+1,1:xStops(n)-xStarts(n)+1)...
                        %                                 = I(yStarts(n):yStops(n),xStarts(n):xStops(n));
                        tile   = I(yStarts(n):yStops(n),xStarts(n):xStops(n));
                        if (size(tile,1) ~= newRes ) & (size(tile,2)~=newRes)
                            tile(newRes,newRes) = 0;
                        end
                        I0t{n,1} = tile;
                        Irowt(n,1) = nrs(n);
                        Icolt(n,1) = ncs(n);
                        imwrite(uint8(tile),[newDir newFilenames{n}]);
                    end
                    I0 = cat(1,I0,I0t);
                    Irow = cat(1,Irow,Irowt);
                    Icol = cat(1,Icol,Icolt);
                    
                end
                
                
            end % if new files are possible
            
        
            
            
            
        end % run all images in section
        toc
        
        
        
        
        
        mipDir = {'1', '2', '3', '4', '5', '6', '7', '8'};

        
        for m = 1:length(mipDir)
            
                newDir = [TPN secDir{s} '\' mipDir{m} '\'];
                    
                newRow = floor(Irow/2);
                newCol = floor(Icol/2);
                newRow2 = mod(Irow,2);
                newCol2 = mod(Icol,2);
                
                for i = 1 :length(I0)
                    
                    
                    
                    
                end
               
                            
            
            
            
        end
        
        
        
        
        
        
        
    end
end



%{
               
               for y = 1:fac
                   for x = 1:fac
                       yStart = (y-1)* newRes + 1;
                       xStart = (x-1) * newRes + 1;
                       yStop = min(y * newRes,ys);
                       xStop = min(x * newRes, xs);
                       
                       if (yStart<ys) & ( xStart<xs)
                           yStop = min(yStop,ys);
                                                      
                           sI = I(yStart:yStop,xStart:xStop);
                                                                                
                           nr = (r)*fac + y-1;
                           nc = (c)*fac + x-1;
                           
                           newFilename = sprintf('%d_%d.png',nr,nc)
                           imwrite(uint8(sI),[newDir newFilename]);
                           
                       end
                       
                   end
               end
               
               
               
               
           end
           
       end
           
   end
    
end


%}





