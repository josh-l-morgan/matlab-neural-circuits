
clear all



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

for s = 302 %1:length(secDir)
    
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
        Ir = {};
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
                        
                        
                        nrs(cnt) = (rs(i))*fac + y-1;
                        ncs(cnt) = (cs(i))*fac + x-1;
                        
                        newFileNames{cnt} = sprintf('%d_%d.png',nrs(cnt),ncs(cnt));
                        
                    end
                    
                end
            end
            
            
            
            if cnt>0
                if 1%~exist([newDir newFileNames{end}],'file')
                    %I = imread([oldDir nam]);
                    
                    I0t = {};
                    Irowt = zeros(cnt,1);
                    Icolt = Irowt;
                    parfor n = 1:cnt
                        countI0 = countI0 + 1;
                        %                             tile(1:yStops(n)-yStarts(n)+1,1:xStops(n)-xStarts(n)+1)...
                        %                                 = I(yStarts(n):yStops(n),xStarts(n):xStops(n));
                        tile = I(yStarts(n):yStops(n),xStarts(n):xStops(n));
                        if (size(tile,1) ~= newRes ) & (size(tile,2)~=newRes)
                            tile(newRes,newRes) = 0;
                        end
                        I0t{n,1} = tile;
                        Irowt(n,1) = nrs(n);
                        Icolt(n,1) = ncs(n);
                        imwrite(uint8(tile),[newDir newFileNames{n}]);
                    end
                    I0 = cat(1,I0,I0t);
                    Irow = cat(1,Irow,Irowt);
                    Icol = cat(1,Icol,Icolt);
                    %                     scatter(Irowt,Icolt,'o','markerfacecolor',[rand rand rand])
                    %                     hold on
                    pause(.1)
                end
                
                
            end % if new files are possible
            
            
            
            
            
        end % run all images in section
        toc
        
        
        
        
        hR = newRes/2;
        mipDir = {'1', '2', '3', '4', '5', '6', '7', '8'};
        
        
        Iold = I0;
        oldRow = Irow;
        oldCol = Icol;
        clear Inew;
        
        for m = 1:length(mipDir)
            clear Inew
            newDir = [TPN secDir{s} '\' mipDir{m} '\']
            if ~exist(newDir,'dir'), mkdir(newDir); end
            
            newRow = ceil(oldRow/2);
            newCol = ceil(oldCol/2);
            newRow2 = mod((oldRow-1),2);
            newCol2 = mod((oldCol-1),2);
            newInd = sub2ind([max(newRow) max(newCol)],newRow,newCol);
            
            
            %transform tiles to downsampled supertiels
            
            Inew{max(newInd)} = [];
            for i = 1 :length(Iold)
                if ~isempty(Iold{i})
                    Inew{newInd(i)}(newRow2(i)*hR+1:newRow2(i)*hR+hR,...
                        newCol2(i)*hR+1:newCol2(i)*hR+hR) = imresize(Iold{i},.5);
                end
            end
            for i = 1:length(Inew)  % Fix sizes to new res.
                if ~isempty(Inew{i})
                    
                    [tys txs] = size(Inew{i});
                    if ~((tys==newRes)&(txs==newRes))
                        Inew{i}(newRes,newRes) = 0;
                    end
                end
            end
            
            
            
            %
            if 1
                clear Itest
                
                
                for i = 1 :length(Iold)
                    if ~isempty(Iold{i})
                    disp(sprintf('show %d of %d',i,length(Iold)))
                    startY = (newRow(i)-1) * newRes + newRow2(i) *hR + 1;
                    startX = (newCol(i)-1) * newRes + newCol2(i)*hR + 1;
                    stopY = (newRow(i)-1) * newRes + newRow2(i) *hR + hR;
                    stopX = (newCol(i)-1) * newRes + newCol2(i)*hR + hR;
                    Itest(startY:stopY, startX:stopX) = imresize(Iold{i},.5);
                    else
                    i
                    end
                    
                end
                image(Itest)
                pause(.01)
            end
           pause
            
            if 1
                clear Itest
                for i = 1 :length(newInd)
                    
                    if ~isempty(Inew{newInd(i)})
                    disp(sprintf('show %d of %d',i,length(newInd)))
                    startY = (newRow(i)-1) * newRes +  1;
                    startX = (newCol(i)-1) * newRes  + 1;
                    stopY = (newRow(i)) * newRes ;
                    stopX = (newCol(i)) * newRes ;
                    Itest(startY:stopY, startX:stopX) = Inew{newInd(i)};
                    end
                end
                image(Itest)
                pause(.01)
            end
            m
            pause
            
            % make file names
            clear newFileNames
            for i = 1:length(newInd)
                newFileNames{i} = sprintf('%d_%d.png',newRow(i),newCol(i));
            end
            
            % Write files
            for i = 1:length(newInd)
               
                    if ~isempty(Inew{newInd(i)})
                        imwrite(uint8(Inew{newInd(i)}),[newDir newFileNames{i}]);
                    end
                
            end
            
            Iold = Inew;
            oldRow = newRow;
            oldCol = newCol;
            
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





