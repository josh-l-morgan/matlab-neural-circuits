
clear all


SPN = 'Z:\Active\morganLab\DATA\LGN\LGNs1\2023transfer\mip0_batch2\'
TPN = 'Z:\Active\morganLab\DATA\LGN\LGNs1\2023transfer\LGNs1_vsvi\';

oldRes = 2048;
newRes = 512;
fac = oldRes/newRes;
tile = zeros(newRes);

dSPN = dir(SPN)
secDir = {dSPN(3:end).name};

clear secNum
for i = 1:length(secDir)
    nam = secDir{i};
    und = regexp(nam,'_w');
    secNum(i) = str2num(nam(1:und(1)-1));
end
[sortNum ia] = sort(secNum,'ascend');
secDir = secDir(ia);


mipDir = {'0', '1', '2', '3', '4', '5', '6', '7', '8'};

for s =1:length(secDir)
    
    disp(sprintf('redicing %d of %d',s,length(secDir)))
    
    oldDir = [SPN secDir{s} '\'];
    newDir = sprintf('%s%4.0f\\0\\',TPN,secNum(s));
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
            u1 = regexp(nam,'_tr');
            u2 = regexp(nam,'-tc');
            d = regexp(nam,'.png');
            rs(i) = str2num(nam(u1+3:u2-1));
            cs(i) = str2num(nam(u2+3:d-1));
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
                        
                        nrs(cnt) = (rs(i))*fac + y;
                        ncs(cnt) = (cs(i))*fac + x;
                        
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
        oldInd = 1:length(Icol);

        clear Inew;
        
        for m = 1:length(mipDir)
            clear Inew
            newDir = sprintf('%s%4.0f\\%s\\',TPN,secNum(s),mipDir{m});
            if ~exist(newDir,'dir'), mkdir(newDir); end
            
            newRow = ceil(oldRow/2);
            newCol = ceil(oldCol/2);
            newRow2 = mod((oldRow-1),2);
            newCol2 = mod((oldCol-1),2);
            newInd = sub2ind([max(newRow) max(newCol)],newRow,newCol);
            uNewInd = unique(newInd);
            [uNewRow uNewCol] = ind2sub([max(newRow) max(newCol)],uNewInd);
            
            %transform tiles to downsampled supertiels
            Inew{max(newInd)} = [];
            for i = 1 :length(oldInd)
                if ~isempty(Iold{oldInd(i)})
                    Inew{newInd(i)}(newRow2(i)*hR+1:newRow2(i)*hR+hR,...
                        newCol2(i)*hR+1:newCol2(i)*hR+hR) = imresize(Iold{oldInd(i)},.5);
                end
            end
            
            % Fix sizes to new res.
            for i = 1:length(Inew)  
                if ~isempty(Inew{i})
                    [tys txs] = size(Inew{i});
                    if ~((tys==newRes)&(txs==newRes))
                        Inew{i}(newRes,newRes) = 0;
                    end
                end
            end
                                    
           % if 0
           %      clear Itest
           % 
           %      for i = 1 :length(oldInd)
           %          if ~isempty(Iold{oldInd(i)})
           %              disp(sprintf('show %d of %d',i,length(Iold)))
           %              startY = (newRow(i)-1) * newRes + newRow2(i) *hR + 1;
           %              startX = (newCol(i)-1) * newRes + newCol2(i)*hR + 1;
           %              stopY = (newRow(i)-1) * newRes + newRow2(i) *hR + hR;
           %              stopX = (newCol(i)-1) * newRes + newCol2(i)*hR + hR;
           % 
           % 
           %              startY = (newRow(i)-1) * newRes +  1;
           %              startX = (newCol(i)-1) * newRes +  1;
           %              stopY = (newRow(i)-1) * newRes +  newRes;
           %              stopX = (newCol(i)-1) * newRes + newRes;
           % 
           % 
           %              %Itest(startY:stopY, startX:stopX) = imresize(Iold{oldInd(i)},1);
           %              Itest(startY:stopY, startX:stopX) = Iold{oldInd(i)};
           % 
           %          else
           %              i
           %          end
           % 
           %      end
           %      image(Itest)
           %              pause(.01)
           %  end
                       
            % 
            % if 0
            %     clear Itest
            %     for i = 1 :length(oldInd)
            % 
            %         if ~isempty(Iold{oldInd(i)})
            %             disp(sprintf('show %d of %d',i,length(oldInd)))
            %             startY = (oldRow(i)-1) * newRes +  1;
            %             startX = (oldCol(i)-1) * newRes  + 1;
            %             stopY = (oldRow(i)) * newRes ;
            %             stopX = (oldCol(i)) * newRes ;
            %             Itest(startY:stopY, startX:stopX) = Iold{oldInd(i)};
            %         end
            %     end
            %     image(Itest)
            %     pause(.01)
            % end
            
           
           % if 0
           %      clear Itest
           %      for i = 1 :length(uNewInd)
           % 
           %          if ~isempty(Inew{uNewInd(i)})
           %          disp(sprintf('show %d of %d',i,length(uNewInd)))
           %          startY = (uNewRow(i)-1) * newRes +  1;
           %          startX = (uNewCol(i)-1) * newRes  + 1;
           %          stopY = (uNewRow(i)) * newRes ;
           %          stopX = (uNewCol(i)) * newRes ;
           %          Itest(startY:stopY, startX:stopX) = Inew{uNewInd(i)};
           %          end
           %      end
           %      image(Itest)
           %      pause(.01)
           %  end
          
            
            % make file names
            clear newFileNames
            for i = 1:length(uNewInd)
                newFileNames{i} = sprintf('%d_%d.png',uNewRow(i),uNewCol(i));
            end
            
            % Write files
            parfor i = 1:length(uNewInd)
               
                    if ~isempty(Inew{uNewInd(i)})
                        try
                        imwrite(uint8(Inew{uNewInd(i)}),[newDir newFileNames{i}]);
                        catch err
                            pause(1)
                            try 
                                                        imwrite(uint8(Inew{uNewInd(i)}),[newDir newFileNames{i}]);
                            catch err
                                disp(sprintf('could not write %s',[newDir newFileNames{i}]))
                            end
                        end
                    end
                
            end
            
            Iold = Inew;
            oldRow = uNewRow;
            oldCol = uNewCol;
            oldInd = uNewInd;
            
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





