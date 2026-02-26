



SPN = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2\'
TPN = 'E:\IxQ_KarlsRetinaVG3_2019\VAST\scm_ixQ_hires_pyr2_reDice3\'

oldRes = 4096;
newRes = 512;
fac = oldRes/newRes;


dSPN = dir(SPN)
secDir = {dSPN(3:end).name};

clear secNum
for i = 1:length(secDir)
 secNum(i) = str2num(secDir{i});
end
[sortNum ia] = sort(secNum,'ascend');
secDir = secDir(ia);


mipDir = {'0'};%, '1', '2', '3', '4', '5', '6', '7', '8'};


for s =760:length(secDir)
   
    disp(sprintf('redicing %d of %d',s,length(secDir)))
   for m = 1:length(mipDir)
       
        oldDir = [SPN secDir{s} '\' mipDir{m} '\'];
        newDir = [TPN secDir{s} '\' mipDir{m} '\'];
        if ~exist(newDir,'dir')
            mkdir(newDir)
        end
       if exist(oldDir,'dir')
           
           dI = dir([oldDir '*.png']);
           inams = {dI.name};
           for i = 1:length(inams)
               
               nam = inams{i};
               u = regexp(nam,'_');
               d = regexp(nam,'.png');
               r = str2num(nam(1:u-1));
               c = str2num(nam(u+1:d-1));
               
               %info = imfinfo([oldDir nam]);
%                ys = info.Height;
%                xs = info.Width;     

               I = imread([oldDir nam]);
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
                           
                                                                                
                           nr = (r)*fac + y-1;
                           nc = (c)*fac + x-1;
                           
                           newFileNames{cnt} = sprintf('%d_%d.png',nr,nc);
                           
                       end
                       
                   end
               end
               
               if cnt>0
                   if ~exist([newDir newFileNames{end}],'file')
                       %I = imread([oldDir nam]);
                       for n = 1:cnt
                           sI = I(yStarts(n):yStops(n),xStarts(n):xStops(n));
                           imwrite(uint8(sI),[newDir newFilenames{n}]);
                       end
                   end
               end % if new files are possible
               
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





