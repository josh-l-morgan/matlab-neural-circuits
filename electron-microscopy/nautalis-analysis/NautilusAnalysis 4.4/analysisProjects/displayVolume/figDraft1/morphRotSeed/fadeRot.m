clear all

%% map directories
rotDir = 'D:\LGNs1\Analysis\movies\morphBout\rot180c\'
mixDir = [rotDir 'mix11\'];


targDir = [rotDir 'targ\'];
boutDir = [rotDir 'bout\'];
axDir = [rotDir 'ax\'];
tcDir = [rotDir 'tcr\'];
if ~exist(mixDir,'dir'), mkdir(mixDir),end;

dirOrd = {targDir boutDir axDir tcDir};


for i = 1:length(dirOrd);
    
    iNams{i} = dir([dirOrd{i} '*.png']);
    ds = zeros(length(iNams{i}),1);
    for f = 1:length(iNams{i})
       nam = iNams{i}(f).name;
       
       deg = regexp(nam,'deg');
       dotpng = regexp(nam,'.png');
       d = str2num(nam(deg(1)+3:dotpng-1));
       ds(f) = d;
        
    end
    iDegs{i} = ds;
end

%% chose turn ons
mixDeg = 0:1:359;
blank = (0:1000) * 0;


degID = 0:1000;
mix{1} = blank;
mix{1}(1:90) = 1;
mix{1}(91:110) = (20:-1:1)/20;
mix{1}(341:360) = (1:20)/20;


mix{2} = blank;
mix{2}(91:180) = 1;
mix{2}(91:110) = (1:20)/20;
mix{2}(181:200) = (20:-1:1)/20;


mix{3} = blank;
mix{3}(181:270) = 1;
mix{3}(181:200) = (1:20)/20;
mix{3}(271:290) = (20:-1:1)/20;



mix{4} = blank;
mix{4}(271:360) = 1;
mix{4}(271:290) = (1:20)/20;
mix{4}(341:360) = (20:-1:1)/20;



%% mix iimages

for i = 1:length(mixDeg);
    d = mixDeg(i);
    disp(sprintf('mixing %d of %d',i,length(mixDeg)))
    clear Imix
    mixTarg = (find(degID==d,1));
    for m = 1:4
        if mix{m}(mixTarg)>0
           targ = find(iDegs{m}==d,1);
           if ~isempty(targ)
              
               I = double(imread([dirOrd{m} iNams{m}(targ).name]));
               if exist('Imix','var')
                   Imix = Imix + I * mix{m}(mixTarg);
               else
                   Imix = I * mix{m}(mixTarg);
               end
               
           end
            
        end
                  
    end
    
    if exist('Imix','var')
        mixName = sprintf('mix_%03.0f.png',d);
        imwrite(uint8(Imix),[mixDir mixName]); 
    end
    
end


  
  
  
  
  
  
  
  
  
  