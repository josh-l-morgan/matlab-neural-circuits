clear all

%% map directories
rotDir = 'D:\LGNs1\Analysis\movies\morphBout\rot180c\'
mixDir = [rotDir 'mix30_0b\'];


targDir = [rotDir 'targ\'];
boutDir = [rotDir 'bout\'];
axDir = [rotDir 'ax\'];
tcDir = [rotDir 'tcr\'];
obDir = [rotDir 'otherBout\'];
if ~exist(mixDir,'dir'), mkdir(mixDir),end;

dirOrd = {targDir boutDir axDir obDir tcDir};


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
mixDeg = 0:359;
blank = (0:1000) * 0;
fade = 20;
change = [1 72 144 216 288 360]
change = round([0:5] .* 360/5 - fade);
change(change<1) = 1;

degID = 0:1000;
mix{1} = blank;
mix{1}(change(1):change(2)) = 1;
mix{1}(change(2)+1:change(2) + fade) = (fade:-1:1)/fade;
%mix{1}(change(6)-fade+1:change(6)) = (1:fade)/fade;
mix{1}(change(6):change(6)+fade-1) = (1:fade)/fade;


mix{2} = blank;
mix{2}(change(2):change(3)) = 1;
mix{2}(change(3)+1:change(3) + fade) = (fade:-1:1)/fade;
mix{2}(change(2):change(2)+ fade-1) = (1:fade)/fade;



mix{3} = blank;
mix{3}(change(3):change(4)) = 1;
mix{3}(change(4)+1:change(4) + fade) = (fade:-1:1)/fade;
mix{3}(change(3):change(3)+ fade-1) = (1:fade)/fade;



mix{4} = blank;
mix{4}(change(4):change(5)) = 1;
mix{4}(change(5)+1:change(5) + fade) = (fade:-1:1)/fade;
mix{4}(change(4):change(4)+ fade-1) = (1:fade)/fade;



mix{5} = blank;
 mix{5}(change(5):change(6)) = 1;
mix{5}(change(6)+1:change(6)+fade) = (fade:-1:1)/fade;
mix{5}(change(5):change(5)+fade-1) = (1:fade)/fade;

%% mix iimages

for i = 1:length(mixDeg);
    d = mixDeg(i);
    disp(sprintf('mixing %d of %d',i,length(mixDeg)))
    clear Imix
    mixTarg = (find(degID==d,1));
    for m = 1:length(mix)
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


  
  
  
  
  
  
  
  
  
  