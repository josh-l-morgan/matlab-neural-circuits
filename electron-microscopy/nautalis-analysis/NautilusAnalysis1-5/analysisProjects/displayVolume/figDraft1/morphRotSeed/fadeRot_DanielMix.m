clear all

%% map directories
rotDir = 'D:\LGNs1\Analysis\movies\Daniel\mix1\'
mixDir = [rotDir 'mix25\'];

dir1 = [rotDir 'render\'];
dir2 = [rotDir 'render2\'];
dir3 = [rotDir 'render3\'];

if ~exist(mixDir,'dir'), mkdir(mixDir),end;

dirOrd = {dir1  dir3};


for i = 1:length(dirOrd);
    
    iNams{i} = dir([dirOrd{i} '*.png']);
    ds = zeros(length(iNams{i}),1);
    for f = 1:length(iNams{i})
       nam = iNams{i}(f).name;
       dotpng = regexp(nam,'.png');
       d = str2num(nam(dotpng(1)-4:dotpng-1));
       ds(f) = d;
        
    end
    iDegs{i} = ds;
end

%% chose turn ons
mixDeg = 0:1:359;
blank = (0:1000) * 0;
fade = 20;
change = [1 60 180 360]
%change = round([0:3] .* 360/5 - fade);
change(change<1) = 1;

degID = 0:1000;
mix{1} = blank + 1;
mix{1}(change(1):change(2)) = 1;
mix{1}(change(2)+1:change(3)) = 0;
mix{1}(change(2)+1:change(2) + fade) = (fade:-1:1)/fade;
%mix{1}(change(6)-fade+1:change(6)) = (1:fade)/fade;
mix{1}(change(3):change(3)+fade-1) = (1:fade)/fade;


mix{2} = blank;
mix{2}(change(2):change(3)) = 1;
mix{2}(change(3)+1:change(3) + fade) = (fade:-1:1)/fade;
mix{2}(change(2):change(2)+ fade-1) = (1:fade)/fade;
% 
% 
% mix{3} = blank;
% mix{3}(change(2):change(3)) = 1;
% mix{3}(change(3)+1:change(3) + fade) = (fade:-1:1)/fade;
% mix{3}(change(2):change(2)+ fade-1) = (1:fade)/fade;

% 
% 
% mix{3} = blank;
% mix{3}(change(3):change(4)) = 1;
% mix{3}(change(4)+1:change(4) + fade) = (fade:-1:1)/fade;
% mix{3}(change(3):change(3)+ fade-1) = (1:fade)/fade;
% 
% 
% 
% mix{4} = blank;
% mix{4}(change(4):change(5)) = 1;
% mix{4}(change(5)+1:change(5) + fade) = (fade:-1:1)/fade;
% mix{4}(change(4):change(4)+ fade-1) = (1:fade)/fade;
% 
% 
% 
% mix{5} = blank;
%  mix{5}(change(5):change(6)) = 1;
% mix{5}(change(6)+1:change(6)+fade) = (fade:-1:1)/fade;
% mix{5}(change(5):change(5)+fade-1) = (1:fade)/fade;

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
        mixName = sprintf('mix_%03.0f.png',d+360);
        imwrite(uint8(Imix),[mixDir mixName]); 
        image(uint8(Imix)),pause(.01)
    end
    
end


  
  
  
  
  
  
  
  
  
  