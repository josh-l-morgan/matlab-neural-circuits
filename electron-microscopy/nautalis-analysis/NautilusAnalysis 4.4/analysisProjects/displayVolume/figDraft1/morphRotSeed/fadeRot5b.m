clear all

%% map directories
rotDir = 'D:\LGNs1\Analysis\movies\morphBout\rot180c\'
mixDir = [rotDir 'mix360_x5k\'];


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
       if ~isempty(deg)
       dotpng = regexp(nam,'.png');
       d = str2num(nam(deg(1)+3:dotpng-1));
       ds(f) = d;
       else
           ds(f) = nan;
       end
        
    end
    iDegs{i} = ds;
end

%% chose turn ons
mixDeg = repmat(0:1:180,[1,5]);
blank = mixDeg * 0;
fade = 60;
change = round([0:5] .* 280 - fade);
change(change<1) = 1;

degID = 0:10000;
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

Itemp =   double(imread([dirOrd{1} iNams{1}(1).name]));
Itemp = Itemp * 0;

for i = 1:1:length(mixDeg);
    disp(sprintf('mixing %d of %d',i,length(mixDeg))),pause(.01)
    
    crop = [100 1369;600 2000];
    [Imix] = funWriteMix(i,mix,iDegs,mixDeg,dirOrd,iNams,Itemp,crop);
    mixName = sprintf('mix_%03.0f.png',i);
    imwrite(uint8(Imix),[mixDir mixName]);
end





  
  
  
  
  
  
  