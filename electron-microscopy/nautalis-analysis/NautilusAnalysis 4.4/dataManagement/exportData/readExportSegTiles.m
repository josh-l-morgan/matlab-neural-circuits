

SPN = 'E:\IxQ_KarlsRetinaVG3_2019\seg2png\joshm_aw2d\';

dSPN = dir([SPN '*.png']);
inams = {dSPN.name};

info = imfinfo([SPN inams{1}]);
width = info.Width;
height = info.Height;

clear r c s
for i = 1:length(inams)
    nam = inams{i};
    ur = regexp(nam,'_r');
    uc = regexp(nam,'_c');
    us = regexp(nam,'_s');
    dot = regexp(nam,'.png');
    r(i) = str2num(nam(ur(1)+2:uc(1)-1));
    c(i) = str2num(nam(uc(1)+2:us(1)-1));
    s(i) = str2num(nam(us(1)+2:dot(1)-1));
    
end

secs = unique(s);

ds = 10;

widthDS = ceil(width/ds);
heightDS = ceil(height/ds);
se = strel('disk',3);
fullWidth = max(c) * widthDS;
fullHeight = max(r) * heightDS;

volDS = zeros(fullHeight,fullWidth,length(secs),...
    'uint16');

for sn = 1:length(secs)
    sprintf('building section %d of %d',secs(sn),length(secs))
    issec = find(s==secs(sn));
    secR = r(issec);
    secC = c(issec);
    %scatter(secC,secR,'.')
    secN = inams(issec);
    Isec = zeros(fullHeight,fullWidth,'uint16');
    w = 10;
    for t = 1:length(issec)
       I = imread([SPN secN{t}]);
       Ids = imresize(I, 1/ds,'nearest');
        [Iy, Ix] = size(Ids);
      
       
       wX = min(w,Ix-1);
       wY = min(w,Iy-1);
       
       Ids(:,1:wX) = 100000;
       Ids(:,end-wX:end) = 100000;
       Ids(1:wY,:) = 1000000;
       Ids(end-wY:end,:) = 100000;
       
       
        Isec((secR(t)-1)*heightDS+1:(secR(t)-1)*heightDS + Iy,...
            (secC(t)-1)*widthDS+1:(secC(t)-1)*widthDS + Ix) = Ids;
       
%         Isec((secR(t)-1)*width+1:secR(t)*width,...
%         (secC(t)-1)*height+1:secC(t)*height) = I;
    end
   
    image(Isec*1000),pause(.01)
    pause(.01)
    volDS(:,:,sn) = Isec;
    
end


for i = 1:size(volDS,3)
   
    image(volDS(:,:,i)*1000)
    pause
    
end




























