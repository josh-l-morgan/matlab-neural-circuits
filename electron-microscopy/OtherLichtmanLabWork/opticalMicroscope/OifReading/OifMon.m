%%Read OIFs and assemble into single matrix that can be written as 2D
%%color tifs
clear all

MPN=GetMyDir


%% Get Log file and read
dMPN = dir(MPN);
for i = 1:length(dMPN)
    nam = dMPN(i).name;
    if regexp(nam,'.log')
        logFile = nam;
        break
    end
end

logX = xmlread([MPN logFile]);
[xtree RootName DOMnode] = xml_read(logX);
Idat = xtree.Mosaic.ImageInfo;
tileNum = length(Idat);

%% Match files

for i = 1:tileNum
   fileName = Idat(i).Filename;
   preNam = fileName(1:end-4);
   for n = 1:length(dMPN)
      if dMPN(n).isdir
         matched = regexp(dMPN(n).name,preNam); 
         if matched == 1
            foldName{i} = dMPN(n).name; 
         end
      end
   end
end

%% Run files
myChan = [2 1 3 4 5 6 7 8 9 10];
for i = 1 : tileNum
   foldNam = foldName{i};
   inam = getTifs([MPN foldNam]);
   mif.t(i).foldNam = foldNam;
   maxChan = 0;
maxZ = 0;
   for p = 1:length(inam)
       nam = inam{p};
       mif.t(i).p(p).fileNam = nam;
       col = regexp(nam,'_C');
       if ~isempty(col)
           channel = str2num(nam(col+2:col+4));
           z = str2num(nam(col+6:col+8));
           
           mif.t(i).p(p).sub = [myChan(channel) z]; 
           maxChan = max(channel,maxChan);
           maxZ = max(z,maxZ);
       end
   end
   mif.t(i).planes = maxZ;
   mif.t(i).channels = maxChan;
end
%%
Info = imfinfo([MPN mif.t(1).foldNam '\' mif.t(1).p(1).fileNam]);

for i = 1 : tileNum
   Istack = zeros(Info.Height,Info.Width,3,mif.t(i).planes);
   sprintf('running tile %d of %d.',i,tileNum)
   for p = 1:length(mif.t(i).p) 
       if ~isempty(col)
           I = imread([MPN mif.t(i).foldNam '\' mif.t(i).p(p).fileNam]);           
           Istack(:,:,mif.t(i).p(p).sub(1), mif.t(i).p(p).sub(2)) = I;
       end
   end
   Imax = max(Istack,[],4);
   Imax = uint8(Imax/16);
   imagesc(Imax),pause(.01) 
   maxStack{i} = Imax; 
end

mif.Idat = Idat;
mif.xtree = xtree;
save([MPN 'maxStack.mat'],'maxStack');


%% Build montage
res = Info.XResolution;
res = .1027;
pix = Info.Height;
fov = pix*res;
for i = 1:length(Idat);
    X(i) = Idat(i).XPos;
    Y(i) = Idat(i).YPos;
end

X = X - min(X) + fov;
Y = Y - min(Y) + fov;
X = round(X /res);
Y = round(Y /res);


mon = zeros(fix(max(Y))+10,fix(max(X))+10,3,'uint8');

for i =length(Idat):-1:1
    mon(Y(i)-pix+1:Y(i),X(i)-pix+1:X(i),:) = flipdim(maxStack{i},2);
end
   image(mon)
   
imwrite(mon,[MPN 'mon.tif'],'Compression','none')



MPN




