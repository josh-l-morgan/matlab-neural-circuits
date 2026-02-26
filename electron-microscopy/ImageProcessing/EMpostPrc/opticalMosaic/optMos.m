clear all
TPN = GetMyDir;
tifNames = getPics(TPN)
saveTPN = [TPN(1:end-1) 'maxesMed'];
if ~exist(saveTPN,'dir'), mkdir(saveTPN); end

for i = 1: length(tifNames)
   nam = tifNames{i};
   zs = regexp(nam,'_z');
   ms = regexp(nam,'_m');
   dots = regexp(nam,'.tif');
   

   z(i) = str2num(nam(zs(end)+2:ms(end)-1))+1;
   if isempty(ms),
       m(i) = 1;
   else
       m(i) = str2num(nam(ms(end)+2:dots(end)-1))+1;
   end
end
mosNum = sort(unique(m),'ascend');

for i = 1:length(mosNum)
    mos = find(m == mosNum(i));
    sprintf('running stack %d of %d',i,length(mosNum))

    info = imfinfo([TPN tifNames{mos(1)}]); 
    I = zeros(info.Height,info.Width,3,length(mos),'uint8');
    if info.BitDepth== 48
        scale = (2^8)/(2^12);
    else 
        scale = 1;
    end
    for p = 1:length(mos)   
        pic = imread([TPN tifNames{mos(p)}]);
        I(:,:,:,p) = double(pic) * scale;
    end  
    [newI shifted] = findSlide(I);
    clear I
    shifted
    for i = 1:size(newI,4)
        for c = 1:size(newI,3)
            pic = newI(:,:,c,i);
            pic = medfilt2(pic,[3 3]);
            newI(:,:,c,i) = pic;
        end
    end
    maxI = max(newI,[],4);
    imwrite(maxI,[saveTPN '\' num2str(mosNum(i)) '.tif'],'Compression','none')
    image(maxI),pause(.1)
    
end

clear newI
