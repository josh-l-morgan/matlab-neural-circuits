TPN = GetMyDir;

%% parse names
nams = GetPics(TPN);

for i = 1:length(nams)
    nam = nams{i};
    zs = regexp(nam,'_z');
    ms = regexp(nam,'_m');
    z(i) = str2num(nam(zs(1)+2:zs(1)+3))+1;
    m(i) = str2num(nam(ms(1)+2:ms(1)+3))+1;
end

%% Read tiles
tileIds = sort(unique(m),'ascend');

for t = 1:length(tileIds)
    planes = find(m==tileIds(t));
    Idat = imfinfo([TPN nams{planes(1)}]);
    xs = Idat.Width;
    ys = Idat.Height;
    tileSum = zeros(ys,xs,Idat.SamplesPerPixel);
    tileMax = tileSum;
    numPlanes = length(planes);
    for p = 1:numPlanes
        sprintf('reading plane %d of %d from tile %d of %d.',...
            p,numPlanes,t,length(tileIds))
        nam = nams{planes(p)};
        I = double(imread([TPN nam]));
        %image(uint8(I*3) ),pause(.01)
        tileSum = tileSum + I;
        tileMax = max(tileMax,I);
        tileStd(p) = std(I(:));
        tileBrightness(p) = sum(I(:));
    end
    tileMean = uint8(tileSum/numPlanes);
    maxStd = find(tileStd == max(tileStd),1);
    tileMaxDif = imread([TPN nams{planes(maxStd)}]);
    brightest = find(tileBrightness == max(tileBrightness),1);
    tileBright = imread([TPN nams{planes(brightest)}]);
    Tiles(t).tileMean = uint8(tileMean);
    Tiles(t).tileMax = uint8(tileMax);
    Tiles(t).tileMaxDif = uint8(tileMaxDif);
    Tiles(tn).tileBright = uint8(tileBright);
    image(tileMean),pause(.01)
end

save([TPN 'Tiles.mat'],'Tiles')
%% Display tiles

for i = 1:length(Tiles)
    i
   I = Tiles(i).tileMean;
   image(I),pause
end


%% Arrange tiles
rows = 4;
cols = 15;
overlap = .05;
xLeft = fix(xs/2);
yLeft = fix(ys/2);
xOver = xs * overlap;
yOver = ys * overlap;

xpos = round((1:cols) * (xs - xOver)-xs/2+xOver);
ypos = round((1:rows) * (ys - yOver)-ys/2+yOver);
xstart = xpos-xLeft+1; xstop = xpos+xs - xLeft;
ystart = ypos-yLeft+1; ystop = ypos+ys - yLeft;

ysMos = max(ystop);
xsMos = max(xstop);
%%
Mos = zeros(ysMos,xsMos,3,'uint8');
t = 0;
for y = 1: rows;
    for x = 1: cols;
        t = t + 1;
        I = Tiles(t).tileMaxDif; 
        %I1 = I(:,:,1);
        Mos(ystart(y):ystop(y),xstart(x):xstop(x),:) = I(yLeft : end - yLeft,;
        image(Mos),pause(.01)
    end
end

%% 
 
mkdir([TPN 'Mos'])
imwrite([TPN 'Mos\tileMaxDif.tif'],Mos,'Compression','none')

    
    
    

