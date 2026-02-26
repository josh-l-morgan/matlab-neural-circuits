function[qual] = checkQuality(TPN)


if ~exist('TPN','var')
    TPN = GetMyDir;
end

dsamp = 100;  %number of samples in each dimension

%% Get Tiles
picNams = GetPics(TPN);

ct = 0;
for i = 1:length(picNams)
    nam = picNams{i}
    rs = regexp(nam,'_r');
    us= regexp(nam,'_');
    cs = regexp(nam,'-c');
    if ~isempty(rs) & ~isempty(cs)
        r = str2num(nam(rs(1)+2:cs(1)-1))
        c = str2num(nam(cs(1)+2:us(2)-1))
        ct = ct+1;
        row(ct) = r;
        col(ct) = c;
        nams{ct} = nam;
    end
end

if isempty(nams)
    nams = picNams;
end

%% Define variables

yshift = [ 0 0 0 1 1 1 2 2 2]; %shifts for extracting 3X3 samples
xshift = [ 0 1 2 0 1 2 0 1 2];

%1 %2 %3
%4 %5 %6
%7 %8 %9


sur = {[2 4 6 8], [1 2 3 6], [2 3 6 9]};
cent = {[1 3 7 9], [9 8 7 4], [8 7 4 1]};
%%Other useful patterns
%         sur ={[1 2 3 4 6 7 8 9],[1 4 7  3 6 9],[1 2 3 7 8 9],...
%             [1 2 4 6 8 9], [2 3 6 4 7 8],[1 3 7 9],[1 3 8]};
%         cent ={[5], [2 5 8],  [4 5 6],[3 5 7],[1 5 9],[2 4 6 9],[7 2 9]};


colormap gray(256)

%% test Focus
for t = 1:length(nams)
    
    sprintf('Running image %d of %d.',t,length(nams))
    
    %% Read Kernals
    checkFile = [TPN nams{t}];
    
    info = imfinfo(checkFile);
    xs = info.Width;
    ys = info.Height;
    dsampy = fix(xs/dsamp)+1; %down sample image reading
    dsampx = fix(ys/dsamp)+1;
    
    I1 = double(imread(checkFile,'PixelRegion',{[ 1 dsampy ys],[1 dsampx xs]}));
    Id = zeros(size(I1,1),size(I1,2),length(yshift));
    
    %grab samples
    for i = 1:length(yshift)
        Is = double(imread(checkFile,'PixelRegion',{[ 1+yshift(i) dsampy ys],[1 + xshift(i) dsampx xs]}));
        Is = Is(1:size(I1,1),1:size(I1,2));
        Id(:,:,i) = Is;
    end
    
    %%   Saturation
    L = length(Id(:));
    tooHigh = sum(Id(:)  == 255);
    tooLow = sum(Id(:) == 0);
    percentSat = (tooHigh + tooLow)/L * 100;
    
    %% Find signal
    meanId = mean(Id,3);
    Sats = Id==255;
    sumSats = sum(Sats,3);
    vals = meanId(sumSats<2);
    %image(meanId)
    sortMean = sort(vals,'descend');
    sigThresh = .20;
    threshVal = sortMean(fix(length(sortMean)*sigThresh));
    useSig = find((meanId>=threshVal) & (sumSats<2));
    image(fitH(meanId>=threshVal))
    pause(1)
    
    

    
    %%     Find contrasts
    for f = 1: length(cent)
        dif = mean(Id(:,:,cent{f}),3)-mean(Id(:,:,sur{f}),3);
        dif = abs(dif);
        difs(:,:,f) = dif;
        useDifs{f} = dif(useSig);
    end
    negCon = [useDifs{1}];
    posCon = [useDifs{2} useDifs{3}];
    maxDif = max(difs(:,:,2:3),[],3);
    image(maxDif*3),pause(.1)
    
    %% Find Focus Target
    %     dome = fspecial('gaussian',size(meanId,1),size(meanId,1));
    %     dome = dome/max(dome(:));
    %     SE = strel('disk',round(size(dome,1)/5),4);
    %     mesa = imdilate(dome,SE);
    %     image(fitH(mesa))
    %     plot(mesa(round(size(mesa,1)/2),:))
    %     ylim([0 max(mesa(:))])
    %     image(fitH(mesa))
    %
    %
    %     bKern = ones(2,4);
    %     freqReg = fastCon(difMap,bKern);
    %     freqReg = freqReg - min(freqReg(:));
    %     freqREg = freqReg/max(freqReg(:));
    %     freqReg = freqReg.*mesa;
    %     image(fitH(freqReg))
    %     [h w] = find(freqReg == max(freqReg(:)),1);
    %     focTarg = [h/size(freqReg,1) w/size(freqReg,2)];
    %
    
    
    
    %% global contrast
    glob = meanId;
    gY = abs(glob(2:end,:)-glob(1:end-1,:)) ;
    gX = abs(glob(:,2:end)-glob(:,1:end-1)) ;
    spacer = zeros(size(glob,1),1);
    gAll = cat(3,[spacer' ; gY], [gY*-1 ;spacer'],[spacer gX], [gX*-1  spacer]);
    maxG = max(gAll,[],3);
    
    
    
    %% calculate quality
    randContrast = mean(negCon);
    kernContrast = mean(posCon(:));
    maxKContrast = mean(max(posCon,[],2));
    scaledDifs = mean(posCon,1)-randContrast;
    meanGlob = mean(maxG(useSig));
    globalContrast = meanGlob/randContrast;
    featureContrast = maxKContrast/randContrast;
    high2low = maxKContrast/meanGlob;
    quality = featureContrast;
    
    %% Record data
    tile(t).quality = quality;
    tile(t).percentSaturation = percentSat;
    tile(t).dat.globalContrast= globalContrast;
    tile(t).dat.featureContrast = featureContrast;
    tile(t).dat.high2low = high2low;
    tile(t).dat.maxKContrast = maxKContrast;
    
end


%% Organize data

allQual = [tile.quality];
[sortQual ranks] = sort(allQual);%sort(tile.use.quality,'ascend')

if ct>0
    
    qual.row = row;
    qual.col = col;
for i = 1:length(allQual)
    qualMos(row(i),col(i)) = allQual(i);
    satMos(row(i),col(i)) = tile(i).percentSaturation;
    globMos(row(i),col(i)) = tile(i).dat.globalContrast;
end
subplot(1,1,1)
image(fitH(qualMos))

qual.mos.qualMos = qualMos;
qual.mos.globMos = globMos;
qual.mos.satMos = satMos;

end

qual.tile = tile;
qual.info = info;
qual.nams = nams;


save([TPN 'qual.mat'],'qual')







