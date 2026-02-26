%function[qual] = checkFileQuality(checkFile)

if ~exist('checkFile','var')
    [TFN TPN] = GetMyFile;
    checkFile = [TPN TFN];
end

dsamp = 100;  %number of samples in each dimension

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

    
    info = imfinfo(checkFile);
    xs = info.Width;
    ys = info.Height;
    
    dsampy = fix(xs/dsamp)+1; %down sample image reading
    dsampx = fix(ys/dsamp)+1;
    maxx = fix((xs-max(xshift))/dsampx);
    maxy = fix((ys-max(yshift))/dsampy);
        
    %I1 = double(imread(checkFile,'PixelRegion',{[ 1 dsampy ys],[1 dsampx xs]}));
    Id = zeros(maxy,maxx,length(yshift));
    
    %grab samples
    for i = 1:length(yshift)
        Is = double(imread(checkFile,'PixelRegion',{[ 1+yshift(i) dsampy ys],[1 + xshift(i) dsampx xs]}));
        Is = Is(1:maxy,1:maxx);
        Id(:,:,i) = Is;
    end
    
    %% Anylize samples
    
    Im = mean(Id,3);
%     for i = 1:size(Id,3)
%         Idif(:,:,i) = Id(:,:,i)-Im;
%     end
%     Idev = sqrt(sum(Idif.^2,3));
    Idev = std(Id,1,3);
    devVals = sort(Idev(:),'ascend');
    lowDev = mean(devVals((1:100)));
    highDev = mean(devVals(end-10:end));
    Dev = highDev-lowDev;
    
    
    Idev = sqrt(sum((Id-Im)^2,3));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    tile.quality = quality
    tile.percentSaturation = percentSat;
    tile.dat.globalContrast= globalContrast;
    tile.dat.featureContrast = featureContrast;
    tile.dat.high2low = high2low;
    tile.dat.maxKContrast = maxKContrast;
    






