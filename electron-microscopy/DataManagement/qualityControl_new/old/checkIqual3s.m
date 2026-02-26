function[tile] = checkIQuality(I)

if ~exist('I','var');
    [TFN TPN] = uigetfile;
    I = imread([TPN TFN]);
end

dsamp = 3;
yshift = [ 0 0 0 1 1 1 2 2 2]; %shifts for extracting 3X3 samples
xshift = [ 0 1 2 0 1 2 0 1 2];

%1 %2 %3
%4 %5 %6
%7 %8 %9

sur = {[2 4 6 ], [1 2 3 ], [4 8 6] [4 5 6],...
    [1 5 7], [1 4 7],[2 6 8], [ 2 5 8]};
cent = {[1 5 3 ], [4 5 6], [7 5 9] [7 8 9],...
    [2 4 8], [ 2 5 8], [3 5 9], [3 6 9]};
%%Other useful patterns
%         sur ={[1 2 3 4 6 7 8 9],[1 4 7  3 6 9],[1 2 3 7 8 9],...
%             [1 2 4 6 8 9], [2 3 6 4 7 8],[1 3 7 9],[1 3 8]};
%         cent ={[5], [2 5 8],  [4 5 6],[3 5 7],[1 5 9],[2 4 6 9],[7 2 9]};


colormap gray(256)

    
    I = max(I,[],3);
    [ys xs] = size(I);
    
    dsampy = fix(ys/dsamp); %down sample image reading
    dsampx = fix(xs/dsamp);
    
    %I1 = double(imread(checkFile,'PixelRegion',{[ 1 dsampy ys],[1 dsampx xs]}));
    Id = zeros(dsampy,dsampx,length(yshift));
    
    %grab samples
    for i = 1:length(yshift)
        %Is = double(imread(checkFile,'PixelRegion',{[ 1+yshift(i) dsampy ys],[1 + xshift(i) dsampx xs]}));
        Is = I([1 + yshift(i) : 3 :ys],[1 + xshift(i) :3:xs]);
        Is = Is(1:dsampy,1:dsampx);
        Id(:,:,i) = 255-Is;
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
    
    %% Find tops
    top = .006;
    groups = {[1 3],[2 4],[5 7],[6 8]};
    for f = 1:length(groups)
       vals = [useDifs{[groups{f}]}];
       sortVals = sort(vals(:),'descend');
       thresh = sortVals(round(length(sortVals)*top));
       topCon(f) = mean(vals(vals>=thresh));
    end
    
    
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
    
   
    vertQual = topCon(2)/topCon(1);
    horzQual = topCon(4)/topCon(3);
    quality = min(vertQual,horzQual)-2;
%     
%     meanGlob = mean(maxG(useSig));
%     globalContrast = meanGlob/negCon;
%     featureContrast = mean(posCon(:));
%     high2low = featureContrast/meanGlob;
%     
    %% Record data
    tile.quality = quality;
    tile.percentSaturation = percentSat;
%     tile.dat.globalContrast= globalContrast;
%     tile.dat.featureContrast = posCon;
%     tile.dat.high2low = high2low;
%     tile.dat.randContrast = negCon;
%     
    tile.horzQual = horzQual;
    tile.vertQual = vertQual;






