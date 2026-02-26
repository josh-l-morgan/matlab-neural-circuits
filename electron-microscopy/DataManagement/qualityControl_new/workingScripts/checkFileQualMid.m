function[tile] = checkFileQual(checkFile)

if ~exist('checkFile','var');
    [TFN TPN] = GetMyFile;
    checkFile = [TPN TFN];
end

dsamp = 100;
sigThresh = .01;

yshift = [ 0 0 0 1 1 1 2 2 2]; %shifts for extracting 3X3 samples
xshift = [ 0 1 2 0 1 2 0 1 2];

%1 %2 %3
%4 %5 %6
%7 %8 %9

sur = {[1 4 7 3 6 9 ], [1 2 3 7 8 9], [2 3 6 4 7 8], [4 1 2 8 9 6],...
    [2 4 5 6 7 9], [ 1 3 4 5 6 8], [ 1 2 5 6 7 8], [ 2 3 4 5 8 9]};
cent = {[2 5 8],[4 5 6],[ 1 5 9], [7 5 3 ], [1 8 3] , [7 2 9],...
    [3 4 9], [1 6 7]};

colormap gray(256)

info = imfinfo(checkFile);
xs = info.Width;
ys = info.Height;

dsampy = fix((ys-2)/dsamp); %down sample image reading
dsampx = fix((xs-2)/dsamp);

Id = zeros(dsamp,dsamp,length(yshift));

%grab samples
for i = 1:length(yshift)
    Is = double(imread(checkFile,'PixelRegion',{[ 1+yshift(i) dsampy ys],[1 + xshift(i) dsampx xs]}));
    %Is = I([1 + yshift(i) : 3 :ys],[1 + xshift(i) :3:xs]);
    Is = Is(1:dsamp,1:dsamp);
    Id(:,:,i) = 255-Is;
end

maxId = max(Id,[],3);
subplot(1,2,1)
image(maxId),pause(.01)
devI = std(Id,1,3);


%% Find signal
Sats = (Id==255) | (Id == 0);
percentSat = sum(Sats(:))/length(Sats(:))*100;
sumSats = sum(Sats,3);
vals = devI(sumSats<2);
sortMean = sort(vals,'ascend');
threshVal = sortMean(end + 1 - (fix(length(sortMean)*sigThresh)));
useSig = find((devI>=threshVal) & (sumSats<2));
subplot(1,2,2)
image((devI>=threshVal)*1000)
pause(.1)

%%     Find contrasts
for f = 1: length(cent)
    dif = mean(Id(:,:,cent{f}),3)-mean(Id(:,:,sur{f}),3);
    dif = abs(dif);
    difs(:,:,f) = dif;
    useDifs{f} = dif(useSig);
end

%% group difs
top = 1; % select top X of sample
groups = {1,2,3,4,5,6,7,8};
for f = 1:length(groups)
    vals = [useDifs{[groups{f}]}];
    sortVals = sort(vals(:),'ascend');
    thresh = sortVals(end+1 - (round(length(sortVals)*top)));
    topCon(f) = mean(vals(vals>=thresh));
end

%% global contrast
  glob = mean(Id,3);
% gY = abs(glob(2:end,:)-glob(1:end-1,:)) ;
% gX = abs(glob(:,2:end)-glob(:,1:end-1)) ;
% spacer = zeros(size(glob,1),1);
% gAll = cat(3,[spacer' ; gY], [gY*-1 ;spacer'],[spacer gX], [gX*-1  spacer]);
% maxG = max(gAll,[],3);

globVals = sort(glob(:),'ascend');
backGround = median(globVals(1:dsamp));
topGlobM = mean(glob(useSig));
range = topGlobM-backGround;

%% calculate quality
nonCon = mean(topCon(5:8));
ori = topCon(1:4)-nonCon;
quality = min(ori)/range*100;

%% Record data
tile.quality = quality;
tile.percentSaturation = percentSat;
tile.range = range;
tile.ori = ori;
tile.topCon = topCon;
tile.stdHigh = mean(devI(useSig))/range;
tile.std = std(Id(:));





