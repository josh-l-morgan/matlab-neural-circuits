

SPN = 'D:\LGNs1\testS8segmentation\s8Sample\'
FPN = 'D:\LGNs1\testS8segmentation\s8Sample_mean2maxDS\'

downSampXY = 4;
downSampZ = 4;

%% Check file
if ~exist(FPN,'dir'),mkdir(FPN);end

%Get  tifs
dSPN = dir([SPN '*.tif']);
inams = {dSPN.name};
inams = sort(inams');

%% plan reshape

imageInfo = imfinfo([SPN inams{1}]);
W = imageInfo.Width;
H = imageInfo.Height;
xs = fix(W/downSampXY);
ys = fix(H/downSampXY);
zs = fix(length(inams)/downSampZ);

%% Read data
tic
for nz = 1:D;
    disp(' ')
    if ~exist([FPN newName],'file')
        disp(sprintf('running %d of %d',nz,D))
        
        %disp('read data')
        smallI = zeros(ys,xs,downSampZ,'uint8');
        parfor stepZ = 1:downSampZ
            getZ = (nz-1)*downSampZ + stepZ;
            rawI = imread([SPN inams{getZ}]);
            smallI(:,:,stepZ) = imresize(rawI,1/downSampXY,'method','bicubic');
        end
        
        %disp('ditch outliers')
        filtI = smallI;
        
        minI = min(smallI,[],3);
        maxI = max(smallI,[],3);
        
        isLow = smallI == repmat(minI,[1 1 downSampZ]);
        isHigh = smallI == repmat(maxI,[1 1 downSampZ]);
        isMid = ~isLow & ~isHigh;
        sumMid = sum(isMid,3);
        keepAnyway = repmat(sumMid<1,[1 1 downSampZ]);
        useVox = keepAnyway | isMid;
        scaleVox = useVox./repmat(sum(useVox,3),[1 1 downSampZ]);
        
        meanI =  sum(double(smallI) .* scaleVox,3);
        
        %disp('write image')
        newName = sprintf('filt%05.0f.tif',nz);
        imwrite(uint8(meanI),[FPN newName])
    end
end
toc
disp('final time')






