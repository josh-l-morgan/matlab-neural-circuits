

fSize = 300;
TPN = 'Z:\Active\morganLab\DATA\Mahsa\ArtificialTrainingData\withDistractor\';

featureNum = 3;
featureRadius = 20;
featureType = 'circle';

distractorNum = 10;
distractorLength = 20;
distractorWidth = 4;
distractorType = 'ring';



imageNumber = 10000;

avoidEdges = 1;
gaussianBlur = 1;
guassianSigma = 10;
addNoise = .3;
showImage = 1;


%% make directories

imageDir = [TPN 'imageDir\'];
segDir = [TPN 'segDir\'];
if ~exist(imageDir,'dir'),mkdir(imageDir),end
if ~exist(segDir,'dir'),mkdir(segDir),end

%% make gaussian
gKern = fspecial('gaussian',guassianSigma,guassianSigma * 5);

%% make feature kernal

if strcmp(lower(featureType),'circle')

    kern = zeros(featureRadius*2+1);
    [y x] = find(kern+1);
    dists = sqrt((y-featureRadius-1).^2 + (x-featureRadius-1).^2);
    kern(:) = double(dists<=featureRadius);
    image(kern*1000)

end

%% make distractor kernal

dSource = zeros(fSize*4);
[y x] = find(dSource+1);
dists = sqrt((y-fSize*2-1).^2 + (x-fSize*2-1).^2);
for r = 1:distractorWidth:fSize
    dSource(dists>r) = r;
end
image(dSource * 255/fSize)
distractorVals = unique(dSource);
distractorVals = distractorVals((distractorVals<(fSize * .9)));

for i = 1:imageNumber

    D = zeros(fSize);
    dGrab = floor(fSize/2);
    for t = 1:distractorNum
        p = ceil(rand*length(distractorVals));
        d = distractorVals(p);
        [y x] = find(dSource==d);
        p = ceil(rand*length(y));
        y = round(y(p) + rand*fSize-fSize/2);
        x = round(x(p) + rand*fSize-fSize/2);
        kSamp = dSource(y-dGrab:y+dGrab-1,x-dGrab:x+dGrab-1);
        kSamp = kSamp==d;
        D = D+kSamp;
    end
    D = double(D>0);

    I = zeros(fSize);
    if avoidEdges
        for n = 1 : featureNum
            x = ceil(rand * (fSize - featureRadius * 2) + featureRadius);
            y = ceil(rand * (fSize - featureRadius * 2) + featureRadius);
            I(y,x) = 1;
        end
    else
        for n = 1 : featureNum
            x = ceil(rand * fSize);
            y = ceil(rand * fSize);
            I(y,x) = 1;
        end
    end

    Icon = conv2(I,kern,'same');
    Icon = double(Icon>0);
    Iseg = uint8(Icon * 255);
    If = Icon + D;
    If = double(If>0);

    if gaussianBlur
        If = imfilter(If,gKern,'same');
    end
    if addNoise>0
        nField = randn(fSize);
        If = If + nField * addNoise;
    end

    If = uint8(If * 100 + 50);
 
    if ~(mod(i,10)-1) 
        disp(sprintf('writing %d of %d',i,imageNumber))
        if showImage
            Ic = cat(3,If,Iseg * .2,If);
            image(Ic)
            drawnow
        end
    end

    %% Write image and seg
    iName = sprintf('i_%05.0f.tif',i);
    sName = sprintf('s_%05.0f.tif',i);

    imwrite(If,[imageDir iName]);
    imwrite(Iseg,[segDir sName]);



end










