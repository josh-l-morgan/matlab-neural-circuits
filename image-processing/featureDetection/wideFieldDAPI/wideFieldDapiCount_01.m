[SFN SPN] = uigetfile('*.*')


[segFN segPN] = uigetfile('*.tif')

%% read image stack
iInfo = imfinfo([SPN SFN])
zs = length(iInfo);

clear I 
for i = 1:zs
    I(:,:,i) = imread([SPN  SFN],i);
end
I = double(I);


%% read segementation file
iInfoSeg = imfinfo([segPN segFN])
zsSeg = length(iInfo);

clear Sf 
for i = 1:zsSeg
    Sf(:,:,i) = imread([segPN  segFN],i);
end
Sf = double(Sf);
Si = Sf* 0;
Sc = Sf*0;

for i = 1:zsSeg
   
    Icol = zeros(size(S,1),size(S,2),3,'uint8');
    Icol(:,:,1) = S(:,:,i) * 1000;
    Icol(:,:,2) = I(:,:,i) * 2;
    image(Icol)
    drawnow
    
end

%% Normalize brightness
maxI = max(I(:));
In = I; %Normalized stack
for i = 1:zs
    Is = I(:,:,i);
    vals = Is(Is>0);
    sVals = sort(vals,'ascend');
    high5 = sVals(round(0.95 * length(vals)));
    v = var(vals);
    m = mode(vals);
    Is = Is-m;
    Is = Is * 100/high5;
    Is = Is + 100;
    image(Is)
    drawnow
    In(:,:,i) = Is;
end

%% Band pass filter
smallLimit = 5;
bigLimit = 15;

kernS = fspecial('gaussian',smallLimit*3,smallLimit);
kernB = fspecial('gaussian',bigLimit*3,bigLimit);
Ibp = In; %bandpass stack
for i = 1:zs
   
    Is = In(:,:,i);
    
    Ismall = imfilter(Is,kernS);
    Ibig = imfilter(Is,kernB);
    Idif = Ismall-Ibig ;
    Ibp(:,:,i) = Idif;
    image(Idif *100)
    drawnow
end

%% Watershed

for i = 1:zs
    Is = Ibp(:,:,i);
    w = watershed(-Is);
    
    Ic = zeros(size(Is,1),size(Is,2),3,'uint8');
    Ic(:,:,2) = Is*2;
    Ic(:,:,1) = w*1000;
    image(Ic)
    drawnow
    Iw(:,:,i) = w;
    
end


%% Feature extraction
areas = [];
maxVals = [];
centroids = [];
propZ = [];
for i = 1:zs
    Is = Iw(:,:,i);
    Iv = Ibp(:,:,i);
    
    props = regionprops(Is,Iv,'Area','MaxIntensity','WeightedCentroid','Centroid');
    
    L = length(areas);
    propNum = length(props);
    
    areas(L+1:L+propNum,:) = [props.Area]';
    maxVals(L+1:L+propNum,:) = [props.MaxIntensity]';
    centroids(L+1:L+propNum,:) = cat(1,props.Centroid);
    propZ(L+1:L+propNum,:) = ones(propNum,1)*i;
    
end

%% grab segmentation valuse
centInd = sub2ind(size(I),round(centroids(:,1)),round(centroids(:,2))...
    ,round(propZ));

isFull = Sf(centInd);
isContra = Sc(centInd);
isIpsi = Si(centInd);



%% Choose cell bodies


col = hsv(zs);

clf
subplot(1,2,1)
hold on
for i = 1:zs
    isZ = find(propZ==i);
    scatter(areas(isZ),maxVals(isZ),5,'filled','markerfacecolor',col(i,:),'markerfacealpha',.1)
    xlim([0 1500])
    
end

%%Draw box around
disp('position crosshair to bottom left of nucleus cluster')
ch = drawcrosshair(gca)
threshArea = ch.Position(1);
threshIntensity = ch.Position(2);
pass = (areas>= threshArea) & ( maxVals >= threshIntensity);
passCent = centroids(pass,:);

%% Show results
for i = 1:zs
    Is = Ibp(:,:,i);
    Icol = zeros(size(Is,1),size(Is,2),3,'uint8');
    Icol(:,:,2) = Is * 2 + 100;
    
    clf
    image(Is*3+50)
    hold on
    showFull = (propZ == i) & pass& (isFull>0);
    scatter(centroids(showFull,1),centroids(showFull,2),'r','o')
    pause
end
%%
