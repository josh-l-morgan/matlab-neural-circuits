clear all
SPN = 'E:\SeanMcCracken\testCaSegmentation\test1\'

%% Parse Files
dSPN = dir([SPN '*.tif'])
iNams = {dSPN.name};

for i = 1:length(iNams)
   nam = iNams{i};
   zStr = regexp(nam,'_z');
   cStr = regexp(nam,'_c');
   dot = regexp(nam,'.tif');
   
   fz(i) = str2num(nam(zStr(1)+2:cStr(1)-1));
   fChan(i) =  str2num(nam(cStr(1)+2:dot(1)-1));
  
    
end

chans = unique(fChan);
planes = unique(fz);
cNum = length(chans);
zs = length(planes);

%% Read Data
Itest = imread([SPN iNams{1}]);
[ys xs] = size(Itest);

Iraw{1} = zeros(ys,xs,zs);
Iraw{2} = Iraw{1};

for i = 1:length(iNams)
    It = imread([SPN iNams{i}]);
    Iraw{fChan(i)}(:,:,fz(i)) = It;
    
end


%% Filter data
kSize = [5 10];
colormap gray(255)
for c = 1:cNum
    I = double(Iraw{c});
    image(uint8(I))
    
    Is = squeeze(max(I,[],1));
    Is = Is-mode(Is(:));
    Is = Is * 255/max(Is(:));
    image(Is);
    
    
    I1 = imgaussfilt3(I,[5 5 10]);
    I2 = imgaussfilt3(I,[7 7 12]);
    I3 = I1-I2;
    
  
    If{c} = I3;
    
    
    
end

%% Normalize

In = If;
for c = 1:cNum
   cMode(c) = mode(In{c}(:));
   In{c} = In{c}-cMode(c);
   cMax(c) = max(In{c}(:));
   In{c} = In{c} * 255/cMax(c);
end


%% Show

Is = In;
for z = 1:zs
    
    Ic(:,:,1) = Is{1}(:,:,z);
    Ic(:,:,2) = Is{2}(:,:,z);
    Ic(:,:,3) = Is{1}(:,:,z);
    image(uint8(Ic+100))
    pause(0.01)
    
    
end


%% Watershed

Im = In{1} + In{2};
Ii = max(Im(:))-Im;
Iw = watershed(Ii,26);
props = regionprops3(Iw,Im,'VoxelIdxList','MeanIntensity','VoxelValues','MaxIntensity');

Ip = Im*0;
for p = 1:length(props.VoxelIdxList);
    ind = props.VoxelIdxList{p};
    Ip(ind) = props.MaxIntensity(p);
    
end

Ip(Ip<1) = 1;
Ir = Im./Ip;
Ir(Ir<0.5) = 0;
Ir(Iw==0) = 0;
Imask = Ir>0;

Imasked = Im;
Imasked(~Imask) = 0;

Ilab = bwlabeln(Imask,26);
props = regionprops3(Ilab,Imasked,'VoxelIdxList','MeanIntensity',...
    'VoxelValues','MaxIntensity','Solidity');



colormap gray(256)
for i = 1:size(Iw,3)
   subplot(1,2,1)
   image(Im(:,:,i));
   
   subplot(1,2,2)
   image(Imasked(:,:,i)*.4)
   pause(.1)
    
end



%% itterative
%{
if 0


Im = In{1} + In{2};

vals = sort(Im(:));
low = vals(round(length(vals)*.01));
high = vals(round(length(vals)*.99));
step = (low-high)/20;


%Im = Im(55:90,360:375,:);
maxI = max(Im(:));
threshes = maxI:step:low;

Itrack = Im * 0;

cmap = hsv(256);
cmap(1,:) = 0;
colormap(cmap)
Itrack = Im*0;
growing = [];
voxNumThresh = 100; % minimum number of voxels for filtering
minorLength = 15;
%colormap gray(255)
for r = 1:length(threshes)
    
    
    It = Im >= threshes(r);
    Ilab = bwlabeln(It,26);
    props = regionprops3(Ilab,'VoxelIdxList');

    Ilab(Ilab>0) = Ilab(Ilab>0) + max(Itrack(:));
    
    if    length(props.VoxelIdxList)>10
        return
    end
        

    
    for p = 1:length(props.VoxelIdxList)
        
        ind = props.VoxelIdxList{p};
        shouldReplace = 1;
        if length(ind) >voxNumThresh;  %analyze shape
            [y3 x3 z3] = ind2sub([ys xs zs],ind);
            y3 = y3 - min(y3)+1;
            x3 = x3 - min(x3) + 1;
            maxX3 = max(x3);
            maxY3 = max(y3);
            clear I2D;
            I2D(y3,x3) = 1;
            
            props2D = regionprops(I2D,'MinorAxisLength');
            if props2D.MinorAxisLength > minorLength
                Im(ind) = 0;
                shoulReplace = 0;
            end
        end
        
        
         %Replace label with tracked ID
        if shouldReplace
         vals = Itrack(ind);
            minVal = min(vals(vals>0));
            if ~isempty(minVal)
                Itrack(ind) = minVal;
            end
        end
    end
%     
%     for g = 1:length(growing)
%        n = Ilab(Itrack == growing(g));
%        change = Ilab==n(1);
%        Ilab(change) = 0;
%        Itrack(change) = growing(g);
%
%     end

%%Shift all remaining labels to Itrack
Itrack(Itrack==0) = Ilab(Itrack ==0);

Ishow = Itrack;
Ishow(Ishow==0) = inf;
Ishow = 256 - Ishow;

subplot(1,2,1)
image(squeeze(max(Ishow,[],3)));
%image(squeeze(sum(Itrack,3)*10/max(Itrack(:))));

    subplot(1,2,2)
    image(max(Im,[],3));
drawnow

%     maxInd = find(Im == max(Im(:)),1);
%     [y x z] = ind2sub(size(Im),maxInd);
%
    
end
%}


end

