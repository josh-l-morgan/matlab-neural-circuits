

clear all
TPN = GetMyDir
look = 5;
steps = 20;

%% read I

dTPN = dir(TPN); dTPN = dTPN(3:length(dTPN));
clear I
for i = length(dTPN):-1:1
    nam = dTPN(i).name;
    if length(nam) > 4
        if strcmp(nam(length(nam)-3:length(nam)),'.tif')
            I(:,:,i,:) = imread([TPN nam]);
        end
    end
end
Imax = squeeze(max(I,[],3));
image(Imax), pause(.01)

bRaw = I(:,:,:,1);
bF = bRaw;
for i = 1: size(bRaw,3)
    H = fspecial('gaussian',look*2,look);
    bF(:,:,i) = imfilter(bRaw(:,:,i),H,'same');
    bF(:,:,i) = imfilter(bF(:,:,i),H,'same');
end
bF = buffI(bF,look);
image(max(bF,[],3)),pause(.01)

%% get all bipMasks
slashes = find(TPN == '\');
upOne = TPN(1:slashes(length(slashes)-1));
subFolders = findFolders(upOne);
%%Find image folders
ImageFolders={}; CellFolders = {};
ImageDir = ImageFolders;
for i = 1:length(subFolders)
    Name=subFolders{i};
    slashes = find(Name == '\');
    lastFold = Name(slashes(length(slashes))+1:length(Name));
    cellFold = Name(1:slashes(length(slashes)));
    if strcmp(lastFold(1:min(7,length(lastFold))),'bipMask')
            ImageFolders = [ImageFolders ; [Name '\']];            
    end
end

%% Run all bip mask

for i = 1: length(ImageFolders)
    %uLine = find(
    
    
        bipDir = dir(ImageFolders{i}); bipDir = bipDir(3:length(bipDir));
        TiffNames={};
        for b = 1: size(bipDir,1)
            siz=length(bipDir(b).name);
            if bipDir(b).name(siz-2:siz)== 'tif'
                if bipDir(b).name(siz-8:siz-7)~='-R'
                    TiffNames(length(TiffNames)+1,1)={bipDir(b).name};
                end
            end
        end

        clear bipMask
        for b = length(TiffNames):-1:1
            bipMask(:,:,b) = sum(imread([ImageFolders{i} TiffNames{b}]),3);
        end
        bipMask = smooth3(bipMask,'box',[5,5,1]);%,'gaussian',[5 5 3]);
        bipMask = bwperim(bipMask>.3,26);
        bipMask = buffI(bipMask,look);
        %I = double(I);
        maxBipMask= sum(bipMask,3);
        image(maxBipMask*10),pause(.01)
        

        %% Climb
        voxList = find(bipMask);
        [ys xs zs] = size(bipMask);
        colorBip = zeros(ys, xs,3,zs);
        idBip = zeros(ys,xs,zs);
        [ys xs zs] = size(bipMask); siz = [ys xs zs];
        %% Make surround
        surDim = look * 2 + 1;
        [sury surx surz] = ind2sub([surDim surDim surDim], find(ones(surDim,surDim,surDim)));
        dists = dist([sury surx surz],[look+1 look+1 look +1]);
        sury = sury(dists<=look)-look -1; 
        surx = surx(dists<=look)-look - 1; 
        surz = surz(dists<=look)-look -1;
        totV = length(voxList);
        for v = 1: totV  %run all voxes
            if ~mod(v,1000),percentDone = v/totV * 100, end
            [y x z] = ind2sub(siz, voxList(v));
            ny = y; nx = x; nz = z;
            for m = 1: steps
                neary = sury + ny; nearx = surx + nx; nearz = surz + nz;
                nearInd = sub2ind(siz,neary,nearx,nearz);
                vals = bF(nearInd);
                [ny nx nz] = ind2sub(siz,nearInd(find(vals == max(vals),1,'first')));
%                 showBip = maxBipMask;
%                 showBip(ny, nx) = 1000;
%                 image(showBip*10),pause(.01)
                
            end
           
            colorBip(y, x, :,z) = [ny nx nz];
        end
            
           
            image(uint8(squeeze(max(colorBip,[],4))))
            Iwrite([upOne 'colors2\col'],uint8(colorBip))
            allY = squeeze(colorBip(:,:,1,:)); allX=squeeze(colorBip(:,:,2,:));
            allZ = squeeze(colorBip(:,:,3,:));  
            allInd = sub2ind(siz,allY(voxList),allX(voxList),allZ(voxList));
            uniqueInd = unique(allInd);
            idBip = zeros(ys, xs,zs);
            for u = 1: length(uniqueInd)
                idBip(voxList(allInd==uniqueInd(u)))= rand * 254 + 1;
            end
            Iwrite([upOne 'IDs'],uint8(idBip))
end
            
















