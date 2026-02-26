%
imageDir{1} = GetMyDir
imageDir{2} = GetMyDir
resultDir = GetMyDir

for d = 1:length(imageDir)
    
    foundImages = dir([imageDir{d} '*.png']);
    imageNames = {foundImages.name};
    
    for i = 1:length(imageNames)
        
        nam = imageNames{i};
        
        sPos = regexp(nam,'_s');
        yPos = regexp(nam,'_Y');
        xPos = regexp(nam,'_X');
        dotPos = regexp(nam,'.png');
        
        if isempty(yPos)
            parseTiles(i).name = nam;
            parseTiles(i).z = str2double(nam(sPos+2:dotPos-1));
            parseTiles(i).col = 0;
            parseTiles(i).row = 0;
            
        else
            parseTiles(i).name = nam;
            parseTiles(i).z = str2double(nam(sPos+2:yPos-1));
            parseTiles(i).col = str2double(nam(xPos+2:dotPos-1));
            parseTiles(i).row = str2double(nam(yPos+2:xPos-1));
        end
        
    end
    allDir(d).parseTiles = parseTiles;
end





%% Compare to first selected directory

nam1 = [imageDir{1} allDir(1).parseTiles(1).name];
        rawImage = double(imread(nam1));
        [ys xs zs] = size(rawImage);
        recVol1 = zeros(ys,xs,zs,'uint8');
        recVol2 = zeros(ys,xs,zs,'uint8');


tilePos2 = cat(2,[allDir(2).parseTiles.row]',[allDir(2).parseTiles.col]',[allDir(2).parseTiles.z]');
parfor i = 1:length(allDir(1).parseTiles)

    disp(sprintf('comparing image %d of %d',i,length(allDir(1).parseTiles)))
        if ~exist([resultDir allDir(1).parseTiles(i).name]);

    targ = find((tilePos2(:,1) == allDir(2).parseTiles(i).row) & (tilePos2(:,2) == allDir(2).parseTiles(i).col) & ...
        (tilePos2(:,3) == allDir(2).parseTiles(i).z));
    if ~isempty(targ)
        nam1 = [imageDir{1} allDir(1).parseTiles(i).name];
        rawImage = double(imread(nam1));
        image1 = rawImage(:,:,1)*256^2 + rawImage(:,:,2)*256 + rawImage(:,:,3);
        
        nam2 = [imageDir{2} allDir(2).parseTiles(targ).name];
        rawImage = double(imread(nam2));
        image2 = rawImage(:,:,1)*256^2 + rawImage(:,:,2)*256 + rawImage(:,:,3);
        
        difImage = abs((image2>0)-(image1>0));
        sumDif = sum(difImage(:)>0);
      
        if sumDif
        tempI = difImage*0;
        tempI(logical(difImage)) = image1(logical(difImage));
       
        resName  = [resultDir  allDir(1).parseTiles(i).name];
        imwrite(tempI,resName)
        end
        
        
        %recVol1(:,:,i) = tempI;
        
%          tempI = difImage*0;
%         tempI(logical(difImage)) = image2(logical(difImage));
%         recVol2(:,:,i) = tempI;
        
    end
    end
end



%% get point list

clear pointList
pointList(1).ids = cat(1,sparseTiles.id1);
pointList(2).ids = cat(1,sparseTiles.id2);
I_info=imfinfo([imageDir{1} allDir(1).parseTiles(i).name]);
Height = I_info.Height;
Width = I_info.Width;

numVox = length(pointList(1).ids);
if numVox
    pointList(1).Z = zeros(numVox,1);
    pointList(1).X = zeros(numVox,1);
    pointList(1).Y = zeros(numVox,1);
    
    c = 0;
    useTiles = 1:length(sparseTiles);
    for i = 1:length(sparseTiles)
        p = useTiles(i);
        [y x] = ind2sub([Height Width],sparseTiles(i).ind);
        newCount = c + length(y);
        pointList(1).X(c+1:newCount) = Width * (parseTiles(p).col) + x;
        pointList(1).Y(c+1:newCount) = Height * (parseTiles(p).row) + y;
        pointList(1).Z(c+1:newCount) = parseTiles(p).z;
        c = newCount;
    end
    
    
    [pointList(2).ids idx] = sort(pointList(2).ids,'ascend');
    pointList(2).X = pointList(1).X(idx);
    pointList(2).Y = pointList(1).Y(idx);
    pointList(2).Z = pointList(1).Z(idx);
    
    [pointList(1).ids idx] = sort(pointList(1).ids,'ascend');
    pointList(1).X = pointList(1).X(idx);
    pointList(1).Y = pointList(1).Y(idx);
    pointList(1).Z = pointList(1).Z(idx);
%     
%     
%     range.minX = min(range.minX,min(pointList(1).X));
%     range.minY = min(range.minY,min(pointList(1).Y));
%     range.minZ = min(range.minZ,min(pointList(1).Z));
%     range.maxX = max(range.maxX,max(pointList(1).X));
%     range.maxY = max(range.maxY,max(pointList(1).Y));
%     range.maxZ = max(range.maxZ,max(pointList(1).Z));
    
    
    
    %% get vastOb
    
    for r = 1:2
        trackIds = [];
        clear subs vastOb
        uniqueIds = unique(pointList(r).ids);
        uniqueIds = uniqueIds(uniqueIds>0);
        if length(uniqueIds)>1
            countIds = hist(pointList(r).ids,uniqueIds);
        else
            countIds = length(pointList(r).ids);
        end
        %trackIds(uniqueIds) = trackIds(uniqueIds)+ countIds;
        idStop = cumsum(countIds);
        idStart = [0 idStop(1:end-1)]+1;
        
        pointList = pointList; % lets parfor know pontList is a variable
        
        parfor i = 1:length(uniqueIds);
            
            subs{i} = cat(2,pointList(r).Y(idStart(i):idStop(i)),...
                pointList(r).X(idStart(i):idStop(i)), ...
                pointList(r).Z(idStart(i):idStop(i)));
        end
        
        vastOb(r).uniqueIds = uniqueIds;
        vastOb(r).countIds = countIds;
        vastOb(r).size = [max(pointList(r).Y) max(pointList(r).X) max(pointList(r).Z)];
        vastOb(r).subs = subs;
        
        %save([PPN newName],'vastOb')
        stopTime = datenum(clock);
        %hours = ((stopTime-startTime)*24)*(numPlanes-plane);
        
        %disp(sprintf('read plane %d of %d, %5.2f hours left',plane,numPlanes,hours))
    end
    
    
    
    
%     tileInfo.range = range;
%     tileInfo.trackIds = trackIds;
    %tileInfo.PPN = PPN;
    
    %save([TPN 'tileInfo.mat'],'tileInfo')
    %%
    
    fsize = vastOb(r).size;
    dims = [1 2];
    I = zeros(fsize(dims));
    lookUpIds = unique(cat(1,vastOb.uniqueIds));
    
    colMap = hsv(256);
    col = colMap(ceil((1:length(lookUpIds))*256/length(lookUpIds)),:);
    col = col(randperm(size(col,1)),:);
    
    for r = 1:2
        IcSum = zeros([fsize(dims) 3]);
        
        for o = 1:length(vastOb(r).subs)
            sub = vastOb(r).subs{o};
            
            
            inds = sub2ind(fsize(dims),sub(:,dims(1)),sub(:,dims(2)));
            uinds = unique(inds);
            
            if length(uinds)>1
                hinds = hist(inds,uinds);
            else
                hinds = length(inds);
            end
            
            I(uinds) = I(uinds) + hinds';
            contrastFactor = 20;
            minInt = 50;
            I = I * contrastFactor;
            I((I>0)) = I(I>0) + minInt;
            I(I>255) = 255;
            I((I>0) & (I<minInt)) = minInt;
                
            
            
            
            IcSum(:,:,1) = IcSum(:,:,1) + I*col(o,1);
            IcSum(:,:,2) = IcSum(:,:,2) + I*col(o,2);
            IcSum(:,:,3) = IcSum(:,:,3) + I*col(o,3);
            
        end
        IcAll{r} = IcSum;
        image(uint8(IcSum ))
        pause(.1)
    end
    subplot(2,1,1)
    image(uint8(IcAll{1}))
    subplot(2,1,2)
    image(uint8(IcAll{2}))
else
    disp('no disparate voxels found')
end
    sum(IcAll{1}(:))
        sum(IcAll{2}(:))

    
    