function[] = tif2point(SPN)

%SPN = 'D:\Joshm\S8\export_TestRead\';
TPN = [SPN(1:end-1) '_mat\'];
if ~exist(TPN),mkdir(TPN),end

PPN = [TPN 'obPlanes\'];
if ~exist(PPN),mkdir(PPN),end


%% Parse name and build structure

if 1% ~exist('parseTiles','var')
    
    if 1 %~exist([TPN 'parseTiles.mat'],'file')
        
        dSPN = dir([SPN '*.png']);
        imageNames = {dSPN.name}';
        
        parfor i = 1:length(imageNames)
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
        
        save([TPN 'parseTiles.mat'],'parseTiles');
        
    else
        
        load([TPN 'parseTiles.mat']);
        
    end
end
%% read in images

allZ = [parseTiles.z];
uniqueZ = unique(allZ);
trackIds = zeros(1,100000);

numPlanes = length(uniqueZ)'

range.minX = 10000000; range.minY = 10000000; range.minZ = 10000000;
range.maxX = 0; range.maxY = 0; range.maxZ = 0;
imageInfo = imfinfo([SPN parseTiles(1).name]);
Height = imageInfo.Height;
Width = imageInfo.Width;



for plane = 1:numPlanes
    
        z = uniqueZ(plane);

        newName = sprintf('%05.0f.mat',z);
        disp(newName)
        
if 1%~exist([PPN newName],'file')
    
    
    startTime = datenum(clock);
    %% get sparse tile list
  
    useTiles = find(allZ == z);
    clear sparseTiles
    for t = 1:length(useTiles) %% Could be parfor
        p = useTiles(t);
        Iraw = double(imread([SPN parseTiles(p).name]));
        if size(Iraw,3)==3;
            I = Iraw(:,:,1) * 256^2 + Iraw(:,:,2)*256 + Iraw(:,:,3);
        else
            I = Iraw(:,:,1);
        end
        ind = find(I>0);
        sparseTiles(t).ind = single(ind);
        sparseTiles(t).id = single(I(ind));
    end
    
    
    %% get point list
   
    clear pointList
    pointList.ids = cat(1,sparseTiles.id);
    numVox = length(pointList.ids);
    if numVox
    pointList.Z = zeros(numVox,1);
    pointList.X = zeros(numVox,1);
    pointList.Y = zeros(numVox,1);
    
    c = 0;
    for i = 1:length(sparseTiles)
        p = useTiles(i);
        [y x] = ind2sub([Height Width],sparseTiles(i).ind);
        newCount = c + length(y);
        pointList.X(c+1:newCount) = Width * (parseTiles(p).col) + x;
        pointList.Y(c+1:newCount) = Height * (parseTiles(p).row) + y;
        pointList.Z(c+1:newCount) = parseTiles(p).z;
        c = newCount;
    end
    
    
    [pointList.ids idx] = sort(pointList.ids,'ascend');
    pointList.X = pointList.X(idx);
    pointList.Y = pointList.Y(idx);
    pointList.Z = pointList.Z(idx);
    
    
    range.minX = min(range.minX,min(pointList.X));
    range.minY = min(range.minY,min(pointList.Y));
    range.minZ = min(range.minZ,min(pointList.Z));
    range.maxX = max(range.maxX,max(pointList.X));
    range.maxY = max(range.maxY,max(pointList.Y));
    range.maxZ = max(range.maxZ,max(pointList.Z));

        
    
  %% get vastOb
  
    clear subs vastOb
    uniqueIds = unique(pointList.ids);
    if length(uniqueIds)>1
        countIds = hist(pointList.ids,uniqueIds);
    else
        countIds = length(pointList.ids);
    end
    trackIds(uniqueIds) = trackIds(uniqueIds)+ countIds;
    idStop = cumsum(countIds);
    idStart = [0 idStop(1:end-1)]+1;

    pointList = pointList; % lets parfor know pontList is a variable

    parfor i = 1:length(uniqueIds);

        subs{i} = cat(2,pointList.Y(idStart(i):idStop(i)),...
            pointList.X(idStart(i):idStop(i)), ...
            pointList.Z(idStart(i):idStop(i)));
    end

    vastOb.uniqueIds = uniqueIds;
    vastOb.countIds = countIds;
    vastOb.size = [max(pointList.Y) max(pointList.X) min(pointList.Z)];    
    vastOb.subs = subs;
    
    save([PPN newName],'vastOb')
    stopTime = datenum(clock);
    hours = ((stopTime-startTime)*24)*(numPlanes-plane);
    
    disp(sprintf('read plane %d of %d, %5.2f hours left',plane,numPlanes,hours))
    end
end
   
        
end

tileInfo.range = range;
tileInfo.trackIds = trackIds;
tileInfo.PPN = PPN;

save([TPN 'tileInfo.mat'],'tileInfo')



