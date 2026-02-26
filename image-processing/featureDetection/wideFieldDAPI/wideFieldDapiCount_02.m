
SPN = 'H:\Victoria\AMIRA\Tiff\';
%SPN = uigetfolder('H:\Victoria\AMIRA\Tiff\')

dSPN = dir([SPN '*.tif']);
ifShow =  0;

iNams = {dSPN.name}';
iNum = length(iNams);

%% Parse image names
lgnId = zeros(iNum,1);
typeId = zeros(iNum,1);
exNam = cell(iNum,1);
for i = 1:length(iNams)
    nam = iNams{i};
    lgnS =  regexp(nam,'_LGN');
    under = regexp(nam,'_');
    dot = regexp(nam,'.tif');
    exNam{i} = nam(1:lgnS-1);
    lgnId(i) = str2num(nam(lgnS+4));
    
    typeS = nam(under(2)+1:dot(1)-1);
    switch lower(typeS)
        case 'dapi'
            typeId(i) = 1;
        case 'contra'
            typeId(i) = 2;
        case 'ipsi'
            typeId(i) = 3;
        case 'full_vol'
            typeId(i) = 4;
    end
    
end

c = 0;
clear expId
for i = 1:length(exNam)
    nam1 = exNam{i};
    for n = 1:length(exNam)
        if strcmp(exNam{i}, exNam{n})
            expId(n) = i;
        end
    end
end

uExpId = unique(expId);
lookUpId = zeros(1,length(expId));
lookUpId(uExpId) = 1:length(uExpId);
expId = lookUpId(expId)';
exps = unique(expId);

id = expId + lgnId * 1000;
uId = unique(id);
lookUpId = zeros(1,length(uId));
lookUpId(uId) = 1:length(uId);
id = lookUpId(id)';
uId = unique(id);


for ix = 1:length(uId)
    e = exps(ix);
    l = lgnId(ix);
    hasImage = zeros(1,4);
    for f = 1:4
        targ = find((expId == e) & (lgnId == l) & (typeId == f));
        if ~isempty(targ) && (length(targ) == 1)
            fileNames{f} = iNams{targ};
            hasImage(f) = 1;
        end
    end
    
    
    if hasImage(1) & hasImage(4) & sum(hasImage(2:3))>0 % then process
        
        
        
        
        %% read image stack
        iInfo = imfinfo([SPN fileNames{1}]);
        zs = length(iInfo);
        Width = iInfo.Width;
        Height = iInfo.Height;
        
        I = zeros(Height,Width,zs,'double');
        for i = 1:zs
            I(:,:,i) = imread([SPN fileNames{1}],i);
        end
        
        %% read segementation file
        if hasImage(2)
            Sc = logical(I * 0);
            for i = 1:zs
                Sc(:,:,i) = imread([SPN  fileNames{2}],i);
            end
        end
        if hasImage(3)
            Si = logical(I * 0);
            for i = 1:zs
                Si(:,:,i) = imread([SPN  fileNames{3}],i);
            end
        end
        
        Sf = logical(I * 0);
        for i = 1:zs
            Sf(:,:,i) = imread([SPN  fileNames{4}],i);
        end
        
        if hasImage(2) & ~hasImage(3)
            Si = Sf & ~Sc;
        elseif ~hasImage(2) & hasImage(3)
            Sc = Sf & ~Si;
        end
        
        
        if ifShow
            for i = 1:zs
                
                Icol = zeros(size(S,1),size(S,2),3,'uint8');
                Icol(:,:,1) = Si(:,:,i) * 1000;
                Icol(:,:,2) = I(:,:,i) * 2;
                Icol(:,:,3) = Sc(:,:,i) * 1000;
                
                image(Icol)
                drawnow
                
            end
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
            if ifShow
                image(Is)
                drawnow
            end
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
            if ifShow
                image(Idif *10)
                drawnow
            end
        end
        
        %% Watershed
        
        for i = 1:zs
            Is = Ibp(:,:,i);
            w = watershed(-Is);
            
            Ic = zeros(size(Is,1),size(Is,2),3,'uint8');
            Ic(:,:,2) = Is*2;
            Ic(:,:,1) = w*1000;
            if ifShow
                image(Ic)
                drawnow
            end
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
        
        %% grab segmentation values
        centInd = sub2ind(size(I),round(centroids(:,2)),round(centroids(:,1))...
            ,round(propZ));
        
        isFull = Sf(centInd)>0;
        isContra = Sc(centInd)>0;
        isIpsi = Si(centInd)>0;
        
        
        
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
        
        
        
        %%Find maxVal Threshold
        clf
        maxValRange = -100:200;
         maxValBinR = 1;
         maxValCount =  maxValRange * 0;
        for a = 1:length( maxValRange)
            A =  maxValRange(a);
             maxValCount(a) = sum(( maxVals>=(A- maxValBinR)) & (maxVals<(A+ maxValBinR)));
        end
        plot( maxValRange, maxValCount)
        hold on
        valPeak = find(maxValCount==max(maxValCount));
        
        countUp = 0;
        needCount = 3;
        startCount = 0;
        for p = valPeak:max(maxValRange)-1
           if maxValCount(p)<maxValCount(p+1)
            countUp = countUp + 1;
            if ~startCount
                startCount = p;
            end
            if countUp>=needCount
                break
            end
               
           else
              countUp = 0;
              startCount = 0;
           end
        end
        threshIntensity = maxValRange(startCount);
        plot([threshIntensity threshIntensity],[0 max(maxValCount)])
        passMaxVals = maxVals>=threshIntensity;
        pause(.1)
        
        %%Find area threshold
        clf
        areaRange = 0:2000;
        areaBinR = 60;
        areaCount = areaRange * 0;
        for a = 1:length(areaRange)
            A = areaRange(a);
            areaCount(a) = sum((areas(passMaxValThresh)>=(A-areaBinR)) & (areas(passMaxValThresh)<(A+areaBinR)));
        end
        plot(areaRange,areaCount)
        hold on
        areaPeak = find(areaCount==max(areaCount));

        countUp = 0;
        needCount = 5;
        startCount = 0;
        for p = areaPeak:-1:2
           if areaCount(p)>areaCount(p+1)
            countUp = countUp + 1;
            if ~startCount
                startCount = p;
            end
            if countUp>=needCount
                break
            end
               
           else
              countUp = 0;
              startCount = 0;
           end
        end
        
        threshArea = areaRange(startCount);
        plot([threshArea threshArea],[0 max(areaCount)])
        passArea = areas>=threshArea;
        pause(.1)
%         
%         %%Draw box around
%         disp('position crosshair to bottom left of nucleus cluster')
%         ch = drawcrosshair(gca)
%         threshArea = ch.Position(1);
%         threshIntensity = ch.Position(2);
        pass = (passArea) & (passMaxVals);
        passCent = centroids(pass,:);
        
        %% Show results
        if 1%ifShow
            for i = 1:zs
                Is = Ibp(:,:,i);
                Icol = zeros(size(Is,1),size(Is,2),3,'uint8');
                Icol(:,:,2) = Is * 2 + 100;
                
                clf
                image(Is*3+50)
                hold on
                showFull = (propZ == i) & pass & (isFull>0);
                scatter(centroids(showFull,1),centroids(showFull,2),'b','+')
                
                showFull = (propZ == i) & pass & (isIpsi>0);
                scatter(centroids(showFull,1),centroids(showFull,2),'r','o')
                showFull = (propZ == i) & pass & (isContra>0);
                scatter(centroids(showFull,1),centroids(showFull,2),'g','o')
                drawnow
            end
        end
        %%
        
        dat(i).
        
        
        
    end
end
