


%% Set variables
SPN = 'H:\Victoria\AMIRA\Tiff\';
SPN = 'F:\Victoria\AMIRA\Tiff\';

%SPN = uigetfolder('H:\Victoria\AMIRA\Tiff\')
autoThreshold = 1; %run auto threshold 1 or manual 0
ifShow =  1; % show each step
checkPrevious = 0; %check to see if files have been done before

useZ =[8:11];




%% Read images
dSPN = dir([SPN '*.tif']);


iNams = {dSPN.name}';
iNum = length(iNams);

if checkPrevious
    if exist([SPN 'dat.mat'],'file')
        previousDat = 1;
        load([SPN 'dat.mat'])

        clear prevExp prevLGN
        for i = 1:length(dat)
            prevExp{i} = dat(i).f.expNam;
            prevLgn(i) = dat(i).f.lgnId;
        end
    else
        previousDat = 0;
    end
else
    previousDat = 0;
end


%% Parse image names
disp('parsing file names')
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

%%Assign ID to each experiment name
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
[exps ex] = unique(expId);
expNams = exNam(ex);

id = expId + lgnId * 1000;
uId = unique(id);
tLgn = floor(uId/1000);
tExp = mod(uId,1000);
tExpNams = expNams(tExp);

lookUpId = zeros(1,length(uId));
lookUpId(uId) = 1:length(uId);
id = lookUpId(id)';
uId = unique(id);

%% Run each experiment
for ix = 1:length(uId)

    clf
    e = tExp(ix);
    l = tLgn(ix);
    eN = tExpNams{ix};
    sprintf('runing experiment %s, LGN %d',eN,l)

    %% check if previously processed
    notDone = 1;
    if previousDat
        for n = 1:length(prevExp)
            if strcmp(prevExp{n},eN) & (l == prevLgn(n));
                notDone = 0;
            end
        end
    end


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

        if isempty(useZ) %use subset of planes or not
            useZ = 1:zs;
        end

        us = length(useZ);

        I = zeros(Height,Width,us,'double');
        for i = 1:us
            I(:,:,i) = imread([SPN fileNames{1}],useZ(i));
        end

        %% read segementation file
        if hasImage(2)
            Sc = logical(I * 0);
            for i = 1:us
                Sc(:,:,i) = imread([SPN  fileNames{2}],useZ(i));
            end
        end
        if hasImage(3)
            Si = logical(I * 0);
            for i = 1:us
                Si(:,:,i) = imread([SPN  fileNames{3}],useZ(i));
            end
        end

        Sf = logical(I * 0);
        for i = 1:us
            Sf(:,:,i) = imread([SPN  fileNames{4}],useZ(i));
        end

        if hasImage(2) & ~hasImage(3)
            Si = Sf & ~Sc;
        elseif ~hasImage(2) & hasImage(3)
            Sc = Sf & ~Si;
        end

        clear areaF areaC areaI
        for i = 1:us
            areaF(i) = sum(sum(Sf(:,:,i)));
            areaC(i) = sum(sum(Sc(:,:,i)));
            areaI(i) = sum(sum(Si(:,:,i)));
        end

        if ifShow
            for i = 1:us

                Icol = zeros(size(Si,1),size(Si,2),3,'uint8');
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
        for i = 1:us
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
        for i = 1:us

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

        for i = 1:us
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
        for i = 1:us
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




        %% Choose cell bodies

        col = hsv(us);

        clf
        subplot(1,2,1)
        hold on
        for i = 1:us
            isZ = find(propZ==i);
            scatter(areas(isZ),maxVals(isZ),5,'filled','markerfacecolor',col(i,:),'markerfacealpha',.1)
            xlim([0 1500])
        end

        if autoThreshold
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
            %valPeak = find(maxValCount==max(maxValCount));
            valPeak = find(maxValRange==0);


            countUp = 0;
            needCount = 3;
            startCount = valPeak;
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
                    startCount = valPeak;
                end
            end
            threshIntensity = maxValRange(startCount);
            plot([threshIntensity threshIntensity],[0 max(maxValCount)])
            passMaxVals = maxVals>=threshIntensity;
            pause(.1)

            %% Find area threshold
            disp('finding area thresholds')
            clf
            areaRange = 0:2000;
            areaBinR = 50;
            areaCount = areaRange * 0;
            for a = 1:length(areaRange)
                A = areaRange(a);
                areaCount(a) = sum((areas(passMaxVals)>=(A-areaBinR)) & (areas(passMaxVals)<(A+areaBinR)));
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
            pass = (passArea);

        else

            %%Draw box around
            disp('position crosshair to bottom left of nucleus cluster')
            ch = drawcrosshair(gca);
            threshArea = ch.Position(1);
            threshIntensity = ch.Position(2);
            pass = (areas>=threshArea) & (maxVals>=threshIntensity);
        end
        passCent = centroids(pass,:);


        %% count by Z
        passPropZ = propZ(pass);
        rangeZ = 1:us;
        countDapiZ = hist(passPropZ,rangeZ);

        %% grab segmentation values
        centInd = sub2ind(size(I),round(passCent(:,2)),round(passCent(:,1))...
            ,round(propZ(pass)));

        isFull = Sf(centInd)>0;
        isContra = Sc(centInd)>0;
        isIpsi = Si(centInd)>0;


        %% Show results
        clear dapiContraZ dapiIpsiZ dapiFullZ
        for i = 1:us
            Is = Ibp(:,:,i);
            Icol = zeros(size(Is,1),size(Is,2),3,'uint8');
            Icol(:,:,2) = Is * 2 + 100;


            clf
            image(Is*3+50)
            hold on
            showFull = (passPropZ == i) & (isFull>0);
            scatter(passCent(showFull,1),passCent(showFull,2),'b','+')
            dapiFullZ(i) = sum(showFull);
            showFull = (passPropZ == i) & (isIpsi>0);
            dapiIpsiZ(i) = sum(showFull);
            scatter(passCent(showFull,1),passCent(showFull,2),'r','o')
            showFull = (passPropZ == i) &  (isContra>0);
            dapiContraZ(i) = sum(showFull);

            scatter(passCent(showFull,1),passCent(showFull,2),'g','o')
            drawnow
        end
    end

    pause(.1)

    %% show hist
    clf
    plot(dapiFullZ./areaF,'b')
    hold on
    plot(dapiContraZ./areaC,'r')
    plot(dapiIpsiZ./areaI,'g')
    pause(1)


    %% save to dat

    dat(i).f.expNam = eN;
    dat(i).f.id = ix;
    dat(i).f.expId = e;
    dat(i).f.lgnId = l;

    dat(i).hasImage = hasImage;
    dat(i).YXZ = cat(2,passCent(:,[2 1]),passPropZ);
    dat(i).areaFCI = [sum(areaF) sum(areaC) sum(areaI)];
    dat(i).countFCI = [length(isFull) length(isContra) length(isIpsi)];

    dat(i).p.maxIntensity = maxVals(pass);
    dat(i).p.area = areas(pass);
    dat(i).p.isFull = find(isFull);
    dat(i).p.isContra = find(isContra);
    dat(i).p.isIpsi = find(isIpsi);
    dat(i).p.unfilteredCentroids = centroids;

    dat(i).z.areaF = areaF;
    dat(i).z.areaC = areaC;
    dat(i).z.areaI = areaI;
    dat(i).z.dapiFullZ = dapiFullZ;
    dat(i).z.dapiContraZ = dapiContraZ;
    dat(i).z.dapiIpsiZ = dapiIpsiZ;


    save([SPN 'dat.mat'],'dat')




end
