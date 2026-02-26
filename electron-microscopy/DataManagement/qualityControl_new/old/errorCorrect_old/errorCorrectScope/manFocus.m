%function[] = manFocus(wif);

if ~exist('wif','var')
    wif = GetMyWafer;
end
fig = gcf;
sid = wif.secID;
colormap gray(256)
xs = wif.imfinfo(1).Width;
ys = wif.imfinfo(1).Height;
q = 0;
readI = 1; readS = 1;
maxRes = 500;

if exist([wif.dir 'focSec.mat'],'file')
    load([wif.dir 'focSec.mat'])
else
    'Checking image quality'
    focSec = checkFocus(wif)
end


brightness = -90;
contrast = 1.5;
button = 1;
t = 1; tNew = 1;
isize = 400;

maxSize = wif.imfinfo(1).Height;
subReg = [round(xs/2 - isize/2) round(xs/2 - isize/2) + isize-1];
subRegX = subReg; subRegY = subReg;
waffer = 1;
strip = 1;
section = 1;
if size(sid,2) == 3
    s = find((sid(:,1) == waffer)& (sid(:,2) == strip) & (sid(:,3)== section));
else
    s = 1;
end
secNam = wif.secNam{s};
maxSec = length(wif.secNam);

numTile = length(wif.sec(s).tile);
good = zeros(1,numTile); rate = good;
I = 255-imread(wif.sec(s).tile{t},'PixelRegion',{subRegY,subRegX});

downSamp = imread([wif.sec(s).ov]);  %Get overview image
%downSamp = imread([wif.dir(1:end-1) 'shaped\downsamp\' wif.secNam{s} '.tif']);

%% Get images of mosic quality
load([wif.dir(1:end-1) 'shaped\quality\qual.mat']);
showMos = mosA;
datDims = size(mosA);
rNum = datDims(1);
cNum = datDims(2);
sNum = datDims(3);
tileID = (sNum-1)*cNum*rNum + ((rNum-1)*cNum)+cNum;
rL = sNum * cNum * rNum;
for c = 1:3; %scale colors
    vals = showMos(:,:,:,c);
    showMos(:,:,:,c) = 255 * (showMos(:,:,:,c)-min(vals(:)))/(max(vals(:))/min(vals(:)));
end
showMos = uint8(showMos);
qualSamp = showMos(:,:,s,:);
qualSamp = squeeze(qualSamp);

%% get previous quality values
if exist([wif.dir 'userEval.mat'],'file')
    load([wif.dir 'userEval.mat'])
    if isfield(userEval,'qualVals')
        qualVals = userEval.qualVals;
    else
        qualVals = zeros(cNum,rNum,sNum);
    end
else
    qualVals = zeros(cNum,rNum,sNum);
end
[rates fails thresh] = findThresh(mosA,qualVals);

%qualVals = zeros(cNum,rNum,sNum);
%Display
downSampAx = subplot(2,3,1);
image(downSamp)
mosAx = subplot(2,3,4);
image(qualSamp), hold on
scatter(wif.sec(s).rc(t,2),wif.sec(s).rc(t,1));
failed = fails{s};
rated = rates{s};
scatter(rated(:,2),rated(:,1),'*','g')
scatter(failed(:,2),failed(:,1),'x','r')
%scatter(rL(rL == 1,1)+.4,rL(rL==2,2),'.','g')
imageAx = subplot(2,3,[2 3 5 6]);
image((I + brightness) * contrast)
rc = wif.sec(s).rc;

%% Start
while ~q  % While not quitting
    figure(fig)
    mouseOK = 1;

    %pause(.03)
    [pressed Ax] = getMyIn;
    
    figure(fig)


    if isstruct(pressed) %was there a keypress?
        key = pressed.Key;
        changeC = 0;
        if strcmp(key,'w'),      subRegY = subRegY - size(I,1)/5; readI = 1;
        elseif strcmp(key,'s'),  subRegY = subRegY + size(I,1)/5; readI = 1;
        elseif strcmp(key,'a'),  subRegX = subRegX - size(I,2)/5; readI = 1;
        elseif strcmp(key,'d'),  subRegX = subRegX + size(I,2)/5; readI = 1;
        elseif strcmp(key,'q') | strcmp(key,'leftbracket'), %zoom out
            mY = mean(subRegY); dif =diff(subRegY);
            subRegY = round([mY-dif mY+dif]);
            mX = mean(subRegX); dif =diff(subRegX);
            subRegX = round([mX-dif mX+dif]);
            readI = 1;
        elseif strcmp(key,'e') | strcmp(key,'rightbracket'), %zoom out
            mY = mean(subRegY); dif =diff(subRegY);
            subRegY = round([mY-dif/2.5 mY+dif/4]);
            mX = mean(subRegX); dif =diff(subRegX);
            subRegX = round([mX-dif/2.5 mX+dif/4]);
            readI = 1;
        elseif strcmp(key,'c'),subRegX = subReg; subRegY = subReg; readI = 1; %reset zoom
        elseif strcmp(key,'f'),readI = 1; %point to quality mosaic
            [c r] = ginput(1);
            t = find((rc(:,1) ==round(r))&(rc(:,2)==round(c)));
        elseif strcmp(key,'downarrow'), brightness = brightness - 10;
        elseif strcmp(key,'uparrow'), brightness = brightness + 10;
        elseif strcmp(key,'rightarrow'), contrast = contrast + .25;
        elseif strcmp(key,'leftarrow'), contrast = contrast - .25;
        elseif strcmp(key,'z'),s = s-1; %move back one section
            if s<1;s = 1; end;   readI = 1; readS = 1;
            rc = wif.sec(s).rc;
        elseif strcmp(key,'x'),s = s+1; %move forward one section
            if s>maxSec, s = maxSec; end; readI = 1; readS = 1;
            rc = wif.sec(s).rc;
        elseif strcmp(key,'space'), good(t) = 1;'good', t = t+1; readI = 1;
            userEval.sec(s).tile(t).retake = 0;
            rL(sub2ind([sNum rNum cNum],s,rc(t,1),rc(t,2))) = 0;
        elseif strcmp(key,'shift'), good(t) = 2;'bad',t = t+1; readI = 1;
            userEval.sec(s).tile(t).retake = 1;
            rL(sub2ind([sNum rNum cNum],s,rc(t,1),rc(t,2))) = 1;
        elseif ~isempty(str2num(key)),  readS = 1;
            rating = str2double(key);
            sprintf('you rated the image quality a %d out of 4',rating)
            qualVals(rc(t,1),rc(t,2),s) = rating;
            userEval.qualVals = qualVals;
            [rates fails thresh] = findThresh(mosA,qualVals);
            userEval.thresh = thresh;
            userEval.rates = rates;
            userEval.fails = fails;
            save([wif.dir 'userEval.mat'],'userEval');
        elseif strcmp(key,'tab'),
            escape = verify;
            if escape
                break
            end

        elseif strcmp(key,'backslash'), readS
            qualVals = zeros(cNum,rNum,sNum)
        elseif strcmp(key,'capslock') | strcmp(key,'control') ,'saving'

        else

        end

        if readS
            [rates fails thresh] = findThresh(mosA,qualVals);
            qualSamp = squeeze(showMos(:,:,s,:));
            subplot(2,3,4),image(qualSamp),hold on
            scatter(wif.sec(s).rc(t,2),wif.sec(s).rc(t,1),'w');
            scatter(wif.sec(s).rc(t,2),wif.sec(s).rc(t,1),'.','k');
            failed = fails{s};
            rated = rates{s};
            scatter(rated(:,2),rated(:,1),'*','g')
            scatter(failed(:,2),failed(:,1),'x','r')
            pause(.01)
            hold off

        end


        if  readI

            subplot(2,3,[2 3 5 6])
            if t<1,t=1;end
            if t> numTile, t = numTile;end
            subRegY = round(subRegY);subRegX = round(subRegX);
            if subRegY(1)<0; subRegY = subRegY-subRegY(1)+1;   end
            if subRegX(1)<0; subRegX = subRegX-subRegX(1)+1;   end
            if subRegY(2)>maxSize; subRegY = subRegY+(maxSize - subRegY(2));   end
            if subRegX(2)>maxSize; subRegX = subRegX+(maxSize - subRegX(2));   end
            if diff(subRegY)>maxSize; subRegY(1) = 1; subRegX(1)=1; end


            rsize = diff(subRegY);
            if rsize>maxRes
                useRegY = [subRegY(1) fix(rsize/maxRes)+1 subRegY(2)];
                useRegX = [subRegX(1) fix(rsize/maxRes)+1 subRegX(2)];
            else
                useRegY = subRegY; useRegX = subRegX;
            end
            I = 255-imread(wif.sec(s).tile{t},'PixelRegion',{useRegY,useRegX});
            readI = 0;

        end

        subplot(2,3,[2 3 5 6])
        image((I + brightness) * contrast)
        pause(.01)

    else
        
        if length(pressed)>1
            inP = pressed(1,:);
        if Ax == mosAx %if clicked on mosaic
            newt = find((rc(:,1) ==round(inP(2)))&(rc(:,2)==round(inP(1))));
        elseif Ax == downSampAx
            c = fix(cNum * inP(1)/size(downSamp,2))+1;
            r = fix(rNum * inP(2)/size(downSamp,1))+1;
            newt = find((rc(:,1) ==r)&(rc(:,2)==c));

        else
            button = get(gcf, 'SelectionType');
            if ~strcmp(button,'open')
                butt = button;
            end
            if strcmp(butt,'normal') %change colors
                newt = t-1;
            elseif strcmp(butt,'alt') %edit labels
                newt = t+1;
            end
        end%which axis clicked
        end
        if (newt>0) & (newt<numTile)
            t = newt;
        end

        rsize = diff(subRegY);
        if rsize>maxRes
            useRegY = [subRegY(1) fix(rsize/maxRes)+1 subRegY(2)];
            useRegX = [subRegX(1) fix(rsize/maxRes)+1 subRegX(2)];
        else
            useRegY = subRegY; useRegX = subRegX;
        end

        qualSamp = squeeze(showMos(:,:,s,:));
        subplot(2,3,4),image(qualSamp),hold on
        scatter(wif.sec(s).rc(t,2),wif.sec(s).rc(t,1),'w');
        scatter(wif.sec(s).rc(t,2),wif.sec(s).rc(t,1),'.','k');
        failed = fails{s};
        rated = rates{s};
        scatter(rated(:,2),rated(:,1),'*','g')
        scatter(failed(:,2),failed(:,1),'x','r')
        pause(.01)
        hold off

        subplot(2,3,[2 3 5 6])
        I = 255-imread(wif.sec(s).tile{t},'PixelRegion',{useRegY,useRegX});
        image((I + brightness) * contrast ) %.3 sec
        use = focSec(s).tile(t).use
        pause(.01)

        %}

    end
end





