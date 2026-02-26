%% Target Dir
TPN = GetMyDir;


dsamp = 100;
nams2 = GetPics(TPN);

%% Define variables

yshift = [ 0 0 0 1 1 1 2 2 2]; %shifts for extracting 3X3 samples
xshift = [ 0 1 2 0 1 2 0 1 2];


%% test Focus
colormap gray(256)
Isamps = zeros(201,201,length(nams2),'uint8');
clear difs

ct = 0;
for i = 1:length(nams2)
    nam = nams2{i}
    rs = regexp(nam,'_r');
    us= regexp(nam,'_');
    cs = regexp(nam,'-c');
    
    
    if ~isempty(rs) & ~isempty(cs)
        r = str2num(nam(rs(1)+2:cs(1)-1))
        c = str2num(nam(cs(1)+2:us(2)-1))
        ct = ct+1;
        row(ct) = r;
        col(ct) = c;
        nams{ct} = nam;
        
    else
        i
    end
end


for t = 1:length(nams)
    
    checkFile = [TPN nams{t}];
    

    
    info = imfinfo(checkFile);
    xs = info.Width;
    ys = info.Height;
    dsampy = fix(xs/dsamp)+1; %down sample image reading
    dsampx = fix(ys/dsamp)+1;



    I1 = double(imread(checkFile,'PixelRegion',{[ 1 dsampy ys],[1 dsampx xs]}));
    Id = zeros(size(I1,1),size(I1,2),length(yshift));
    Isamp(:,:,t) = imread(checkFile,'PixelRegion',{[ ys/2-100 ys/2+100],[xs/2-100 xs/2+100]});
    
    
    %grab samples
    for i = 1:length(yshift)
        Is = double(imread(checkFile,'PixelRegion',{[ 1+yshift(i) dsampy ys],[1 + xshift(i) dsampx xs]}));
        Is = 255-Is(1:size(I1,1),1:size(I1,2));
        Id(:,:,i) = Is;
    end
    
    
    %% Find signal
    meanId = mean(Id,3);
    image(meanId)
    sortMean = sort(meanId(:),'ascend');
    sigThresh = .20;
    threshVal = sortMean(fix(length(sortMean)*sigThresh));
    useSig = find(meanId<=threshVal);
    
    
    %%   Saturation
    L = length(Id(:));
    tooHigh = sum(Id(:)  == 255);
    tooLow = sum(Id(:) == 0);
    percentSat = (tooHigh + tooLow)/L * 100;



    %%     Find contrasts
    %1 %2 %3
    %4 %5 %6
    %7 %8 %9
    %
    %         I = double(imread(wif.sec(s).tile{t},'PixelRegion',...
    %             {[ round(ys/2) - 50 round(ys/2)+50],[round(xs/2)- 100 round(xs/2)+100]}));
    %         subplot(2,1,1)
    %         image(I)

    %%Other useful patterns
    %         sur ={[1 2 3 4 6 7 8 9],[1 4 7  3 6 9],[1 2 3 7 8 9],...
    %             [1 2 4 6 8 9], [2 3 6 4 7 8],[1 3 7 9],[1 3 8]};
    %         cent ={[5], [2 5 8],  [4 5 6],[3 5 7],[1 5 9],[2 4 6 9],[7 2 9]};

    sur = {[2 4 6 8], [1 2 3 6], [2 3 6 9]};
    cent = {[1 3 7 9], [9 8 7 4], [8 7 4 1]};
    box = [1 3 7 9];
    cols = ['k' 'r' 'b' 'm' 'c' 'g' '.'];
    subplot(2,1,2)
    clear difs
    for f = 1: length(cent)
        dif = mean(Id(:,:,cent{f}),3)-mean(Id(:,:,sur{f}),3);
        dif = abs(dif);
        meanDifs(f) = (mean(dif(useSig)));
        difs(:,:,f) = dif;
        %             range = min(dif(:)):30:max(dif(:));
        %             range = 0:100:2000;
        %             H(:,f) = hist(dif(:),range);%,-10:1:1000);
        %             plot(range,H(:,f),cols(f))
        %             hold on
        %             ylim([0 300])
        %             xlim([0 2000])

    end

    difMap = max(difs(:,:,2),difs(:,:,3));




    %% global contrast
    glob = mean(Id(:,:,box),3);
    gY = abs(glob(1:end-1,:) - glob(2:end,:));
    gX = abs(glob(:,1:end-1) - glob(:,2:end));
    image(fitH(gX))
    meanGlob = (mean([gX(:);gY(:)]))-meanDifs(1);

    scaledDifs = meanDifs-meanDifs(1);
    dmap = [scaledDifs(2) scaledDifs(3)];
%    mos(rc(t,1),rc(t,2),:) = [dmap(1) dmap(2) meanGlob];
%     %% Find Focus Target
%     dome = fspecial('gaussian',size(meanId,1),size(meanId,1));
%     dome = dome/max(dome(:));
%     SE = strel('disk',round(size(dome,1)/5),4);
%     mesa = imdilate(dome,SE);
%     image(fitH(mesa))
%     plot(mesa(round(size(mesa,1)/2),:))
%     ylim([0 max(mesa(:))])
%     image(fitH(mesa))
% 
% 
%     bKern = ones(2,4);
%     freqReg = fastCon(meanId,bKern);
%     freqReg = freqReg - min(freqReg(:));
%     freqREg = freqReg/max(freqReg(:));
%     freqReg = freqReg.*mesa;
%     image(fitH(freqReg))
%     [h w] = find(freqReg == max(freqReg(:)),1);
%     focTarg = [h/size(freqReg,1) w/size(freqReg,2)];


  
    pause(.01)
    %% Record data
      quality = scaledDifs(2) * scaledDifs(3);
    tile.quality(t) = quality;
    tile.percentSaturation(t) = percentSat;
    tile.meanGlob(t) = meanGlob;
    tile.meanDifs(t,:) = meanDifs;
    tile.scaledDifs(t,:) = scaledDifs;
    tile.sdOverGlob(t,:) = scaledDifs/meanGlob;
    %tile.focusTarget(t,:) = focTarg;
    %%
    tile.use.quality(t) = quality;
    tile.use.percentSaturation(t) = percentSat;
    tile.use.globalContrast(t) = meanGlob;

end


%% Display data


[sortQual ranks] = sort(tile.use.quality,'ascend')
allqual = [tile.use.quality]

for i = 1:length(allqual)
    qualMos(row(i),col(i)) = allqual(i);
end
subplot(1,1,1)
image(fitH(qualMos))

r = 1; %define starting rank
while 1
    [pressed Ax] = getMyIn;
    %figure(fig)


    if isstruct(pressed) %was there a keypress?
        key = pressed.Key
end
         button = get(gcf, 'SelectionType');
            if ~strcmp(button,'open')
                butt = button;
            end
            if strcmp(butt,'normal')
                r = r-1;
            elseif strcmp(butt,'alt') 
                r = r+1;
            end
      
    r = max(1,r);
    r = min(length(ranks),r);
    subplot(2,1,1)
    plot(sortQual),
    hold on
    scatter(r,sortQual(r),'r')
    hold off
    pause(.1)
    subplot(2,1,2)
    image(256-Isamp(:,:,ranks(r))),pause(.01)
end




% %% Save
% TPNsav = [TPN '\quality\'];
% if ~exist(TPNsav),mkdir(TPNsav);end
% save([TPNsav 'qual.mat'],'mosA')
% save([wif.dir 'focSec.mat'],'focSec')







