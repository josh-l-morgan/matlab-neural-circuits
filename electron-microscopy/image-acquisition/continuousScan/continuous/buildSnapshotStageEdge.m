%function[] = builSnapshot(TPN)
%%Take fast image of entire stage area

if ~exist('TPN','var')
    TPN = GetMyDir; %% Get directory to place all images
end
rawDir = [TPN 'raw\'];
isoDir = [TPN 'isoDir\'];
BEDir = [TPN 'BEDir\'];

if ~exist(isoDir),mkdir(isoDir),end
colormap gray(255)

load([ TPN 'dat\fibSettings.mat']);

%% Collect images

nams = getPics(BEDir);
BE = {};
for i = 1:length(nams)
    nam = nams{i};
    if regexp(nam,'BEstop')
        id = str2num(nam(7:end-4));
        flip = 2-mod(id,2);
        BE{flip,id} = imread([BEDir nam]);
    elseif regexp(nam,'BEstart')
        id = str2num(nam(8:end-4));
        flip = 2-mod(id+1,2);
        BE{flip,id} = imread([BEDir nam]);
    end
end

inams = getPics(isoDir);
for i = 1:length(inams)
    nam = inams{i};
    id = str2num(nam(1:end-4));
    flip = 2-mod(id,2);
    if flip-1
        Is{i} = imread([isoDir inams{i}]);
    else
        Is{i} = flipud(imread([isoDir inams{i}]));
    end
end

%% Fix with book ends
FOV = fibSet.BEstop(1).FOV;
[ys xs] = size(BE{1,1});
res = FOV/ys;

%%Get Positions
for i = 1:length(fibSet.BEstart)
    flip = 2-mod(i+1,2);
    fibSet.BEstart(i).X
    X(flip,i) = fibSet.BEstart(i).X;
    Y(flip,i) = fibSet.BEstart(i).Y;
    flip = 2-mod(i,2);
    X(flip,i) = fibSet.BEstop(i).X;
    Y(flip,i) = fibSet.BEstop(i).Y;
    
end
Xp = X * 1000000/res;
Yp = Y * 1000000/res;

xshift = min(Xp(:))-xs/2;
yshift = min(Yp(:))-ys/2;

Xp = Xp - xshift;
Yp = Yp - yshift;

edgeI = zeros(fix(max(Yp(:))+ys/2)+1,fix(max(Xp(:))+xs/2)+1,'uint8');

for i = 1:length(Is)
    for s = 1:2
        I = BE{s,i};
        yrange = round(Yp(s,i) +1 - ys/2:Yp(s,i)+.5+ys/2);
        xrange = round(Xp(s,i) +1 - xs/2:Xp(s,i)+.5+xs/2);
        edgeI(yrange,xrange) = fliplr(I);
    end
end
image(edgeI), pause(.01)
subplot(1,1,1)
save([TPN 'dat\ring.mat'],'edgeI')
clear edgeI


%% Correct flips

maxdim = [0 0];

%shrinkFlip = Ymove/(Ymove + FOV * .75);

shrinkFlip = 1/1.08532

shrinkIt(1) = 1.086;;
shrinkIt(2) = 1;

clipEdge = 11;
%shrinkFlip = 1/shrinkFlip;
for i = 1: length(Is)
    flip = 2 - mod(i,2);
    I = Is{i};
    I = flipud(I);
    I = imresize(I,[size(I,1)*shrinkIt(flip),size(I,2)],'nearest');
    maxdim = max(maxdim,[size(I,1) size(I,2)]);
    %image(I),pause(.01)
    Is{i} = I;
end


%% align strips to book ends

for i = 1:length(Is)
    
    I = Is{i};
    top = BE{2,i};
    top = top(:,1:size(I,2));
    bottom = BE{1,i};
    bottom = bottom(:,1:size(I,2));

    [bYs bXs] = size(top);
    
    sW = 100; %define scan width
    topSamp = single(top(end-sW:end,:));
    scaleDiff = sum(abs(topSamp(:)-mean(topSamp(:))));
    clear topDif
    startSearch =  1:10:size(I,1)/2-sW;
    for ns = 1:length(startSearch)
        s = startSearch(ns);
        Isamp = single(I(s:s+sW,:));
        topDif(ns)= sum(abs(Isamp(:) - topSamp(:)))/scaleDiff;
    end
    plot(topDif)
    pause(.01)
    lowest = startSearch(find( topDif == min(topDif),1));
    newrange = lowest - 100 : lowest + 100;
    clear topDif2
    for ns = 1:length(newrange)
        s = newrange(ns);
        Isamp = single(I(s:s+sW,:));
        topDif2(ns)= sum(abs(Isamp(:) - topSamp(:)))/scaleDiff;
    end
    plot(topDif2),pause(.01)
    toplow(i) = newrange(find(topDif2 == min(topDif2),1));
    
    %%Run bottom
    bottomSamp = single(bottom(1:sW+1,:));
    scaleDiff = sum(abs(bottomSamp(:)-mean(bottomSamp(:))));
    clear bottomDif
    startSearch = 1: 10 :size(I,1)/2-sW;
    for ns = 1:length(startSearch)
        s = startSearch(ns);
        Isamp = single(I(end - s - sW: end - s,:));
        bottomDif(ns)= sum(abs(Isamp(:) - bottomSamp(:)))/scaleDiff;
        
    end
    
    plot(bottomDif,'r')
    pause(.01)
    lowest = startSearch(find( bottomDif == min(bottomDif),1));
    newrange = lowest - 100 : lowest + 100;
    clear bottomDif2
    for ns = 1:length(newrange)
        s = newrange(ns);
        Isamp = single(I(end - s - sW: end - s,:));
        bottomDif2(ns)= sum(abs(Isamp(:) - bottomSamp(:)))/scaleDiff;

    end
    plot(bottomDif2,'r'),pause(.01)
    bottomlow(i) = newrange(find(bottomDif2 == min(bottomDif2),1));

end

plot(toplow,'b')
hold on
plot(bottomlow,'r')
hold off

%% Apply shifts

stripI = zeros((fix(max(Yp(:))+ys/2)+1 ) * 2,fix(max(Xp(:))+xs/2)+1,'uint8');
[ys xs] = size(BE{1,1});
topAlign0 = ys/2 - sW;
for i = 1:length(It)
        I = Is{i};
        
        Xc = Xp(1,i); % X center
        Yref = Yp(2,i)+topAlign0-toplow(i);
        yrange = round(Yref: Yref + size(I,1)-1);
        xrange = round(Xc +1 - xs/2:Xc+.5+xs/2);
        stripI(yrange+10000,xrange) = fliplr(I);
end
stripI = stripI(10001:10000+ size(buildI,1),:);
image(stripI), pause(.01)
subplot(1,1,1)

imwrite(stripI,[TPN 'build\' 'stripI.tif'])

%% combine images
combI = buildI;
combI(buildI==0) = stripI(buildI == 0);
image(combI),pause(.01)
imwrite(combI,[TPN 'build\' 'combI.tif']);










