%function[] = builSnapshot(TPN)
%%Take fast image of entire stage area
%clear all

if ~exist('TPN','var') 
    TPN = GetMyDir; %% Get directory to place all images
end
rawDir = [TPN 'raw\'];
isoDir = [TPN 'isoDir\'];
BEDir = [TPN 'BEDir\'];
buildDir = [TPN 'build\'];

if ~exist(isoDir),mkdir(isoDir),end
colormap gray(255)

if ~exist(buildDir),mkdir(buildDir),end
colormap gray(255)

load([ TPN 'dat\fibSettings.mat']);

%% Collect images
'Collecting Images'
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
   if mod(id,2)
        Is{i} = imread([isoDir inams{i}]);
    else
        Is{i} = flipud(imread([isoDir inams{i}]));
    end
end

%% Make book end image
FOV = fibSet.BEstop(1).FOV;
[ys xs] = size(BE{1,1});
res = FOV/ys;

%%Get Positions
for i = 1:length(fibSet.BEstart)
    flip = 2-mod(i+1,2);
    X(flip,i) = fibSet.BEstart(i).X;
    Y(flip,i) = fibSet.BEstart(i).Y;
    flip = 2-mod(i,2);
    X(flip,i) = fibSet.BEstop(i).X;
    Y(flip,i) = fibSet.BEstop(i).Y;
    
end
Xp = (.110/ res) - X * 1000000/res;
Yp = Y * 1000000/res;

xshift = min(Xp(:))-xs/2;
yshift = min(Yp(:))-ys/2;

Xp = Xp - xshift;
Yp = Yp - yshift;

H = (fix(max(Yp(:))+ys/2)+1);
W = fix(max(Xp(:))+xs/2)+1;
edgeI = zeros(H,W,'uint8');

for i = 1:length(Is)
    for s = 1:2
        I = BE{s,i};
        yrange = round(Yp(s,i) +1 - ys/2:Yp(s,i)+.5+ys/2);
        xrange = round(Xp(s,i) +1 - xs/2:Xp(s,i)+.5+xs/2);
        edgeI(yrange,xrange) =  edgeI(yrange,xrange) + fliplr(I);
    end
end
image(edgeI), pause(.01)
subplot(1,1,1)
imwrite(edgeI,[TPN 'build\' 'edgeI.tif'])
save([TPN 'dat\edgeI.mat'],'edgeI')



%% Correct flips
maxdim = [0 0];

shrinkIt(1) = 1.2 * 0.9791;
shrinkIt(2) =  1.2 * .9104 * .9798;

for i = 1: length(Is)
    flip = 2 - mod(i,2);
    I = Is{i};
    I = imresize(I,[size(I,1)*shrinkIt(flip),size(I,2)],'nearest');
    maxdim = max(maxdim,[size(I,1) size(I,2)]);
    %image(I),pause(.01)
    Is{i} = I;
end


%% align strips to book ends

for i = 1:length(Is)
    
    I = Is{i};
    gKern = gaus3d([5,5,1],2);
    
    top = BE{2,i};
 
    top = fastCon(top,gKern);
     top = top(:,1:size(I,2));
%     bottom = BE{1,i};
%     bottom = fastCon(bottom,gKern);
%     bottom = bottom(:,1:size(I,2));

    [bYs bXs] = size(top);
    
    sW = round(bYs/2); %define scan width
    topSamp = single(top(end-sW:end,:));
       image(topSamp)
    scaleDiff = sum(abs(topSamp(:)-mean(topSamp(:))));
    clear topDif
    startSearch =  1:20:size(I,1)-sW;
    for ns = 1:length(startSearch)
        s = startSearch(ns);
         Isamp = single(I(s:s+sW,:));
%         subplot(1,2,1)
%         image(Isamp)
%         subplot(1,2,2)
%         image(topSamp)
%         pause(.1)
        topDif(ns)= sum(abs(Isamp(:) - topSamp(:)))/scaleDiff;
    end
    lowest = startSearch(find( topDif == min(topDif),1));
    newrange = lowest - 20 : lowest + 20;
    newrange = max(newrange,1);
    newrange = min(newrange,size(I,1)-sW);
    clear topDif2
    for ns = 1:length(newrange)
        s = newrange(ns);
        Isamp = single(I(s:s+sW,:));
        topDif2(ns)= sum(abs(Isamp(:) - topSamp(:)))/scaleDiff;
    end
    plot(topDif2),ylim([0 3]),pause(.01)
    toplow(i) = newrange(find(topDif2 == min(topDif2),1));
    
    %%Run bottom
%     bottomSamp = single(bottom(1:sW+1,:));
%     scaleDiff = sum(abs(bottomSamp(:)-mean(bottomSamp(:))));
%     clear bottomDif
%     startSearch = 1: 20 :size(I,1)-sW;
%     for ns = 1:length(startSearch)
%         s = startSearch(ns);
%         Isamp = single(I(end - s - sW: end - s,:));
%         bottomDif(ns)= sum(abs(Isamp(:) - bottomSamp(:)))/scaleDiff;
%         
%     end
%     lowest = startSearch(find( bottomDif == min(bottomDif),1));
%     newrange = lowest - 20 : lowest + 20;
%     clear bottomDif2
%     for ns = 1:length(newrange)
%         s = newrange(ns);
%         Isamp = single(I(end - s - sW: end - s,:));
%         bottomDif2(ns)= sum(abs(Isamp(:) - bottomSamp(:)))/scaleDiff;
% 
% %   end
%     plot(bottomDif2,'r'),Ylim([0 3]),pause(.01)
%     bottomlow(i) = newrange(find(bottomDif2 == min(bottomDif2),1));

end

plot(toplow,'b')
% hold on
% plot(bottomlow,'r')
% hold off

%% find Stretch
% 
% [ys xs] = size(BE{1,1});
% [Iys Ixs] = size(Is{i});
% topAlign0 = ys/2 - sW;
% bottomAlign0 = ys/-2 + sW; %I end point
% for i = 1:length(Is)
%         Xc = Xp(1,i); % X center
%         topYref = Yp(2,i)+topAlign0-toplow(i);
%         bottomYref = Yp(1,i) - bottomAlign0 + bottomlow(i);
%         refDif = Yp(1,i)+ bottomAlign0 - Yp(2,i)+topAlign0;
%         alignDif = Iys - bottomlow(i) - toplow(i);
%         stretchIt(i) =refDif/alignDif;
% end
% 
% plot(stretchIt)
% oddStretch = mean(stretchIt(1))
% evenStretch = mean(stretchIt(2))

%% Apply shifts

[ys xs] = size(BE{1,1});
topAlign0 = ys/2 - sW;
BufY = H + ys *2;

stripI = zeros(BufY * 3,W,'uint8');
for i = 1:length(Is)
        I = Is{i};
        
        Xc = Xp(1,i); % X center
        Yref = Yp(2,i)+topAlign0-toplow(i);
        yrange = round(Yref: Yref + size(I,1)-1);
        xrange = round(Xc +1 - size(I,2)/2:Xc+.5+size(I,2)/2);
        stripI(yrange+BufY,xrange) = stripI(yrange+BufY,xrange) + fliplr(I);
end
stripI = stripI(BufY+1:BufY + H,:);
image(stripI), pause(.01)
subplot(1,1,1)

imwrite(stripI,[TPN 'build\' 'stripI.tif'])


%% combine images

load([TPN 'dat\edgeI.mat'])
combI = edgeI;
% combI(:,:,2) = stripI;
% combI(:,:,3) = 0;
combI(edgeI==0) = stripI(edgeI == 0);
combI = fliplr(combI);
image(combI),pause(.01)
imwrite(combI,[TPN 'build\' 'combI.tif']);
[TPN ' is finished']









