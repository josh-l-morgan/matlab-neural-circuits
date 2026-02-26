

colormap gray(256)
SPN = 'D:\KarlsRetina\HxQ\test60cf\'
aspect = 3;

%% Get image
SPNd = dir([SPN '*.tif']);
I1 = imread([SPN SPNd(1).name]);
[ys xs cs] = size(I1);
zs = length(SPNd);
Iraw = zeros(ys, xs, zs, 'double');
for i = 1:zs
    disp(sprintf('reading %d of %d',i,zs))
   Iraw(:,:,i) = imread([SPN SPNd(i).name]); 
   image(Iraw(:,:,i)/100);
   pause(.01)
end

image(mean(Iraw,3)/100);
image(squeeze(mean(Iraw,2))/200);
    
for i = 1:zs
    image(Iraw(:,:,i)/100);
    pause(.01)
end
I = Iraw(1:100,1:100,45:95);
Isum = sum(I,3);
Isum = Isum*256/max(Isum(:));
image(Isum)

%% resample
[ys xs zs] = size(I);
I = imresize3(I,[ys xs round(zs * aspect)],'lanczos3');
I = I * 256/max(I(:));

%% Start Filter

oc = 8;
minG = .2;

h = floor(sqrt(oc+1));
w = ceil((oc+1)/h);

subplot(h,w,1)
image(mean(I,3)/100)

sigs = (minG * 2.^(1:oc));
%Is = imgaussfilt3(I,minG);
for i = 1:oc
  %Ic = Is;
  Ig{i} = imgaussfilt3(I,sigs(i));
  %Id = Ic-Is;
  %subplot(oc+1,1,i+1)
  %image(mean(Id,3)/3)
  %pause(1)
end


for i = 1:oc-2
  Id{i} = Ig{i}-Ig{i+2};
  subplot(h,w,i+1)
  image(mean(Id{i},3)*100)
  pause(1)
end



%% Something similar

vol = Id{2};
[y1 x1 z1] = find3(vol == max(vol(:)),1);


%% Crawler
buf = 2;
gRaw = gaus3d([3 3 3],1,[1 1 1]);
[y x z] = find3(gRaw>0);
shifty = y - 2;
shiftx = x -2;
shiftz = z - 2;
shiftv = gRaw(gRaw>0);

vol = Id{2};
[ys xs zs] = size(vol);
volB = zeros(ys + buf*2, xs + buf*2, zs +buf*2);
volB(buf+1:buf+ys,buf+1:buf+xs,buf+1:buf+zs) = vol;

[y1 x1 z1] = find3(volB == max(volB(:)),1);
[bys bxs bzs] = size(volB);
s = 1; %segment number
for i = 1:length(vol(:))
    val(i) = volB(y1,x1,z1);
    pos(i,:) = [y1 x1 z1];
    volB(y1,x1,z1) = 0;
    
    Y = y1 + shifty;
    X = x1 + shiftx;
    Z = z1 + shiftz;
    
    ind = sub2ind([bys bxs bzs],Y,X,Z);
    V = volB(ind);
    
    maxV = max(V);
    if maxV>0
    targ = find(V==max(V));
    y1 = Y(targ);
    x1 = X(targ);
    z1 = Z(targ);
    else
        s = s + 1;
        seg(i) = s;
        [y1 x1 z1] = find3(volB == max(volB(:)),1);
    end
    image(sum(volB,3)*10);
    hold on
    scatter(x1,y1,300,'r','.');
    hold off
    pause(.01)
end








%% 
anaBulbs4(Id{2});
%anaBulbs2(Id{1});
%anaBulbs3(Id{1},waterMorphVar);
%Ig = gpuArray(I);

%% Get Local Maxima

subplot(1,3,1)
image(mean(Id{1},3)*.3)

If2 = imgaussfilt3(Id{1},sigs(1));

subplot(1,3,2)
image(mean(If2,3)*.3)

subplot(1,3,3)

Im = imregionalmax(If2);
%Im = Im .* Id{1};%(Id{1}>0)
image(sum(Im,3)*1000)

%% Convolve with barrel

gr = 8;
gs = gr*2+1;
gRaw = gaus3d([gs gs gs] * 2,2,[1 1 .01]);
g2 = gaus3d([gs gs gs] * 2,6,[1 1 .01]);
gc = gaus3d([gs gs gs],1,[1 1 1]);
gc = gc>=gc(gr+1,gr+1,1);
% 
% gRaw = gRaw/max(gRaw(:));
% g2 = g2/max(g2(:));
gDif = gRaw * 1.2 - g2;%(gRaw-g2);

as = [0 :10 :170];
Imax = I * 0;
Isum = I * 0;
Itrack = I * 0;

aNum = length(as);
redA = [0:2/aNum:1 (1-2/aNum):-2/aNum:0];
greenA = [1:-2/aNum:0 2/aNum:2/aNum:1];
blueA = [ 1:-2/aNum:0 2/aNum:2/aNum:1];

Icol = zeros(size(I,1),size(I,2),size(I,3),3);
angInd = 0; 

for x = 1:length(as)
    for y = 1:length(as)
        angInd = angInd + 1;
        
        ax = as(x);
        ay = as(y);
        
        g = gDif;
        
        w = [1 0 0];
        g = imrotate3(g,ay,w,'cubic');
                
%         [w1 w2 w3] = sph2cart(0, ay, 1)
%         w = [w1 w2 w3];
        w = [ 0 0 1];
        g = imrotate3(g,ax,w,'cubic');
        [ny nx nz] = size(g);
        g = g(round(ny/2)-gr:round(ny/2)+gr, round(nx/2)-gr:round(nx/2)+gr,...
            round(nz/2)-gr:round(nz/2)+gr);
        g = g .* gc;
        %g = g/sum(g(:));
        
        
        gScale = 500/max(g(:));
        subplot(2,2,1)
        image(squeeze(mean(g,1)).*gScale+120);
        subplot(2,2,2)
        image(squeeze(mean(g,2))*gScale+120);
        subplot(2,2,3)
        image(squeeze(mean(g,3))*gScale+120);
        Ic = convn(I,g,'same');
                       
        Isum = Isum + Ic;
        bigger = Ic>(Imax);
        Imax(bigger) = Ic(bigger);
        Itrack(bigger) = angInd; 
        Imax = max(Imax,Ic);
        
        subplot(2,2,1)
        image(max(I,[],3))    
        subplot(2,2,4)
        image(max(Imax,[],3)*1)
        subplot(2,2,4)
        image(max(Ic,[],3)*1)
        
        colTemp = Icol(:,:,:,1);
        colTemp(bigger) = Imax(bigger) * redA(x);
        Icol(:,:,:,1) = colTemp;
        
        colTemp = Icol(:,:,:,2);
        colTemp(bigger) = Imax(bigger) * greenA(x);
        Icol(:,:,:,2) = colTemp;
        
        colTemp = Icol(:,:,:,3);
        colTemp(bigger) = Imax(bigger) * blueA(y);
        Icol(:,:,:,3) = colTemp;
        
        subplot(2,2,2)
        IcolMax = squeeze(mean(Icol,3));
        image(uint8(IcolMax)*5);
        %image(uint8(squeeze(Icol(:,:,10,:)))*3);

        pause(.01)
        
        %gSum = gSum +g;
    end
end

image(max(Imax,[],3)*6)
Imean = Isum/angInd;
anaBulbs4(Imax);

%% Itterative Thresh

colormap gray(256)
for t = 0 :256
    It = Imax>t;
    image(sum(It,3));
    pause(.1)
end
































