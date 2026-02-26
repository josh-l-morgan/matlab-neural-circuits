
colormap gray(256)
fileName = 'C:\Users\jlmorgan\Documents\Projects\DPRPI\test1Image.png'
Iraw = imread(fileName);
image(Iraw)

%% Convert image into field of points
I =double(Iraw(:,:,1)*3 - Iraw(:,:,2) - Iraw(:,:,3));

ind = find(I>0);
[y x] = find(I>0);
v = I(ind);

vRes = 256;
I = I/max(I(:))*vRes;
ys = [];
xs = [];
for i = 1:vRes
    [y x] = find(I>=i);
    ys = [ys; y];
    xs = [xs; x];
end

%% principle component analysis on field of points
[coef,score,latent,tsquared] = pca([xs ys]);
R = atan(coef(2,1)/coef(2,2));
X = score(:,1);
Y = score(:,2);

Xrange = [min(X) max(X)];
Yrange = [min(Y) max(Y)];

Ir = imrotate(I,rad2deg(R),'bicubic');

xSize = round(max(X)-min(X)+10);
ySize = round(max(Y)-min(Y)+10);
xCent = size(Ir,2)/2;
yCent = size(Ir,1)/2;


image(Ir/max(Ir(:))*256)
hold on
scatter(X+xCent,Y+yCent,'r','.')
hold off



%% Image Properties

Ib = Ir>0;
Ip = bwperim(Ib);
props = regionprops(Ib,'Centroid','MajorAxisLength',...
    'MinorAxisLength','Area','Eccentricity', 'Solidity','perimeter');

Is = double(Ir) + cat(3,double(Ip*1000))

[p,S,mu] = polyfit(X,Y,2);

x1 = unique(X);
y1 = polyval(p,x1,S,mu);

ends = [x1(1) x1(end);y1(1) y1(end)];



%%
image(uint8(Ir))
hold on
%plot(ends(1,:)+ xCent,ends(2,:)+yCent,'b','lineWidth',3)
plot(Xrange+xCent,[yCent yCent],'b','lineWidth',3)
plot([xCent xCent],Yrange + yCent,'b','lineWidth',3)
plot(x1+xCent,y1+yCent,'r','linewidth',3)

hold off



















