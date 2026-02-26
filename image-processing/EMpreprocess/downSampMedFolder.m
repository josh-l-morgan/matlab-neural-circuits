function[medI] = downSampMedFolder(fold)

%%Downsample all tiles in a folder.  Must be same size

dFold = dir([fold '\Tile*.tif']);
newFold = [fold '_med'];
if ~exist(newFold,'dir'), mkdir(newFold); end
iNams = {dFold.name};

inf = imfinfo([fold '\' iNams{1}]);
ys = inf.Height;
xs = inf.Width;

%[ys xs] = size(I);
nys = floor(ys/2);
nxs = floor(xs/2);

medI = zeros(nys,nxs);
[ny nx] = find(medI==0);
y = (ny-1)*2; 
x = (nx-1)*2;
shiftInd = [1 2; 3 4];

ySub = cat(1,y+1,y+2,y+1,y+2);
xSub = cat(1,x+1,x+1,x+2,x+2);

nySub = cat(1,ny,ny,ny,ny);
nxSub = cat(1,nx,nx,nx,nx);
nzSub = cat(1,ny*0+1,ny*0+2,ny*0+3,ny*0+4);

takeInd = sub2ind([ys xs],ySub,xSub);
dropInd = sub2ind([nys nxs 4],nySub, nxSub, nzSub);


stackI = zeros(nys,nxs,4);

for i = 1:length(iNams)
    I = imread([fold '\' iNams{i}]);
    stackI(dropInd) = stackI(takeInd);    
    medI = median(stackI,3);
    imwrite(medI,[newFold '\' iNams{i}],'uint8')
end

% 
% for shiftY = 1:2
%         for shiftX = 1:2
%             takeInd = sub2ind([ys xs],y +shiftY , x  + shiftX);
%             dropInd = sub2ind([nys nxs 4],ny,nx,ny*0+shiftInd(shiftY,shiftX));
%             stackI(dropInd) = I(takeInd);
%         end
%     end





