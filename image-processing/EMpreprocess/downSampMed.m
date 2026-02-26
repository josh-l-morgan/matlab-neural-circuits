function[medI] = downSampMed(I)

[ys xs] = size(I);
nys = floor(ys/2);
nxs = floor(xs/2);

medI = zeros(nys,nxs);
[ny nx] = find(medI==0);
y = (ny-1)*2; 
x = (nx-1)*2;


stackI = zeros(nys,nxs,4);
shiftInd = [1 2; 3 4];
for shiftY = 1:2
    for shiftX = 1:2
       takeInd = sub2ind([ys xs],y +shiftY , x  + shiftX);
       dropInd = sub2ind([nys nxs 4],ny,nx,ny*0+shiftInd(shiftY,shiftX));
       stackI(dropInd) = I(takeInd);
    end
end
% 
% tic
% ySub = cat(1,y+1,y+2,y+1,y+2);
% xSub = cat(1,x+1,x+1,x+2,x+2);
% 
% nySub = cat(1,ny,ny,ny,ny);
% nxSub = cat(1,nx,nx,nx,nx);
% nzSub = cat(1,ny*0+1,ny*0+2,ny*0+3,ny*0+4);
% 
% takeInd = sub2ind([ys xs],ySub,xSub);
% dropInd = sub2ind([nys nxs 4],nySub, nxSub, nzSub);
% stackI(dropInd) = stackI(takeInd);
% toc

medI = median(stackI,3);