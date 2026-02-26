function[newI shifted] = findSlide(I,shifted)

maxShift = 100;
[ys xs cs zs ] = size(I);
e = xs-maxShift;
listLines = 2*(1:fix(ys/2));

top = single(I(listLines-1,:,:,:));
bottom = single(I(listLines,:,:,:));
clear I

if nargin == 1
s = 2;
sProd = top(:,s:e+s,:,:) .* bottom(:,1:e+1,:,:);
rProd = mean(sProd(:));
sProd = bottom(:,s:e+s,:,:) .* top(:,1:e+1,:,:);
lProd = mean(sProd(:));
if rProd>lProd
    sDir = 1;
else
    sDir = -1;
    holdLines = top;
    top = bottom;
    bottom = holdLines;
    clear holdLines;
end


for s = 1:maxShift
    s
    sProd = top(:,s:e+s,:,:) .* bottom(:,1:e+1,:,:);
    xProd(s) = mean(sProd(:));
    if xProd(s)<xProd(max(1,s-1)),break,end
end

%plot(xProd,'r')

shifted = find(xProd == max(xProd))*sDir;
end


newI = zeros(ys,xs+abs(shifted),cs,zs,'uint8');
if shifted >0
    newI(listLines-1,1:xs,:,:) = top;
    newI(listLines,shifted:xs+shifted-1,:,:)= bottom;

else
    newI(listLines-1,abs(shifted):xs+abs(shifted)-1,:,:) = bottom;
    newI(listLines,1:xs,:,:)= top;
end

%image(max(newI,[],4))


    
