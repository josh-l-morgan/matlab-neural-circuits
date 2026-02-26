

maxShift = 20;
[ys xs zs ] = size(I);
e = xs-maxShift;
listLines = 2*(1:fix(ys/2));

top = single(I(listLines-1,:,:));
bottom = single(I(listLines,:,:));
clear sProd
for s = 1:maxShift
    sProd = top(:,s:e+s,:) .* bottom(:,1:e+1,:);
    showProd = ( top(:,s:e+s,:) + bottom(:,1:e+1,:));
    image(uint8(showProd(1:200,1:200,:))); pause(.01)
    rProd(s) = mean(sProd(:));
    sProd = bottom(:,s:e+s,:) .* top(:,1:e+1,:);
    lProd(s) = mean(sProd(:));
end

plot(rProd,'r')
hold on
plot(lProd,'b')
hold off
% 
% if shifted == maxShift
%     maxShift = maxShift *2;
% end


rmax = max(rProd);
lmax = max(lProd);
if rmax > lmax
    shifted = find(rProd == rmax,1);

else
    shifted = find(rProd == rmax,1) * -1;
end

newI = zeros(ys,xs+abs(shifted),zs);
if shifted >0

    newI(listLines-1,1:xs,:) = top;
    newI(listLines,shifted:xs+shifted-1,:)= bottom;

else
    newI(listLines-1,abs(shifted):xs+abs(shifted)-1,:) = top;
    newI(listLines,1:xs,:)= bottom;
    
    
end
newI = uint8(newI);
image(newI)

    
