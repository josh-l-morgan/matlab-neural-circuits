colormap gray(256)
clear all
[TFN TPN] = GetMyFile

Iraw = imread([TPN TFN]);

I = Iraw;%(1070:1100,1000:1300);
subplot(1,2,1)
image(I)


%%
sR = 50;
shiftPix = -1*sR : sR;

Ishift = zeros(size(I,1),size(I,2)-length(shiftPix)+1,'uint8');
Islide = I(1,:);
Ishift(1,:) = Islide(sR+1:end-sR);
bestShift = zeros(size(I,1),1);

Iref = double(Islide(sR+1:end-sR));
for y = 2:size(I,1)
    Islide = double(I(y,:));
    for s = 1:length(shiftPix);
        m = shiftPix(s);
        sumProd(s) = mean(Iref.*Islide(sR+m+1:end-sR+m));
    end
    %plot(sumProd),pause(.01)
    bestShift(y) = shiftPix(find(sumProd == max(sumProd),1));
    Iref = Islide(sR+bestShift(y)+1:end-sR+bestShift(y));
    Ishift(y,:) = Iref;
    
    if mod(y,1000)==0
        sprintf('running %d of %d', y, size(I,1))
    end
end


subplot(3,1,1)
image(I(:,sR:end-sR+1))
subplot(3,1,2)
image(Ishift)
subplot(3,1,3)
plot(bestShift)
%%
% Itop = zeros(size(I,1)/2,size(I,2),'uint8');
% 
% for i = 1:size(I,1)
%     Itop(i,:) = I(i*2,:);
%     
% end
% 
% image(Itop)
