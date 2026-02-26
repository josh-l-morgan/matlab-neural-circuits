colormap gray(256)
clear all
[TFN TPN] = GetMyFile

Iraw = imread([TPN TFN]);

I = Iraw(1070:1100,1000:1300);
subplot(1,2,1)
image(I)


%%
sR = 5;
shiftPix = -1*sR : sR;

Ishift = I * 0;
Islide = I(1,:);
Ishift(1,:) = Islide;
bestShift = zeros(size(I,1),1);
totalShift = 0;

Iref = double(Islide);
for y = 2:size(I,1)
    Islide = double(I(y,:));
    Islide = circshift(Islide,[0 totalShift]);
    for s = 1:length(shiftPix);
        m = shiftPix(s);
        sumProd(s) = mean(Iref.*circshift(Islide,[0 m]));
    end
    %plot(sumProd),pause(.01)
    bestShift(y) = shiftPix(find(sumProd == max(sumProd),1));
    totalShift = totalShift + bestShift(y);
    Iref = circshift(Islide,[0 bestShift(y)]);
    Ishift(y,:) = Iref;
    
    if mod(y,1000)==0
        sprintf('running %d of %d', y, size(I,1))
    end
end

%%

subplot(2,1,1)
image(I);%(500:600,4000:4060))
subplot(2,1,2)
image(Ishift);%(500:600,4000:4060))

% 
% subplot(3,1,1)
% image(I(:,sR:end-sR+1))
% subplot(3,1,2)
% image(Ishift)
% subplot(3,1,3)
% hist(bestShift)
%%
% Itop = zeros(size(I,1)/2,size(I,2),'uint8');
% 
% for i = 1:size(I,1)
%     Itop(i,:) = I(i*2,:);
%     
% end
% 
% image(Itop)
